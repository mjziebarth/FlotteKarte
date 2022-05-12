# Plotting routines concerning the axes.

import numpy as np
from .extensions import compute_axes_ticks, invert_proj, gradients_east_north
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection

def generate_axes_grid(ax, xlim, ylim, proj_str, linewidth=0.8,
                       tick_spacing=1.0, tick_bot='lon', tick_top='lon',
                       tick_left='lat', tick_right='lat', proj=None):
    """
    Generates the axes grid.
    """
    # Allocate space *on the Figure axis* for the map axes labels:
    Dx = xlim[1] - xlim[0]
    Dy = ylim[1] - ylim[0]
    margin_ticks = 0.02 * max(Dx,Dy)
    marginx = 0.09*Dx + margin_ticks # very elaborate formula.
    marginy = 0.04*Dx + margin_ticks# very elaborate formula.
    ax.set_xlim(xlim[0]-marginx, xlim[1]+marginx)
    ax.set_ylim(ylim[0]-marginy, ylim[1]+marginy)

    # Get maximum zorder:
    zorder=1
    for h in ax.get_children():
        try:
            zorder = max(h.get_zorder(),zorder)
        except:
            for h1 in h:
                zorder = max(h1.get_zorder(), zorder)

    # Add the axes rect at new maximum zorder:
    ax_rect = Rectangle((xlim[0],ylim[0]), xlim[1]-xlim[0], ylim[1]-ylim[0], facecolor='none',
                        edgecolor='k', linewidth=linewidth, zorder=zorder+10)
    hpatch = ax.add_patch(ax_rect)

    # Set this rect to clip path for everything plotted so far:
    for h in ax.get_children():
        if h == hpatch:
            continue
        try:
            h.set_clip_path(ax_rect)
        except:
            for h1 in h:
                h.set_clip_path(ax_rect)

    if any(t not in ('lon','lat',None) for t in [tick_bot, tick_top, tick_left, tick_right]):
        raise ValueError("Only 'lon', 'lat', or None are valid values for tick_* parameters.")


    # Compute ticks:
    ticks_bot, ticks_top, ticks_left, ticks_right \
       = compute_axes_ticks(proj_str, *xlim, *ylim, tick_spacing)

    # Plot the ticks:
    if proj is None:
        proj = Proj(proj_str)
    tick_x, tick_y = proj(*np.concatenate((ticks_bot, ticks_top, ticks_left, ticks_right)).T)

    tick_id = np.zeros(tick_x.size, dtype=int)
    i0 = len(ticks_bot)
    i1 = i0 + len(ticks_top)
    i2 = i1 + len(ticks_left)
    tick_id[i0:i1] = 1
    tick_id[i1:i2] = 2
    tick_id[i2:] = 3
    tick_lon, tick_lat = invert_proj(proj_str, tick_x, tick_y)
    tick_dir_east, tick_dir_north = gradients_east_north(proj_str, tick_lon, tick_lat)
    # Ensure the right direction (pointing outward of the map) for all ticks:
    tick_dir_north[:i0][tick_dir_north[:i0,1] > 0,:] *= -1
    tick_dir_north[i0:i1][tick_dir_north[i0:i1,1] < 0,:] *= -1
    tick_dir_east[i1:i2][tick_dir_east[i1:i2,0] > 0,:] *= -1
    tick_dir_east[i2:][tick_dir_east[i2:,0] < 0,:] *= -1
    tl = margin_ticks
    tick_dir_east *= tl / np.linalg.norm(tick_dir_east, axis=1)[:,np.newaxis]
    tick_dir_north *= tl / np.linalg.norm(tick_dir_north, axis=1)[:,np.newaxis]
    ticks = [[[x,y],[x+gn[0],y+gn[1]]] if i < 2 else [[x,y],[x+ge[0],y+ge[1]]]
             for x,y,i,ge,gn in zip(tick_x, tick_y, tick_id, tick_dir_east, tick_dir_north)]
    ax.add_collection(LineCollection(ticks, color='k', linewidth=linewidth))

    # Plot the ticklabels:
    for x,(tick,_) in zip(tick_x[:i0],ticks_bot):
        ax.text(x, ylim[0]-margin_ticks, str(tick), ha='center', va='top')
    for x,(tick,_) in zip(tick_x[i0:i1],ticks_top):
        ax.text(x, ylim[1]+margin_ticks, str(tick), ha='center', va='bottom')
    for y,(_,tick) in zip(tick_y[i1:i2],ticks_left):
        ax.text(xlim[0]-margin_ticks, y, str(tick), ha='right', va='center')
    for y,(_,tick) in zip(tick_y[i2:],ticks_right):
        ax.text(xlim[1]+margin_ticks, y, str(tick), ha='left', va='center')

    ax.set_axis_off()
