# Plotting routines concerning the axes.
#
# Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
#
# Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam
#
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by
# the European Commission - subsequent versions of the EUPL (the "Licence");
# You may not use this work except in compliance with the Licence.
# You may obtain a copy of the Licence at:
#
# https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the Licence is distributed on an "AS IS" basis,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and
# limitations under the Licence.

import numpy as np
from .extensions import compute_axes_ticks, invert_proj, gradients_east_north
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection

def tick_text(tick: float, which: str = 'lon'):
    """
    Generate text for a tick.
    """
    # TODO Update this!
    if abs(int(tick) - tick) < 1e-9:
        # Integer tick.
        if which == 'lon':
            suffix = "° E" if tick > 0 else "° W"
            return str(abs(int(tick))) + suffix
        else:
            suffix = "° N" if tick > 0 else "° S"
            return str(abs(int(tick))) + suffix
    return str(tick)


def generate_axes_grid(ax, xlim, ylim, proj_str, linewidth=0.8,
                       tick_spacing=1.0, tick_bot='lon', tick_top='lon',
                       tick_left='lat', tick_right='lat', proj=None,
                       fontsize=8):
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
    ax_rect = Rectangle((xlim[0],ylim[0]), xlim[1]-xlim[0], ylim[1]-ylim[0],
                        facecolor='none', edgecolor='k', linewidth=linewidth,
                        zorder=zorder+10)
    hpatch = ax.add_patch(ax_rect)

    # Set this rect to clip path for everything plotted so far:
    for h in ax.get_children():
        if h == hpatch:
            continue
        try:
            # TODO: Here we have to check whether a clip path is already
            # set! If so, we have to combine both clips:
            h.set_clip_path(ax_rect)
        except:
            for h1 in h:
                h.set_clip_path(ax_rect)

    if any(t not in ('lon','lat',None) for t in [tick_bot, tick_top, tick_left,
                                                 tick_right]):
        raise ValueError("Only 'lon', 'lat', or None are valid values for tick_* parameters.")


    # Compute ticks:
    ticks_bot, ticks_top, ticks_left, ticks_right \
       = compute_axes_ticks(proj_str, *xlim, *ylim, tick_spacing, bot=tick_bot,
                            top=tick_top, left=tick_left, right=tick_right)

    # Plot the ticks:
    if proj is None:
        proj = Proj(proj_str)
    tick_x, tick_y = proj(*np.concatenate((ticks_bot, ticks_top, ticks_left,
                                           ticks_right)).T)

    tick_id = np.zeros(tick_x.size, dtype=int)
    i0 = len(ticks_bot)
    i1 = i0 + len(ticks_top)
    i2 = i1 + len(ticks_left)
    tick_id[i0:i1] = 1
    tick_id[i1:i2] = 2
    tick_id[i2:] = 3
    tick_lon, tick_lat = invert_proj(proj_str, tick_x, tick_y)
    tick_dir_east, tick_dir_north = gradients_east_north(proj_str, tick_lon,
                                                         tick_lat)
    # Ensure the right direction (pointing outward of the map) for all ticks:
    tick_off = np.zeros((tick_x.size,2))
    # Bottom:
    if tick_bot == 'lat':
        tick_off[:i0,:] = tick_dir_east[:i0,:]
    else:
        tick_off[:i0,:] = tick_dir_north[:i0,:]
    tick_off[:i0][tick_off[:i0,1] > 0,:] *= -1

    # Top:
    if tick_top == 'lat':
        tick_off[i0:i1,:] = tick_dir_east[i0:i1,:]
    else:
        tick_off[i0:i1,:] = tick_dir_north[i0:i1,:]
    tick_off[i0:i1][tick_off[i0:i1,1] < 0,:] *= -1

    # Left:
    if tick_left == 'lat':
        tick_off[i1:i2,:] = tick_dir_east[i1:i2,:]
    else:
        tick_off[i1:i2,:] = tick_dir_north[i1:i2,:]
    tick_off[i1:i2][tick_off[i1:i2,0] > 0,:] *= -1

    # Right
    if tick_right == 'lat':
        tick_off[i2:,:] = tick_dir_east[i2:,:]
    else:
        tick_off[i2:,:] = tick_dir_north[i2:,:]
    tick_off[i2:][tick_off[i2:,0] < 0,:] *= -1

    tl = margin_ticks
    tick_off *= tl / np.linalg.norm(tick_off, axis=1)[:,np.newaxis]

    # Plot the ticklabels:
    for x,tick in zip(tick_x[:i0] + tick_off[:i0,0], ticks_bot):
        ax.text(x, ylim[0]-margin_ticks,
                tick_text(tick[int(tick_bot == 'lat')], which=tick_bot),
                ha='center', va='top', fontsize=fontsize)
    for x,tick in zip(tick_x[i0:i1] + tick_off[i0:i1,0], ticks_top):
        ax.text(x, ylim[1]+margin_ticks,
                tick_text(tick[int(tick_top == 'lat')], which=tick_top),
                ha='center', va='bottom', fontsize=fontsize)
    for y,tick in zip(tick_y[i1:i2] + tick_off[i1:i2,1], ticks_left):
        ax.text(xlim[0]-margin_ticks, y,
                tick_text(tick[int(tick_left == 'lat')], which=tick_left),
                ha='right', va='center', fontsize=fontsize)
    for y,tick in zip(tick_y[i2:] + tick_off[i2:,1], ticks_right):
        ax.text(xlim[1]+margin_ticks, y,
                tick_text(tick[int(tick_right == 'lat')], which=tick_right),
                ha='left', va='center', fontsize=fontsize)

    # Select only a set of tick labels that does not overlap:


    ticks = [[[x,y],[x+g[0],y+g[1]]] for x,y,g in zip(tick_x, tick_y, tick_off)]
    ax.add_collection(LineCollection(ticks, color='k', linewidth=linewidth))
    ax.set_axis_off()
