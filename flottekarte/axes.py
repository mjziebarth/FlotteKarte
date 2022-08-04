# Plotting routines concerning the axes.
#
# Authors: Malte J. Ziebarth (ziebarth@gfz-potsdam.de)
#
# Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
#                    Malte J. Ziebarth
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
from matplotlib.patches import Polygon
from matplotlib.collections import LineCollection
from matplotlib.axes import Axes
from typing import Tuple, Optional, Union
from pyproj import Proj

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

def extract_grid_ticks_for_rectangular_axes(grid_cuts, xlim, ylim, tol):
    """
    From the grid cuts extract those ticks that coincide with
    a rectangular axes.
    """
    tick_b = [(c[0],ylim[0],c[2],c[3],180) for c in grid_cuts
              if abs(c[1] - ylim[0]) <= tol]
    tick_t = [(c[0],ylim[1],c[2],c[3],0) for c in grid_cuts
              if abs(c[1] - ylim[1]) <= tol]
    tick_l = [(xlim[0],c[1],c[2],c[3],-90) for c in grid_cuts
              if abs(c[0] - xlim[0]) <= tol]
    tick_r = [(xlim[1],c[1],c[2],c[3],90) for c in grid_cuts
              if abs(c[0] - xlim[1]) <= tol]
    return tick_b, tick_t, tick_l, tick_r


def generate_axes_grid(ax: Axes, boundary: np.ndarray,
                       boundary_angles: np.ndarray,
                       proj_str: str, linewidth: float = 0.8,
                       tick_spacing: float = 1.0,
                       which_ticks: Union[str,list] = 'auto',
                       proj: Optional[Proj] = None, fontsize: int = 8,
                       grid_cuts: Optional[list] = None,
                       rotate_labels: bool = True,
                       ticks_bot: bool = True, ticks_top: bool = True,
                       ticks_left: bool = True, ticks_right: bool = True,
                       remove_overlaps: bool = True):
    """
    Generates the axes grid.
    """
    # Compute the limits:
    xlim = float(boundary[:,0].min()), float(boundary[:,0].max())
    ylim = float(boundary[:,1].min()), float(boundary[:,1].max())

    # Boundary sanity:
    Nseg = boundary.shape[0]
    if boundary_angles.shape != (Nseg,):
        raise RuntimeError("`boundary` and `boundary_angles` do not fit "
                           "together.")

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

    # Add the axes boundary at new maximum zorder:
    ax_boundary = Polygon(boundary, facecolor='none', edgecolor='k',
                          linewidth=linewidth, zorder=zorder+10)
    hpatch = ax.add_patch(ax_boundary)

    # Set this boundary to clip path for everything plotted so far:
    for h in ax.get_children():
        if h == hpatch:
            continue
        try:
            # TODO: Here we have to check whether a clip path is already
            # set! If so, we have to combine both clips:
            h.set_clip_path(ax_boundary)
        except:
            for h1 in h:
                h.set_clip_path(ax_boundary)

    # Sanity check on the tick parameters:
    def tick_param_error(t):
        """
        Return True if there is an error with the tick parameter `t`
        or False if not.
        """
        # `t` can be one of 'auto', 'lon', 'lat', None:
        if t in ['auto','lon','lat',None]:
            return False
        # Otherwise, `t` can be a list like this: [tt0, tt1, ..., ttn],
        # where tti is a tuple(float,s) and `s` indicates the coordinate,
        # that is, `s` is either "lon" or "lat".
        if not isinstance(t, list):
            return True
        for tt in t:
            if not isinstance(tt,tuple) or len(tt) != 2:
                return True
            if not isinstance(tt[0],int) and not isinstance(tt[0],float):
                return True
            if tt[1] not in ["lon","lat"]:
                return True
        return False

    if isinstance(which_ticks,str):
        which_ticks = [which_ticks for i in range(boundary.shape[0])]

    if any(tick_param_error(t) for t in which_ticks):
        raise ValueError("Only 'auto', 'lon', 'lat', or None are valid values "
                         "for which_ticks parameter.")

    # Make sure that we can project:
    if proj is None:
        proj = Proj(proj_str)

    # Compute desired ticks based on angles:
    def desired_tick(angle):
        if (angle <= -45 and angle >= -135) or (angle >= 45 and angle <= 135):
            return 'lat'
        return 'lon'
    desired_ticks = [desired_tick(ang) for ang in boundary_angles]

    # Create a weight key for the ticks. If two ticks collide because space
    # is too small, keep the one with the higher weight. The weight key is
    # (conceptually) composed like this
    #   [is_requested, is_grid, is_desired, -- "good" tick value --]
    # Left is more important than right.
    # We compose this ordering in a binary fashion and represent it
    # as a single integer with the relevant bits set or unset.
    def generate_weight(is_requested, is_grid, is_desired, tick_value):
        weight = 128 * is_requested + 64 * is_grid + 32 * is_desired
        # Some cosmetic judgement:
        if tick_value % 2 == 0:
            weight += 1
        if tick_value % 3 == 0:
            weight += 2
        if tick_value % 5 == 0:
            weight += 4
        if tick_value % 10 == 0:
            weight += 8
        if tick_value == 0:
            weight += 16
        return weight

    def compute_tick_weights(desired, tick_values, which, is_grid):
        if isinstance(desired,list):
            # The desired ticks for this boundary segment are explicitly
            # given. We weight fully if a tick is contained and zero
            # if it is not within `desired`.
            tick_weight = [
                generate_weight(
                   is_requested = any(w == d[1] and # Type correct.
                                      d[0] == tick_spacing * t  # Val
                                      for d in desired),
                   is_grid = is_grid,
                   is_desired = False,
                   tick_value = tick_spacing * t)
                for t,w in zip(tick_values,which)
            ]
        else:
            tick_weight = [
                generate_weight(
                   is_requested = False,
                   is_grid = is_grid,
                   is_desired = (desired == w),
                   tick_value = tick_spacing * t)
                for t,w in zip(tick_values,which)
            ]
        return tick_weight

    # Tick candidates from tick spacing.
    # Here, compute both longitude and latitude ticks for all axes
    # so that we can later select the best-fitting ticks.
    tick_candidates = {}

    # Use C++ backend to compute quickly all ticks along the boundary:
    ticks_axes, which_ticks_axes, segments_ticks_axes \
       = compute_axes_ticks(proj_str, boundary, tick_spacing)

    for i in range(Nseg):
        # Find all ticks belonging to this segment:
        tick_ids = np.argwhere(segments_ticks_axes == i).reshape(-1)
        if tick_ids.size == 0:
            continue

        # Extract the ticks:
        ticks = ticks_axes[tick_ids,:]
        which_i = which_ticks_axes[tick_ids]
        which = [("lon","lat")[j] for j in which_i]

        desired = desired_ticks[i]
        # Shortcut if no ticks desired:
        if desired is None:
            continue

        # Determine the tick values (e.g. lon = 15.0)
        # in multiples of the tick spacing.
        tick_values = [int(round(t[w] / tick_spacing)) for t,w in
                       zip(ticks, which_i)]

        # Weight according to the desired ticks:
        tick_weight = compute_tick_weights(desired, tick_values, which, False)

        # Project:
        ticks_x, ticks_y = proj(*ticks.T)

        # Add the tick candidates:
        tick_candidates.update({
            (i,whch,val) : (weight,(x,y), (float(lola[0]),float(lola[1])),
                            boundary_angles[i])
            for val, weight, x, y, lola, whch in zip(tick_values, tick_weight,
                                                     ticks_x, ticks_y, ticks,
                                                     which)
        })

    # Tick candidates from grid ticks:
    if grid_cuts is not None:
        # Tolerance based on the extents:
        tol = 1e-4 * max(Dx,Dy)

        # Compute those ticks that
        ticks_btlr \
           = extract_grid_ticks_for_rectangular_axes(grid_cuts, xlim, ylim, tol)

        for i,(ticks,desired) in enumerate(zip(ticks_btlr,desired_ticks)):
            # If no ticks desired here, continue:
            if desired is None:
                continue

            # Tick values (as above):
            tick_values = [int(round(t[3] / tick_spacing)) for t in ticks]

            # Weighting:
            tick_weight = compute_tick_weights(desired, tick_values,
                                               [t[2] for t in ticks], True)

            # Lon lat:
            tick_lola = invert_proj(proj_str, [t[0] for t in ticks],
                                    [t[1] for t in ticks])

            # The normal angle of the segment at the tick location.
            # TODO: Hotfix. In the future, we should consider the orientation
            # of a potentially complex boundary path.
            if i == 0:
                angle = 180
            elif i == 1:
                angle = 0
            elif i == 2:
                angle = -90
            elif i == 3:
                angle = 90

            # Add the tick candidates:
            tick_candidates.update({
                (i,which,val) : (weight, (x,y), (lon, lat), angle)
                for val, weight, (x,y,which,_,angle), lon, lat
                                 in zip(tick_values, tick_weight, ticks,
                                        *tick_lola)
            })



    # Create an ordering of the ticks according to their weight:
    ticks_ordered = list(sorted(tick_candidates.items(),
                                key=lambda kv : kv[1][0], reverse=True))

    # Plot the ticks.
    # 1) Compute the tick marker coordinates.
    # 2) Create the labels iteratively, ensuring that no previous one
    #    is overlapped.
    # 3) Plot the remaining labels the the corresponding tick markers.
    ticks_xy = np.array([kv[1][1] for kv in ticks_ordered])
    tick_lola = np.array([kv[1][2] for kv in ticks_ordered])
    tick_normal_angle = np.array([kv[1][3] for kv in ticks_ordered])
    tick_dir_east, tick_dir_north = gradients_east_north(proj_str, *tick_lola.T)

    tick_off = np.zeros(ticks_xy.shape)
    mask_east = np.array([kv[0][1] == 'lat' for kv in ticks_ordered])
    tick_off[mask_east, :] = tick_dir_east[mask_east, :]
    mask_north = np.invert(mask_east,out=mask_east)
    tick_off[mask_north, :] = tick_dir_north[mask_north, :]

    # Ensure the right direction (pointing outward of the map)
    tna_rad = np.deg2rad(tick_normal_angle)
    outward = np.stack((np.sin(tna_rad), np.cos(tna_rad)), axis = 1)
    points_inward = (tick_off * outward).sum(axis=1) < 0
    tick_off[points_inward] *= -1

    tl = margin_ticks
    tick_off *= tl / np.linalg.norm(tick_off, axis=1)[:,np.newaxis]

    # Computing label rotation based on angle:
    def rotation_angle(angle):
        if not rotate_labels:
            return 0
        if angle == -90:
            return 90
        return ((-angle + 90.0) % 180.0) - 90.0

    # Make sure that a draw has been issued so that we can get the
    # window extents of the labels:
    ax.get_figure().canvas.draw_idle()
    transform = ax.transData.inverted()
    dot2data = np.sqrt(np.sum((ax.get_figure().dpi_scale_trans + transform)
                              .get_matrix()[:2,:2]**2)) / 72.

    # Iterate through all of the ticks, plot the labels (temporarily)
    # and decide whether to keep them:
    labels = []

    for kv, off, xy, angle, out in zip(ticks_ordered, tick_off, ticks_xy,
                                       tick_normal_angle, outward):
        # Offset due to the font size height, rotated:
        rot_off = 0.5 * fontsize * dot2data * out

        txt = ax.text(*(xy + off + rot_off),
                      tick_text(kv[0][2]*tick_spacing, which=kv[0][1]),
                      ha='center', va='center', fontsize=fontsize,
                      rotation=rotation_angle(angle))
        labels.append(txt)

    # Make sure that a draw has been issued so that we can get the
    # window extents of the labels:
    ax.get_figure().canvas.draw_idle()

    # Collect all of the window extents:
    bboxes = []
    renderer = ax.get_figure().canvas.get_renderer()
    for txt in labels:
        # Get the extent of the text:
        extent = txt.get_window_extent(renderer).transformed(transform)
        bboxes.append((float(extent.xmin), float(extent.xmax),
                       float(extent.ymin), float(extent.ymax)))

    # Check the overlap and remove one of each overlapping label pair
    # (with higher weights taking precedence):
    mask = []
    def overlaps(bb1, bb2):
        for x,y in [(bb2[0],bb2[2]), (bb2[1],bb2[2]), (bb2[1],bb2[3]),
                    (bb2[0],bb2[3])]:
            if x >= bb1[0] and x <= bb1[1] and y >= bb1[2] and y <= bb1[3]:
                return True
        for x,y in [(bb1[0],bb1[2]), (bb1[1],bb1[2]), (bb1[1],bb1[3]),
                    (bb1[0],bb1[3])]:
            if x >= bb2[0] and x <= bb2[1] and y >= bb2[2] and y <= bb2[3]:
                return True
        return False

    def out_of_bounds(bb):
        if not ticks_right and bb[1] > xlim[1]:
            return True
        if not ticks_left and bb[0] < xlim[0]:
            return True
        if not ticks_bot and bb[2] < ylim[0]:
            return True
        if not ticks_top and bb[3] > ylim[1]:
            return True
        return False

    if remove_overlaps:
        for i,bbi in enumerate(bboxes):
            if out_of_bounds(bbi):
                maski = False
                labels[i].remove()
            else:
                maski = True
                for bbj in bboxes[:i]:
                    if overlaps(bbi,bbj):
                        maski = False
                        labels[i].remove()
                        break
            # Check the overlap:
            mask.append(maski)
    else:
        mask = [True] * len(bboxes)

    # Plot the ticks:
    ticks = [[xy, xy+off] for xy,off,m in zip(ticks_xy, tick_off, mask) if m]
    ax.add_collection(LineCollection(ticks, color='k', linewidth=linewidth))
    ax.set_axis_off()
