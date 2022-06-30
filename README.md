flottekarte
=======

**flottekarte** is a Python library for quick and versatile cartography
based on PROJ4-string syntax and using Matplotlib, NumPy, and PyPROJ under
the hood.

## Installing
The install system is based on the Meson build system and Python setuptools.
A typical install might be performed by executing the following two commands
(of which the second should suffice on its own) within this project directory:
```sh
./compile.sh
pip install --user
```

### Requirements
The following software has to be installed:
 - PROJ
 - OpenMP
 - NumPy
 - Matplotlib
 - PyProj
 - SciPy

The following software will be automatically downloaded during Meson installation:
 - [ProjWrapCpp](https://github.com/mjziebarth/ProjWrapCpp)

## Usage
FlotteKarte is a low-overhead plotting routine. The conceptual idea behind this package
is that a map is fully defined through the 2D cartesian coordinates that result from applying the
map projection to different geographical data. For displaying data on a two-dimensional
canvas, Matplotlib is a powerful tool. Conversion between geographic and projected
coordinates can easily be done using PyProj. The gap between these two powerful tools
and a polished map lies in potential difficulties when translating spherical line topology
to 2D cartesian space, and by introducing typical map decorations such as grids or ticks.
FlotteKarte aims to fill this gap with a simple interface.

FlotteKarte's philosophy is to work completely within the 2D projected coordinates,
that is, very close to the projected data. If projected coordinates of data can be 
obtained, the data can be drawn directly on the underlying Matplotlib Axes. The
`Map` class can then be used to add typical map decoration to that axes using information
that it derives from the numerics of the PROJ projection.

The basic usage pattern is:
```python
# Given:
#   lon:      some longitudes
#   lat:      some latitudes of equal size
#   proj_str: A PROJ4 string

import matplotlib.pyplot as plt
from flottekarte import Map
from pyproj import Proj
fig = plt.figure()
ax = fig.add_subplot(111)

# Project the data:
x,y = Proj(proj_str)(lon,lat)

# Create a map on axis 'ax' from the PROJ string.
# Alternative:
#    mp = Map(proj_str, ax, xlim, ylim)
# where xlim and ylim would be determined through x and y.
mp = Map.for_data(lon, lat, proj_str, ax)

# Scatter the data:
ax.scatter(x,y)
# ... and maybe some other here.

# Plot a grid with 5 degree spacing:
mp.plot_grid(5)

# Complete the plot by this call (ticks every 5 degrees).
mp.plot_axes(5)
```

## License
This software is licensed under the European Public License (EUPL) version 1.2 or later.

## Changelog
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### [0.1.0] - 2022-06-30
#### Added
 - First version
