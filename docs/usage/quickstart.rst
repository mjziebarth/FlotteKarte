Quickstart
==========
Suppose you have a PROJ string that marks a projection (let this be
``proj=eqc`` in this example -- not good, but illustrative) and suppose that you
have a set of longitude and latitude coordinates that you want to plot.
For the sake of this quick introduction, create some random data:

.. code:: python

    import numpy as np
    rng = np.random.default_rng(1234)
    lon = 45.0 + 20.0 * rng.random(10)
    lat = 10.0 + 0.2 * (lon - 45.0) + 5 * rng.random(10)

A very simple scatter plot of these data using *flottekarte* could be

.. code:: python

    from pyproj import Proj
    import matplotlib.pyplot as plt
    from flottekarte import Map

    fig = plt.figure()
    ax = fig.add_subplot(111)
    proj = Proj("proj=eqc")
    mp = Map.for_data(lon, lat, "proj=eqc", ax)
    ax.scatter(*proj(lon, lat))
    mp.plot_axes(5.0)
    mp.plot_grid(5.0)


which results in a map like this:

.. figure:: quickstart-map.svg
   :scale: 13 %
   :alt: Map of the random data