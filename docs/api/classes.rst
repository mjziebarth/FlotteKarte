==========
Python API
==========

The :class:`Map` class
----------------------

FlotteKarte is designed for the use of the :class:`.Map` class, which
provides a simple interface for converting a Matplotlib
:class:`matplotlib.axes.Axes` instance into a map by means of a PROJ
string.

The :class:`.Map` class does not offer a large number of methods but a
small number of cartography-related functionality. All plotting
functionality can be performed normally on the
:class:`matplotlib.axes.Axes` instance in projected coordinates.
A call to :meth:`.Map.plot_axes` adds coordinate ticks and map boundary,
turning the two-dimensional canvas into a map.

.. autoclass:: flottekarte.Map
   :members:



The :class:`GeoJSON` class
--------------------------

This class allows to load data from a GeoJSON and plot it using
:meth:`.Map.add_data`.

.. autoclass:: flottekarte.GeoJSON
   :members: