Installation
============
Currently, installation is a bit hacky. It uses the **Meson Build system** to
compile the C++ binary extension (loaded with `ctypes`). Then, a `pip` build
using `setuptools` `pyproj.toml` can be performed and will copy the compiled
binary to the right location.

In short, perform the following two steps in the project's main directory:

.. code:: Bash

    ./compile.sh
    pip install --user .

That should work (hopefully).