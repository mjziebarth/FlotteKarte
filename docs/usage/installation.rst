Installation
============
Currently, installation is a bit hacky. It uses the **Meson Build system** to
compile the C++ binary extension (loaded with `ctypes`). This command is invoked
whenever the backend binary library fails to be imported at runtime. Hence, the
library is compiled at first run after the (Python-only) install.

In short, perform the following command in the project's main directory:

.. code:: Bash

    pip install --user .

That should work (hopefully). A working compiler and the Meson build system
are required to run the package afterwards. The local install is critical here:
if the Python process loading the library does not have write access to the
package directory, the backend cannot be compiled.