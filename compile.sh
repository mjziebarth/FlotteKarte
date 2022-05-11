#!/bin/bash
# Compilation hack.

meson compile -C builddir

# Link to the python module:
if [ ! -f projplot/extensions/libinverse.so ]; then
    cd projplot/extensions
    ln -s ../../builddir/libinverse.so
fi
