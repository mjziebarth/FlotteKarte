#!/bin/bash
# Compilation hack.

meson compile -C builddir

# Link to the python module:
if [ ! -f flottekarte/extensions/libinverse.so ]; then
    cd flottekarte/extensions
    ln -s ../../builddir/libinverse.so
fi
