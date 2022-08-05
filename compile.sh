#!/bin/bash
# Compilation hack.

# First check if the meson build has already been set up:
if [ ! -d builddir ]; then
    meson setup builddir
fi

# Once setup, can compile:
meson compile -C builddir

# Link to the python module:
if [ ! -f flottekarte/extensions/libflottekarte.so ]; then
    cd flottekarte/extensions
    ln -s ../../builddir/libflottekarte.so
fi
