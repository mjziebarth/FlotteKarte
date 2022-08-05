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
    cd ../..
fi

# Link the C++ source code into the Python package tree
# so as to be able to potentially recompile the extension
# from within the package:
if [ ! -d flottekarte/extensions/include/ ]; then
    cd flottekarte/extensions/
    ln -s ../../include
    cd ../..
fi
if [ ! -d flottekarte/extensions/src/ ]; then
    cd flottekarte/extensions/
    ln -s ../../src
    cd ../..
fi
if [ ! -f flottekarte/extensions/meson.build ]; then
    cd flottekarte/extensions
    ln -s ../../meson.build
    cd ../..
fi