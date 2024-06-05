# FlotteKarte meson build file.
#
# Authors: Malte J. Ziebarth (malte.ziebarth@tum.de)
#
# Copyright (C) 2022 Deutsches GeoForschungsZentrum Potsdam,
#               2024 Technische Universität München
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

project('flottekarte', 'cpp')
add_project_arguments('-O3',  language : 'cpp')

# Use natively tuned code by default - turn off by use of the 'portable'
# option:
if not get_option('portable')
    add_project_arguments(['-march=native','-mtune=native'], language : 'cpp')
endif

# Projwrapper subprocject:
libprojwrap_proj = subproject('libprojwrap')
libprojwrap_dep = libprojwrap_proj.get_variable('libprojwrap_dep')
include = [include_directories('/usr/include'), include_directories('include')]
dependencies = [dependency('proj'),dependency('openmp'),libprojwrap_dep]
sources = ['src/pyapi.cpp', 'src/invert.cpp',
           'src/types.cpp', 'src/linalg.cpp', 'src/griddedinverter.cpp',
           'src/gradient.cpp', 'src/tickfinder.cpp', 'src/grid.cpp',
           'src/paths.cpp', 'src/boundary.cpp', 'src/augmentedproj.cpp',
           'src/interpolate.cpp', 'src/streamlines.cpp', 'src/azimuth.cpp']
library('flottekarte',
        sources : sources,
        include_directories : include,
        dependencies : dependencies,
        override_options : ['cpp_std=c++20'])
