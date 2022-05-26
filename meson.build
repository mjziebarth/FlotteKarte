project('flottekarte', 'cpp')
add_project_arguments('-O3', '-march=native',  language : 'cpp')
# Projwrapper subprocject:
libprojwrap_proj = subproject('libprojwrap')
libprojwrap_dep = libprojwrap_proj.get_variable('libprojwrap_dep')
include = [include_directories('/usr/include'), include_directories('include')]
dependencies = [dependency('proj'),dependency('openmp'),libprojwrap_dep]
sources = ['src/pyapi.cpp', 'src/invert.cpp',
           'src/types.cpp', 'src/linalg.cpp', 'src/griddedinverter.cpp',
           'src/gradient.cpp', 'src/tickfinder.cpp', 'src/grid.cpp',
           'src/paths.cpp']
library('inverse',
        sources : sources,
        include_directories : include,
        dependencies : dependencies,
        override_options : ['cpp_std=c++20'])
