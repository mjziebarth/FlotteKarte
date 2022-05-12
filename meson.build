project('projplot', 'cpp')
add_project_arguments('-O3', '-march=native',  language : 'cpp')
include = [include_directories('/usr/include'), include_directories('include')]
dependencies = [dependency('proj'),dependency('openmp')]
sources = ['src/projwrapper.cpp', 'src/pyapi.cpp', 'src/invert.cpp',
           'src/types.cpp', 'src/linalg.cpp', 'src/griddedinverter.cpp',
           'src/gradient.cpp', 'src/tickfinder.cpp']
library('inverse',
        sources : sources,
        include_directories : include,
        dependencies : dependencies)
