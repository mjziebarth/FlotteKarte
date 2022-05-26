project('libprojwrap', 'cpp')

inc = [include_directories('include')]

dependencies = [dependency('proj')]

sources = ['src/projwrapper.cpp','src/projwrappertypes.cpp']

libprojwrap = static_library('projwrap',
    sources : sources,
    include_directories : inc,
    dependencies : dependencies,
    install : true)

libprojwrap_dep = declare_dependency(include_directories : inc,
    link_with : libprojwrap)