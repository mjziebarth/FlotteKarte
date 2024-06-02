# Test the C++ files.

project('test_flottekarte', 'cpp')
add_project_arguments('-O3', language : 'cpp')

# Projwrapper subprocject:
libprojwrap_proj = subproject('libprojwrap')
libprojwrap_dep = libprojwrap_proj.get_variable('libprojwrap_dep')
include = [
    include_directories('/usr/include'),
    include_directories('../../../flottekarte/extensions/include')
]
dependencies = [dependency('proj'),dependency('openmp'),libprojwrap_dep]

sources_fk = files(
    '../../flottekarte/extensions/src/aabbtree.cpp',
    '../../flottekarte/extensions/src/geometry.cpp'
)

sources = ['src/test_aabb.cpp', 'src/main.cpp']
executable('test_flottekarte',
        sources : [sources_fk, sources],
        include_directories : include,
        dependencies : dependencies,
        override_options : ['cpp_std=c++20'])