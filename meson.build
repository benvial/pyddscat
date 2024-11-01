


project('pyddscatt','cpp', 'fortran',
  meson_version: '>=1.3.2',
  default_options : ['warning_level=2'],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)

subdir('ddscat')
subdir('pyddscat')
