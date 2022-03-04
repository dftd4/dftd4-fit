# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

if get_option('openmp')
  omp_dep = dependency('openmp')
  exe_deps += omp_dep
endif

mctc_dep = dependency(
  'mctc-lib',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['default_library=static'],
)
exe_deps += mctc_dep

dftd4_dep = dependency(
  'dftd4',
  fallback: ['dftd4', 'dftd4_dep'],
  default_options: ['default_library=static'],
)
exe_deps += dftd4_dep

nlopt_dep = dependency(
  'nlopt-f',
  fallback: ['nlopt-f', 'nlopt_dep'],
  default_options: ['default_library=static'],
)
exe_deps += nlopt_dep

minpack_dep = dependency(
  'minpack',
  fallback: ['minpack', 'minpack_dep'],
  default_options: ['default_library=static'],
)
exe_deps += minpack_dep
