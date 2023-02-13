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

project(
  'dftd4-fit',
  'fortran',
  version: '0.2.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.55',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)

# General configuration information
inc_dirs = []
exe_deps = []
subdir('config')

# Add applications
subdir('app')

# Package the license files
dftd4_lic = files(
  'COPYING',
  'COPYING.LESSER',
)

# Distribute the license files in share/licenses/<name>
install_data(
  dftd4_lic,
  install_dir: get_option('datadir')/'licenses'/meson.project_name()
)

asciidoc = find_program('asciidoctor', required: false)
if asciidoc.found()
  install_man(
    configure_file(
      command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
      input: files('man/dftd4-fit.1.adoc'),
      output: '@BASENAME@',
    )
  )
endif
