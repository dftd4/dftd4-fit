! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Definition of the command line interface
module dftd4_cli
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : get_filetype
   use dftd4, only : rational_damping_param, get_dftd4_version
   use dftd4_argument, only : argument_list, len
   implicit none
   private

   public :: cli_config, fit_config, get_arguments
   public :: prog_name, header

   !> The name of the program
   character(len=*), parameter :: prog_name = "dftd4-fit"

   !> Base command line configuration
   type, abstract :: cli_config
   end type cli_config

   !> Configuration for parameter optimization
   type, extends(cli_config) :: fit_config
      !> Input file for data set
      character(len=:), allocatable :: input
      !> Basename of the input file
      character(len=:), allocatable :: basename
      !> Verbosity of the output
      integer :: verbosity = 2
      !> Selected optimizer
      character(len=:), allocatable :: optimizer
      !> Selected algorithm
      character(len=:), allocatable :: algorithm
      !> Directory containing the data set
      character(len=:), allocatable :: directory
      !> Tolerance for convergence of error in parameters
      real(wp), allocatable :: xtol
      !> Tolerance for convergence of error on data set
      real(wp), allocatable :: ftol
      !> Starting parameters
      real(wp), allocatable :: x(:)
   end type fit_config

contains


subroutine get_arguments(config, error)

   !> Configuation data
   class(cli_config), allocatable, intent(out) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(argument_list) :: list
   integer :: iarg

   list = argument_list()
   iarg = 0

   block
      type(fit_config), allocatable :: fit
      allocate(fit)
      call get_fit_arguments(fit, list, iarg, error)
      call move_alloc(fit, config)
   end block
end subroutine get_arguments


!> Read configuration for the fit driver
subroutine get_fit_arguments(config, list, start, error)

   !> Configuation data
   type(fit_config), intent(out) :: config

   !> List of command line arguments
   type(argument_list), intent(in) :: list

   !> First command line argument
   integer, intent(in) :: start

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg

   iarg = start
   getopts = .true.
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
      if (.not.getopts) then
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end if
      select case(arg)
      case("--")
         getopts = .false.
      case("-h", "--help")
         call help(output_unit)
         stop
      case("--version")
         call version(output_unit)
         stop
      case("--license")
         call license(output_unit)
         stop
      case("--citation")
         call citation(output_unit)
         stop
      case("-v", "--verbose")
         config%verbosity = config%verbosity + 1
      case("-s", "--silent")
         config%verbosity = config%verbosity - 1
      case("-f", "--file")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for directory")
            exit
         end if
         call move_alloc(arg, config%basename)
      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         if (arg(1:1) == "-") then
            call fatal_error(error, "Unknown command option '"//arg//"'")
         else
            call fatal_error(error, "Too many positional arguments present")
         end if
         exit
      case("-C", "--directory")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for directory")
            exit
         end if
         call move_alloc(arg, config%directory)
      case("--algorithm")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for algorithm")
            exit
         end if
         call move_alloc(arg, config%algorithm)
      case("--optimizer")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for optimizer")
            exit
         end if
         call move_alloc(arg, config%optimizer)
      case("--ftol")
         allocate(config%ftol)
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%ftol, error)
         if (allocated(error)) exit
      case("--xtol")
         allocate(config%xtol)
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%xtol, error)
         if (allocated(error)) exit
      case("--param")
         allocate(config%x(3))
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%x(1), error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%x(2), error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%x(3), error)
         if (allocated(error)) exit
      end select
   end do
   if (allocated(error)) return

   if (.not.allocated(config%optimizer)) then
      config%optimizer = "MINPACK"
   end if

   if (.not.allocated(config%algorithm)) then
      config%algorithm = "LN_NEWUOA_BOUND"
   end if

   if (.not.allocated(config%x)) then
      config%x = [1.0_wp, 0.4_wp, 5.0_wp]
   end if

   if (.not.allocated(config%basename)) then
      config%basename = "mol.xyz"
   end if

   if (.not.any([allocated(config%xtol), allocated(config%ftol)])) then
      config%xtol = 1.0e-4_wp
   end if

   if (.not.allocated(config%input)) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_fit_arguments


subroutine get_argument_as_real(arg, val, error)

   !> Index of command line argument, range [0:command_argument_count()]
   character(len=:), intent(in), allocatable :: arg

   !> Real value
   real(wp), intent(out) :: val

   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat

   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read real value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read real value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_real


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Optimization driver for damping parameters in DFT-D4 using NLOpt or MINPACK.", &
      "As input a data set is required containing the missing dispersion energies", &
      "for the respective functional in Hartree.", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-v, --verbose", "Show more, can be used multiple times", &
      "-s, --silent", "Show less, use twice to supress all output", &
      "    --optimizer <str>", "Name of the library used for optimization", &
      "", "options: MINPACK (default) and NLOPT", &
      "    --algorithm <str>", "Name of the algorithm used for optimization", &
      "", "available algorithms for NLOPT can be found at https://nlopt.rtfd.io", &
      "-C, --directory <dir>", "Directory containing data set", &
      "    --ftol <real>", "Tolerance for convergence of error on data set", &
      "    --xtol <real>", "Tolerance for convergence of parameters", &
      "    --param <real>...", "Initial parameters (s8, a1, a2)", &
      "    --version", "Print program version and exit", &
      "    --citation", "Print citation information and exit", &
      "    --license", "Print license header and exit", &
      "-h, --help", "Show this help message"

   write(unit, '(a)')

end subroutine help

subroutine header(unit)
   integer, intent(in) :: unit

   write(unit,'(a)') &
      "===================================",&
      "       D F T - D 4   F I T",&
      "===================================", ""
end subroutine header


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_dftd4_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


subroutine citation(unit)
   integer, intent(in) :: unit

   write(unit, '(a)') &
      "Please include the appropriate citations when using DFTD4 in your work.", &
      "", &
      "Original DFTD4 idea:", &
      "Eike Caldeweyher, Christoph Bannwarth and Stefan Grimme,", &
      "J. Chem. Phys., 2017, 147, 034112.", &
      "DOI: 10.1063/1.4993215", &
      "", &
      "DFTD4 model:", &
      "Eike Caldeweyher, Sebastian Ehlert, Andreas Hansen, Hagen Neugebauer,", &
      "Sebastian Spicher, Christoph Bannwarth and Stefan Grimme,", &
      "J. Chem Phys, 2019, 150, 154122.", &
      "DOI: 10.1063/1.5090222", &
      "ChemRxiv: 10.26434/chemrxiv.7430216.v2", &
      "", &
      "Periodic DFTD4 model:", &
      "Eike Caldeweyher, Jan-Michael Mewes, Sebastian Ehlert", &
      "and Stefan Grimme, Phys. Chem. Chem. Phys., 2020, 22, 8499-8512.", &
      "DOI: 10.1039/D0CP00502A", &
      "ChemRxiv: 10.26434/chemrxiv.10299428.v1", &
      ""

end subroutine citation


subroutine license(unit)
   integer, intent(in) :: unit

   write(unit, '(a)') &
      "dftd4 is free software: you can redistribute it and/or modify it under", &
      "the terms of the Lesser GNU General Public License as published by", &
      "the Free Software Foundation, either version 3 of the License, or", &
      "(at your option) any later version.", &
      "", &
      "dftd4 is distributed in the hope that it will be useful,", &
      "but WITHOUT ANY WARRANTY; without even the implied warranty of", &
      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the", &
      "Lesser GNU General Public License for more details.", &
      "", &
      "You should have received a copy of the Lesser GNU General Public License", &
      "along with dftd4.  If not, see <https://www.gnu.org/licenses/>."
end subroutine license


end module dftd4_cli
