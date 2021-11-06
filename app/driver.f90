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

!> Entry point for running single point calculations with dftd4
module dftd4_driver
   use, intrinsic :: iso_fortran_env, only : output_unit, input_unit
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type, read_structure, filetype
   use mctc_io_convert, only : autokcal
   use dftd4, only : get_dispersion, d4_model, new_d4_model, &
      realspace_cutoff, damping_param, rational_damping_param
   use dftd4_utils, only : wrap_to_central_cell
   use dftd4_cli, only : cli_config, fit_config, header
   use nlopt_wrap
   use nlopt_enum
   implicit none
   private

   public :: main

   type :: entry_type
      character(len=:), allocatable :: dir
   end type entry_type

   type :: record_type
      integer, allocatable :: coeffs(:)
      integer, allocatable :: idx(:)
      real(wp) :: reference = 0.0_wp
   end type record_type

   type :: job_type
      type(structure_type) :: mol
      type(d4_model) :: d4
   end type job_type

   type :: dataset_type
      type(record_type), allocatable :: records(:)
      type(job_type), allocatable :: jobs(:)
   end type dataset_type

contains


!> Main entry point for the driver
subroutine main(config, error)

   !> Configuration for this driver
   class(cli_config), intent(in) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call header(output_unit)

   select type(config)
   class default
      call fatal_error(error, "Unknown runtime selected")
   type is(fit_config)
      call fit_main(config, error)
   end select
end subroutine main


!> Entry point for the fit driver
subroutine fit_main(config, error)

   !> Configuration for this driver
   type(fit_config), intent(in) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(dataset_type) :: dataset
   logical :: exist
   integer :: stat
   real(wp) :: obj
   real(wp), allocatable :: x(:)
   type(nlopt_opt) :: opt

   inquire(file=config%input, exist=exist)
   if (.not.exist) error stop

   call read_dataset(config%input, dataset, config%directory, error)
   if (allocated(error)) return

   print '(a, 1x, i0)', "Number of data points", size(dataset%records)
   print '(a, 1x, i0)', "Number of evaluations", size(dataset%jobs)

   associate(ialg => algorithm_from_string(config%algorithm))
      if (ialg <= 0) then
         call fatal_error(error, "Invalid algorithm '"//config%algorithm//"'")
         return
      end if
      call create(opt, ialg, 3)
   end associate

   call opt%set_xtol_rel(1.0e-4_wp, stat)
   call opt%set_lower_bounds1(0.0_wp, stat)
   call opt%set_upper_bounds1(10.0_wp, stat)
   x = [0.6_wp, 0.5_wp, 6.0_wp]

   associate(f => nlopt_func(evaluator, dataset))
      call opt%set_min_objective(f, stat)
      call opt%optimize(x, obj, stat)
   end associate

   print '(a)', "Final parameters"
   print *, x
end subroutine fit_main


function evaluator(x, gradient, func_data) result(f)
   real(wp), intent(in) :: x(:)
   real(wp), intent(inout), optional :: gradient(:)
   class(*), intent(in), optional :: func_data
   real(wp) :: f

   class(damping_param), allocatable :: param
   type(error_type), allocatable :: error
   real(wp), parameter :: step = 1.0e-4_wp
   integer :: i
   real(wp) :: fl, fr
   real(wp), allocatable :: y(:)
   type(dataset_type) :: dataset

   if (.not.present(func_data)) error stop

   param = from_array(x)

   select type(func_data)
   type is(dataset_type)
      call single_eval(func_data, param, f, 1, error)

      if (present(gradient)) then
         y = x
         dataset = func_data
         !$omp parallel do default(shared) schedule(static) &
         !$omp shared(dataset, param, x, gradient) &
         !$omp private(i, fl, fr) firstprivate(y)
         do i = 1, size(x)
            y(i) = x(i) + step
            param = from_array(y)
            call single_eval(dataset, from_array(y), fl, 0, error)
            y(i) = x(i) - 2*step
            param = from_array(y)
            call single_eval(dataset, from_array(y), fr, 0, error)
            gradient(i) = (fl - fr)/(2*step)
            y(i) = x(i) + step
         end do
      end if
   end select
end function evaluator

function from_array(x) result(param)
   real(wp), intent(in) :: x(:)
   type(rational_damping_param) :: param

   param%s6 = 1.0_wp
   param%s8 = x(1)
   param%a1 = x(2)
   param%a2 = x(3)
   param%s9 = 1.0_wp
end function from_array


subroutine single_eval(dataset, param, obj, verbosity, error)
   type(dataset_type), intent(in) :: dataset
   class(damping_param), intent(in) :: param
   real(wp), intent(out) :: obj
   integer, intent(in) :: verbosity
   type(error_type), allocatable, intent(out) :: error

   integer :: ijob, irec
   real(wp), allocatable :: actual(:)
   real(wp), allocatable :: energy(:)

   allocate(energy(size(dataset%jobs)))
   !$omp parallel do default(none) schedule(dynamic) &
   !$omp shared(dataset, energy, param) private(ijob, error)
   do ijob = 1, size(dataset%jobs)
      call run_job(dataset%jobs(ijob), param, energy(ijob))
   end do

   obj = 0.0_wp
   allocate(actual(size(dataset%records)), source=0.0_wp)
   do ijob = 1, size(dataset%records)
      associate(record => dataset%records(ijob))
         do irec = 1, size(record%idx)
            actual(ijob) = actual(ijob) + energy(record%idx(irec)) * record%coeffs(irec)
         end do
      end associate
   end do
   obj = sum((actual - dataset%records%reference)**2)/size(dataset%records)
   if (verbosity > 0) then
      print '(*(1x, a, 1x, f10.6))', &
         & "MD", sum(actual - dataset%records%reference)/size(actual) * autokcal, &
         & "MAD", sum(abs(actual - dataset%records%reference))/size(actual) * autokcal, &
         & "RMSD", sqrt(obj) * autokcal
   end if
end subroutine single_eval


subroutine create_job(job, lentry, directory, error)
   type(job_type), intent(out) :: job
   type(entry_type), intent(in) :: lentry
   !> Working directory
   character(len=*), intent(in), optional :: directory
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: input
   real(wp) :: charge
   integer :: stat, unit
   logical :: exist

   if (present(directory)) then
      input = directory//"/"//lentry%dir//"/mol.xyz"
   else
      input = lentry%dir//"/mol.xyz"
   end if

   call read_structure(job%mol, input, error)
   if (allocated(error)) return

   if (present(directory)) then
      input = directory//"/"//lentry%dir//"/.CHRG"
   else
      input = lentry%dir//"/.CHRG"
   end if

   inquire(file=input, exist=exist)
   if (exist) then
      open(file=input, newunit=unit)
      read(unit, *, iostat=stat) charge
      if (stat == 0) then
         job%mol%charge = charge
      end if
      close(unit)
   end if
   call wrap_to_central_cell(job%mol%xyz, job%mol%lattice, job%mol%periodic)

   call new_d4_model(job%d4, job%mol)

end subroutine create_job


subroutine run_job(job, param, energy)
   type(job_type), intent(in) :: job
   class(damping_param), intent(in) :: param
   real(wp), intent(out) :: energy

   character(len=:), allocatable :: input
   type(structure_type) :: mol
   type(d4_model) :: d4
   real(wp) :: charge
   integer :: stat, unit
   logical :: exist

   call get_dispersion(job%mol, job%d4, param, realspace_cutoff(), energy)

end subroutine run_job


!> Read the dataset from the input file
subroutine read_dataset(filename, dataset, directory, error)
   
   !> File name
   character(len=*), intent(in) :: filename

   !> Data set for fitting
   type(dataset_type), intent(out) :: dataset

   !> Working directory
   character(len=*), intent(in), optional :: directory

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: unit, stat
   character(len=:), allocatable :: line
   type(record_type) :: record
   type(entry_type), allocatable :: entries(:)
   integer :: ijob

   allocate(dataset%records(0))
   
   open(file=filename, newunit=unit, iostat=stat)
   do while (stat == 0)
      call getline(unit, line, stat)
      if (stat /= 0) exit
      call read_record(line, record, entries)
      dataset%records = [dataset%records, record]
   end do
   close(unit)

   allocate(dataset%jobs(size(entries)))
   do ijob = 1, size(entries)
      call create_job(dataset%jobs(ijob), entries(ijob), directory, error)
   end do
end subroutine read_dataset


!> Read record from a line
subroutine read_record(line, record, entries)

   !> Line
   character(len=*), intent(in) :: line

   !> Record
   type(record_type), intent(out) :: record

   !> Entries
   type(entry_type), allocatable, intent(inout) :: entries(:)

   type(entry_type) :: lentry
   integer :: first, last, irec, coeff, stat

   allocate(record%idx(0), record%coeffs(0))
   coeff = -1
   first = 1
   last = 0
   do
      last = index(line(first:), ',') + first - 2
      if (last < first) then
         last = len(line)
         exit
      end if

      !print '(*(a))', line, new_line('a'), repeat(' ', first-1), repeat('=', last-first+1)

      lentry%dir = trim(adjustl(line(first:last)))
      call push_back(entries, lentry)

      record%idx = [record%idx, find(entries, lentry%dir)]
      record%coeffs = [record%coeffs, coeff]

      first = last + 2
      coeff = 1
   end do

   read(line(first:), *, iostat=stat) record%reference
end subroutine read_record


subroutine push_back(entries, lentry)
   !> Entries
   type(entry_type), allocatable, intent(inout) :: entries(:)
   !> Entry
   type(entry_type), intent(in) :: lentry

   integer :: pos

   if (.not.allocated(entries)) allocate(entries(0))

   pos = find(entries, lentry%dir)
   if (pos == 0) entries = [entries, lentry]
end subroutine push_back


function find(entries, dir) result(pos)
   !> Entries
   type(entry_type), allocatable, intent(inout) :: entries(:)
   !> Directory
   character(len=*), intent(in) :: dir
   !> Position
   integer :: pos

   integer :: i

   pos = 0
   do i = 1, size(entries)
      if (entries(i)%dir == dir) then
         pos = i
         exit
      end if
   end do
end function find


!> Consume a whole line from a formatted unit
subroutine getline(unit, line, iostat, iomsg)
   !> Formatted IO unit
   integer, intent(in) :: unit
   !> Line to read
   character(len=:), allocatable, intent(out) :: line
   !> Status of operation
   integer, intent(out) :: iostat
   !> Error message
   character(len=:), allocatable, optional :: iomsg

   integer, parameter :: bufsize = 512
   character(len=bufsize) :: buffer, msg
   integer :: size, stat
   intrinsic :: is_iostat_eor, present, trim

   allocate(character(len=0) :: line)
   do
      read(unit, '(a)', advance='no', iostat=stat, iomsg=msg, size=size) &
         & buffer
      if (stat > 0) exit
      line = line // buffer(:size)
      if (stat < 0) exit
   end do

   if (is_iostat_eor(stat)) stat = 0

   if (stat /= 0) then
      if (present(iomsg)) iomsg = trim(msg)
   end if
   iostat = stat

end subroutine getline


end module dftd4_driver
