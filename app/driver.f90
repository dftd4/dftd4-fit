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

!> Entry point for fitting damping parameters with dftd4
module dftd4_driver
   use, intrinsic :: iso_fortran_env, only: output_unit, input_unit
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type, read_structure, filetype
   use mctc_io_convert, only: autokcal
   use dftd4, only: get_dispersion, d4_model, new_d4_model, &
      realspace_cutoff, damping_param, rational_damping_param
   use dftd4_utils, only: wrap_to_central_cell
   use dftd4_cli, only: cli_config, fit_config, header
   use nlopt_wrap
   use nlopt_enum
   use minpack_module, only: lmder, lmdif
   implicit none
   private

   public :: main

   !> An entry in the data set
   type :: entry_type
      !> Directory containing the structure input
      character(len=:), allocatable :: dir
   end type entry_type

   !> Record in the data set
   type :: record_type
      !> Coefficients for the individual entries
      integer, allocatable :: coeffs(:)
      !> Indices of the entries creatin the record
      integer, allocatable :: idx(:)
      !> Reference to optimize against
      real(wp) :: reference = 0.0_wp
   end type record_type

   !> Working information for evaluating an entry
   type :: job_type
      !> Molecular structure data
      type(structure_type) :: mol
      !> Dispersion model
      type(d4_model) :: d4
   end type job_type

   !> Complete data set for the optimization
   type :: dataset_type
      !> All records forming the data set
      type(record_type), allocatable :: records(:)
      !> Jobs representing the entries in the records
      type(job_type), allocatable :: jobs(:)
      !> Unit for output
      integer :: io = output_unit
   end type dataset_type

contains
   !> Main entry point for the driver
   subroutine main(config, error)
      !> Configuration for this driver
      class(cli_config), intent(in) :: config
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      call header(output_unit)

      select type (config)
      class default
         call fatal_error(error, "Unknown runtime selected")
      type is (fit_config)
         call fit_main(output_unit, config, error)
      end select
   end subroutine main

   !> Entry point for the fit driver
   subroutine fit_main(io, config, error)
      !> Unit for output
      integer, intent(in) :: io
      !> Configuration for this driver
      type(fit_config), intent(in) :: config
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      ! Complete data set for the optimization
      type(dataset_type) :: dataset
      logical :: exist
      real(wp), allocatable :: x(:)

      inquire (file=config%input, exist=exist)
      if (.not. exist) then
         call fatal_error(error, "Data file '"//config%input//" not found")
         return
      end if

      call read_dataset(config%input, config%format, config%basename, dataset, config%directory, error)
      if (allocated(error)) return
      dataset%io = io

      write (io, '(a, 1x, i0)') "Number of data points", size(dataset%records)
      write (io, '(a, 1x, i0)') "Number of evaluations", size(dataset%jobs)

      if (config%optimizer == "MINPACK") then
         call run_minpack(dataset, config, x, error)
      else if (config%optimizer == "NLOPT") then
         call run_nlopt(dataset, config, x, error)
      else
         call fatal_error(error, "Invalid optimizer '"//config%optimizer//"'")
         return
      end if

      write (io, '(/, a)') "Final parameters"
      write (io, *) x
   end subroutine fit_main

   !> Wrapper for MINPACK optimizer
   subroutine run_minpack(dataset, config, x, error)
      !> Configuration for this driver
      type(fit_config), intent(in) :: config
      !> Complete data set for the optimization
      type(dataset_type), intent(in) :: dataset
      !> Independant variable vector
      real(wp), allocatable, intent(out) :: x(:)
      !> Error handling
      type(error_type), allocatable, intent(inout) :: error

      !> Number of datapoints/functions
      integer :: m
      !> Number of variables (n must not exceed m)
      integer :: n
      !> Termination criterium
      real(wp) :: tol
      !> Value of function at x
      real(wp), allocatable :: fvec(:)
      !> Jacobian matrix
      real(wp), allocatable :: fjac(:, :)
      !> Work array for minpack
      real(wp), allocatable :: diag(:), qtf(:), wa1(:), wa2(:), wa3(:), wa4(:)
      !> Integer work array for minpack of length n
      integer, allocatable :: iwa(:)

      ! Unit for output
      integer :: io
      ! Exit parameter
      integer :: info
      integer :: nfev, njev, maxfev, mode
      real(wp), parameter :: factor = 1.0_wp
      character(len=2) tmp

      intrinsic :: norm2

      ! initial guess for parameters
      x = config%x
      ! output unit
      io = dataset%io

      n = size(x)
      m = size(dataset%records)

      ! Set "tol" to the square root of the machine precision. Unless high
      ! precision solutions are required, this is the recommended setting.
      tol = sqrt(epsilon(1._wp))

      allocate (fvec(m))
      allocate (diag(n), qtf(n), wa1(n), wa2(n), wa3(n), wa4(m))
      allocate (iwa(n))
      allocate (fjac(m, n))

      write (io, '(/, a)') 'Optimizing (using MINPACK) ...'
      diag(:) = 1.0_wp
      nfev = 0
      njev = 0

      if (config%algorithm == "lmdif") then
         maxfev = 200 * (n + 1)
         mode = 2
         call lmdif(single_eval, m, n, x, fvec, tol, tol, 0.0_wp, maxfev, 0.0_wp, diag, &
            & mode, factor, config%verbosity, info, nfev, fjac, m, iwa, qtf, &
            & wa1, wa2, wa3, wa4)
      else
         maxfev = 100 * (n + 1)
         mode = 1
         call lmder(fcn_lmder, m, n, x, fvec, fjac, m, tol, tol, 0.0_wp, maxfev,   &
            & diag, mode, factor, config%verbosity, info, nfev, njev, iwa, qtf, &
            & wa1, wa2, wa3, wa4)
      end if

      if (info < 0) then
         write (tmp, '(I1)') info
         call fatal_error(error, "Optimization failed with "//tmp)
      else
         write (tmp, '(I1)') info
         write (io, '(/, a)') "Optimization terminated with status "//tmp
         write (io, '(a, f15.7)') 'Final L2 norm of the residuals:', norm2(fvec)
      end if

   contains
      !> User-supplied subroutine which calculates the functions
      subroutine fcn_lmder(m, n, x, fvec, fjac, ldfjac, iflag)
         !> Number of datapoints/functions
         integer, intent(in) :: m
         !> Number of variables (n must not exceed m)
         integer, intent(in) :: n
         !> leading dimension of the array fjac.
         integer, intent(in) :: ldfjac
         !> Independant variable vector
         real(wp), intent(in) :: x(n)
         !> Value of function at x
         real(wp), intent(inout) :: fvec(m)
         !> jacobian matrix at `x`
         real(wp), intent(inout) :: fjac(ldfjac, n)
         !> Status
         integer, intent(inout) :: iflag

         integer :: ix, info
         real(wp), parameter :: step = 1.0e-3_wp
         real(wp), allocatable :: fl(:), fr(:), y(:)

         if (iflag == 0) then
            write (io, '(*(4x, a, 1x, f7.4))') &
               & "MD", sum(fvec)/size(fvec)*autokcal, &
               & "MAD", sum(abs(fvec))/size(fvec)*autokcal, &
               & "RMSD", sqrt(sum((fvec)**2)/size(fvec))*autokcal, &
               & "SD", sqrt(sum(abs(fvec - sum(fvec)/size(fvec))**2)/size(fvec))*autokcal
            return
         end if

         if (iflag == 1) then
            call single_eval(m, n, x, fvec, iflag)
         else
            info = 1
            y = x
            fl = fvec
            fr = fvec
            do ix = 1, n
               y(ix) = x(ix) + step
               call single_eval(m, n, y, fl, info)
               y(ix) = x(ix) - step
               call single_eval(m, n, y, fr, info)
               y(ix) = x(ix)
               fjac(:, ix) = (fl - fr) / (2 * step)
            end do
         end if
      end subroutine fcn_lmder

      !> User-supplied subroutine which calculates the functions
      subroutine single_eval(m, n, x, fvec, iflag)
         !> Number of datapoints/functions
         integer, intent(in) :: m
         !> Number of variables (n must not exceed m)
         integer, intent(in) :: n
         !> Independant variable vector
         real(wp), intent(in) :: x(n)
         !> Value of function at x
         real(wp), intent(out) :: fvec(m)
         !> Status
         integer, intent(inout) :: iflag

         integer :: ijob, irec
         real(wp), allocatable :: actual(:)
         real(wp), allocatable :: energy(:)

         class(damping_param), allocatable :: param

         if (iflag == 0) then
            write (io, '(*(4x, a, 1x, f7.4))') &
               & "MD", sum(fvec)/size(fvec)*autokcal, &
               & "MAD", sum(abs(fvec))/size(fvec)*autokcal, &
               & "RMSD", sqrt(sum((fvec)**2)/size(fvec))*autokcal, &
               & "SD", sqrt(sum(abs(fvec - sum(fvec)/size(fvec))**2)/size(fvec))*autokcal
            return
         end if

         param = from_array(x)

         allocate (energy(size(dataset%jobs)))
         !$omp parallel do default(none) schedule(dynamic) &
         !$omp shared(dataset, energy, param) private(ijob)
         do ijob = 1, size(dataset%jobs)
            call run_job(dataset%jobs(ijob), param, energy(ijob))
         end do

         allocate (actual(size(dataset%records)), source=0.0_wp)
         do ijob = 1, size(dataset%records)
            associate (record => dataset%records(ijob))
               do irec = 1, size(record%idx)
                  actual(ijob) = actual(ijob) + energy(record%idx(irec))*record%coeffs(irec)
               end do
            end associate
         end do

         fvec(:) = dataset%records%reference - actual
      end subroutine single_eval
   end subroutine run_minpack

   !> Wrapper for NLOPT optimizer
   subroutine run_nlopt(dataset, config, x, error)
      !> Configuration for this driver
      type(fit_config), intent(in) :: config
      !> Complete data set for the optimization
      type(dataset_type), intent(in) :: dataset
      !> Independant variable vector
      real(wp), allocatable, intent(out) :: x(:)
      !> Error handling
      type(error_type), allocatable, intent(inout) :: error

      !> Unit for output
      integer :: io, stat
      real(wp) :: obj
      type(nlopt_opt) :: opt

      ! initial guess for parameters
      x = config%x
      ! output unit
      io = dataset%io

      associate (ialg => algorithm_from_string(config%algorithm))
         if (ialg <= 0) then
            call fatal_error(error, "Invalid algorithm '"//config%algorithm//"'")
            return
         end if
         call create(opt, ialg, size(x))
      end associate

      write (io, '(/, a)') "Optimizing (using NLOPT with '"//config%algorithm//"' algorithm) ..."

      if (allocated(config%xtol)) then
         call opt%set_xtol_rel(config%xtol)
      end if
      if (allocated(config%ftol)) then
         call opt%set_ftol_rel(config%ftol)
      end if

      if (allocated(config%bound_lower)) then
         call opt%set_lower_bounds1(config%bound_lower)
      end if
      if (allocated(config%bound_upper)) then
         call opt%set_upper_bounds1(config%bound_upper)
      end if

      associate (f => nlopt_func(evaluator, dataset))
         call opt%set_min_objective(f)
         call opt%optimize(x, obj, stat)
      end associate
      if (stat < 0) then
         call fatal_error(error, "Optimization failed with "//result_to_string(stat))
      else
         write (io, '(/, a)') "Optimization terminated with status "//result_to_string(stat)
      end if
   end subroutine run_nlopt

   !> Call back for evaulation of the data set
   function evaluator(x, gradient, func_data) result(f)
      !> Current parameters
      real(wp), intent(in) :: x(:)
      !> Gradient of objective function
      real(wp), intent(inout), optional :: gradient(:)
      !> Data set for evaluation
      class(*), intent(in), optional :: func_data
      !> Objective function
      real(wp) :: f

      class(damping_param), allocatable :: param

      if (.not. present(func_data)) &
         error stop "Implementation error: data pointer not passed to callback"

      param = from_array(x)

      select type (func_data)
      type is (dataset_type)
         call single_eval(func_data, param, f, 1)

         if (present(gradient)) then
            call grad_eval(func_data, x, gradient)
         end if
      class default
         error stop "Implementation error: data set not available in callback"
      end select
   end function evaluator

   !> Create damping parameters from parameter array
   function from_array(x) result(param)
      !> Current parameter set
      real(wp), intent(in) :: x(:)
      !> Damping parameter object
      type(rational_damping_param) :: param

      param%s6 = 1.0_wp
      param%s8 = x(1)
      param%a1 = x(2)
      param%a2 = x(3)
      param%s9 = 1.0_wp
   end function from_array

   !> Evaluate objective function
   subroutine single_eval(dataset, param, obj, verbosity)
      !> Data set for optimization
      type(dataset_type), intent(in) :: dataset
      !> Current damping parameters
      class(damping_param), intent(in) :: param
      !> Objective function
      real(wp), intent(out) :: obj
      !> Verbosity of printout
      integer, intent(in) :: verbosity

      integer :: ijob, irec
      real(wp), allocatable :: actual(:)
      real(wp), allocatable :: energy(:)
      real(wp), allocatable :: err(:)

      allocate (energy(size(dataset%jobs)))
      !$omp parallel do default(none) schedule(dynamic) &
      !$omp shared(dataset, energy, param) private(ijob)
      do ijob = 1, size(dataset%jobs)
         call run_job(dataset%jobs(ijob), param, energy(ijob))
      end do

      obj = 0.0_wp
      allocate (actual(size(dataset%records)), source=0.0_wp)
      do ijob = 1, size(dataset%records)
         associate (record => dataset%records(ijob))
            do irec = 1, size(record%idx)
               actual(ijob) = actual(ijob) + energy(record%idx(irec))*record%coeffs(irec)
            end do
         end associate
      end do
      obj = sum((actual - dataset%records%reference)**2)/size(dataset%records)

      if (verbosity > 0) then
         allocate (err(size(dataset%records)))
         err = actual - dataset%records%reference

         write (dataset%io, '(*(4x, a, 1x, f7.4))') &
          & "MD", sum(err)/size(err)*autokcal, &
          & "MAD", sum(abs(err))/size(err)*autokcal, &
          & "RMSD", sqrt(obj)*autokcal, &
          & "SD", sqrt(sum(abs(err - sum(err)/size(err))**2)/size(err))*autokcal
      end if
   end subroutine single_eval

   !> Evaluate gradient of objective function by numerical differentation
   subroutine grad_eval(dataset, x, gradient)
      !> Data set for optimization
      type(dataset_type), intent(in) :: dataset
      !> Current parameter vector
      real(wp), intent(in) :: x(:)
      !> Gradient of objective function
      real(wp), intent(inout) :: gradient(:)

      real(wp), parameter :: step = 1.0e-4_wp
      integer :: i
      real(wp) :: fl, fr
      real(wp), allocatable :: y(:)
      class(damping_param), allocatable :: param

      y = x
      !$omp parallel do default(shared) schedule(static) &
      !$omp shared(dataset, param, x, gradient) &
      !$omp private(i, fl, fr) firstprivate(y)
      do i = 1, size(x)
         y(i) = x(i) + step
         param = from_array(y)
         call single_eval(dataset, from_array(y), fl, 0)
         y(i) = x(i) - step
         param = from_array(y)
         call single_eval(dataset, from_array(y), fr, 0)
         gradient(i) = (fl - fr)/(2*step)
         y(i) = x(i)
      end do
   end subroutine grad_eval

   !> Create job for evaluation of an entry of a data set
   subroutine create_job(job, lentry, basename, directory, error)
      !> Working instructions for entry
      type(job_type), intent(out) :: job
      !> Entry to be expanded
      type(entry_type), intent(in) :: lentry
      !> Name of the geometry input
      character(len=*), intent(in) :: basename
      !> Working directory
      character(len=*), intent(in), optional :: directory
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      character(len=:), allocatable :: input
      real(wp) :: charge
      integer :: stat, unit
      logical :: exist

      if (present(directory)) then
         input = directory//"/"//lentry%dir//"/"//basename
      else
         input = lentry%dir//"/"//basename
      end if

      call read_structure(job%mol, input, error)
      if (allocated(error)) return

      if (present(directory)) then
         input = directory//"/"//lentry%dir//"/.CHRG"
      else
         input = lentry%dir//"/.CHRG"
      end if

      inquire (file=input, exist=exist)
      if (exist) then
         open (file=input, newunit=unit)
         read (unit, *, iostat=stat) charge
         if (stat == 0) then
            job%mol%charge = charge
         end if
         close (unit)
      end if
      call wrap_to_central_cell(job%mol%xyz, job%mol%lattice, job%mol%periodic)

      call new_d4_model(job%d4, job%mol)

   end subroutine create_job

   !> Evaluate entry of a data set
   subroutine run_job(job, param, energy)
      !> Work to be runned
      type(job_type), intent(in) :: job
      !> Current damping parameters
      class(damping_param), intent(in) :: param
      !> Dispersion energy
      real(wp), intent(out) :: energy

      call get_dispersion(job%mol, job%d4, param, realspace_cutoff(), energy)

   end subroutine run_job

   !> Read the dataset from the input file
   subroutine read_dataset(filename, format, basename, dataset, directory, error)
      !> File name
      character(len=*), intent(in) :: filename
      !> Format of the input file
      integer :: format
      !> Basename of the coordinate file
      character(len=*), intent(in) :: basename
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

      allocate (dataset%records(0))

      open (file=filename, newunit=unit, iostat=stat)
      do while (stat == 0)
         call getline(unit, line, stat)
         if (stat /= 0) exit
         call read_record(line, record, entries, format, error)
         if (allocated(error)) return
         dataset%records = [dataset%records, record]
      end do
      close (unit)

      allocate (dataset%jobs(size(entries)))
      do ijob = 1, size(entries)
         call create_job(dataset%jobs(ijob), entries(ijob), basename, directory, error)
         if (allocated(error)) return
      end do
   end subroutine read_dataset

   !> Read record from a line
   subroutine read_record(line, record, entries, format, error)
      !> Line
      character(len=*), intent(in) :: line
      !> Record
      type(record_type), intent(out) :: record
      !> Entries
      type(entry_type), allocatable, intent(inout) :: entries(:)
      !> Format of input file
      integer :: format
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(entry_type) :: lentry
      integer :: first, last, coeff, stat, counter
      character(len=10) :: format_string, coeff_string

      allocate (record%idx(0), record%coeffs(0))
      first = 1
      last = 0

      select case(format)
      case(1)
         ! Stoichiometry factors are not required. First entry is taken as the 
         ! product. All others are educts.
         !
         ! Example:
         ! S22x5/01-0.9, S22x5/01-A, S22x5/01-B, 1.0007611865e-03
         ! S22x5/01-1.0, S22x5/01-A, S22x5/01-B, 1.5228237266e-03
         ! S22x5/01-1.2, S22x5/01-A, S22x5/01-B, 1.6586059147e-03
         ! S22x5/01-1.5, S22x5/01-A, S22x5/01-B, 1.2297590834e-03
         ! S22x5/01-2.0, S22x5/01-A, S22x5/01-B, 6.2420992500e-04
         ! ...
         coeff = -1
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
      case(2)
         ! Stoichiometry factors are explicitly given after the directory.
         !
         ! Example:
         ! S22x5/01-0.9, -1, S22x5/01-A, 1, S22x5/01-B, 1, 1.0007611865e-03
         ! S22x5/01-1.0, -1, S22x5/01-A, 1, S22x5/01-B, 1, 1.5228237266e-03
         ! S22x5/01-1.2, -1, S22x5/01-A, 1, S22x5/01-B, 1, 1.6586059147e-03
         ! S22x5/01-1.5, -1, S22x5/01-A, 1, S22x5/01-B, 1, 1.2297590834e-03
         ! S22x5/01-2.0, -1, S22x5/01-A, 1, S22x5/01-B, 1, 6.2420992500e-04
         ! ...
         counter = 1
         do
            last = index(line(first:), ',') + first - 2
            if (last < first) then
               last = len(line)
               exit
            end if

            ! print '(*(a))', line, new_line('a'), repeat(' ', first-1), repeat('=', last-first+1)
            
            if (modulo(counter, 2) == 0) then
               coeff_string = trim(adjustl(line(first:last)))
               read(coeff_string, "(I4)", iostat=stat) coeff
               if (stat /= 0 ) then
                  call fatal_error(error, "Cannot convert stoichiometry coefficient '"//coeff_string//"' to integer")
                  return
               end if

               record%coeffs = [record%coeffs, coeff]
            else
               lentry%dir = trim(adjustl(line(first:last)))
               call push_back(entries, lentry)
               
               record%idx = [record%idx, find(entries, lentry%dir)]
            end if

            first = last + 2
            counter = counter + 1
         end do
      case default
         write (format_string, "(I4)") format
         call fatal_error(error, "Unknown format option '"//trim(format_string)//"'")
         return
      end select

      read (line(first:), *, iostat=stat) record%reference
      if (stat /= 0 ) then
         call fatal_error(error, "Cannot read reference energy '"//line(first:)//"'")
         return
      end if
   end subroutine read_record

   !> Add new entry to table
   subroutine push_back(entries, lentry)
      !> Entries
      type(entry_type), allocatable, intent(inout) :: entries(:)
      !> Entry
      type(entry_type), intent(in) :: lentry

      integer :: pos

      if (.not. allocated(entries)) allocate (entries(0))

      pos = find(entries, lentry%dir)
      if (pos == 0) entries = [entries, lentry]
   end subroutine push_back

   !> Find existing entry in table
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

      allocate (character(len=0) :: line)
      do
         read (unit, '(a)', advance='no', iostat=stat, iomsg=msg, size=size) &
            & buffer
         if (stat > 0) exit
         line = line//buffer(:size)
         if (stat < 0) exit
      end do

      if (is_iostat_eor(stat)) stat = 0

      if (stat /= 0) then
         if (present(iomsg)) iomsg = trim(msg)
      end if
      iostat = stat

   end subroutine getline

end module dftd4_driver
