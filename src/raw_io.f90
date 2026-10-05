! ===========================================================================
! raw_io.f90
! ===========================================================================
!> Raw binary output for VP_PIC, selected with output_format="raw". Third
!! alternative to the ASCII output of utils.f90 and the HDF5 output of
!! hdf5_io.f90, with the same content as the latter.
!!
!! Frequent snapshots are dominated by the per-object bookkeeping of HDF5
!! rather than by the bytes themselves, so a plain stream of records is
!! cheaper to write. The price is that the format is not self-describing: a
!! reader needs the exact layout below (tools/raw_io.py implements it).
!!
!! One file per run, "<directory>/vlasov_output.raw", opened with
!! access="stream" (a plain byte stream, no Fortran record markers), so the
!! layout below is exact and is what numpy.fromfile expects. All integers and
!! reals have 8 bytes, in the byte order of the machine that wrote the file.
!!
!!   HEADER (once, at the start of the file):
!!     int64             Nr
!!     int64             1 with autointeraction, 0 without
!!     float64(Nr)       r          -- radial grid, fixed for the whole run
!!
!!   RECORD (one per saved snapshot, back to back, in save order):
!!     int64             l          -- time step index
!!     float64           time
!!     float64           kinetic_energy
!!     float64           potential_energy
!!     float64           total_energy
!!     int64             Npart      -- number of particles of this record
!!     float64(Nr)       rho        -- r**2 * rho
!!     float64(Nr)       avg_rho    -- r**2 * avg_rho
!!     float64(Nr)       curr       -- r**2 * curr
!!     float64(Nr)       force      -- only with autointeraction
!!     float64(Nr)       potential  -- only with autointeraction
!!     float64(Npart)    r_part
!!     float64(Npart)    p_part
!!     float64(Npart)    fl         -- l_part*f
!!     float64(Npart)    l_part
!!
!! Npart can change during a run (reduceparticles), so the records do not
!! have a fixed size: a reader walks the file once, reading Npart at each
!! record to know how far to skip. The layout is that of vlasov-poisson_PIC
!! with fl and l_part in place of f.

module raw_io

  use parameters
  use arrays

  implicit none

  integer, private :: raw_unit
  logical, private :: header_written = .false.

contains

  !> Create the raw output file of this run. Call once, before the first
  !! call to save_data_raw.
  subroutine open_raw_file

    integer :: ios

    open(newunit=raw_unit, file=trim(directory)//'/vlasov_output.raw', &
         form='unformatted', access='stream', status='replace', iostat=ios)

    if (ios/=0) then
       print *
       print *, 'ERROR: could not create raw output file in directory ',trim(directory)
       print *, 'Aborting ...'
       print *
       stop 1
    end if

    header_written = .false.

  end subroutine open_raw_file


  !> Close the raw output file. Call once, at the end of the run.
  subroutine close_raw_file

    close(raw_unit)

  end subroutine close_raw_file


  !> Save one snapshot as one record (layout at the top of this file).
  !! Mirrors save_data_hdf5 in hdf5_io.f90.
  subroutine save_data_raw(l)

    implicit none

    integer, intent(in) :: l

    integer(8) :: autoint8

    if (.not. header_written) then
       write(raw_unit) int(Nr,8)
       autoint8 = 0_8
       if (autointeraction) autoint8 = 1_8
       write(raw_unit) autoint8
       write(raw_unit) r(1:Nr)
       header_written = .true.
    end if

    write(raw_unit) int(l,8)
    write(raw_unit) t
    write(raw_unit) kinetic
    write(raw_unit) potential
    write(raw_unit) total_energy
    write(raw_unit) int(Npart,8)

    write(raw_unit) r(1:Nr)**2*rho(1:Nr)
    write(raw_unit) r(1:Nr)**2*avg_rho(1:Nr)
    write(raw_unit) r(1:Nr)**2*curr(1:Nr)

    if (autointeraction) then
       write(raw_unit) force(1:Nr)
       write(raw_unit) pot(1:Nr)
    end if

    write(raw_unit) r_part
    write(raw_unit) p_part
    write(raw_unit) l_part*f
    write(raw_unit) l_part

  end subroutine save_data_raw

end module raw_io
