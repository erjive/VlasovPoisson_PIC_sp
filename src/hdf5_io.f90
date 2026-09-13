! ===========================================================================
! hdf5_io.f90
! ===========================================================================
!> HDF5 output for VP_PIC.  Alternative to the ASCII output in utils.f90
!! (save0Ddata/save1Ddata/save2Ddata_particles/save_data), selected at
!! run time via the "output_format" parameter ("ascii" or "hdf5").
!!
!! On-disk layout: one file per run, "<directory>/vlasov_output.h5".
!! The radial grid r(1:Nr) never changes during a run, so it is written
!! once to "/grid/r".  Every saved snapshot becomes its own group,
!! "/step_<l>" (l = the time step index, zero-padded), holding:
!!
!!   attributes: time, kinetic_energy, potential_energy, total_energy
!!   datasets:   rho, avg_rho, curr           (r**2 * quantity, size Nr)
!!               force, potential             (size Nr, only if autointeraction)
!!               r_part, p_part, fl            (size Npart; fl = l_part*f)
!!
!! This mirrors the block structure of the ASCII files (one block per
!! saved time), just as HDF5 groups instead of text blocks, which keeps
!! this first implementation simple: no resizable/extendible datasets,
!! and it works unchanged if Npart varies at run time (reduce_arrays).

module hdf5_io

  use hdf5
  use parameters
  use arrays

  implicit none

  integer(HID_T), private :: file_id
  logical, private :: grid_written = .false.

contains

  !> Create the HDF5 output file for this run.  Call once, before the
  !! first call to save_data_hdf5.
  subroutine open_hdf5_file

    integer :: error

    call h5open_f(error)
    call h5fcreate_f(trim(directory)//'/vlasov_output.h5',H5F_ACC_TRUNC_F,file_id,error)

    if (error/=0) then
       print *
       print *, 'ERROR: could not create HDF5 output file in directory ',trim(directory)
       print *, 'Aborting ...'
       print *
       stop
    end if

  end subroutine open_hdf5_file


  !> Close the HDF5 output file.  Call once, at the end of the run.
  subroutine close_hdf5_file

    integer :: error

    call h5fclose_f(file_id,error)
    call h5close_f(error)

  end subroutine close_hdf5_file


  !> Write a 1D real(8) dataset "name" of size n under group loc_id.
  subroutine write_dataset_1d(loc_id,name,values,n)

    implicit none

    integer(HID_T), intent(in) :: loc_id
    character(*), intent(in) :: name
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: values

    integer(HID_T) :: dspace_id,dset_id,plist_id
    integer(HSIZE_T) :: dims(1),chunk_dims(1)
    integer :: error

! Minimum size worth chunking+compressing: for small arrays the fixed
! per-chunk/per-dataset metadata overhead outweighs any space saved
! (measured directly: a run with Nr=51 and ~5 saved steps produced an
! uncompressed .h5 file *larger* than the ASCII files it replaced).
! Below the threshold, fall back to a plain contiguous dataset.

    dims(1) = int(n,HSIZE_T)

    call h5screate_simple_f(1,dims,dspace_id,error)

    if (n>=64) then
       chunk_dims(1) = dims(1)
       call h5pcreate_f(H5P_DATASET_CREATE_F,plist_id,error)
       call h5pset_chunk_f(plist_id,1,chunk_dims,error)
       call h5pset_deflate_f(plist_id,4,error)
       call h5dcreate_f(loc_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error,plist_id)
       call h5pclose_f(plist_id,error)
    else
       call h5dcreate_f(loc_id,name,H5T_NATIVE_DOUBLE,dspace_id,dset_id,error)
    end if

    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,values,dims,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)

  end subroutine write_dataset_1d


  !> Write a scalar real(8) attribute "name" on loc_id (a file, group
  !! or dataset identifier).
  subroutine write_attr_real(loc_id,name,value)

    implicit none

    integer(HID_T), intent(in) :: loc_id
    character(*), intent(in) :: name
    real(8), intent(in) :: value

    integer(HID_T) :: aspace_id,attr_id
    integer(HSIZE_T) :: adims(1) = (/1/)
    real(8) :: buf(1)
    integer :: error

    buf(1) = value

    call h5screate_f(H5S_SCALAR_F,aspace_id,error)
    call h5acreate_f(loc_id,name,H5T_NATIVE_DOUBLE,aspace_id,attr_id,error)
    call h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,buf,adims,error)
    call h5aclose_f(attr_id,error)
    call h5sclose_f(aspace_id,error)

  end subroutine write_attr_real


  !> Save one snapshot (grid profiles + particle sample + energy) as a
  !! new group "/step_<l>".  Mirrors save_data in utils.f90.
  subroutine save_data_hdf5(l)

    implicit none

    integer, intent(in) :: l

    character(20) :: groupname
    integer(HID_T) :: group_id,grid_group_id
    integer :: error

    if (.not. grid_written) then
       call h5gcreate_f(file_id,'/grid',grid_group_id,error)
       call write_dataset_1d(grid_group_id,'r',r(1:Nr),Nr)
       call h5gclose_f(grid_group_id,error)
       grid_written = .true.
    end if

    write(groupname,'(A,I10.10)') '/step_',l

    call h5gcreate_f(file_id,trim(groupname),group_id,error)

    call write_attr_real(group_id,'time',t)
    call write_attr_real(group_id,'kinetic_energy',kinetic)
    call write_attr_real(group_id,'potential_energy',potential)
    call write_attr_real(group_id,'total_energy',total_energy)

    call write_dataset_1d(group_id,'rho',    r(1:Nr)**2*rho(1:Nr),    Nr)
    call write_dataset_1d(group_id,'avg_rho',r(1:Nr)**2*avg_rho(1:Nr),Nr)
    call write_dataset_1d(group_id,'curr',   r(1:Nr)**2*curr(1:Nr),   Nr)

    if (autointeraction) then
       call write_dataset_1d(group_id,'force',    force(1:Nr),Nr)
       call write_dataset_1d(group_id,'potential',pot(1:Nr),  Nr)
    end if

    call write_dataset_1d(group_id,'r_part',r_part,  Npart)
    call write_dataset_1d(group_id,'p_part',p_part,  Npart)
    call write_dataset_1d(group_id,'fl',    l_part*f,Npart)

    call h5gclose_f(group_id,error)

  end subroutine save_data_hdf5

end module hdf5_io
