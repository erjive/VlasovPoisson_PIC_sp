! ===========================================================================
! utils.f90
! ===========================================================================
!> Utility file


module utils

  use parameters
  use arrays

 contains

  !> Read the parameters of the program
  subroutine read_initial_param

    read(*,*) dr
    read(*,*) Nrc
    read(*,*) Npc
    read(*,*) Nlc
    read(*,*) courant
    read(*,*) Nt
    read(*,*) rmin
    read(*,*) rmax
    read(*,*) rminc
    read(*,*) rmaxc
    read(*,*) pminc
    read(*,*) pmaxc
    read(*,*) lminc
    read(*,*) lmaxc
    read(*,*) reduceparticles
    read(*,*) Nreduce
    read(*,*) time_output
    read(*,*) spatial_output
    read(*,*) directory
    read(*,*) a0
    read(*,*) r0
    read(*,*) p0
    read(*,*) l0
    read(*,*) sr
    read(*,*) sp
    read(*,*) sl
    read(*,*) state
    read(*,*) cutoff
    read(*,*) bsplineorder
    read(*,*) integrator
    read(*,*) spatialorder
    read(*,*) forcetype
    read(*,*) BGtype
    read(*,*) autointeraction

  end subroutine read_initial_param


! Sanity check.
  subroutine test_consistency

  if (rmin<0.0D0) then
     print *
     print *, 'rmin must be greater than or equal to zero'
     print *, 'Aborting ...'
     print *
     stop
  end if

  if (rmax<=rmin) then
     print *
     print *, 'rmin must be smaller than rmax'
     print *, 'Aborting ...'
     print *
     stop
  end if

!  if (pmax<=0.d0) then
!     print *
!     print *, 'pmax must be greater then zero'
!     print *, 'Aborting ...'
!     print *
!     stop
!  end if

!  if (pmax<=pmin) then
!     print *
!     print *, 'pmin must be smaller than pmax'
!     print *, 'Aborting ...'
!     print *
!     stop
!  end if

  if (lmin<0.D0) then
     print *
     print *, 'lminc must be greater than or equal to zero'
     print *, 'Aborting ...'
     print *
     stop
  end if

  if (lmax<0.D0) then
     print *
     print *, 'lmaxc must be greater than or equal to zero'
     print *, 'Aborting ...'
     print *
     stop
  end if

  if (lmax<lmin) then
     print *
     print *, 'lminc must be smaller than lmaxc'
     print *, 'Aborting ...'
     print *
     stop
  end if

  end subroutine test_consistency



  !> Set all the parameters that are not set in the input file
  subroutine set_grid_size
! ***************************
! ***   FIND GRID SIZES   ***
! ***************************

    ! Find out number of grid points in r direction.
    Nr = int((rmax-rmin)/dr)
    if (rmin==0.0d0) Nr=Nr+1

    ! Find out number of grid points in p direction.
    !Np = 2*int(pmax/dp)

    print *
    print *, 'Number of points in r direction ',Nr
    !print *, 'Number of points in p direction ',Np

    ! Minimum radius when the angular momentum is fix.
    !eps = Lfix/(10.0D0*pmax)
     eps = 0.0D0
  end subroutine set_grid_size

  !> Allocate all memory
  subroutine alloc_mem_set0

! Find the number of ghost zones

  if (rmin>0.d0) then
    ghost = 0
  else if (rmin == 0.d0) then
    if (spatialorder == "two") then
        ghost = 2
    else if (spatialorder == "four") then
        ghost = 3
    end if
  end if

! Find the number of particles

  Npart = Nrc*Npc*Nlc

  print *, "Number of computational particles=",Npart
! **************************************
! ***   ALLOCATE MEMORY FOR ARRAYS   ***
! **************************************

! Position and momentum of particles

  allocate(r_part      (1:Npart))
  allocate(r_part_p    (1:Npart))
  allocate(p_part      (1:Npart))
  allocate(p_part_p    (1:Npart))
  allocate(p_part_h    (1:Npart))
  allocate(l_part      (1:Npart))
  allocate(pot_part    (1:Npart))
  allocate(force_part  (1:Npart))

  r_part     = 0.d0
  r_part_p   = 0.d0
  p_part     = 0.d0
  p_part_p   = 0.d0
  p_part_h   = 0.d0
  l_part     = 0.d0
  pot_part   = 0.d0
  force_part = 0.d0

! Coordinates, force and potential.

  allocate(r(1-ghost:Nr))
  r = 0.0d0

  if (autointeraction) then
    allocate(force  (1-ghost:Nr))
    allocate(pot    (1-ghost:Nr))
    allocate(dev_pot(1-ghost:Nr))
    force   = 0.0d0
    pot     = 0.0d0
    dev_pot = 0.0d0
  end if


! Density function f and sources.

  allocate(f  (1:Npart))
  !allocate(f_p(1:Npart))

  f   = 0.0d0
  !f_p = 0.0d0

! Integrates density, current and continuity equation.

  allocate(rho   (1-ghost:Nr))
  allocate(rho_p (1-ghost:Nr))
  allocate(avg_rho (1-ghost:Nr))
  allocate(curr  (1-ghost:Nr))
  allocate(curr_p(1-ghost:Nr))
  allocate(cont  (1-ghost:Nr))

  rho    = 0.0d0
  rho_p  = 0.0d0
  avg_rho= 0.0d0
  curr   = 0.0d0
  curr_p = 0.0d0
  cont   = 0.0d0

  end subroutine alloc_mem_set0


  !> Free all the memory in the allocated arrays
  subroutine deallocate_mem

  deallocate(r_part)
  deallocate(r_part_p)
  deallocate(p_part)
  deallocate(p_part_p)
  deallocate(p_part_hp)
  deallocate(pot_part)
  deallocate(force_part)

  deallocate(r)
  deallocate(force)
  deallocate(pot)
  deallocate(dev_pot)

! Density function f

  deallocate(f)
!  deallocate(f_p)

! Array for residual evaluation (only for convergence study)
  if (conv_test=="on") then
    deallocate(res)
  end if

  deallocate(rho)
  deallocate(rho_p)
  deallocate(avg_rho)
  deallocate(curr)
  deallocate(curr_p)
  deallocate(cont)

  if (autointeraction) then
    deallocate(pot    )
    deallocate(dev_pot)
  end if
  
  print *, "Memory deallocated"
  
  end subroutine deallocate_mem


subroutine construct_grid
! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

! Position in direction r.  In case that rmin=0
! we make sure to stagger the origin and we add
! two ghost points to the left of the origin for 
! second order integration, and three points for 
! fourth order integration.

  integer i

  if (rmin>0.d0) then
     do i=0,Nr
        r(i) = rmin + dble(i)*dr
     end do
  else

  do i=1-ghost,Nr
    r(i) = (dble(i)-0.5d0)*dr
  end do
      
  end if

! Position in direction p. Notice that the grid in the
! p direction will always cover the region (-pmax,pmax).

!  do j=0,Np
!     p(j) = pmin + dble(j) *dp
!  end do


end subroutine construct_grid

  !> Set the time step.
  !! Here we find the time step using information from the
  !! Courant factor and the maximum value of the momentum.
  subroutine set_timestep

! *********************
! ***   TIME STEP   ***
! *********************

! Here we find the time step using information from the
! Courant factor and the maximum value of the momentum.

! Make sure that the Courant condition is satisfied in
! the r direction.

  integer i
  real(8) dtr,dtp       ! Auxiliary variables.
!  real(8) pmax_aux,Fmax_aux
  
!  pmax_aux = MAXVAL(abs(p_part))
!  dtr = courant*dr/pmax_aux
  dtr = courant*dr/pmax
! Now make sure that it is also satisfied in the p direction.
! For this we need to find the maximum value of the force.

  if (BGtype /= "null") then
    Fmax = 0.0d0

!    Fmax=maxval(force_part)
    do i=1,Npart
       Fmax = max(Fmax,abs(force_part(i)))
    end do

    dtp = courant*dpc/Fmax

! Fix the time step as the smallest bvalue between dtr and dtp.

    dt = min(dtr,dtp)
  else 
    dt = dtr
  end if
    dt = dtr
  end subroutine set_timestep



  !> Save all the data to the corresponding files
  subroutine save_data

    character(100) :: filename     !< Name of output file

!   Save distribution function f.
    filename = 'vlasov_fdist'
    !call save3Ddata_particles(directory,filename,Npart,t,r_part,p_part,l_part,f)
    call save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,l_part*f)

!   Save rho, curr and cont (multiplied by r**2).
    filename = 'vlasov_density'
    call save1Ddata(directory,filename,Nr,t,r,r**2*rho)
    filename = 'vlasov_avg_density'
    call save1Ddata(directory,filename,Nr,t,r,r**2*avg_rho)
    filename = 'vlasov_curr'
    call save1Ddata(directory,filename,Nr,t,r,r**2*curr)
    filename = 'vlasov_energy'
    call save0Ddata(directory,filename,t,total_energy)
    filename = 'vlasov_k_phi_e'
    call save_energy(directory,filename,t,kinetic,potential,total_energy)
!    filename = 'vlasov_cont'
!    call save1Ddata(directory,filename,Nr,Np,t,r,r**2*cont)


!   Save force and potential.
    if (autointeraction) then

       filename = 'vlasov_force'
       call save1Ddata(directory,filename,Nr,t,r,force)

       filename = 'vlasov_potential'
       call save1Ddata(directory,filename,Nr,t,r,pot)

       filename = 'vlasov_potential_r0'
       call save0Ddata(directory,filename,t,pot(1))

       filename = 'vlasov_force_r0'
       call save0Ddata(directory,filename,t,force(1))
    end if

  end subroutine save_data


  subroutine save0Ddata(directory,filename,t,var)

! **********************
! ***   SAVE0DDATA   ***
! **********************

! This subroutine saves 1D data to files.

  implicit none

  real(8) t,var

  character(20) directory,filename,filestatus


! **************************
! ***   OPEN DATA FILE   ***
! **************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open file.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim(filename)//'.tl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.tl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(2ES16.8)") t,var


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save0Ddata



  subroutine save1Ddata(directory,filename,Nr,t,r,var)

! **********************
! ***   SAVE1DDATA   ***
! **********************

! This subroutine saves 1D data to files.

  implicit none

  integer i
  integer Nr

  real(8) t

  real(8), dimension(1:Nr) :: r,var

  character(20) directory,filename,filestatus


! **************************
! ***   OPEN DATA FILE   ***
! **************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open file.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=1,Nr
     if (dabs(var(i)).gt.1.0D-50) then
        write(101,"(2ES16.8)") r(i),var(i)

     else
        write(101,"(2ES16.8)") r(i),0.0d0
     end if
  end do

! Leave two blank spaces before next time.
! The reason to leave two spaces is that 'gnuplot' asks
! for two spaces to distinguish different records.

  write (101,*)
  write (101,*)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save1Ddata

subroutine save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,var)

! **********************
! ***   SAVE2DDATA   ***
! **********************

! This subroutine saves 2D data to files.

  implicit none

  integer i
  integer Npart

  real(8) t

  real(8), dimension(1:Npart) :: r_part
  real(8), dimension(1:Npart) :: p_part

  real(8), dimension(1:Npart) :: var

  character(20) directory,filename,filestatus


! ***************************
! ***   OPEN DATA FILES   ***
! ***************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open files.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim(filename)//'.2D',form='formatted',status=filestatus)

  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.2D',form='formatted',status=filestatus,position='append')
  end if


! ************************
! ***   SAVE 2D DATA   ***
! ************************

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=1,Npart
    if (dabs(var(i)).gt.1.0D-50) then
      write(101,"(3ES16.8)") r_part(i),p_part(i),var(i)
    else
      write(101,"(3ES16.8)") r_part(i),p_part(i),0.0D0
    end if
  end do

  write (101,*)
  write (101,*)


! ****************************
! ***   CLOSE DATA FILES   ***
! ****************************

  close(101)

! ***************
! ***   END   ***
! ***************

  end subroutine save2Ddata_particles

  subroutine save3Ddata_particles(directory,filename,Npart,t,r_part,p_part,l_part,var)

! **********************
! ***   SAVE2DDATA   ***
! **********************

! This subroutine saves 2D data to files.

  implicit none

  integer i
  integer Npart

  real(8) t

  real(8), dimension(1:Npart) :: r_part
  real(8), dimension(1:Npart) :: p_part
  real(8), dimension(1:Npart) :: l_part

  real(8), dimension(1:Npart) :: var

  character(20) directory,filename,filestatus


! ***************************
! ***   OPEN DATA FILES   ***
! ***************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open files.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim(filename)//'.2D',form='formatted',status=filestatus)

  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.2D',form='formatted',status=filestatus,position='append')
  end if


! ************************
! ***   SAVE 2D DATA   ***
! ************************

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=1,Npart
    if (dabs(var(i)).gt.1.0D-50) then
      write(101,"(4ES16.8)") r_part(i),p_part(i),l_part(i),var(i)
    else
      write(101,"(4ES16.8)") r_part(i),p_part(i),l_part(i),0.0D0
    end if
  end do

  write (101,*)
  write (101,*)


! ****************************
! ***   CLOSE DATA FILES   ***
! ****************************

  close(101)

! ***************
! ***   END   ***
! ***************

  end subroutine save3Ddata_particles


  subroutine save_energy(directory,filename,t,kinetic,potential,energy)

! **********************
! ***   SAVE ENERGY  ***
! **********************

! This subroutine saves the energy (kinectic,potential,virial) data to files.

  implicit none

  real(8) t,kinetic,potential,energy

  character(20) directory,filename,filestatus


! **************************
! ***   OPEN DATA FILE   ***
! **************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open file.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim(filename)//'.tl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.tl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  write(101,"(5ES16.8)") t,kinetic,potential,energy


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save_energy



  subroutine save_dens_current(directory,filename,Nr,Np,t,r,density,current,error)

! ************************************************
! ***   SAVE DENSITY, CURRENT AND ERROR DATA   ***
! ************************************************

! This subroutine saves the density, current and error in the
! continuity equation.

  implicit none

  integer i
  integer Nr,Np

  real(8) t

  real(8), dimension(0:Nr) :: r,density,current,error

  character(20) directory,filename,filestatus


! **************************
! ***   OPEN DATA FILE   ***
! **************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open file.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  density = density + 1.0D-50
  current = current + 1.0D-50
  error   = error   + 1.0D-50

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(101,"(4ES16.8)") r(i),density(i),current(i),error(i)
  end do

  write (101,*)
  write (101,*)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save_dens_current


  subroutine save_force_pot(directory,filename,Nr,Np,t,r,force,pot)

! *****************************************
! ***   SAVE FORCE AND POTENTIAL DATA   ***
! *****************************************

! This subroutine saves the force and potential data for the
! self gravitating case

  implicit none

  integer i
  integer Nr,Np

  real(8) t

  real(8), dimension(0:Nr) :: r,force,pot

  character(20) directory,filename,filestatus


! **************************
! ***   OPEN DATA FILE   ***
! **************************

! Is this the first time step?

  if (t==0) then
     filestatus = 'replace'
  else
     filestatus = 'old'
  end if

! Open file.

  if (filestatus=='replace') then
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus)
  else
     open(101,file=trim(directory)//'/'//trim(filename)//'.rl',form='formatted',status=filestatus,position='append')
  end if


! *********************
! ***   SAVE DATA   ***
! *********************

  force = force + 1.0D-50
  pot   = pot   + 1.0D-50

  write(101,"(A8,ES14.6)") '#Time = ',t

  do i=0,Nr
     write(101,"(3ES16.8)") r(i),force(i),pot(i)
  end do

  write (101,*)
  write (101,*)


! ***************************
! ***   CLOSE DATA FILE   ***
! ***************************

  close(101)


! ***************
! ***   END   ***
! ***************

  end subroutine save_force_pot


subroutine reduce_arrays
  use parameters
  use arrays

 
  integer i,j
  integer counter,Npart_aux
  real(8), allocatable, dimension (:) :: r_aux,p_aux,l_aux,f_aux
!  real(8), dimension (1:Npart) :: r_aux,p_aux,f_aux


! Count How many particles are still in the grid

  counter = 0

  do i=1,Npart
    if (r_part(i)<= rmax) then
      counter = counter + 1
    end if 
  end do

! Reduce the size of arrays if we have less than 95% 
! of the original number of particles.

  if (dble(counter) < 0.95D0* Npart) then 
  allocate(r_aux(1:Npart))
  allocate(p_aux(1:Npart))
  allocate(l_aux(1:Npart))
  allocate(f_aux(1:Npart))

! Copy the position and momentum of particles

  r_aux = r_part
  p_aux = p_part
  l_aux = l_part
  f_aux = f

! Deallocate arrays

  deallocate(r_part)
  deallocate(r_part_p)
  deallocate(p_part)
  deallocate(p_part_p)
  deallocate(p_part_h)
  deallocate(l_part)
  deallocate(pot_part)
  deallocate(force_part)
  deallocate(f)
!  deallocate(f_p)

! Allocate arrays 

  Npart_aux = Npart
  Npart     = counter

  allocate(r_part(1:Npart))
  allocate(r_part_p(1:Npart))
  allocate(p_part(1:Npart))
  allocate(p_part_p(1:Npart))
  allocate(p_part_h(1:Npart))
  allocate(l_part(1:Npart))
  allocate(pot_part(1:Npart))
  allocate(force_part(1:Npart))
  allocate(f(1:Npart))
!  allocate(f_p(1:Npart))

  j = 1
  do i=1,Npart_aux
    if (r_aux(i)<rmax) then
      r_part(j) = r_aux(i) 
      p_part(j) = p_aux(i)
      l_part(j) = l_aux(i)
      f(j)      = f_aux(i)
      j = j+1
    end if

  end do
  print *, "In the grid are still", Npart, "computational particles"
  end if
end subroutine reduce_arrays

end module utils

