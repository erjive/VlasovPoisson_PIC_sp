! ===========================================================================
! utils.f90
! ===========================================================================
!> Utility file


module utils

  use parameters
  use arrays

 contains

  ! The run parameters are read by read_parameters (paramfile.f90), which
  ! also checks them; the positional reader and test_consistency are gone.



  !> Set all the parameters that are not set in the input file
  subroutine set_grid_size
! ***************************
! ***   FIND GRID SIZES   ***
! ***************************

    ! Find out number of grid points in r direction.
    ! Staggered grid, r(i) = rmin + (i-1/2) dr, covering [rmin,rmax].
    Nr = int((rmax-rmin)/dr) + 1

    print *
    print *, 'Number of points in r direction ',Nr

    ! No softening of the centrifugal term: the angle-action maps (analysish,
    ! aa_quad, analytic) use the unsoftened L**2/(2 r**2), and a softened
    ! potential would not match them.
    eps = 0.0D0
  end subroutine set_grid_size

  !> Allocate all memory
  subroutine alloc_mem_set0

! Find the number of ghost zones

! Two ghost points at the origin reach the support of the widest weight
! function (order 3, two grid spacings on each side).

  if (rmin>0.d0) then
    ghost = 0
  else
    ghost = 2
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
  allocate(potself_part(1:Npart))
  allocate(force_part  (1:Npart))

  r_part     = 0.d0
  r_part_p   = 0.d0
  p_part     = 0.d0
  p_part_p   = 0.d0
  p_part_h   = 0.d0
  l_part     = 0.d0
  pot_part   = 0.d0
  potself_part = 0.d0
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

  f   = 0.0d0

! Integrates density, current and continuity equation.

  allocate(rho   (1-ghost:Nr))
  allocate(avg_rho (1-ghost:Nr))
  allocate(curr  (1-ghost:Nr))

  rho    = 0.0d0
  avg_rho= 0.0d0
  curr   = 0.0d0

  end subroutine alloc_mem_set0


  !> Free all the memory in the allocated arrays.
  !!
  !! Must mirror alloc_mem_set0 exactly: deallocate every array
  !! allocated there (and only those), each guarded by the same
  !! condition used to allocate it.  Deallocating an array that was
  !! never allocated is a runtime error, not a no-op.
  subroutine deallocate_mem

! q0_part/j0_part exist only for integrator="analytic".
  if (allocated(q0_part)) deallocate(q0_part)
  if (allocated(j0_part)) deallocate(j0_part)

  deallocate(r_part)
  deallocate(r_part_p)
  deallocate(p_part)
  deallocate(p_part_p)
  deallocate(p_part_h)
  deallocate(l_part)
  deallocate(pot_part)
  deallocate(potself_part)
  deallocate(force_part)

  deallocate(r)

! force, pot and dev_pot only exist with self-gravity (alloc_mem_set0).

  if (autointeraction) then
    deallocate(force)
    deallocate(pot)
    deallocate(dev_pot)
  end if

! Density function f

  deallocate(f)

  deallocate(rho)
  deallocate(avg_rho)
  deallocate(curr)

  print *, "Memory deallocated"

  end subroutine deallocate_mem


subroutine construct_grid
! *************************************
! ***   FIND GRID POINT POSITIONS   ***
! *************************************

! Staggered grid: r(i) = rmin + (i-1/2) dr. With rmin = 0 no point sits at
! the origin, and the ghost points 1-k sit at -r(k), the mirror images used by
! the symmetry f(r,p) = f(-r,-p).

  integer i

  do i=1-ghost,Nr
    r(i) = rmin + (dble(i)-0.5d0)*dr
  end do

end subroutine construct_grid

  !> Build a cell list for the current particles.
  !!
  !! Groups particle indices (1:Npart) by the grid point nearest to them, so
  !! that the deposit on grid point i only scans the particles of a few
  !! nearby cells instead of all of them.
  !!
  !! On output, the particles assigned to grid cell c (1<=c<=Nr) are
  !! particle_order(cell_start(c):cell_start(c+1)-1). Both arrays
  !! are allocated here; the caller must deallocate them.
  !!
  !! The grid is uniform, r(k) = r(1) + (k-1)*dr, so the nearest grid index
  !! to a position x is nint((x-r(1))/dr) + 1. Particles beyond either end of
  !! the grid (including negative radii within a time step) are filed in the
  !! end cell; the deposit still applies the exact support of the weight.
  subroutine build_cell_list(cell_start,particle_order)

    implicit none

    integer, allocatable, intent(out) :: cell_start(:)
    integer, allocatable, intent(out) :: particle_order(:)

    integer :: j,c
    integer, allocatable :: cell_count(:),cursor(:),ic(:)

    allocate(cell_start(1:Nr+1))
    allocate(particle_order(1:Npart))
    allocate(cell_count(1:Nr))
    allocate(cursor(1:Nr))
    allocate(ic(1:Npart))

    cell_count = 0

    do j=1,Npart
      ic(j) = nint((r_part(j)-r(1))/dr) + 1
      ic(j) = max(1,min(Nr,ic(j)))
      cell_count(ic(j)) = cell_count(ic(j)) + 1
    end do

    cell_start(1) = 1
    do c=1,Nr
      cell_start(c+1) = cell_start(c) + cell_count(c)
    end do

    cursor(1:Nr) = cell_start(1:Nr)

    do j=1,Npart
      c = ic(j)
      particle_order(cursor(c)) = j
      cursor(c) = cursor(c) + 1
    end do

    deallocate(cell_count,cursor,ic)

  end subroutine build_cell_list

  !> Set the time step.
  !! Here we find the time step using information from the
  !! Courant factor, the maximum value of the momentum, the
  !! maximum value of the force (acceleration) and the passage
  !! of the particles through their pericentres.
  subroutine set_timestep

! *********************
! ***   TIME STEP   ***
! *********************

! Make sure that the Courant condition is satisfied in
! the r direction: a particle should not cross more than
! a fraction "courant" of a grid cell dr per time step.

  integer i
  real(8) dtr,dtp,vmax  ! Auxiliary variables.
  real(8) dtl
  real(8), save :: rpmin = 0.0d0, tpmin = -1.0d0   ! Pericentre bound; negative: not computed yet.
  logical first

! pmax > 0 is the velocity scale chosen by the user; pmax <= 0 takes the
! largest |p| of the particles now. Neither criterion knows the orbital
! frequency: with a deep background (nfw, test case of the audit) the step
! from the particles' own |p| made a run diverge. The step has to be checked
! with a convergence test in any case.

  if (pmax > 0.0d0) then
     vmax = pmax
  else
     vmax = maxval(abs(p_part))
  end if

  if (vmax > 0.0d0) then
     dtr = courant*dr/vmax
  else
     dtr = huge(1.0d0)
  end if

! The step is also bounded with the acceleration criterion of symplectic
! N-body integration (e.g. Gadget-2, Springel 2005): the time to move one grid
! spacing dr, the scale on which the force is resolved, from rest under the
! largest force on any particle (background, self-gravity and centrifugal
! term),
!
!   dtp = courant * sqrt(2*dr/Fmax).
!
! By default this runs once, before the main loop (dt_switch="fix").

  Fmax = 0.0d0
  do i=1,Npart
     Fmax = max(Fmax,abs(force_part(i)))
  end do

  if (Fmax>0.0d0) then
     dtp = courant*sqrt(2.0d0*dr/Fmax)
  else
     dtp = huge(1.0d0)
  end if

  dt = min(dtr,dtp)

  if (dt >= huge(1.0d0)) then
     print *
     print *, 'No momentum and no force: the time step is undefined; set pmax > 0.'
     print *, 'Aborting ...'
     print *
     stop 1
  end if

! Third bound, from the pericentre. Near it the length scale of the orbit is
! not dr but r_p, and the particle moves at v_p = L/r_p, so the Courant
! condition with that scale is
!
!   dtl = courant * min_j r_p,j/v_p,j = courant * min_j r_p,j**2/L_j,
!
! with r_p,j the pericentre of particle j. Neither bound above sees it: with
! L = 1e-3 and dt = 0.01 an isochrone orbit had |dE/E| = 3.5 (leapfrog) and
! 200 (yoshida4). The energy error peaks at each pericentre, at about
! 0.03 (Omega_p dt)**2 with leapfrog and 0.13 (Omega_p dt)**4 with yoshida4,
! Omega_p = L/r_p**2, so courant sets the accuracy and the bound only makes
! sure that the passage is resolved: away from the pericentre that orbit
! keeps its energy to 9e-4 with courant = 0.5 and to 5e-9 with 0.25
! (yoshida4; AUDITORIA_FISICA_2026-09-21.md, N1; same bound as in
! vlasov-poisson_PIC). For L of order one it is not the smallest of the
! three. A particle with L = 0 has no pericentre: it goes through the
! centre, where the field is regular, and is left out.
!
! The pericentres are computed once, at the first call, in the field at
! t = 0, also with dt_switch="var": one bisection per particle at every step
! made a self-gravitating run ten times slower. With self-gravity that field
! can change: in a cold collapse the potential at the centre deepens and the
! pericentres shrink, and a step estimated at t = 0 does not resolve them.
!
! The smallest pericentre and the largest Omega_p*dt are reported then,
! whichever bound sets the step.

  first = (tpmin < 0.0d0)
  if (first) call pericentre_min(rpmin,tpmin)

  if (tpmin < huge(1.0d0)) then
     dtl = courant*tpmin
     if (dtl < dt) then
        dt = dtl
        if (first) then
           print *
           print *, 'Time step set by the pericentre, courant*min(r_p**2/L).'
        end if
     end if
     if (first) then
        print *
        print *, 'Smallest pericentre r_p = ',rpmin
        print *, 'Largest Omega_p*dt = L*dt/r_p**2 = ',dt/tpmin
        print *, '(energy error at each pericentre: about 0.03 (Omega_p*dt)**2 with leapfrog,'
        print *, ' 0.13 (Omega_p*dt)**4 with yoshida4)'
     end if
  end if

  end subroutine set_timestep


  !> Smallest pericentre of the particles, rpmin, and smallest time of passage
  !! through it, tpmin = min_j r_p,j**2/L_j, in the field at the time of the
  !! call (the one set by the last grav_force), from
  !!
  !!   p**2/2 + L**2/(2 r**2) + Phi(r) = E,   r_p = smallest root in (0,|r|].
  !!
  !! Phi is the background (closed form) plus, with self-gravity, the self
  !! potential of the grid: pot minus the background, linear between points,
  !! constant inside r(1) (it is even and smooth at the origin) and
  !! Phi(r_Nr) r_Nr/r beyond the grid. The energy uses the same Phi, so the
  !! root always lies in (0,|r|]. The effective potential of a potential that
  !! grows with r has a single minimum, so the root is unique and a bisection
  !! (in log r, down to 1e-12 |r|) finds it; the lower end of the bracket is
  !! returned, so the estimate errs on the safe side. The softening eps is
  !! included as in grav_force. Particles with L = 0 are skipped; if every
  !! particle has L = 0 both results are huge(1.0d0).
  subroutine pericentre_min(rpmin,tpmin)

    real(8), intent(out) :: rpmin,tpmin

    real(8), allocatable :: ps(:)
    real(8) :: x,E,lo,hi,mid
    integer :: j,it
    logical :: sg

    sg = autointeraction

    if (sg) then
       allocate(ps(1:Nr))
       ps = pot(1:Nr) - bgpot(r(1:Nr))
    end if

    rpmin = huge(1.0d0)
    tpmin = huge(1.0d0)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,E,lo,hi,mid,it) REDUCTION(min:rpmin,tpmin)
    do j=1,Npart
       if (l_part(j) == 0.0d0) cycle
       x  = abs(r_part(j))
       E  = 0.5d0*p_part(j)**2 + veff(x,l_part(j))
       lo = 1.0d-12*x
       hi = x
       if (veff(lo,l_part(j)) > E) then
          do it=1,64
             mid = sqrt(lo*hi)
             if (veff(mid,l_part(j)) > E) then
                lo = mid
             else
                hi = mid
             end if
          end do
       end if
       rpmin = min(rpmin,lo)
       tpmin = min(tpmin,lo**2/l_part(j))
    end do
    !$OMP END PARALLEL DO

    if (sg) deallocate(ps)

  contains

    real(8) function veff(y,am)
      real(8), intent(in) :: y,am
      integer :: k
      real(8) :: w
      veff = 0.5d0*am**2/(y**2 + eps**2) + bgpot(y)
      if (sg) then
         if (y <= r(1)) then
            veff = veff + ps(1)
         else if (y >= r(Nr)) then
            veff = veff + ps(Nr)*r(Nr)/y
         else
            k = min(Nr-1,int(y/dr + 0.5d0))
            w = (y - r(k))/dr
            veff = veff + (1.0d0-w)*ps(k) + w*ps(k+1)
         end if
      end if
    end function veff

  end subroutine pericentre_min


! Backgrounds given by closed formulas for the potential and the force
! (force = -dpot/dr, checked symbolically for every one), in units of the
! mass and the scale of the background. grav_force has its own fused loops
! for the isochrone and the point mass; their cases here are for
! set_timestep.
!
! The closed forms are evaluated at |r|, and the force carries the sign of r.
! A particle may step to r < 0 inside a time step (main.f90 reflects it only
! afterwards) and the grid has ghost points at negative radii, so the
! background has to be even in r and its force odd. Written in terms of r
! itself, "nfw" and "burkert" are neither (nfw even reverses the sign of the
! force, pushing a particle that crosses the origin further out), "sphere"
! takes its inner branch for every r < 1 including r < -1, and "iso" is not
! defined at all for r < 0 (AUDITORIA_2026-09-20.md, point 4).

  elemental real(8) function bgpot(x)

    real(8), intent(in) :: x
    real(8) :: a

    a = abs(x)

    select case (BGtype)
    case ("Isochrone")
       bgpot = -1.0d0/(1.0d0+sqrt(1.0d0+a**2))
    case ("Central")
       bgpot = -1.0d0/a
    case ("sphere")
!      Constant density star of mass 1 and radius 1.
       if (a<1.d0) then
          bgpot = 0.5d0*(a**2 - 3.d0)
       else
          bgpot = - 1.d0/a
       end if
    case ("iso")
       bgpot = 3.0d0*log(a)
    case ("isotrun")
       bgpot = (10.0d0/6.0d0)*( 2.0d0*atan(a)/a + log(a**2+1) )
    case ("nfw")
       bgpot = -16.0d0*log( 1.0d0+a )/a
    case ("burkert")
       bgpot = ( 10.0d0/(3.0d0*a) )*( 2.0d0*(1.0d0+a)*atan(a) -2.0d0*(1.0d0+a)*log(1.0d0+a) &
               -(1.0d0-a)*log(1.0d0+a**2) )
    case default
       bgpot = 0.0d0
    end select

  end function bgpot

  elemental real(8) function bgforce(x)

    real(8), intent(in) :: x
    real(8) :: a

    a = abs(x)

    select case (BGtype)
    case ("Isochrone")
       bgforce = -a/(sqrt(1.0d0+a**2)*(1.0d0+sqrt(1.0d0+a**2))**2)
    case ("Central")
       bgforce = -1.0d0/a**2
    case ("sphere")
       if (a<1.d0) then
          bgforce = - a
       else
          bgforce = - 1.d0/a**2
       end if
    case ("iso")
       bgforce = -3.0d0/a
    case ("isotrun")
       bgforce = -(10.0d0/3.0d0)*( a-atan(a) )/a**2
    case ("nfw")
       bgforce = -16.0d0*( log(1.0d0+a)-a/(1.0d0+a) )/a**2
    case ("burkert")
       bgforce = -( 10.0d0/(3.0d0*a*a) )*( log( (1.0d0+a**2)*(1.0d0+a)**2 ) - 2.0d0*atan(a) )
    case default
       bgforce = 0.0d0
    end select

!   Odd extension to r < 0.

    if (x < 0.0d0) bgforce = - bgforce

  end function bgforce



  !> Save all the data to the corresponding files
  !> Time series (energies, potential and force at the centre), written every
  !! spatial_output steps, when they are computed, independently of the
  !! snapshot cadence field_output (ASCII output; HDF5 keeps them as
  !! attributes of each snapshot).
  subroutine save_series

    character(100) :: filename

    filename = 'vlasov_energy'
    call save0Ddata(directory,filename,t,total_energy)
    filename = 'vlasov_k_phi_e'
    call save_energy(directory,filename,t,kinetic,potential,total_energy)

    if (autointeraction) then
       filename = 'vlasov_potential_r0'
       call save0Ddata(directory,filename,t,pot(1))
       filename = 'vlasov_force_r0'
       call save0Ddata(directory,filename,t,force(1))
    end if

  end subroutine save_series


  subroutine save_data

    character(100) :: filename     !< Name of output file

!   Save distribution function f.
    filename = 'vlasov_fdist'
    call save2Ddata_particles(directory,filename,Npart,t,r_part,p_part,l_part*f)

!   Save rho, curr and cont (multiplied by r**2).
!
!   NOTE: r, rho, avg_rho and curr are allocated with ghost points
!   (bounds 1-ghost:Nr, see alloc_mem_set0), but save1Ddata's dummy
!   arguments are explicit-shape (1:Nr).  Passing the whole arrays (or
!   an expression built from them) here associates them by sequence,
!   silently shifting the data by "ghost" points: the output would
!   start with the unphysical ghost points and drop the last "ghost"
!   physical points near r=rmax.  Slicing to (1:Nr) selects exactly
!   the physical points and matches the dummy's shape.
    filename = 'vlasov_density'
    call save1Ddata(directory,filename,Nr,t,r(1:Nr),r(1:Nr)**2*rho(1:Nr))
    filename = 'vlasov_avg_density'
    call save1Ddata(directory,filename,Nr,t,r(1:Nr),r(1:Nr)**2*avg_rho(1:Nr))
    filename = 'vlasov_curr'
    call save1Ddata(directory,filename,Nr,t,r(1:Nr),r(1:Nr)**2*curr(1:Nr))

!   Save force and potential.
    if (autointeraction) then

       filename = 'vlasov_force'
       call save1Ddata(directory,filename,Nr,t,r(1:Nr),force(1:Nr))

       filename = 'vlasov_potential'
       call save1Ddata(directory,filename,Nr,t,r(1:Nr),pot(1:Nr))

    end if

  end subroutine save_data


  subroutine save0Ddata(directory,filename,t,var)

! **********************
! ***   SAVE0DDATA   ***
! **********************

! This subroutine saves 1D data to files.

  implicit none

  real(8) t,var

  character(*) :: directory,filename
  character(20) :: filestatus


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

  character(*) :: directory,filename
  character(20) :: filestatus


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
     write(101,"(2ES16.8)") r(i),var(i)
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

  character(*) :: directory,filename
  character(20) :: filestatus


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
    write(101,"(3ES16.8)") r_part(i),p_part(i),var(i)
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



  subroutine save_energy(directory,filename,t,kinetic,potential,energy)

! **********************
! ***   SAVE ENERGY  ***
! **********************

! This subroutine saves the energy (kinectic,potential,virial) data to files.

  implicit none

  real(8) t,kinetic,potential,energy

  character(*) :: directory,filename
  character(20) :: filestatus


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







  !> Keep only the particles with keep(j) true, in their order, resizing every
  !! per-particle array. The potential and force of the kept particles are no
  !! longer valid afterwards: the caller must recompute them before they are
  !! used (main calls grav_force right after).
  subroutine remove_particles(keep)

    implicit none

    logical, intent(in) :: keep(:)
    integer :: n
    real(8), allocatable :: tmp(:)

    n = count(keep)
    if (n == Npart) return

    call pack_array(r_part)
    call pack_array(p_part)
    call pack_array(l_part)
    call pack_array(f)

    deallocate(r_part_p,p_part_p,p_part_h,pot_part,potself_part,force_part)
    allocate(r_part_p(n),p_part_p(n),p_part_h(n),pot_part(n),potself_part(n),force_part(n))
    r_part_p = 0.0d0
    p_part_p = 0.0d0
    p_part_h = 0.0d0
    pot_part = 0.0d0
    potself_part = 0.0d0
    force_part = 0.0d0

    print *, "Particles removed: ",Npart-n,", remaining: ",n
    Npart = n

  contains

    subroutine pack_array(a)
      real(8), allocatable, intent(inout) :: a(:)
      tmp = pack(a,keep)
      call move_alloc(tmp,a)
    end subroutine pack_array

  end subroutine remove_particles


  !> Discard the particles beyond rmax (reduceparticles=.true.).
  subroutine reduce_arrays

    implicit none

    call remove_particles(r_part <= rmax)

  end subroutine reduce_arrays


  !> Solve the Kepler-like equation Q = eta - ecc*sin(eta) for eta, with
  !! 0 <= Q < 2 pi and 0 <= ecc <= 1.
  !!
  !! First the plain Newton-Raphson from eta = Q that the code always used
  !! (the same operations, so the same result wherever it worked). It failed
  !! for ecc >~ 0.98: in 190 of the 10000 cases of test U6, with residuals up
  !! to 1e25, and it returned a plausible but wrong eta without notice
  !! (AUDITORIA_CIENTIFICA_2026-09-21.md, D3). So the result is checked,
  !! without a new evaluation of the residual (the analytic integrator inverts
  !! every particle at every step): it is accepted if the iteration stopped
  !! because the residual g at eta_k fell below tol, the last step |g/g'| was
  !! at most sqrt(tol), so that the residual at eta_k+1 is below about
  !! |g/g'|**2/2 <= tol/2, and eta lies in [Q - ecc, Q + ecc], where the root
  !! must be (|eta - Q| = ecc |sin(eta)|). Otherwise the equation is solved
  !! again with a safeguarded Newton: the left side minus Q grows
  !! monotonically with eta, the bracket shrinks with the sign of the residual
  !! at every iterate, and a step that would leave it is replaced by a
  !! bisection. Same routine as in vlasov-poisson_PIC.
  elemental real(8) function kepler_eta(Q,ecc,tol) result(eta)

    implicit none

    real(8), intent(in) :: Q,ecc,tol

    real(8) :: g,gp,lo,hi,eta_new
    integer :: it
    logical :: ok

    ok  = .false.
    eta = Q
    do it=1,50
      g  = eta - ecc*sin(eta) - Q
      gp = 1.d0 - ecc*cos(eta)
      eta = eta - g/gp
      if (abs(g) < tol) then
        ok = (abs(g/gp) <= sqrt(tol))
        exit
      end if
    end do

    if (ok .and. eta >= Q - ecc .and. eta <= Q + ecc) return

    lo  = Q - ecc
    hi  = Q + ecc
    eta = Q
    do it=1,200
      g  = eta - ecc*sin(eta) - Q
      if (g == 0.0d0) exit
      if (g < 0.0d0) then
        lo = max(lo,eta)
      else
        hi = min(hi,eta)
      end if
      gp = 1.d0 - ecc*cos(eta)
      eta_new = eta - g/gp
      if (.not. (eta_new >= lo .and. eta_new <= hi)) eta_new = 0.5d0*(lo+hi)
      eta = eta_new
      if (abs(g) < tol .or. hi - lo <= 4.0d0*epsilon(1.0d0)*max(abs(eta),1.0d0)) exit
    end do

  end function kepler_eta


  !> Invert the angle-action pair (Q,J) of a particle with angular momentum
  !! L back to (r,p_r), in the isochrone of unit mass and scale.
  !!
  !! J and L fix the energy, E = -1/(2 (J+c)**2) with c = (L+sqrt(L**2+4))/2,
  !! and the energy fixes the turning points. With s = 1+sqrt(1+r**2), the
  !! radial motion is s = (s1+s2)/2 - (s2-s1)/2 cos(eta), and the angle obeys
  !! the Kepler-like equation Q = eta - ecc sin(eta), solved by kepler_eta.
  !! p_r >= 0 on the way out, eta in [0,pi]. This is the inverse of the
  !! forward map in analysish.f90.
  subroutine invert_QJ_to_rp(Qv,Jv,Lv,rv,pv)

    implicit none

    real(8), intent(in)  :: Qv,Jv,Lv
    real(8), intent(out) :: rv,pv

    real(8) :: Eg,er1,er2,s1,s2,ecc,eta,sg,pv2,smallpi

    smallpi = acos(-1.0d0)

    Eg = -1.d0/(2.d0*(Jv+0.5d0*(Lv+sqrt(Lv**2+4.d0)))**2)

!   Every radicand is clamped at zero, as in analysish: the discriminant
!   vanishes on circular orbits and the inner turning point for L = 0, so
!   rounding alone can turn them negative. Without the clamp the circular
!   orbit J = 0 with L = 0.25 came out at r = 0, and with L = 0 the round
!   trip of test U6 returned a wrong angle for 650 of 1600 orbits (D4 of the
!   same audit).
    er1 = dsqrt(max((1.d0+Eg*(2.d0+Lv**2)-dsqrt(max(1.d0+2.d0*Eg*(2.d0+2.d0*Eg+Lv**2),0.0d0)))/(2.d0*Eg**2),0.0d0))
    er2 = dsqrt(max((1.d0+Eg*(2.d0+Lv**2)+dsqrt(max(1.d0+2.d0*Eg*(2.d0+2.d0*Eg+Lv**2),0.0d0)))/(2.d0*Eg**2),0.0d0))
    s1 = 1.d0 + sqrt(1.d0+er1**2)
    s2 = 1.d0 + sqrt(1.d0+er2**2)

!   Same eccentricity as in the forward map; it vanishes on circular orbits,
!   where rounding can make the radicand slightly negative.
    ecc = sqrt((-2.d0*Eg)**3)*sqrt(max(-Lv**2-2.d0*Eg-2.d0-0.5D0/Eg,0.0d0))/(-2.d0*Eg)

    eta = kepler_eta(modulo(Qv,2.0d0*smallpi),ecc,1.0d-14)

    sg = (s1+s2-cos(eta)*(s2-s1))/2.0d0
    rv = sqrt(max((sg-1.d0)**2-1.d0,0.0d0))

    pv2 = 2.d0*(Eg + 1.d0/(1.d0+dsqrt(1.d0+rv**2)) - 0.5d0*Lv**2/max(rv**2,1.0d-300))
    pv2 = sqrt(max(pv2,0.0d0))
    if (modulo(eta,2.0d0*smallpi) > smallpi) then
      pv = -pv2
    else
      pv =  pv2
    end if

  end subroutine invert_QJ_to_rp

  !> Store each particle's initial angle and radial action, for
  !! integrator="analytic". Uses the same forward map as analysish.f90, so
  !! it works from any initial state; unbound particles have no actions and
  !! stop the run.
  subroutine init_action_angle

    implicit none

    integer :: i,nunbound
    real(8) :: en,disc,er1,er2,s1,s2,ss,argaux,eta,smallpi

    smallpi = acos(-1.0d0)

    allocate(q0_part(1:Npart))
    allocate(j0_part(1:Npart))

    nunbound = 0

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(en,disc,er1,er2,s1,s2,ss,argaux,eta) REDUCTION(+:nunbound)
    do i=1,Npart
      en = -1.0d0/(1.0D0+dsqrt(1.0D0+r_part(i)**2)) + 0.5d0*l_part(i)**2/(r_part(i)**2) + 0.5D0*p_part(i)**2
      if (en >= 0.0d0) then
        nunbound = nunbound + 1
        q0_part(i) = 0.0d0
        j0_part(i) = 0.0d0
        cycle
      end if
!     Radicands clamped at zero and the phase of a circular orbit (s1 = s2)
!     set to a fixed value, as in analysish and in state "aa".
      disc = dsqrt(max(1.d0+2.d0*en*(2.d0+2.d0*en+l_part(i)**2),0.0d0))
      er1 = dsqrt(max((1.d0+en*(2.d0+l_part(i)**2)-disc)/(2.d0*en**2),0.0d0))
      er2 = dsqrt(max((1.d0+en*(2.d0+l_part(i)**2)+disc)/(2.d0*en**2),0.0d0))
      s1 = 1.d0 + dsqrt(1.d0+er1**2)
      s2 = 1.d0 + dsqrt(1.d0+er2**2)
      ss = 1.d0 + dsqrt(1.d0+r_part(i)**2)
      if (s2 > s1) then
        argaux = (s1+s2-2.0d0*ss)/(s2-s1)
      else
        argaux = 0.0d0
      end if
      if (p_part(i)>=0.d0) then
        eta = dacos(sign(min(abs(argaux),1.0D0),argaux))
      else
        eta = dacos(-sign(min(abs(argaux),1.0D0),argaux))+smallpi
      end if
      q0_part(i) = eta - dsqrt((-2.d0*en)**3)*sqrt(max(-l_part(i)**2-2.d0*en-2.d0-0.5D0/en,0.0d0))/(-2.d0*en)*dsin(eta)
      j0_part(i) = 1.d0/dsqrt(-2.d0*en)-0.5d0*(l_part(i)+dsqrt(l_part(i)**2+4.d0))
    end do
    !$OMP END PARALLEL DO

    if (nunbound > 0) then
      print *
      print *, 'integrator="analytic": ',nunbound,' particles are unbound (E >= 0)'
      print *, 'and have no angle-action variables.'
      print *, 'Aborting ...'
      print *
      stop 1
    end if

  end subroutine init_action_angle


  !> Advance every particle analytically to the absolute time "tnow".
  !!
  !! Without self-gravity, in the isochrone, J and L are conserved and the
  !! angle advances linearly, Q(t) = Q(0) + omega(J,L) t, with
  !! omega = 1/(J+c)**3 and c = (L+sqrt(L**2+4))/2. The state is recovered
  !! by inverting (Q,J) back to (r,p_r). Exact at any t, so it carries no
  !! phase error and is the reference for the symplectic integrators.
  subroutine advance_analytic(tnow)

    implicit none

    real(8), intent(in) :: tnow
    integer :: i
    real(8) :: om,Qt

    !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(om,Qt)
    do i=1,Npart
      om = 1.0d0/(j0_part(i)+0.5d0*(l_part(i)+sqrt(l_part(i)**2+4.d0)))**3
      Qt = q0_part(i) + om*tnow
      call invert_QJ_to_rp(Qt,j0_part(i),l_part(i),r_part(i),p_part(i))
    end do
    !$OMP END PARALLEL DO

  end subroutine advance_analytic

end module utils

