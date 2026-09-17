! ===========================================================================
! main.f90
! ===========================================================================
!> Particle-in-cell evolution of the spherically symmetric Vlasov-Poisson
!! system with a distribution of angular momentum: each particle carries
!! (r, p_r, L) and a phase-space weight f, moves in a background potential
!! plus, optionally, its own gravity, and L is conserved.
!!
!!     ./VP_PIC [file] [name=value] ...      (see paramfile.f90)

program VP_PIC

  use parameters
  use paramfile
  use arrays
  use utils
  use hdf5_io

  implicit none

  integer i,j,k,l       ! Counters
  logical :: sync       ! Does this step need the force at its final positions?

! Coefficients of the 4th-order symplectic integrator "yoshida4" (Yoshida,
! Phys. Lett. A 150, 262 (1990)): three drift-kick sub-steps whose weights
! cancel the dt^3 local error of a single leapfrog step, leaving a dt^5 local
! (dt^4 global) error.

  real(8), parameter :: yg_cbrt2 = 1.2599210498948732d0            ! 2**(1/3)
  real(8), parameter :: yg_w1    = 1.d0/(2.d0-yg_cbrt2)
  real(8), parameter :: yg_w0    = -yg_cbrt2/(2.d0-yg_cbrt2)
  real(8), parameter :: yg_c1 = yg_w1*0.5d0,      yg_c4 = yg_c1     ! drift weights
  real(8), parameter :: yg_c2 = (yg_w0+yg_w1)*0.5d0, yg_c3 = yg_c2
  real(8), parameter :: yg_d1 = yg_w1,  yg_d2 = yg_w0,  yg_d3 = yg_w1 ! kick weights


  call read_parameters()

  call set_grid_size()

  call alloc_mem_set0()

  call construct_grid()

  call initial_data()

! integrator="analytic" advances each particle with the closed-form solution
! Q(t) = Q(0) + omega(J,L) t, valid only where J and L are exact constants
! of motion: static isochrone background, no self-gravity, and eps = 0,
! since the angle-action map is built from the unsoftened L**2/(2 r**2).

  if (integrator == 'analytic') then

     if (autointeraction .or. BGtype /= "Isochrone" .or. eps /= 0.0d0) then
        print *
        print *, 'integrator="analytic" requires an integrable setup:'
        print *, '  autointeraction = .false.,  BGtype = "Isochrone",  eps = 0.'
        print *, 'Aborting ...'
        print *
        stop 1
     end if

!    reduce_arrays resizes the particle arrays but not q0_part/j0_part.
     if (reduceparticles) then
        print *
        print *, 'integrator="analytic" is incompatible with reduceparticles=.true.'
        print *, 'Aborting ...'
        print *
        stop 1
     end if

     call init_action_angle()

  end if

! ***************************
! ***   OUTPUT DIRECTORY  ***
! ***************************

! Create output directory, copy the parameter file to it and write the
! complete configuration in force (file plus command line overrides) to
! params_usados.par, which is itself a valid input file.

  call system('mkdir -p '//trim(directory))
  call system('cp '//trim(parameter_file)//' '//trim(directory))
  call dump_parameters()

  if (output_format=="hdf5") call open_hdf5_file()

  t = 0.0d0

! *********************************************
! ***   INITIAL DENSITY, FORCE AND ENERGY   ***
! *********************************************

  call density

  call grav_force()

  call energy

  call analysish

! The time step is set once, from the initial force (see set_timestep).

  call set_timestep()

  print *
  print *, 'Time step fixed at size: ',dt

! The Courant bound in set_timestep uses the parameter pmax, not the momenta
! the particles actually have, and the angle-action variables of analysish
! and of the states aa and aa_quad are those of the isochrone. Say so when
! the run steps outside those assumptions.

  if (maxval(abs(p_part)) > pmax) then
     print *
     print *, 'WARNING: particles reach |p| = ',maxval(abs(p_part)),' > pmax = ',pmax
     print *, '         The Courant bound dt <= courant*dr/pmax does not hold for them.'
  end if

  if (BGtype /= "Isochrone") then
     print *
     print *, 'WARNING: BGtype = ',trim(BGtype),'. h_k (analysish) and the states aa and'
     print *, '         aa_quad use the angle-action variables of the isochrone.'
  end if

  print *
  print *,'------------------------------'
  print *,'|  Time step  |     Time     |'
  print *,'------------------------------'

  write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t,'  | '

  if (output_format=="hdf5") then
     call save_data_hdf5(0)
  else
     call save_data()
     call save_series()
  end if


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do l=1,Nt

     t = t + dt

!    yoshida4 and analytic never use the force at the end of a step to move
!    the particles (the next step starts with a drift), only the potential
!    and force of the diagnostics written every spatial_output steps, and
!    set_timestep when dt_switch="var". Skipping that evaluation otherwise
!    saves one of the four force evaluations per yoshida4 step.

     sync = (mod(l,spatial_output) == 0) .or. (dt_switch == "var")

!    Only euler and leapfrog read the values at the start of the step;
!    yoshida4 updates in place.

     if (integrator=='euler' .or. integrator=='leapfrog') then
       r_part_p = r_part
       p_part_p = p_part
     end if

     if (integrator=='euler') then

!      Forward differencing in time, first order.

       r_part = r_part_p + p_part*dt
       p_part = p_part_p + force_part*dt

       call grav_force()

     else if (integrator == 'leapfrog') then

!      Second order, kick-drift-kick.

       p_part_h = p_part_p + force_part*dt*0.5D0
       r_part   = r_part_p + p_part_h * dt

       call  grav_force()

       p_part   = p_part_h + force_part*dt*0.5D0

     else if (integrator == 'yoshida4') then

!      Fourth order symplectic (coefficients above): drift, then three
!      kick-drift pairs. Each kick is fused with the drift that follows it
!      into one parallel pass over the particles. Three force evaluations per
!      step, plus one after the final drift on the steps that need it (sync).

       !$OMP PARALLEL DO SCHEDULE(STATIC)
       do i=1,Npart
         r_part(i) = r_part(i) + yg_c1*dt*p_part(i)
       end do
       !$OMP END PARALLEL DO
       call grav_force()

       !$OMP PARALLEL DO SCHEDULE(STATIC)
       do i=1,Npart
         p_part(i) = p_part(i) + yg_d1*dt*force_part(i)
         r_part(i) = r_part(i) + yg_c2*dt*p_part(i)
       end do
       !$OMP END PARALLEL DO
       call grav_force()

       !$OMP PARALLEL DO SCHEDULE(STATIC)
       do i=1,Npart
         p_part(i) = p_part(i) + yg_d2*dt*force_part(i)
         r_part(i) = r_part(i) + yg_c3*dt*p_part(i)
       end do
       !$OMP END PARALLEL DO
       call grav_force()

       !$OMP PARALLEL DO SCHEDULE(STATIC)
       do i=1,Npart
         p_part(i) = p_part(i) + yg_d3*dt*force_part(i)
         r_part(i) = r_part(i) + yg_c4*dt*p_part(i)
       end do
       !$OMP END PARALLEL DO
       if (sync) call grav_force()

     else if (integrator == 'analytic') then

!      Exact advance in angle-action variables (no self-gravity only).

       call advance_analytic(t)
       if (sync) call grav_force()

     end if

!    At the origin impose the symmetry condition f(r,p) = f(-r,-p).

     if (rmin == 0) then
       !$OMP PARALLEL DO SCHEDULE(STATIC)
       do i=1,Npart
         if (r_part(i)<0.d0) then
           r_part(i) = -r_part(i)
           p_part(i) = -p_part(i)
         end if
       end do
       !$OMP END PARALLEL DO
     end if

!    Density and energies for the output. With self-gravity the averaged
!    density that Poisson needs is recomputed inside grav_force.

     if (mod(l,spatial_output).eq.0) then
        call density
        call energy
     end if

!    By default (dt_switch="fix") the time step computed from the initial
!    force is kept for the whole run: a step that changes from one step to
!    the next breaks the symplectic character of leapfrog and yoshida4, and
!    the energy error stops being bounded. Choose dt with a convergence test
!    instead. dt_switch="var" recomputes dt from the force at the new
!    positions, for runs where the force grows a lot (e.g. a collapsing
!    cloud); the new dt advances the next step.

     if (dt_switch == "var") then
       call set_timestep()
     end if

!    field_output gates the particle and grid snapshot, the bulk of the
!    disk footprint, so h_k (analysish, every spatial_output steps) can be
!    sampled finely without an equally frequent snapshot.

     if (mod(l,field_output).eq.0) then

       if (output_format=="hdf5") then
          call save_data_hdf5(l)
       else
          call save_data()
       end if

     end if

     if (mod(l,spatial_output).eq.0) then
        call analysish
        if (output_format/="hdf5") call save_series()
     end if

!    Discard particles beyond rmax.

     if (reduceparticles .and. (mod(l,Nreduce).eq.0)) then
       call reduce_arrays
     end if

     if (mod(l,time_output).eq.0) then
        write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',l,'   | ',t,'  | '
     end if

  end do

  print *,'------------------------------'


! ***************
! ***   END   ***
! ***************

  print *, 'Minimum radii of particles = ', minval(r_part)
  print *, 'Maximum radii of particles = ', maxval(r_part)
  print *, 'Minimum momenta of particles = ', minval(p_part)
  print *, 'Maximum momenta of particles = ', maxval(p_part)
  print *, 'Minimum ang. momenta of particles = ', minval(l_part)
  print *, 'Maximum ang. momenta of particles = ', maxval(l_part)

  print *
  print *, 'PROGRAM HAS FINISHED'
  print *
  print *, 'Have a nice day!'
  print *
  print *
  print *

  if (output_format=="hdf5") call close_hdf5_file()

  call deallocate_mem()

end program VP_PIC
