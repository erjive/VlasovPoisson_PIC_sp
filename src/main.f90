! ************************
! ***   MAIN PROGRAM   ***
! ************************

program VP_PIC

! Include modules

  use parameters
  use paramfile
  use arrays
  use utils
  use hdf5_io


! Declare variables.

  implicit none

  integer i,j,k,l       ! Counters

! Coefficients for the 4th-order symplectic integrator "yoshida4"
! (Yoshida, Phys. Lett. A 150, 262 (1990)): composes three
! leapfrog-like drift/kick sub-steps, with weights chosen so that the
! dt^3 local error term of a single 2nd-order (leapfrog) step cancels
! between them, leaving a dt^5 local (dt^4 global) error instead.

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

! Initialize time.

  t = 0.0d0

! **************************************
! ***   FIND DENSITY AND FLUX IN r   ***
! **************************************
  call density

! **************************************************
! ***   FIND GRAVITATIONAL POTENTIAL AND FORCE   ***
! **************************************************
  call grav_force()

! **************************************************************
! ***   FIND TOTAL NUMBER OF PARTICLES                       ***
! ***   AVERAGE KINETIC ENERGY, POTENTIAL AND TOTAL ENERGY   ***
! **************************************************************

!  call integrate
   call energy

! ************************
! *** INITIAL ANALYSIS ***
! ************************

   call analysish

!  ***************************
!  ***  INITIAL TIME STEP  ***
!  ***************************
  call set_timestep()

  print *
  print *, 'Time step fixed at size: ',dt

! *****************************
! ***   OUTUPUT TO SCREEN   ***
! *****************************

  print *
  print *,'------------------------------'
  print *,'|  Time step  |     Time     |'
  print *,'------------------------------'

  write(*,"(A5,I7,A5,ES11.4,A4)") ' |   ',0,'   | ',t,'  | '


! *********************************
! ***   SAVE THE INITIAL DATA   ***
! *********************************
   if (output_format=="hdf5") then
      call save_data_hdf5(0)
   else
      call save_data()
   end if


! *************************************
! ***   START MAIN EVOLUTION LOOP   ***
! *************************************

  do l=1,Nt

!    Time.

     t = t + dt

!    Save old time step.

      r_part_p = r_part

      p_part_p = p_part

 
!    Euler method (forward differencing in time, first order).

     if (integrator=='euler') then

       r_part = r_part_p + p_part*dt
       p_part = p_part_p + force_part*dt

       call grav_force()
!   Second order leapfrog method

    else if (integrator == 'leapfrog') then

!     Leapfrog integration 'kick-drift-kick' form

      p_part_h = p_part_p + force_part*dt*0.5D0
      r_part   = r_part_p + p_part_h * dt

      call  grav_force()

      p_part   = p_part_h + force_part*dt*0.5D0

!   Fourth order symplectic integrator (Yoshida 1990), composed of
!   three leapfrog-like drift-kick sub-steps -- see the coefficient
!   definitions near the top of this program.  Costs 4 force
!   evaluations per step (vs 1 for leapfrog: 3 for the sub-steps,
!   plus one more after the final drift so that force_part stays
!   synced with r_part on exit, matching what euler/leapfrog leave
!   behind), but should tolerate a larger dt for the same
!   energy-conservation accuracy.

    else if (integrator == 'yoshida4') then

      r_part = r_part_p + yg_c1*dt*p_part_p
      call grav_force()
      p_part = p_part_p + yg_d1*dt*force_part

      r_part = r_part + yg_c2*dt*p_part
      call grav_force()
      p_part = p_part + yg_d2*dt*force_part

      r_part = r_part + yg_c3*dt*p_part
      call grav_force()
      p_part = p_part + yg_d3*dt*force_part

      r_part = r_part + yg_c4*dt*p_part
      call grav_force()

!    Fourth order Runge-Kutta.

     else if (integrator=='rk4') then

        print *
        print *, 'Fourth order Runge-Kutta not yet implemented.'
        print *, 'Aborting ...'
        print *
        stop 1

!    Unknown integration method.

     else

        print *, 'Unknown integration method.'
        print *, 'Aborting ...'
        print *
        stop 1

     end if

!   At the origin impose symmetry condition f(r,p) = f(-r,-p)

    do i=1,Npart

      if (rmin == 0 .and. r_part(i)<0.d0) then

        r_part(i) = -r_part(i)
        p_part(i) = -p_part(i)

      end if

    end do

!    **************************************
!    ***   FIND DENSITY AND FLUX IN r   ***
!    **************************************

     if (forcetype=="self") then

!    Save old value of rho and curr.
        rho_p  = rho
!        curr_p = curr

!    Integrate over phase space.
        !call density

     else

        if (mod(l+1,spatial_output).eq.0) then

!    Integrate over phase space.
!           call density

!    Save old value of rho and curr in order to calculate the continuity equation
!           rho_p  = rho
!           curr_p = curr

        end if

        if (mod(l,spatial_output).eq.0) then

!    Integrate over phase space.
           call density
           call energy
        end if

     end if


! ******************************************************************
! ***   FIND GRAVITATIONAL POTENTIAL AND FORCE ON THE PARTICLES  ***
! ******************************************************************
     !call grav_force()


!    **********************************
!    ***   SOLVE POISSON EQUATION   ***
!    **********************************

!    For the self gravitating case solve
!    the Poisson equation again.

!     if (forcetype=="self") then
!        !call poisson
!     end if



!    *****************************************
!    ***   CALCULATE CONTINUITY EQUATION   ***
!    *****************************************

!    The continuity equation has the form:
!
!    cont  =  0  =  d(rho)/dt + div(curr)  = d(rho)/dt + (1/r**2) d(r**2 curr)/dr
!
!                =  d(rho)/dt + d(curr)/dr + 2 curr / r
!
!    Notice that this should converge to zero.  The expression below is only
!    second order accurate.

!     if (mod(l,spatial_output).eq.0) then
!        do i=1,Nr-1
!           cont(i) = (rho(i) - rho_p(i))/dt &
!                + 0.25d0*(curr(i+1) + curr_p(i+1) - curr(i-1) - curr_p(i-1))/dr &
!                + (curr(i) + curr_p(i))/r(i)
!        end do
!     end if


!    ***************************
!    ***   ADAPT TIME STEP   ***
!    ***************************

!    For the self-gravitating case the force can change with time
!    (e.g. it grows as the cloud collapses, see the "compactness"
!    runs in the article), so the time step needs to adapt -- it was
!    otherwise only ever computed once, from the *initial* force,
!    before the main loop even starts.  Notice that the time step can
!    go up and down in response to the size of the force.  This uses
!    force_part right after grav_force() was called for the new
!    r_part above, so Fmax reflects the force at the position the
!    particles were just moved to, and the resulting dt is the one
!    used to advance the *next* step.
!
!    This was originally guarded by forcetype=="self", but forcetype
!    never actually gets set to "self" anywhere meaningful (see the
!    "forcetype=self is a no-op" item in BUGS_TODO.md) -- the flag
!    that actually controls whether the force can change over time is
!    autointeraction.

     if (autointeraction) then
       call set_timestep()
     end if


!    *****************************
!    ***   SAVE DATA TO FILE   ***
!    *****************************

     if (mod(l,spatial_output).eq.0) then

       if (output_format=="hdf5") then
          call save_data_hdf5(l)
       else
          call save_data()
       end if

     end if

     if (mod(l,spatial_output).eq.0) then

        call analysish

     end if

!    *************************************************
!    ***   IF POSSIBLE REDUCE SIZE OF THE ARRAYS   ***
!    *************************************************

     if (reduceparticles .and. (mod(l,Nreduce).eq.0)) then
!     if (mod(l,time_output).eq.0) then
       call reduce_arrays

     end if

!    ***********************************
!    ***   END MAIN EVOLUTION LOOP   ***
!    ***********************************

!    Time step information to screen.

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
