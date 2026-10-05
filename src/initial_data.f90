! ===========================================================================
! initial_data.f90
! ===========================================================================
!> Initial particles: position r, radial momentum p, angular momentum L and
!! phase-space weight f, on a regular grid (states gaussian1, aa, aa_quad),
!! on a grid with displaced nodes (aa_halton), drawn at random (aa_random)
!! or read from a file (checkpoint). The mass is normalized to a0.


  subroutine initial_data

    use parameters
    use arrays
    use utils
    use distribution

    implicit none

    integer :: i,j,k,indx,nunbound
    real(8) :: smallpi,f_max
    real(8) :: raux,paux,laux
    real(8) :: gaussian
    real(8) :: energy

!   Auxiliary variables for a distribution function
!   that depends on angle-action variables
    real(8) :: J3, Q3, s, s1, s2, er1, er2, eta, argaux
!   Midpoint grid of state aa_quad
    real(8) :: djc, dqc
!   State checkpoint
    integer :: unit_ic, ios
!   States aa_halton and aa_random
    logical :: bound
    real(8) :: rand4(4), fl_max
    integer(8) :: ntry

    smallpi = acos(-1.0d0)

! Cells of the support in (r,p,L), with the nodes at the cell midpoints. The
! gaussian1 profile adds the mirror image at (-r0,-p0) so that
! f(-r,-p) = f(r,p) holds at the origin.

! Find the size of the cell

    drc = (rmaxc-rminc)/dble(Nrc)
    dpc = (pmaxc-pminc)/dble(Npc)
    dlc = (lmaxc-lminc)/dble(Nlc)

    print *, "(drc,dpc,dlc)=",drc,dpc,dlc
    if(state.eq."gaussian1") then

!     Regular grid in (r,p,L): cell (k,i,j) is particle
!     indx = (k-1)*Nrc*Npc + (i-1)*Npc + j, computed directly so that every
!     thread writes only its own entries.

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,j,k,indx,raux,paux,laux) COLLAPSE(3)
      do k=1,Nlc
      do i=1,Nrc
        do j=1,Npc

            indx = (k-1)*Nrc*Npc + (i-1)*Npc + j

            raux = rminc+(dble(i)-0.5d0)*drc
            paux = pminc+(dble(j)-0.5d0)*dpc
            laux = lminc+(dble(k)-0.5d0)*dlc
            r_part(indx)  = raux
            p_part(indx)  = paux
            l_part(indx)  = laux
            f(indx)       = gaussian(1.0D0,r0,p0,l0,raux,paux,laux,sr,sp,sl)

          end do
        end do
      end do
      !$OMP END PARALLEL DO

      f_max = maxval(f)

!     Cutoff filtering only depends on the particle's own indx, not
!     on i/j/k individually, so it's just a flat loop over 1:Npart.

!     Nodes where f is not above cutoff*f_max carry no mass and are removed.
      call remove_particles(f > cutoff*f_max)



      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc
    
    else if(state .eq."aa") then

!     Regular grid in (r,p,L), indexed as in gaussian1, with f evaluated at
!     the isochrone angle-action variables of each node.

      nunbound = 0

      !$OMP PARALLEL DO SCHEDULE(GUIDED) &
      !$OMP PRIVATE(i,j,k,indx,raux,paux,laux,energy,er1,er2,s1,s2,s,argaux,eta,Q3,J3) &
      !$OMP REDUCTION(+:nunbound) COLLAPSE(3)
      do k=1,Nlc
        do i=1,Nrc
          do j=1,Npc

            indx = (k-1)*Nrc*Npc + (i-1)*Npc + j

            raux = rminc+(dble(i)-0.5d0)*drc
            paux = pminc+(dble(j)-0.5d0)*dpc
            laux = lminc+(dble(k)-0.5d0)*dlc

            r_part(indx) = raux
            p_part(indx) = paux
            l_part(indx) = laux

            energy = -1.0d0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*laux**2/(raux**2) + 0.5D0*paux**2

!           An unbound node has no angle-action variables and carries no mass
!           in this distribution. It used to be caught further down by a NaN
!           test, which also swallowed the rounding NaN of a near-circular
!           orbit -- a perfectly valid node that was being given zero mass in
!           silence. Now the two are separated and the count is reported
!           (AUDITORIA_2026-09-20.md, point 5).

            if (energy >= 0.0d0) then

              f(indx) = 0.0d0
              nunbound = nunbound + 1

            else

!           Every radicand is clamped at zero: all of them vanish on a
!           circular orbit, so rounding alone can turn them negative.

            er1 = dsqrt(max(1.d0+2.d0*energy*(2.d0+2.d0*energy+laux**2),0.0d0))
            er2 = dsqrt(max((1.d0+energy*(2.d0+laux**2)+er1)/(2.d0*energy**2),0.0d0))
            er1 = dsqrt(max((1.d0+energy*(2.d0+laux**2)-er1)/(2.d0*energy**2),0.0d0))
            s1 = 1.d0 + sqrt(1.d0+er1**2)
            s2 = 1.d0 + sqrt(1.d0+er2**2)
            s  = 1.d0 + sqrt(1.d0+raux**2)
            if (s2 > s1) then
              argaux = (s1+s2-2.0d0*s)/(s2-s1)
            else
              argaux = 0.0d0
            end if

            if (paux>=0.d0) then
              eta = dacos(sign(min(abs(argaux),1.0d0),argaux))

            else
              eta = dacos(-sign(min(abs(argaux),1.0d0),argaux))+smallpi
            end if

            Q3 = eta - sqrt((-2.d0*energy)**3) &
                 *sqrt(max(-laux**2-2.d0*energy-2.d0-0.5D0/energy,0.0d0))/(-2.d0*energy)*sin(eta)

            J3 = 1.d0/sqrt(-2.d0*energy)-0.5d0*(laux+sqrt(laux**2+4.d0))

            f(indx) = exp(-sin(0.5d0*Q3)**2/sp**2)*exp(-J3**2/sr**2)*J3**2*exp(-(laux-l0)**2/sl**2)

!           Last resort: with the clamps above no NaN should reach this point.
            if (f(indx) /= f(indx)) f(indx) = 0.D0

            end if

          end do
        end do
      end do
      !$OMP END PARALLEL DO

      if (nunbound > 0) then
        print *
        print *, 'state="aa": ',nunbound,' of ',Npart,' nodes are unbound (E >= 0)'
        print *, '            and were given zero mass.'
      end if

      f_max = maxval(f)

!     Nodes where f is not above cutoff*f_max carry no mass and are removed.
      call remove_particles(f > cutoff*f_max)




      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc


    else if(state .eq."aa_halton") then

!     The grid of "aa", with each node displaced inside its own cell by a
!     three-dimensional Halton sequence (bases 2, 3 and 5 for r, p and L),
!     and F0 = df0(Q,J,L) of dftype evaluated at the displaced node. The
!     displacement breaks the regular lattice in the actions, which makes all
!     particles rephase at the same time (the recurrence of h_k), while the
!     sample stays more uniform than with random displacements. It also
!     gives up the convergence of a quadrature: the error of h_k goes back
!     to the Monte Carlo rate. Kept for comparisons; aa_quad is the state to
!     use (same state as in vlasov-poisson_PIC, there in two dimensions).

      call df0_report

      nunbound = 0

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,j,k,indx,raux,paux,laux,Q3,J3,bound) &
      !$OMP REDUCTION(+:nunbound) COLLAPSE(3)
      do k=1,Nlc
        do i=1,Nrc
          do j=1,Npc

            indx = (k-1)*Nrc*Npc + (i-1)*Npc + j

            raux = rminc+(dble(i)-0.5d0)*drc + (halton(indx,2)-0.5d0)*drc
            paux = pminc+(dble(j)-0.5d0)*dpc + (halton(indx,3)-0.5d0)*dpc
            laux = lminc+(dble(k)-0.5d0)*dlc + (halton(indx,5)-0.5d0)*dlc

            r_part(indx) = raux
            p_part(indx) = paux
            l_part(indx) = laux

            call rp_to_QJ(raux,paux,laux,Q3,J3,bound)

            if (bound) then
              f(indx) = df0(Q3,J3,laux)
            else
              f(indx) = 0.0d0
              nunbound = nunbound + 1
            end if

          end do
        end do
      end do
      !$OMP END PARALLEL DO

      if (nunbound > 0) then
        print *
        print *, 'state="aa_halton": ',nunbound,' of ',Npart,' nodes are unbound (E >= 0)'
        print *, '                   and were given zero mass.'
      end if

      f_max = maxval(f)

!     Nodes where f is not above cutoff*f_max carry no mass and are removed.
      call remove_particles(f > cutoff*f_max)

      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc

    else if(state .eq."aa_random") then

!     Monte Carlo sample of F0 = df0(Q,J,L) of dftype: Npart = Nrc*Npc*Nlc
!     particles of equal mass, drawn by acceptance-rejection in the box
!     [rminc,rmaxc] x [pminc,pmaxc] x [lminc,lmaxc].
!
!     The mass in a cell is 8 pi**2 F L dr dp dL, so the particles are drawn
!     with number density proportional to F L: a point of the box, uniform,
!     is accepted with probability F L/(F_max lmaxc), with F_max a bound of
!     F (df0_max). Equal masses then mean f L constant, that is f = 1/L up
!     to the normalization below; where the particles land carries the shape
!     of the distribution, and a weight F on top of it would represent F**2.
!     Unbound points (E >= 0) have no actions and F = 0 there.
!
!     The error of any diagnostic goes as 1/sqrt(Npart). The sample is
!     repeatable: it only depends on "seed" (see init_rng), and the loop is
!     serial so that it does not depend on the number of threads either.

      call df0_report
      call init_rng

      fl_max = df0_max()*lmaxc

      ntry = 0
      i = 1
      do while (i <= Npart)

        call random_number(rand4)
        ntry = ntry + 1

        raux = rminc + (rmaxc-rminc)*rand4(1)
        paux = pminc + (pmaxc-pminc)*rand4(2)
        laux = lminc + (lmaxc-lminc)*rand4(3)

        if (raux > 0.0d0 .and. laux > 0.0d0) then
          call rp_to_QJ(raux,paux,laux,Q3,J3,bound)
          if (bound) then
            if (rand4(4)*fl_max <= df0(Q3,J3,laux)*laux) then
              r_part(i) = raux
              p_part(i) = paux
              l_part(i) = laux
              f(i)      = 1.0d0/laux
              i = i + 1
            end if
          end if
        end if

!       Where F0 is negligible in the whole box the loop would never end.
        if (ntry == 10000000_8 .and. i <= 10) then
          print *
          print *, 'state="aa_random": fewer than 10 points accepted in 1e7 draws.'
          print *, 'The box [rminc,rmaxc] x [pminc,pmaxc] x [lminc,lmaxc] holds almost'
          print *, 'none of the distribution.'
          print *, 'Aborting ...'
          print *
          stop 1
        end if

      end do

      print *, 'state="aa_random": seed = ',seed,', draws = ',ntry, &
               ', accepted fraction = ',dble(Npart)/dble(ntry)

      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc

    else if(state .eq."aa_quad") then

!     Tensor-product quadrature directly in the angle-action variables:
!     midpoints of a regular grid in J (Nrc cells in [jminc,jmaxc]),
!     Q (Npc cells in [0,2pi)) and L (Nlc cells in [lminc,lmaxc]), each
!     node mapped back to (r,p_r) with the isochrone map. No jitter: the
!     ordered nodes are what makes this a quadrature, whose error in h_k
!     stays at the level of the integrator, while any displacement leaves
!     the 1/sqrt(N) Monte Carlo rate (vlasov-poisson_PIC, commit c4df269).
!
!     Without self-gravity J and L are conserved and Q = Q0 + omega(J,L) t,
!     so h_k(t) is a quadrature of an increasingly oscillatory integral.
!     The grid resolves it while, in each direction,
!        N_J >~ 4 k Delta_omega_J t_max / 2pi,  N_L >~ 4 k Delta_omega_L t_max / 2pi,
!     with Delta_omega the spread of omega = 1/(J+c)**3 across the support.
!     In Q the rule is the periodic trapezoid and converges spectrally.
!
!     Since dQ dJ = dr dp_r, every node carries the same weight dQ dJ dL;
!     the constant is absorbed by the mass normalization below, which uses
!     drc*dpc*dlc like every diagnostic, so f holds the raw distribution.

      call df0_report

      djc = (jmaxc-jminc)/dble(Nrc)
      dqc = 2.0d0*smallpi/dble(Npc)
      print *, "aa_quad: (dJ,dQ,dL)=",djc,dqc,dlc

      !$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(i,j,k,indx,raux,paux,laux,Q3,J3) COLLAPSE(3)
      do k=1,Nlc
        do i=1,Nrc
          do j=1,Npc

            indx = (k-1)*Nrc*Npc + (i-1)*Npc + j

            J3   = jminc + (dble(i)-0.5d0)*djc
            Q3   = (dble(j)-0.5d0)*dqc
            laux = lminc + (dble(k)-0.5d0)*dlc

            call invert_QJ_to_rp(Q3,J3,laux,raux,paux)

            r_part(indx) = raux
            p_part(indx) = paux
            l_part(indx) = laux
            f(indx) = df0(Q3,J3,laux)

          end do
        end do
      end do
      !$OMP END PARALLEL DO

      f_max = maxval(f)

!     Nodes where f is not above cutoff*f_max carry no mass and are removed.
      call remove_particles(f > cutoff*f_max)


      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc

    else if(state .eq."checkpoint") then

!     Particles read from a file, one line "r p_r L f" per particle, exactly
!     Npart = Nrc*Npc*Nlc lines. This is how a state the code cannot build
!     enters a run: a self-consistent equilibrium plus a perturbation, or a
!     snapshot of another run (HDF5 output: r_part, p_part, l_part, fl/l_part).
!
!     f holds raw values of the distribution function on cells of equal
!     weight, as in aa_quad, so the same normalization to the mass a0
!     applies; the diagnostics multiply by the same drc*dpc*dlc, which
!     cancels.

      open(newunit=unit_ic,file=trim(CheckPointfile),status='old',action='read',iostat=ios)
      if (ios /= 0) then
         print *
         print *, 'state="checkpoint": cannot open ',trim(CheckPointfile)
         print *, 'Aborting ...'
         print *
         stop 1
      end if
      do i=1,Npart
         read(unit_ic,*,iostat=ios) r_part(i),p_part(i),l_part(i),f(i)
         if (ios /= 0) then
            print *
            print *, 'state="checkpoint": ',trim(CheckPointfile), &
                     ' has fewer than Npart=Nrc*Npc*Nlc =',Npart,' valid lines'
            print *, 'Aborting ...'
            print *
            stop 1
         end if
      end do
      read(unit_ic,*,iostat=ios) raux
      if (ios == 0) then
         print *
         print *, 'state="checkpoint": ',trim(CheckPointfile), &
                  ' has more than Npart=Nrc*Npc*Nlc =',Npart,' lines'
         print *, 'Aborting ...'
         print *
         stop 1
      end if
      close(unit_ic)

      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Read",Npart," particles from ",trim(CheckPointfile)
      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc

    else

!     read_parameters only accepts the states above; this guards against a
!     new option added there but not here.
      print *
      print *, 'Unknown initial state: ',trim(state)
      print *, 'Aborting ...'
      print *
      stop 1

    endif

  end subroutine initial_data

  function gaussian(a0,r0,p0,l0,x,y,z,sr,sp,sl)

  implicit none

  real(8) a0
  real(8) gaussian
  real(8) r0,p0,l0,x,y,z,sr,sp,sl
  real(8) smallpi

  smallpi = acos(-1.0d0)

  gaussian = dble(a0)*&
             (dexp(-(x-r0)**2/sr**2)*dexp(-(y-p0)**2/sp**2)+ &
              dexp(-(x+r0)**2/sr**2)*dexp(-(y+p0)**2/sp**2))*&
              dexp(-(z-l0)**2/sl**2)

  end function gaussian


