! ===========================================================================
! initial_data.f90
! ===========================================================================
!> Initial particles: position r, radial momentum p, angular momentum L and
!! phase-space weight f, on a regular grid (states gaussian1, aa, aa_quad) or
!! read from a file (checkpoint). The mass is normalized to a0.


  subroutine initial_data

    use parameters
    use arrays
    use utils
    use distribution

    implicit none

    integer :: i,j,k,indx
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

    smallpi = acos(-1.0d0)

! Cells of the support in (r,p,L). The gaussian1 profile adds the mirror
! image at (-r0,-p0) so that f(-r,-p) = f(r,p) holds at the origin.

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

            raux = rminc+(dble(i))*drc
            paux = pminc+dble(j)*dpc
            laux = lminc+dble(k)*dlc
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

      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do indx=1,Npart
         if (f(indx)<= cutoff*f_max) then
            f(indx)=0.0D0
            r_part(indx) = 1000000.D0
         end if
      end do
      !$OMP END PARALLEL DO

      call reduce_arrays


      f = a0/(8.D0*smallpi**2*drc*dpc*dlc*sum(f*l_part))*f

      print *, "Initial total mass=", sum(f*l_part)*8.0*smallpi**2*drc*dpc*dlc
    
    else if(state .eq."aa") then

!     Regular grid in (r,p,L), indexed as in gaussian1, with f evaluated at
!     the isochrone angle-action variables of each node.

      !$OMP PARALLEL DO SCHEDULE(GUIDED) &
      !$OMP PRIVATE(i,j,k,indx,raux,paux,laux,energy,er1,er2,s1,s2,s,argaux,eta,Q3,J3) &
      !$OMP COLLAPSE(3)
      do k=1,Nlc
        do i=1,Nrc
          do j=1,Npc

            indx = (k-1)*Nrc*Npc + (i-1)*Npc + j

            raux = rminc+dble(i)*drc
            paux = pminc+dble(j)*dpc
            laux = lminc+dble(k)*dlc

            r_part(indx) = raux
            p_part(indx) = paux
            l_part(indx) = laux

            energy = -1.0/(1.0D0+dsqrt(1.0D0+raux**2)) + 0.5d0*laux**2/(raux**2) + 0.5D0*paux**2
            er1 = dsqrt((1.d0+energy*(2.d0+laux**2)-dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+laux**2)))/(2.d0*energy**2))
            er2 = dsqrt((1.d0+energy*(2.d0+laux**2)+dsqrt(1.d0+2.d0*energy*(2.d0+2.d0*energy+laux**2)))/(2.d0*energy**2))
            s1 = 1.d0 + sqrt(1.d0+er1**2)
            s2 = 1.d0 + sqrt(1.d0+er2**2)
            s  = 1.d0 + sqrt(1.d0+raux**2)
            argaux = (s1+s2-2.0*s)/(s2-s1)

            if (paux>=0.d0) then
              eta = dacos(sign(min(abs(argaux),1.0),argaux))

            else
              eta = dacos(-sign(min(abs(argaux),1.0),argaux))+smallpi
            end if

            Q3 = eta - sqrt((-2.d0*energy)**3)*sqrt(-laux**2-2.d0*energy-2.d0-0.5D0/energy)/(-2.d0*energy)*sin(eta)

            J3 = 1.d0/sqrt(-2.d0*energy)-0.5d0*(laux+sqrt(laux**2+4.d0))

            f(indx) = exp(-sin(0.5d0*Q3)**2/sp**2)*exp(-J3**2/sr**2)*J3**2*exp(-(laux-l0)**2/sl**2)

            if (f(indx) /= f(indx) ) then
              f(indx) = 0.D0
              r_part(indx) = 10000.D0
            end if
          end do
        end do
      end do
      !$OMP END PARALLEL DO

      f_max = maxval(f)

      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do indx=1,Npart
         if (f(indx)<= cutoff*f_max) then
            f(indx)=0.0D0
            r_part(indx) = 1000000.D0
         end if
      end do
      !$OMP END PARALLEL DO


      call reduce_arrays


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

      !$OMP PARALLEL DO SCHEDULE(GUIDED)
      do indx=1,Npart
         if (f(indx)<= cutoff*f_max) then
            f(indx)=0.0D0
            r_part(indx) = 1000000.D0
         end if
      end do
      !$OMP END PARALLEL DO

      call reduce_arrays

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


