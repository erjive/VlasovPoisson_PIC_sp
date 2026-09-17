! ===========================================================================
! initial_data.f90
! ===========================================================================
!> Here are initialized all the functions defined on the grid.


  subroutine initial_data

    use parameters
    use arrays
    use utils

    implicit none

    logical :: accepted
    integer :: i,j,k,indx
    real(8) :: smallpi,f_max
    real(8) :: raux,paux,laux
    real(8) :: rand3(3)
    real(8) :: gaussian
    real(8) :: w,x,y,z
    real(8) :: energy 

!   Auxiliary variables for a distribution function 
!   that depends on action-angle varialbes
    real(8) :: J3, Q3, w3, s, s1, s2, er1, er2, eta, argaux
!   Midpoint grid of state aa_quad
    real(8) :: djc, dqc

    smallpi = acos(-1.0d0)

! 
! For a fixed value of L, We generate particles for an 
! arbitrary distribution function f(r,pr,L) via an acceptance-rejection method.
! Let fmax the maximum value of f. We generate arbitrary (x,y,z) numbers 
! in the range of (rmin,rmax), (pmin,pmax), (0,fmax) respectively. 
! Then evaluate W=f(x,y,L), if z<=W, accept the point, 
! otherwise, repeat until the condition in fulfilled.


! Initial data for the density function. Notice that
! we add two copies of the function in order to guarantee 
! that the ! boundary condition f(-r,-p) = f(r,p) is satisfied.

! Find the size of the cell

    drc = (rmaxc-rminc)/dble(Nrc)
    dpc = (pmaxc-pminc)/dble(Npc)
    dlc = (lmaxc-lminc)/dble(Nlc)

    print *, "(drc,dpc,dlc)=",drc,dpc,dlc
    if(state.eq."gaussian1") then

!     indx used to be incremented by hand each iteration, which is a
!     loop-carried dependency that blocks parallelizing the triple
!     loop below (hence the disabled "!!$OMP" that used to be here).
!     Since (k,i,j) -> indx = (k-1)*Nrc*Npc + (i-1)*Npc + j is exactly
!     the index this nesting order produces one at a time, computing
!     it directly removes that dependency: every thread writes only
!     to its own indx, so the loop is safe to collapse fully.

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

!     Same indx dependency removed as in the gaussian1 branch above
!     (the "aa" branch's own !!$OMP was disabled for the same reason).
!     This loop body is far more expensive per iteration (several
!     sqrt/trig evaluations), so it's the branch that benefits most
!     from actually running in parallel.

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

            !w3 = 8.0D0/(laux+2.0D0*J3+sqrt(laux**2+4.0D0))**3

            Q3 = Q3

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
            f(indx) = exp(-sin(0.5d0*Q3)**2/sp**2)*exp(-J3**2/sr**2)*J3**2*exp(-(laux-l0)**2/sl**2)

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


