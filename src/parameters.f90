! ===========================================================================
! parameters.f90
! ===========================================================================
!> Module with the global parameters for the physical system.


module parameters

  implicit none

  !Grid

   integer   :: Nr = 1000                !< Number of points in the grid r
   integer   :: Nrc= 100                 !< Number of cells used in r in the support of f
   integer   :: Npc= 100                 !< Number of cells used in p in the support of f
   integer   :: Nlc= 20                  !< Number of cells used in l in the support of f

   real(8)   :: dr = 0.1D0               !< Grid spacing in the r direction
   real(8)   :: dp = 0.1D0               !< Grid spacing in the p direction
   real(8)   :: dl = 0.1D0               !< Grid spacing in the l direction

   real(8)   :: drc = 0.1D0              !< Size of the cell in r direction
   real(8)   :: dpc = 0.1D0              !< Size of the cell in p direction
   real(8)   :: dlc = 0.1D0              !< Size of the cell in l direction
   integer   :: Npart = 0                !< Number of computational particles

   real(8)   :: rmin = 0.0D0             !< Minimum radius of the grid
   real(8)   :: rmax = 20.0D0            !< Maximum radius of the grid
   real(8)   :: pmax = 2.0D0             !< Maximum momentum (to stimate the CFL condition)

   real(8)   :: rminc= 0.0D0             !< Minimum radius in the support of f.
   real(8)   :: rmaxc= 4.0D0             !< Maximum radius in the support of f.

   real(8)   :: pminc = -2.0D0           !< Minimum momentum in the support of f.
   real(8)   :: pmaxc =  2.0D0           !< Maximum momentum in the support of f.

   real(8)   :: lminc = -2.0D0           !< Minimum angular momentum in the support of f.
   real(8)   :: lmaxc =  2.0D0           !< Maximum angular momentum in the support of f.

   real(8)   :: Fmax = 0.0D0             !< Maximum absolute value of the force.
!   real(8)   :: Lfix = 1.0D0             !< Angular momentum
   real(8)   :: eps  = 0.0D0             !< Softening length for angular momentum.
   integer   :: ghost = 0                !< Number of ghost zones.
   !Time
   real(8)   :: courant = 0.5D0           !< Courant factor
   integer   :: Nt = 1000000              !< Total number of time steps
   real(8)   :: dt = 0.1D0                !< Time step
   real(8)   :: t = 0.0D0                 !< Time
   character(20) :: dt_switch = "fix"     !< Time step fixed from the initial force, or recomputed every step (fix,var)

   !Gravitational potentials
!    real(8)   :: rhoc = 1.0D0              !< Density parameter for gravitational potential
!    real(8)   :: rc = 1.0D0                !< Radius parameter for gravitational potential

    !Parameters for the Initial States


   real(8)   :: m0 = 1.0D0                !< Mass of the particles
   character(10) :: state = "gaussian1"   !< Initial distribution (gaussian1,aa)

   !Gaussian distribution

   real(8)   :: a0 = 1.0D0                !< Amplitude of the gaussian
   real(8)   :: r0  = 3.0D0               !< Center of the gaussian in r
   real(8)   :: sr  = 1.0D0               !< Width of the gaussian in r
   real(8)   :: p0  = -0.5D0              !< Center of the gaussian in p
   real(8)   :: sp  = 0.5D0               !< Width of the gaussian in p
   real(8)   :: l0  = 2.0D0               !< Center of the gaussian in l
   real(8)   :: sl  = 0.5D0               !< Width of the gaussian in l
   real(8)   :: cutoff = 0.0D0            !< Cutoff percentage

   !CheckPoint
   character(100) :: CheckPointfile = "input_file.2D"  !< Initial state file

   !Output
   character(100) :: directory = "test"     !< Output directory (may include subdirectories)
   character(20)  :: output_format = "ascii" !< Output format (ascii,hdf5)
   logical        :: reduceparticles = .false. !< Do we discard particles outside the domain? 
   integer        :: Nreduce = 100000       !< How often do we do discard particles outside domain?
   integer        :: spatial_output = 10000 !< Spatial output
   integer        :: time_output = 10000    !< Time output
   integer        :: field_output = 0       !< How often to save particles and grid fields (<=0: every spatial_output)
   integer        :: time_reduce_arr = 10000!< Reduce array size every
   character(20)  :: conv_test = "off"      !< Convergence test switch (on,off)

   !Methods
   integer        :: bsplineorder = 1     !< B-spline order
   character(20)  :: integrator = "euler"   !< Time integrator method (euler,icn,rk4)
   character(20)  :: spatialorder = "two"        !< Spatial order of discretization.
   character(20)  :: forcetype = "bg"     !< Background force (bg). Self-gravity: autointeraction; no background: BGtype=null
   character(20)  :: BGtype = "sphere"    !< Type of background. When the gravitational force is fix (sphere,iso,isotrun,nfw,burkert)
   logical        :: autointeraction = .false. !< Self interaction of particles
   !Energy variables

   real(8)        :: total_energy   = 0.D0
   real(8)        :: kinetic  = 0.D0
   real(8)        :: potential= 0.D0

   !Test functions of h_k (analysish.f90). A negative width or L center
   !means "take the value of the initial distribution": sq=sp, sj=sr,
   !lt=l0, slt=sl. Resolved when the parameters are read.

   real(8)        :: j1 = 0.0D0, j2 = 0.0D0       !< Center in J
   real(8)        :: sj1 = -1.0D0, sj2 = -1.0D0   !< Width in J
   real(8)        :: sq1 = -1.0D0, sq2 = -1.0D0   !< Width in Q
   real(8)        :: lt1 = -1.0D0, lt2 = -1.0D0   !< Center in L
   real(8)        :: slt1 = -1.0D0, slt2 = -1.0D0 !< Width in L

   !Averaging window

   real(8)        :: r1 = 2.0D0
   real(8)        :: r2 = 3.0D0

end module parameters
