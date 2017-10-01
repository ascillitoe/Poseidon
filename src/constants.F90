MODULE OPS_CONSTANTS
! -----------------------------------------------------------------------------
! constants.inc - COMMON blocks of constants used by all routines
!
! History:
! 16/09/17 - First build.                                                - ashS
!
! -----------------------------------------------------------------------------

#ifdef OPS_WITH_CUDAFOR
    use cudafor
#endif

!   ---- Runge-Kutta ---------------------------------------------------------
#ifdef OPS_WITH_CUDAFOR
    real(8),constant  :: cfl, totaltime
    real(8),constant  :: alfas(4) = (/ 970286171893.0_8   / 4311952581923.0_8,  &
                                       6584761158862.0_8  / 12103376702013.0_8, &
                                       2251764453980.0_8  / 15575788980749.0_8, &
                                       26877169314380.0_8 / 34165994151039.0_8 /) 
    real(8),constant  :: betas(5) = (/ 1153189308089.0_8  / 22510343858157.0_8, &
                                       1772645290293.0_8  / 4653164025191.0_8,  &
                                       1672844663538.0_8  / 4480602732383.0_8,  &
                                       2114624349019.0_8  / 3568978502595.0_8,  &
                                       5198255086312.0_8  / 14908931495163.0_8 /)

    integer, constant :: ntmax
#else
    real(8) :: cfl, totaltime
    real(8) :: alfas(4) = (/ 970286171893.0_8   / 4311952581923.0_8,  &
                             6584761158862.0_8  / 12103376702013.0_8, &
                             2251764453980.0_8  / 15575788980749.0_8, &
                             26877169314380.0_8 / 34165994151039.0_8 /) 
    real(8) :: betas(5) = (/ 1153189308089.0_8  / 22510343858157.0_8, &
                             1772645290293.0_8  / 4653164025191.0_8,  &
                             1672844663538.0_8  / 4480602732383.0_8,  &
                             2114624349019.0_8  / 3568978502595.0_8,  &
                             5198255086312.0_8  / 14908931495163.0_8 /)
    integer ntmax
#endif

    real(8) simtime

!   ---- Real-Gas ------------------------------------------------------------
#ifdef OPS_WITH_CUDAFOR
    real(8),constant :: gamma = 1.4_8, R=287.058_8
#else
    real(8) :: gamma = 1.4_8, R=287.058_8
#endif


!   ---- Temporary -----------------------------------------------------------
    real(8) :: dx, dy, dz
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax, phi
    integer :: xhalo, yhalo, zhalo


! Description of constants
! ---- Runge-Kutta -----------------------------------------------------------
! alfas     -- alpha coefficients for RK4(3)5[2R+]C scheme from Kennedy (2000) 
! betas     -- beta coefficients for RK4(3)5[2R+]C scheme from Kennedy (2000) 
! cfl       -- CFL number
! totaltime -- Total simulation time
! ntmax     -- Max no. of timesteps
! simtime   -- Simulation time

! ---- Real-gas --------------------------------------------------------------
! gamma     -- Ratio of specific heats
! R         -- Specific gas constant

! ---- Temporary -------------------------------------------------------------
! dx        -- Grid spacing in x-direction
! dy        -- Grid spacing in y-direction
! dz        -- Grid spacing in z-direction

END MODULE OPS_CONSTANTS

