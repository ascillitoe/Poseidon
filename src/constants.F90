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

!   ---- Temporary -----------------------------------------------------------
    real(8),constant :: dx, dy, dz
    real(8),constant :: xmin,xmax,ymin,ymax,zmin,zmax
    integer,constant :: xhalo, yhalo, zhalo

!   ---- Runge-Kutta ---------------------------------------------------------
    real(8),constant :: alfas, betas, cfl, totaltime
    real(8) simtime

    integer,constant :: nrk
    integer nt

#else

!   ---- Temporary -----------------------------------------------------------
    real(8) :: dx, dy, dz
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
    integer :: xhalo, yhalo, zhalo

!   ---- Runge-Kutta ---------------------------------------------------------
    real(8) :: alfas, betas, cfl, totaltime, simtime
    integer :: nrk, nt

#endif


! Description of constants
! ---- Temporary -------------------------------------------------------------
! dx        -- Grid spacing in x-direction
! dy        -- Grid spacing in y-direction
! dz        -- Grid spacing in z-direction

! ---- Runge-Kutta -----------------------------------------------------------
! alfas     -- RK coeffs for inviscid fluxes
! betas     -- RK coeffs for viscous fluxes
! cfl       -- CFL number
! nrk       -- Number of RK steps
! totaltime -- Total simulation time

END MODULE OPS_CONSTANTS

