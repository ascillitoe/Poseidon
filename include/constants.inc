! -----------------------------------------------------------------------------
! constants.inc - COMMON blocks of constants used by all routines
!
! History:
! 16/09/17 - First build.                                                - ashS
!
! -----------------------------------------------------------------------------
!
! ---- Version info -----------------------------------------------------------

  common / ver_info / version, verdate
  character           version*16, verdate*11

! ---- Files ------------------------------------------------------------------

  common / files / inputfile, gridfile, flowfile
  character*256 inputfile, gridfile, flowfile

! inputfile -- Formatted input file 
! gridfile  -- HDF grid file
! flowfile  -- HDF flow file

! ---- Logicals ---------------------------------------------------------------


! ---- Runge-Kutta ------------------------------------------------------------

  common / RK / alfas(5), betas(5), cfl, nrk, totaltime

  real(8) :: alfas, betas, cfl, totaltime
  integer nrk

! alfas     -- RK coeffs for inviscid fluxes
! betas     -- RK coeffs for viscous fluxes
! cfl       -- CFL number
! nrk       -- Number of RK steps
! totaltime -- Total simulation time



