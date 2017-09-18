! -----------------------------------------------------------------------------
! CFD4turbo  - Compact Finite Differences for turbomachinery
!
! Author(s) - Ashley Scillitoe
!
! This is a multi-block structured CFD code using 6th order finite
! differences and low storage RK to solve the compressible Navier-Stokes 
! equations. 
!
! History:
! 16/09/17  - First build.                                             - as2341
!
! -----------------------------------------------------------------------------

program main

  use OPS_Fortran_Reference
  use OPS_FORTRAN_HDF5_DECLARATIONS
  use OPS_CONSTANTS

  use, intrinsic :: ISO_C_BINDING

  implicit none


! -----------------------------------------------------------------------------
! Read in hdf5 grid and flow files
! -----------------------------------------------------------------------------
  !call input TODO (eventually do ops block definition etc etc in here)

! ops blocks TODO (single block for now, eventually multiblock)
  type(ops_block) :: grid
  !type(ops_block) :: blocks(nbk) see poisson app for e.g. of multiblock
 
! Vars for stencils
  integer S3D_000_array(3) /0,0,0/
  type(ops_stencil) :: S3D_000
  integer S3D_000_0M1_0P1_array(21) / 0,0,0,         &
                                      1,0,0, -1,0,0, &
                                      0,1,0, 0,-1,0, &
                                      0,0,1, 0,0,-1  /
  type(ops_stencil) :: S3D_000_0M1_0P1

! ops datasets
  type(ops_dat) :: x
  type(ops_dat) :: rho_s, rhou_s, rhov_s, rhow_s, Et_s
  type(ops_dat) :: rho, rhou, rhov, rhow, Et
  type(ops_dat) :: mul, mut

! ops_reduction
  type(ops_reduction) :: dtmin

! iteration ranges
  integer iter_range(6)

! vars for halo_depths
  integer d_p(3) /2,2,2/   !max halo depths for the dat in the possitive direction
  integer d_m(3) /-2,-2,-2/  !max halo depths for the dat in the negative direction

! base
  integer base(3) /1,1,1/ ! this is in fortran indexing

! size of block (and data)
  integer size(3) /128,128,3/

! null array 
  real(kind=c_double), dimension(:), allocatable :: temp

! profiling
  real(kind=c_double) :: startTime = 0
  real(kind=c_double) :: endTime = 0

! Initialise constants
  xhalo = 2
  yhalo = 2
  zhalo = 2
  xmin = -5.0_8
  ymin =  0.0_8
  zmin =  0.0_8
  xmax =  5.0_8
  ymax =  0.5_8
  zmax =  0.1_8

  dx = (xmax-xmin)/(size(1)-(1.0_8 + 2.0_8*xhalo))
  dy = (ymax-ymin)/(size(2)-(1.0_8 + 2.0_8*yhalo))
  dz = (zmax-zmin)/(size(3)-(1.0_8 + 2.0_8*zhalo))

! ---- Initialisation --------------------------------------------------------
! OPS initialisation
  call ops_init(1)
         
! ---- OPS Declarations-------------------------------------------------------
! declare block
  call ops_decl_block(3, grid, "grid")

! declare stencils
  call ops_decl_stencil( 3,  1, S3D_000_array, S3D_000, "000")
  call ops_decl_stencil( 3,  7, S3D_000_0M1_0P1_array, S3D_000_0M1_0P1, "000:100:-100:010:0-10:001:00-1")

  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp,     x, "real(8)",     "x")

  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  rho_s, "real(8)",  "rho_s")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp, rhou_s, "real(8)", "rhou_s")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp, rhov_s, "real(8)", "rhov_s")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp, rhow_s, "real(8)", "rhow_s")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,   Et_s, "real(8)",   "Et_s")

  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  rho, "real(8)",  "rho")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp, rhou, "real(8)", "rhou")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp, rhov, "real(8)", "rhov")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp, rhow, "real(8)", "rhow")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,   Et, "real(8)",   "Et")

  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  mul, "real(8)",  "mul")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  mut, "real(8)",  "mut")

 !reduction handle for rms variable 
  call ops_decl_reduction_handle(8, dtmin, "real(8)", "dtmin")

 !---- Partition the mesh ----------------------------------------------------
  call ops_partition("3D_block_decompose")

! -----------------------------------------------------------------------------
! Initialise mesh and flowfield (for now...)
! -----------------------------------------------------------------------------
  iter_range(1) = 1
  iter_range(2) = size(1)
  iter_range(3) = 1
  iter_range(4) = size(2)
  iter_range(5) = 1
  iter_range(6) = size(3)
  call ops_par_loop(initialise_kernel, "initialise_kernel", grid, 3, iter_range, &
                     ops_arg_dat(   x, 3, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat( rho, 1, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat(rhou, 1, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat(rhov, 1, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat(rhow, 1, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat(  Et, 1, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_idx())

! -----------------------------------------------------------------------------
! Main computational loop
! -----------------------------------------------------------------------------
  ! start timer
  call ops_timers(startTime)

! -----------------------------------------------------------------------------
! Output flow and grid to hdf
! -----------------------------------------------------------------------------
  call ops_fetch_block_hdf5_file(grid, "test.h5")
  call ops_fetch_dat_hdf5_file(   rho, "test.h5")
  call ops_fetch_dat_hdf5_file(  rhou, "test.h5")
  call ops_fetch_dat_hdf5_file(  rhov, "test.h5")
  call ops_fetch_dat_hdf5_file(  rhow, "test.h5")
  call ops_fetch_dat_hdf5_file(    Et, "test.h5")


! -----------------------------------------------------------------------------
! Output runtime and exit
! -----------------------------------------------------------------------------
  call ops_timers(endTime)

  if (ops_is_root() .eq. 1) then
    write (*,*) 'Max total runtime =', endTime - startTime,'seconds'
  endif

  call ops_exit( )

end program
