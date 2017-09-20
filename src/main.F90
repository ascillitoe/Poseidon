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
  integer S1D_000_array(3) /0,0,0/                       ! nodal
  integer S3D_000_0M1_0P1_array(21) / 0,0,0,           & ! 2nd order 3D
                                      1,0,0, -1, 0, 0, &
                                      0,1,0,  0,-1, 0, &
                                      0,0,1,  0, 0,-1 /
  integer S3D_000_0P1_array(12) / 0,0,0,  & ! 2nd order i-i+1
                                  1,0,0,  &
                                  0,1,0,  &
                                  0,0,1  /
  integer S1D_iper1_array(6) /  0,0,0,   5, 0, 0 /   ! periodic at i=0
  integer S1D_iper2_array(6) /  0,0,0,  -5, 0, 0 /   ! periodic at i=I+1
  integer S1D_jper1_array(6) /  0,0,0,   0, 5, 0 /   ! periodic at j=0
  integer S1D_jper2_array(6) /  0,0,0,   0,-5, 0 /   ! periodic at j=I+1
  integer S1D_kper1_array(6) /  0,0,0,   0, 0, 5 /   ! periodic at k=0
  integer S1D_kper2_array(6) /  0,0,0,   0, 0,-5 /   ! periodic at k=I+1

  type(ops_stencil) :: S1D_000
  type(ops_stencil) :: S3D_000_0M1_0P1
  type(ops_stencil) :: S3D_000_0P1
  type(ops_stencil) :: S1D_iper1
  type(ops_stencil) :: S1D_iper2
  type(ops_stencil) :: S1D_jper1
  type(ops_stencil) :: S1D_jper2
  type(ops_stencil) :: S1D_kper1
  type(ops_stencil) :: S1D_kper2

! ops datasets
  type(ops_dat) :: x
  type(ops_dat) :: rho_s, rhou_s, Et_s
  type(ops_dat) :: rho, rhou, Et
  type(ops_dat) :: rho_res, rhou_res, Et_res
  type(ops_dat) :: mul, mut

! ops_reduction
  type(ops_reduction) :: dtmin

! iteration ranges
  integer iter_range(6)

! vars for halo_depths
  integer d_p(3) / 1, 1, 1/     !max halo depths for the dat in the possitive direction
  integer d_m(3) /-1,-1,-1/  !max halo depths for the dat in the negative direction

! Mesh size
  integer size(3) / 6,6,6 /

! base
  integer base(3) /1,1,1/ ! this is in fortran indexing

! null array 
  real(kind=c_double), dimension(:), allocatable :: temp

! profiling
  real(kind=c_double) :: startTime = 0
  real(kind=c_double) :: endTime = 0

! Local variables
  real(8) dt
  character*132 line

! Initialise constants
  xmin = -1.0_8
  ymin = -1.0_8
  zmin = -1.0_8
  xmax =  1.0_8
  ymax =  1.0_8
  zmax =  1.0_8

  dx = (xmax-xmin)/(size(1)-1.0_8)
  dy = (ymax-ymin)/(size(2)-1.0_8)
  dz = (zmax-zmin)/(size(3)-1.0_8)

  totaltime = 10.0_8
  cfl       = 1.0_8

! ---- Initialisation --------------------------------------------------------
! OPS initialisation
  call ops_init(1)
         
! ---- OPS Declarations-------------------------------------------------------
! declare block
  call ops_decl_block(3, grid, "grid")

! declare stencils
  call ops_decl_stencil( 3,  1, S1D_000_array        , S1D_000        , "000")
  call ops_decl_stencil( 3,  7, S3D_000_0M1_0P1_array, S3D_000_0M1_0P1, "000:100:-100:010:0-10:001:00-1")
  call ops_decl_stencil( 3,  4, S3D_000_0P1_array, S3D_000_0P1, "000:100:010:001")
  call ops_decl_stencil( 3,  2, S1D_iper1_array      , S1D_iper1      , "000:I00")
  call ops_decl_stencil( 3,  2, S1D_iper2_array      , S1D_iper2      , "000:-I00")
  call ops_decl_stencil( 3,  2, S1D_jper1_array      , S1D_jper1      , "000:0J0")
  call ops_decl_stencil( 3,  2, S1D_jper2_array      , S1D_jper2      , "000:0-J0")
  call ops_decl_stencil( 3,  2, S1D_kper1_array      , S1D_kper1      , "000:00K")
  call ops_decl_stencil( 3,  2, S1D_kper2_array      , S1D_kper2      , "000:00-K")

  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp,     x, "real(8)",     "x")

  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  rho_s, "real(8)",  "rho_s")
  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp, rhou_s, "real(8)", "rhou_s")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,   Et_s, "real(8)",   "Et_s")

  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  rho, "real(8)",  "rho")
  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp, rhou, "real(8)", "rhou")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,   Et, "real(8)",   "Et")

  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  rho_res, "real(8)", "rho_res")
  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp, rhou_res, "real(8)", "rhou_res")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,   Et_res, "real(8)", "Et_res")

  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  mul, "real(8)",  "mul")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  mut, "real(8)",  "mut")

 !reduction handle for dtmin
  call ops_decl_reduction_handle(8, dtmin, "real(8)", "dtmin")

 !---- Partition the mesh ----------------------------------------------------
  call ops_partition("3D_block_decompose")

! -----------------------------------------------------------------------------
! Initialise mesh and flowfield (for now...)
! -----------------------------------------------------------------------------
  iter_range(1) = 0
  iter_range(2) = size(1)+1
  iter_range(3) = 0
  iter_range(4) = size(2)+1
  iter_range(5) = 0
  iter_range(6) = size(3)+1
  call ops_par_loop(initialise_kernel, "initialise_kernel", grid, 3, iter_range, &
                     ops_arg_dat(   x, 3, S1D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat( rho, 1, S1D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat(rhou, 3, S1D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat(  Et, 1, S1D_000, "real(8)", OPS_WRITE), &
                     ops_arg_idx())


! -----------------------------------------------------------------------------
! Main computational loop
! -----------------------------------------------------------------------------
  ! start timer
  call ops_timers(startTime)

  iter_range(1) = 1
  iter_range(2) = size(1)
  iter_range(3) = 1
  iter_range(4) = size(2)
  iter_range(5) = 1
  iter_range(6) = size(3)

  nt = 1
  do while (simtime.le.totaltime)
    write(line,'(a,i8,2g18.9)') NEW_LINE('A') // 'nt, simtime, dt =', nt, simtime, dt
    call ops_printf(line)
   
!   Find timestep based on CFL condition
!   ------------------------------------
    call ops_par_loop(calc_dt_kernel, "calc_dt_kernel", grid, 3, iter_range, &
                         ops_arg_dat(    x, 3, S3D_000_0P1, "real(8)",OPS_READ), &
                         ops_arg_dat(  rho, 1, S3D_000_0P1, "real(8)",OPS_READ), &
                         ops_arg_dat( rhou, 3, S3D_000_0P1, "real(8)",OPS_READ), &
                         ops_arg_reduce(dtmin, 1,          "real(8)", OPS_MIN)  )
    call ops_reduction_result(dtmin, dt)

!   Save variables to start values before RK loop
!   ---------------------------------------------
    call ops_par_loop(save_kernel, "save_kernel", grid, 3, iter_range, &
            & ops_arg_dat( rho_s, 1, S1D_000, "real(8)", OPS_WRITE), &
            & ops_arg_dat(rhou_s, 3, S1D_000, "real(8)", OPS_WRITE), &
            & ops_arg_dat(  Et_s, 1, S1D_000, "real(8)", OPS_WRITE), &
            & ops_arg_dat(   rho, 1, S1D_000, "real(8)", OPS_READ ), &
            & ops_arg_dat(  rhou, 3, S1D_000, "real(8)", OPS_READ ), &
            & ops_arg_dat(    Et, 1, S1D_000, "real(8)", OPS_READ ))

!   Start RK loop
!   -------------
    do nrk = 1,3

!     Apply periodic boundary conditions
!     ----------------------------------
      call periodics(grid,size,S1D_iper1,S1D_iper2,S1D_jper1,  &
                     S1D_jper2,S1D_kper1,S1D_kper2,rho,rhou,Et )

!     Calculate residuals
!     -------------------
      ! continuity eqn: d(rho*u_i)/dx_i term
      call ops_par_loop(conres_kernel, "conres_kernel", grid, 3, iter_range, &
                        ops_arg_dat(      x, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(   rhou, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(rho_res, 1, S1D_000        , "real(8)", OPS_WRITE)  )

      ! momentum eqns: d(rho*u_i*u_j - P_ij)/dx_i terms (TODO - need to calculate P_ij stress tensor first)
      call ops_par_loop(momres_kernel, "momres_kernel", grid, 3, iter_range, &
                        ops_arg_dat(       x, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(     rho, 1, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(    rhou, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(rhou_res, 3, S1D_000        , "real(8)", OPS_WRITE)  )

      ! energy eqn: d(u_i*E - u_j*P_ij + q_i)/dx_i term
      call ops_par_loop(energyres_kernel, "energyres_kernel", grid, 3, iter_range, &
                        ops_arg_dat(     x, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(   rho, 1, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(  rhou, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(    Et, 1, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(Et_res, 1, S1D_000        , "real(8)", OPS_WRITE)  )

!     Apply filter (here or after update?)
!     ------------

!     RK update
!     ---------
    enddo

!   Update to next timestep (RK loop eventually)
!   ----------------------- 
    simtime = simtime + 0.1_8 !dt
    nt      = nt + 1
  enddo

! -----------------------------------------------------------------------------
! Output flow and grid to hdf
! -----------------------------------------------------------------------------
  call ops_print_dat_to_txtfile(       x,        "x.dat")
  call ops_print_dat_to_txtfile(     rho,      "rho.dat")
  call ops_print_dat_to_txtfile( rho_res,  "rho_res.dat")
  call ops_print_dat_to_txtfile(rhou_res, "rhou_res.dat")
  call ops_print_dat_to_txtfile(  Et_res,   "Et_res.dat")

  call ops_fetch_block_hdf5_file(  grid, "test.h5")
  call ops_fetch_dat_hdf5_file(     rho, "test.h5")
  call ops_fetch_dat_hdf5_file(    rhou, "test.h5")
  call ops_fetch_dat_hdf5_file(      Et, "test.h5")
  call ops_fetch_dat_hdf5_file( rho_res, "test.h5")
  call ops_fetch_dat_hdf5_file(rhou_res, "test.h5")
  call ops_fetch_dat_hdf5_file(  Et_res, "test.h5")

! -----------------------------------------------------------------------------
! Output runtime and exit
! -----------------------------------------------------------------------------
  call ops_timers(endTime)

  write (line,*) NEW_LINE('A') // 'Max total runtime =', endTime - startTime,'seconds'
  call ops_printf(line)

  call ops_exit( )

end program
