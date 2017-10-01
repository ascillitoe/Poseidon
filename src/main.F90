! -----------------------------------------------------------------------------
! Poseidon - A high-order multi-block structured CFD code for executation on
!            clusters of GPU's, many-core CPU's and accelerators.
!
! Author(s) - Ashley Scillitoe
!
! This is a multi-block structured CFD code using 6th order finite
! differences and low storage RK to solve the compressible Navier-Stokes 
! equations. The code is made parallel using the open source OPS 
! framework (http://www.oerc.ox.ac.uk/projects/ops).
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

  intrinsic :: sqrt, real

! -----------------------------------------------------------------------------
! Read in hdf5 grid and flow files
! -----------------------------------------------------------------------------
  !call input TODO (eventually do ops block definition etc etc in here)

! ops blocks TODO (single block for now, eventually multiblock)
  type(ops_block) :: grid
  !type(ops_block) :: blocks(nbk) see poisson app for e.g. of multiblock
 
! Vars for stencils
  integer S3D_000_array(3) /0,0,0/                       ! nodal
  integer S3D_000_0M1_0P1_array(21) / 0,0,0,           & ! 2nd order 3D
                                      1,0,0, -1, 0, 0, &
                                      0,1,0,  0,-1, 0, &
                                      0,0,1,  0, 0,-1 /
  integer stride1D_x(3) /1,0,0/ ! 1D profile in i-dir

  type(ops_stencil) :: S3D_000
  type(ops_stencil) :: S3D_000_0M1_0P1
  type(ops_stencil) :: S1D_000_STRIDEx

! ops datasets
  type(ops_dat) :: x
  type(ops_dat) :: rho_s, rhou_s, Et_s
  type(ops_dat) :: rho, rhou, Et
  type(ops_dat) :: rho_res, rhou_res, Et_res
  type(ops_dat) :: tau
  type(ops_dat) :: mul, mut
  type(ops_dat) :: work
  type(ops_dat) :: xprof

! ops_reduction
  type(ops_reduction) :: dtmin

! iteration ranges
  integer iter_range(6)

! vars for halo_depths (for inter-block halos)
  integer d_p(3) / 1, 1, 1/     !max halo depths for the dat in the possitive direction
  integer d_m(3) /-1,-1,-1/  !max halo depths for the dat in the negative direction

! Mesh size
  integer size(3) / 200,7,7 /
  integer size2(3) 

! base
  integer base(3) /1,1,1/ ! this is in fortran indexing

! halos
  type(ops_halo), DIMENSION(18) :: halo_array
  type(ops_halo_group) :: halogrp_per
  integer halo_iter(3), base_from(3), base_to(3), dir_from(3), dir_to(3)

! null array 
  real(kind=c_double), dimension(:), allocatable :: temp

! profiling
  real(kind=c_double) :: startTime = 0
  real(kind=c_double) :: endTime = 0

! Local variables
  real(8) :: a1(3), a2(3)   ! RK constants
  real(8) :: dt
  integer nt, nrk

! Initialise constants
  xmin = 0.0_8
  ymin = 0.0_8
  zmin = 0.0_8
  xmax = 1.0_8
  ymax = 1.0_8
  zmax = 1.0_8

  dx = (xmax-xmin)/(size(1)-1.0_8)
  dy = dx!(ymax-ymin)/(size(2)-1.0_8)
  dz = dx!(zmax-zmin)/(size(3)-1.0_8)

  totaltime = 100.0_8
  ntmax     = 1000
  phi       = 0.05
  cfl       = 0.1_8
  a1 = (/ 2.0_8/3.0_8, 5.0_8/12.0_8, 3.0_8/5.0_8 /)
  a2 = (/ 1.0_8/4.0_8, 3.0_8/20.0_8, 3.0_8/5.0_8 /)

! ---- Initialisation --------------------------------------------------------
! OPS initialisation
  call ops_init(1)

! ---- OPS Declarations-------------------------------------------------------
! declare block
  if (ops_is_root() .eq. 1) print*, 'Declaring block'
  call ops_decl_block(3, grid, "grid")

! declare stencils
  if (ops_is_root() .eq. 1) print*, 'Declaring stencils'
  call ops_decl_stencil( 3,  1, S3D_000_array        , S3D_000        , "000")
  call ops_decl_stencil( 3,  7, S3D_000_0M1_0P1_array, S3D_000_0M1_0P1, "000:100:-100:010:0-10:001:00-1")
  call ops_decl_strided_stencil(3,1,S3D_000_array, stride1D_x, S1D_000_STRIDEx,  "000_STRIDEx")

  if (ops_is_root() .eq. 1) print*, 'Declaring datasets'
  ! coordinates
  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp,     x, "real(8)",     "x")

  ! flow var at start of RK loop
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  rho_s, "real(8)",  "rho_s")
  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp, rhou_s, "real(8)", "rhou_s")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,   Et_s, "real(8)",   "Et_s")

  ! flow var
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  rho, "real(8)",  "rho")
  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp, rhou, "real(8)", "rhou")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,   Et, "real(8)",   "Et")

  ! residuals
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  rho_res, "real(8)", "rho_res")
  call ops_decl_dat(grid, 3, size, base, d_m, d_p, temp, rhou_res, "real(8)", "rhou_res")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,   Et_res, "real(8)", "Et_res")

  ! Viscous stress tensor
  call ops_decl_dat(grid, 6, size, base, d_m, d_p, temp,   tau, "real(8)", "tau")

  ! Viscosity
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  mul, "real(8)",  "mul")
  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  mut, "real(8)",  "mut")

  call ops_decl_dat(grid, 1, size, base, d_m, d_p, temp,  work, "real(8)",  "work")

  ! 1D profile
  size2 = (/ size(1), 1, 1 /)
  d_m   = (/ 0, 0, 0 /)
  d_p   = (/ 0, 0, 0 /)
  call ops_decl_dat(grid, 1, size2, base, d_m, d_p, temp,  xprof, "real(8)",  "xprof")


! declare halos for periodics
  if (ops_is_root() .eq. 1) print*, 'Declaring halos'

! i-dir periodics: i=I-2 to i=1  and  i=3 to i=I
  halo_iter(1:3)  = (/      1,size(2)-2,size(3)-2 /)
  base_from(1:3)  = (/      3,        2,        2 /)
  base_to(1:3)    = (/size(1),        2,        2 /)
  dir_from(1:3)   = (/      1,        2,        3 /)
  dir_to(1:3)     = (/      1,        2,        3 /)
  call ops_decl_halo( rho,  rho, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(1))
  call ops_decl_halo(rhou, rhou, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(2))
  call ops_decl_halo(  Et,   Et, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(3))
  base_from(1) = size(1)-2
  base_to(1)   = 1
  call ops_decl_halo( rho,  rho, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(4))
  call ops_decl_halo(rhou, rhou, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(5))
  call ops_decl_halo(  Et,   Et, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(6))

! j-dir periodics: j=J-2 to j=1  and  j=3 to j=J
  halo_iter(1:3)  = (/ size(1)-2,      1,size(3)-2 /)
  base_from(1:3)  = (/         2,      3,        2 /)
  base_to(1:3)    = (/         2,size(2),        2 /)
  dir_from(1:3)   = (/         1,      2,        3 /)
  dir_to(1:3)     = (/         1,      2,        3 /)
  call ops_decl_halo( rho,  rho, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(7))
  call ops_decl_halo(rhou, rhou, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(8))
  call ops_decl_halo(  Et,   Et, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(9))
  base_from(2) = size(2)-2
  base_to(2)   = 1
  call ops_decl_halo( rho,  rho, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(10))
  call ops_decl_halo(rhou, rhou, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(11))
  call ops_decl_halo(  Et,   Et, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(12))

! k-dir periodics: k=K-2 to k=1  and  k=3 to k=K
  halo_iter(1:3)  = (/size(1)-2,size(2)-2,      1 /)
  base_from(1:3)  = (/        2,        2,      3 /)
  base_to(1:3)    = (/        2,        2,size(3) /)
  dir_from(1:3)   = (/        1,        2,      3 /)
  dir_to(1:3)     = (/        1,        2,      3 /)
  call ops_decl_halo( rho,  rho, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(13))
  call ops_decl_halo(rhou, rhou, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(14))
  call ops_decl_halo(  Et,   Et, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(15))
  base_from(3) = size(3)-2
  base_to(3)   = 1
  call ops_decl_halo( rho,  rho, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(16))
  call ops_decl_halo(rhou, rhou, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(17))
  call ops_decl_halo(  Et,   Et, halo_iter, base_from, base_to, dir_from, dir_to, halo_array(18))

  call ops_decl_halo_group(18,halo_array, halogrp_per)

 !reduction handle for dtmin
  if (ops_is_root() .eq. 1) print*, 'Declaring reduction handles'
  call ops_decl_reduction_handle(8, dtmin, "real(8)", "dtmin")

 !---- Partition the mesh ----------------------------------------------------
  if (ops_is_root() .eq. 1) print*, 'Partitioning the mesh'

  call ops_partition("3D_block_decompose")
  call ops_diagnostic_output()
  
! -----------------------------------------------------------------------------
! Initialise mesh and flowfield (for now...)
! -----------------------------------------------------------------------------
  iter_range(1:6) = (/ 1,size(1), 1,size(2), 1,size(3) /)
  call ops_par_loop(initialise_kernel, "initialise_kernel", grid, 3, iter_range, &
                     ops_arg_dat(   x, 3, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat( rho, 1, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat(rhou, 3, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_dat(  Et, 1, S3D_000, "real(8)", OPS_WRITE), &
                     ops_arg_idx())
  call ops_halo_transfer(halogrp_per)

! -----------------------------------------------------------------------------
! Main computational loop
! -----------------------------------------------------------------------------
  ! start timer
  call ops_timers(startTime)

  nt = 1
  iter_range(1:6) = (/ 2,size(1)-1, 2,size(2)-1, 2,size(3)-1 /)
  do while (simtime.le.totaltime .and. nt.le.ntmax)
   
!   Find timestep based on CFL condition
!   ------------------------------------
    call ops_par_loop(calc_dt_kernel, "calc_dt_kernel", grid, 3, iter_range, &
                         ops_arg_dat(    x, 3, S3D_000_0M1_0P1, "real(8)",OPS_READ), &
                         ops_arg_dat(  rho, 1, S3D_000_0M1_0P1, "real(8)",OPS_READ), &
                         ops_arg_dat( rhou, 3, S3D_000_0M1_0P1, "real(8)",OPS_READ), &
                         ops_arg_dat(   Et, 1, S3D_000_0M1_0P1, "real(8)",OPS_READ), &
                         ops_arg_reduce(dtmin, 1,          "real(8)", OPS_MIN) )

    call ops_reduction_result(dtmin, dt)

!   Save variables to start values before RK loop
!   ---------------------------------------------
    call ops_par_loop(save_kernel, "save_kernel", grid, 3, iter_range, &
            & ops_arg_dat( rho_s, 1, S3D_000, "real(8)", OPS_WRITE), &
            & ops_arg_dat(rhou_s, 3, S3D_000, "real(8)", OPS_WRITE), &
            & ops_arg_dat(  Et_s, 1, S3D_000, "real(8)", OPS_WRITE), &
            & ops_arg_dat(   rho, 1, S3D_000, "real(8)", OPS_READ ), &
            & ops_arg_dat(  rhou, 3, S3D_000, "real(8)", OPS_READ ), &
            & ops_arg_dat(    Et, 1, S3D_000, "real(8)", OPS_READ ))


    if (ops_is_root() .eq. 1 .and. mod(nt,1).eq.0) &
      write(*,'(a,i8,3g18.9)')  'nt, simtime, dt =', nt, simtime, dt

!   Start RK loop
!   -------------
    do nrk = 1,3
      
!     Apply periodic boundary conditions
!     ----------------------------------
      call ops_halo_transfer(halogrp_per)

!     Calculate laminar viscosity
!     ---------------------------
      call ops_par_loop(viscosity_kernel, "viscosity_kernel", grid, 3, iter_range, &
                        ops_arg_dat(    rho, 1, S3D_000, "real(8)",  OPS_READ), &
                        ops_arg_dat(     Et, 1, S3D_000, "real(8)",  OPS_READ), &
                        ops_arg_dat(    mul, 1, S3D_000, "real(8)", OPS_WRITE)  )

!     Calculate viscous stress tensor
!     -------------------------------
      call ops_par_loop(stress_kernel, "stress_kernel", grid, 3, iter_range, &
                        ops_arg_dat(      x, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(    rho, 1, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(   rhou, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(    mul, 1, S3D_000        , "real(8)",  OPS_READ), &
                        ops_arg_dat(    tau, 6, S3D_000        , "real(8)", OPS_WRITE)  )

!     Calculate residuals
!     -------------------
      ! continuity eqn: d(rho*u_i)/dx_i term
      call ops_par_loop(conres_kernel, "conres_kernel", grid, 3, iter_range, &
                        ops_arg_dat(      x, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(   rhou, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(rho_res, 1, S3D_000        , "real(8)", OPS_WRITE)  )

      ! momentum eqns: d(rho*u_i*u_j - P_ij)/dx_i terms (TODO - need to calculate tau_ij stress tensor first)
      call ops_par_loop(momres_kernel, "momres_kernel", grid, 3, iter_range, &
                        ops_arg_dat(       x, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(     rho, 1, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(    rhou, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(      Et, 1, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(     tau, 6, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(rhou_res, 3, S3D_000        , "real(8)", OPS_WRITE)  )

      ! energy eqn: d(u_i*E - u_j*P_ij + q_i)/dx_i term
      call ops_par_loop(energyres_kernel, "energyres_kernel", grid, 3, iter_range, &
                        ops_arg_dat(     x, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(   rho, 1, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(  rhou, 3, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(    Et, 1, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(   tau, 6, S3D_000_0M1_0P1, "real(8)",  OPS_READ), &
                        ops_arg_dat(Et_res, 1, S3D_000        , "real(8)", OPS_WRITE)  )

!     RK update
!     ---------
      call ops_par_loop(updateRK_kernel, "updateRK_kernel", grid, 3, iter_range, &
                        ops_arg_dat(     rho, 1, S3D_000, "real(8)", OPS_RW),    &
                        ops_arg_dat(    rhou, 3, S3D_000, "real(8)", OPS_RW),    &
                        ops_arg_dat(      Et, 1, S3D_000, "real(8)", OPS_RW),    &
                        ops_arg_dat(   rho_s, 1, S3D_000, "real(8)", OPS_RW),    &
                        ops_arg_dat(  rhou_s, 3, S3D_000, "real(8)", OPS_RW),    &
                        ops_arg_dat(    Et_s, 1, S3D_000, "real(8)", OPS_RW),    &
                        ops_arg_dat( rho_res, 1, S3D_000, "real(8)", OPS_READ),  &
                        ops_arg_dat(rhou_res, 3, S3D_000, "real(8)", OPS_READ),  &
                        ops_arg_dat(  Et_res, 1, S3D_000, "real(8)", OPS_READ),  &
                        ops_arg_gbl( a1(nrk), 1,          "real(8)", OPS_READ),  &
                        ops_arg_gbl( a2(nrk), 1,          "real(8)", OPS_READ),  &
                        ops_arg_gbl(      dt, 1,          "real(8)", OPS_READ)   )

    enddo                    

!   Apply filter - simple 7 point average for now...
!   ------------
    call ops_par_loop(filter_kernel, "filter_kernel", grid, 3, iter_range,                &
                      ops_arg_dat(       rho, 1  , S3D_000_0M1_0P1, "real(8)",  OPS_RW),  &
                      ops_arg_dat(      rhou, 3  , S3D_000_0M1_0P1, "real(8)",  OPS_RW),  &
                      ops_arg_dat(        Et, 1  , S3D_000_0M1_0P1, "real(8)",  OPS_RW)   )


!   Update to next timestep
!   ----------------------- 
    call ops_par_loop(update_kernel, "update_kernel", grid, 3, iter_range,    &
                      ops_arg_dat(     rho, 1, S3D_000, "real(8)", OPS_RW),   &
                      ops_arg_dat(    rhou, 3, S3D_000, "real(8)", OPS_RW),   &
                      ops_arg_dat(      Et, 1, S3D_000, "real(8)", OPS_RW),   &
                      ops_arg_dat( rho_res, 1, S3D_000, "real(8)", OPS_READ), &
                      ops_arg_dat(rhou_res, 3, S3D_000, "real(8)", OPS_READ), &
                      ops_arg_dat(  Et_res, 1, S3D_000, "real(8)", OPS_READ), &
                      ops_arg_gbl(      dt, 1,          "real(8)", OPS_READ)  )

    simtime = simtime + dt
    nt      = nt + 1
  enddo


  iter_range(1:6) = (/ 1,size(1), 3, 3, 3, 3 /)
  call ops_par_loop(profile_kernel, "profile_kernel", grid, 3, iter_range,        &
                    ops_arg_dat(         x, 3  , S3D_000, "real(8)",  OPS_READ),  &
                    ops_arg_dat(       rho, 1  , S3D_000, "real(8)",  OPS_READ),  &
                    ops_arg_dat(      rhou, 3  , S3D_000, "real(8)",  OPS_READ),  &
                    ops_arg_dat(        Et, 1  , S3D_000, "real(8)",  OPS_READ),  &
                    ops_arg_dat(     xprof, 1  , S1D_000_STRIDEx, "real(8)",  OPS_INC), &
                    ops_arg_idx())
  call ops_print_dat_to_txtfile(       xprof,        "xprof.dat")

!! -----------------------------------------------------------------------------
!! Output flow and grid to hdf
!! -----------------------------------------------------------------------------
  iter_range(1:6) = (/ 1,size(1), 1,size(2), 1,size(3) /)
  call ops_par_loop(copy_kernel, "copy_kernel", grid, 3, iter_range,         &
                    ops_arg_dat(       x, 3, S3D_000, "real(8)", OPS_READ),  &
                    ops_arg_dat(    work, 1, S3D_000, "real(8)", OPS_WRITE), &
                    ops_arg_gbl(       1, 1,       "integer(4)", OPS_READ)  )
  call ops_fetch_block_hdf5_file(   grid, "x.h5")
  call ops_fetch_dat_hdf5_file(     work, "x.h5")

  call ops_par_loop(copy_kernel, "copy_kernel", grid, 3, iter_range,         &
                    ops_arg_dat(       x, 3, S3D_000, "real(8)", OPS_READ),  &
                    ops_arg_dat(    work, 1, S3D_000, "real(8)", OPS_WRITE), &
                    ops_arg_gbl(       2, 1,       "integer(4)", OPS_READ)  )
  call ops_fetch_block_hdf5_file(   grid, "y.h5")
  call ops_fetch_dat_hdf5_file(     work, "y.h5")

  call ops_par_loop(copy_kernel, "copy_kernel", grid, 3, iter_range,         &
                    ops_arg_dat(       x, 3, S3D_000, "real(8)", OPS_READ),  &
                    ops_arg_dat(    work, 1, S3D_000, "real(8)", OPS_WRITE), &
                    ops_arg_gbl(       3, 1,       "integer(4)", OPS_READ)  )
  call ops_fetch_block_hdf5_file(   grid, "z.h5")
  call ops_fetch_dat_hdf5_file(     work, "z.h5")

  call ops_par_loop(copy_kernel, "copy_kernel", grid, 3, iter_range,         &
                    ops_arg_dat(    rhou, 3, S3D_000, "real(8)", OPS_READ),  &
                    ops_arg_dat(    work, 1, S3D_000, "real(8)", OPS_WRITE), &
                    ops_arg_gbl(       1, 1,       "integer(4)", OPS_READ)  )
  call ops_fetch_block_hdf5_file(   grid, "rhou.h5")
  call ops_fetch_dat_hdf5_file(     work, "rhou.h5")

! -----------------------------------------------------------------------------
! Output runtime and exit
! -----------------------------------------------------------------------------
  call ops_timers(endTime)

  if (ops_is_root() .eq. 1) write (*,*)  'Max total runtime =', endTime - startTime,'seconds'


  call ops_exit( )

end program
