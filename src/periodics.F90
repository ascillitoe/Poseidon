subroutine periodics(grid,size,S1D_iper1,S1D_iper2,S1D_jper1,       &
                     S1D_jper2,S1D_kper1,S1D_kper2,rho,rhou,Et      )

  use OPS_Fortran_Reference
  use OPS_CONSTANTS

  use, intrinsic :: ISO_C_BINDING

  implicit none

  integer iter_range(6)
  integer size(3)

  type(ops_block) :: grid
  type(ops_stencil) :: S1D_iper1
  type(ops_stencil) :: S1D_iper2
  type(ops_stencil) :: S1D_jper1
  type(ops_stencil) :: S1D_jper2
  type(ops_stencil) :: S1D_kper1
  type(ops_stencil) :: S1D_kper2

  type(ops_dat) :: x
  type(ops_dat) :: rho, rhou, Et

  iter_range(1) = 0
  iter_range(2) = 0
  iter_range(3) = 1
  iter_range(4) = size(2)
  iter_range(5) = 1
  iter_range(6) = size(3)
! copies i = I-1 variables to i = 0 ghost node for i-dir periodics
  call ops_par_loop(iper1_kernel, "iper1_kernel", grid, 3, iter_range,    &
                    ops_arg_dat( rho, 1, S1D_iper1, "real(8)",  OPS_RW), &
                    ops_arg_dat(rhou, 3, S1D_iper1, "real(8)",  OPS_RW), &
                    ops_arg_dat(  Et, 1, S1D_iper1, "real(8)",  OPS_RW)  )

  iter_range(1) = size(1)+1
  iter_range(2) = size(1)+1
  iter_range(3) = 1
  iter_range(4) = size(2)
  iter_range(5) = 1
  iter_range(6) = size(3)
! copies i = 1 variables to i = I+1 ghost node for i-dir periodics
  call ops_par_loop(iper2_kernel, "iper2_kernel", grid, 3, iter_range,    &
                    ops_arg_dat( rho, 1, S1D_iper2, "real(8)",  OPS_RW), &
                    ops_arg_dat(rhou, 3, S1D_iper2, "real(8)",  OPS_RW), &
                    ops_arg_dat(  Et, 1, S1D_iper2, "real(8)",  OPS_RW)  )

  iter_range(1) = 1
  iter_range(2) = size(1)
  iter_range(3) = 0
  iter_range(4) = 0
  iter_range(5) = 1
  iter_range(6) = size(3)
! copies j = J-1 variables to j = 0 ghost node for j-dir periodics
  call ops_par_loop(jper1_kernel, "jper1_kernel", grid, 3, iter_range,    &
                    ops_arg_dat( rho, 1, S1D_iper1, "real(8)",  OPS_RW), &
                    ops_arg_dat(rhou, 3, S1D_iper1, "real(8)",  OPS_RW), &
                    ops_arg_dat(  Et, 1, S1D_iper1, "real(8)",  OPS_RW)  )

  iter_range(1) = 1
  iter_range(2) = size(1)
  iter_range(3) = size(2)+1
  iter_range(4) = size(2)+1
  iter_range(5) = 1
  iter_range(6) = size(3)
! copies j = 1 variables to j = J+1 ghost node for j-dir periodics
  call ops_par_loop(jper2_kernel, "jper2_kernel", grid, 3, iter_range,    &
                    ops_arg_dat( rho, 1, S1D_iper2, "real(8)",  OPS_RW), &
                    ops_arg_dat(rhou, 3, S1D_iper2, "real(8)",  OPS_RW), &
                    ops_arg_dat(  Et, 1, S1D_iper2, "real(8)",  OPS_RW)  )

  iter_range(1) = 1
  iter_range(2) = size(1)
  iter_range(3) = 1
  iter_range(4) = size(2)
  iter_range(5) = 0
  iter_range(6) = 0
! copies k = K-1 variables to k = 0 ghost node for k-dir periodics
  call ops_par_loop(kper1_kernel, "kper1_kernel", grid, 3, iter_range,    &
                    ops_arg_dat( rho, 1, S1D_iper1, "real(8)",  OPS_RW), &
                    ops_arg_dat(rhou, 3, S1D_iper1, "real(8)",  OPS_RW), &
                    ops_arg_dat(  Et, 1, S1D_iper1, "real(8)",  OPS_RW)  )

  iter_range(1) = 1
  iter_range(2) = size(1)
  iter_range(3) = 1
  iter_range(4) = size(2)
  iter_range(5) = size(3)+1
  iter_range(6) = size(3)+1
! copies k = 1 variables to k = K+1 ghost node for k-dir periodics
  call ops_par_loop(kper2_kernel, "kper2_kernel", grid, 3, iter_range,    &
                    ops_arg_dat( rho, 1, S1D_iper2, "real(8)",  OPS_RW), &
                    ops_arg_dat(rhou, 3, S1D_iper2, "real(8)",  OPS_RW), &
                    ops_arg_dat(  Et, 1, S1D_iper2, "real(8)",  OPS_RW)  )

end subroutine
