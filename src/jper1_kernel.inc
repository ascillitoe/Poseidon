subroutine jper1_kernel(rho, rhou, Et)
! copies j = J-1 variables to j = 0 ghost node for j-dir periodics

  implicit none
  real (kind=8), DIMENSION(3) :: rhou
  real (kind=8), DIMENSION(1) :: rho, Et
  integer f

  rho(OPS_ACC1(0,0,0)) = rho(OPS_ACC1(0,5,0))
  do f = 1,3
    rhou(OPS_ACC_MD2(f,0,0,0)) = rhou(OPS_ACC_MD2(f,0,5,0))
  enddo
  Et(OPS_ACC3(0,0,0)) = Et(OPS_ACC3(0,5,0))

end subroutine
