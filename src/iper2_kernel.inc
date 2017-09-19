subroutine iper2_kernel(rho, rhou, Et)
! copies i = 1 variables to i = I+1 ghost node for i-dir periodics

  implicit none
  real (kind=8), DIMENSION(3) :: rhou
  real (kind=8), DIMENSION(1) :: rho, Et
  integer f

  rho(OPS_ACC1( 0,0,0)) = rho(OPS_ACC1(-5,0,0))
  do f = 1,3
    rhou(OPS_ACC_MD2(f,0,0,0)) = rhou(OPS_ACC_MD2(f,-5,0,0))
  enddo
  Et(OPS_ACC3(0,0,0)) = Et(OPS_ACC3(-5,0,0))

end subroutine
