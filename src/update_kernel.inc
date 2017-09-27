subroutine update_kernel( rho, rhou, Et, rho_res, rhou_res, Et_res, dt)

  implicit none

  real (kind=8), DIMENSION(1), INTENT(IN) :: rho_res, Et_res
  real (kind=8), DIMENSION(3), INTENT(IN) :: rhou_res
  real (kind=8), DIMENSION(1)             :: rho, Et
  real (kind=8), DIMENSION(3)             :: rhou
  real (kind=8), INTENT(IN) :: dt

  integer d

  ! Update new variables
  rho(OPS_ACC1(0,0,0)) = rho(OPS_ACC1(0,0,0)) - dt*rho_res(OPS_ACC4(0,0,0))
  do d = 1,3
    rhou(OPS_ACC_MD2(d,0,0,0)) = rhou(OPS_ACC_MD2(d,0,0,0)) - &
                                 dt*rhou_res(OPS_ACC_MD5(d,0,0,0))
  enddo
  Et(OPS_ACC3(0,0,0))  = Et(OPS_ACC3(0,0,0)) - dt*Et_res(OPS_ACC6(0,0,0))

end subroutine
