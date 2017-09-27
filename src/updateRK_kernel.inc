subroutine updateRK_kernel( rho, rhou, Et, rho_s, rhou_s, Et_s, &
                            rho_res, rhou_res, Et_res, a1, a2, dt )

  implicit none

  real (kind=8), DIMENSION(1), INTENT(IN) :: rho_res, Et_res
  real (kind=8), DIMENSION(3), INTENT(IN) :: rhou_res
  real (kind=8), DIMENSION(1)             :: rho, Et, rho_s, Et_s
  real (kind=8), DIMENSION(3)             :: rhou, rhou_s
  real (kind=8) :: a1, a2, dt

  integer d
 
  ! Update new variables
  rho(OPS_ACC1(0,0,0)) = rho_s(OPS_ACC4(0,0,0)) + a1*(-rho_res(OPS_ACC7(0,0,0)))*dt 
  do d = 1,3
    rhou(OPS_ACC_MD2(d,0,0,0)) = rhou_s(OPS_ACC_MD5(d,0,0,0)) + &
                                 a1*(-rhou_res(OPS_ACC_MD8(d,0,0,0)))*dt 
  enddo
  Et(OPS_ACC3(0,0,0))  = Et_s(OPS_ACC6(0,0,0)) + a1*(-Et_res(OPS_ACC9(0,0,0)))*dt 

  ! Update old variables
  rho_s(OPS_ACC4(0,0,0)) = rho_s(OPS_ACC4(0,0,0)) + a2*(-rho_res(OPS_ACC7(0,0,0)))*dt 
  do d = 1,3
    rhou_s(OPS_ACC_MD5(d,0,0,0)) = rhou_s(OPS_ACC_MD5(d,0,0,0)) + &
                                   a2*(-rhou_res(OPS_ACC_MD8(d,0,0,0)))*dt
  enddo
  Et_s(OPS_ACC6(0,0,0)) = Et_s(OPS_ACC6(0,0,0)) + a2*(-Et_res(OPS_ACC9(0,0,0)))*dt 

end subroutine
