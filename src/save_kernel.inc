subroutine save_kernel(rho_s, rhou_s, Et_s, rho, rhou, Et)

  implicit none

  real (kind=8) , DIMENSION(1) :: rho_s, Et_s
  real (kind=8) , DIMENSION(3) :: rhou_s
  real (kind=8) , INTENT(IN), DIMENSION(1) :: rho, Et
  real (kind=8) , INTENT(IN), DIMENSION(3) :: rhou
 
  integer d

  rho_s(OPS_ACC1(0,0,0))=rho(OPS_ACC4(0,0,0))
  do d = 1,3
    rhou_s(OPS_ACC_MD2(d,0,0,0))=rhou(OPS_ACC_MD5(d,0,0,0))
  enddo
  Et_s(OPS_ACC3(0,0,0))=Et(OPS_ACC6(0,0,0))

end subroutine
