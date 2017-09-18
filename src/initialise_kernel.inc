subroutine initialise_kernel( x, rho, rhou, rhov, rhow, Et, idx)
  implicit none
  integer (kind=4), DIMENSION(3), INTENT(IN) :: idx
  real (kind=8), DIMENSION(3) :: x
  real (kind=8), DIMENSION(1) :: rho, rhou, rhov, rhow, Et

  x(OPS_ACC_MD1(1,0,0,0)) = dx * (idx(1)-1.0_8)
  x(OPS_ACC_MD1(2,0,0,0)) = dy * (idx(2)-1.0_8)
  x(OPS_ACC_MD1(3,0,0,0)) = dz * (idx(3)-1.0_8)

  rho(OPS_ACC2(0,0,0))  = 1.177_8
  rhou(OPS_ACC3(0,0,0)) = x(OPS_ACC_MD1(1,0,0,0)) 
  rhov(OPS_ACC4(0,0,0)) = x(OPS_ACC_MD1(2,0,0,0)) 
  rhoW(OPS_ACC5(0,0,0)) = x(OPS_ACC_MD1(3,0,0,0)) 
  Et(OPS_ACC6(0,0,0))   = 0.0_8

end subroutine
