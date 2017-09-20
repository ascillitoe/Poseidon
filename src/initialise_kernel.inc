subroutine initialise_kernel( x, rho, rhou, Et, idx)
  implicit none
  integer (kind=4), DIMENSION(3), INTENT(IN) :: idx
  real (kind=8), DIMENSION(3) :: x, rhou
  real (kind=8), DIMENSION(1) :: rho, Et

  ! x,y,z Coords
  x(OPS_ACC_MD1(1,0,0,0)) = xmin + dx * (idx(1)-1.0_8)
  x(OPS_ACC_MD1(2,0,0,0)) = ymin + dy * (idx(2)-1.0_8)
  x(OPS_ACC_MD1(3,0,0,0)) = zmin + dz * (idx(3)-1.0_8)
  
  ! rho
  rho(OPS_ACC2(0,0,0))  = 2.0_8

  ! rho*u
  rhou(OPS_ACC_MD3(1,0,0,0)) = 1.0_8 !x(OPS_ACC_MD1(1,0,0,0))

  ! rho*v
  rhou(OPS_ACC_MD3(2,0,0,0)) = 0.0_8

  ! rho*w
  rhou(OPS_ACC_MD3(3,0,0,0)) = 0.0_8

  ! Et
  Et(OPS_ACC4(0,0,0))   = 10.0_8

end subroutine
