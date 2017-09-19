subroutine energyres_kernel( x, rho, rhou, Et, Et_res)
  implicit none

  ! Declare arguments 
  real (kind=8), INTENT(in), DIMENSION(3) :: x,rhou
  real (kind=8), INTENT(in), DIMENSION(1) :: rho,Et
  real (kind=8), DIMENSION(3) :: Et_res

  ! Local variables
  real (kind=8) f_ip1, f_im1, &
                f_jp1, f_jm1, &
                f_kp1, f_km1, &
                deriv(3)

! Build terms at stencil points
  f_ip1 = rhou(OPS_ACC_MD3(1,1,0,0))/rho(OPS_ACC2(1,0,0))* &
          Et(OPS_ACC4(1,0,0)) ! TODO - add ujPij and qi terms
  f_im1 = rhou(OPS_ACC_MD3(1,-1,0,0))/rho(OPS_ACC2(-1,0,0))* &
          Et(OPS_ACC4(-1,0,0)) ! TODO - add ujPij and qi terms
  f_jp1 = rhou(OPS_ACC_MD3(2,0,1,0))/rho(OPS_ACC2(0,1,0))* &
          Et(OPS_ACC4(0,1,0)) ! TODO - add ujPij and qi terms
  f_jm1 = rhou(OPS_ACC_MD3(2,0,-1,0))/rho(OPS_ACC2(0,-1,0))* &
          Et(OPS_ACC4(0,-1,0)) ! TODO - add ujPij and qi terms
  f_kp1 = rhou(OPS_ACC_MD3(3,0,0,1))/rho(OPS_ACC2(0,0,1))* &
          Et(OPS_ACC4(0,0,1)) ! TODO - add ujPij and qi terms
  f_km1 = rhou(OPS_ACC_MD3(3,0,0,-1))/rho(OPS_ACC2(0,0,-1))* &
          Et(OPS_ACC4(0,0,-1)) ! TODO - add ujPij and qi terms

  ! Take derivatives
  deriv(1) = 0.5*(f_ip1 - f_im1)/(x(OPS_ACC_MD1(1,1,0,0))-x(OPS_ACC_MD1(1,-1, 0, 0)))
  deriv(2) = 0.5*(f_jp1 - f_jm1)/(x(OPS_ACC_MD1(2,0,1,0))-x(OPS_ACC_MD1(2, 0,-1, 0)))
  deriv(3) = 0.5*(f_kp1 - f_km1)/(x(OPS_ACC_MD1(3,0,0,1))-x(OPS_ACC_MD1(3, 0, 0,-1)))

  ! Sum derivatives to get residual
  Et_res(OPS_ACC5(0,0,0)) = deriv(1) + deriv(2) + deriv(3)

end subroutine
