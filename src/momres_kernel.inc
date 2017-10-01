subroutine momres_kernel( x, rho, rhou, Et, tau, rhou_res)

  implicit none

  ! Declare arguments 
  real (kind=8), INTENT(in), DIMENSION(3) :: x,rhou
  real (kind=8), INTENT(in), DIMENSION(1) :: rho, Et
  real (kind=8), INTENT(in), DIMENSION(6) :: tau
  real (kind=8), DIMENSION(3) :: rhou_res

  ! Local variables
  real (kind=8) f_ip1, f_im1, &
                f_jp1, f_jm1, &
                f_kp1, f_km1, &
                deriv(3)

  integer i, ij_map(3,3)

  ij_map = reshape((/ 1,4,5, 4,2,6, 5,6,3 /), shape(ij_map))

! loop through x,y,z momentum eqns and calculate the 3 terms to get gradients off
  do i = 1,3

! Build rho*u_i*u_j terms at stencil points
    f_ip1 = rhou(OPS_ACC_MD3(i,1,0,0)) *                     &
            rhou(OPS_ACC_MD3(1,1,0,0))/rho(OPS_ACC2(1,0,0))  &
            - tau(OPS_ACC_MD5(ij_map(1,i),1,0,0))

    f_im1 = rhou(OPS_ACC_MD3(i,-1,0,0)) *                     &
            rhou(OPS_ACC_MD3(1,-1,0,0))/rho(OPS_ACC2(-1,0,0)) &
            - tau(OPS_ACC_MD5(ij_map(1,i),-1,0,0))

    f_jp1 = rhou(OPS_ACC_MD3(i,0,1,0)) *                     &
            rhou(OPS_ACC_MD3(2,0,1,0))/rho(OPS_ACC2(0,1,0))  &
            - tau(OPS_ACC_MD5(ij_map(2,i),0,1,0))

    f_jm1 = rhou(OPS_ACC_MD3(i,0,-1,0)) *                     &
            rhou(OPS_ACC_MD3(2,0,-1,0))/rho(OPS_ACC2(0,-1,0)) &
            - tau(OPS_ACC_MD5(ij_map(2,i),0,-1,0))

    f_kp1 = rhou(OPS_ACC_MD3(i,0,0,1)) *                     &
            rhou(OPS_ACC_MD3(3,0,0,1))/rho(OPS_ACC2(0,0,1))  &
            - tau(OPS_ACC_MD5(ij_map(3,i),0,0,1))

    f_km1 = rhou(OPS_ACC_MD3(i,0,0,-1)) *                     &
            rhou(OPS_ACC_MD3(3,0,0,-1))/rho(OPS_ACC2(0,0,-1)) &
            - tau(OPS_ACC_MD5(ij_map(3,i),0,0,-1))

    ! Add static pressure p=gamma*(gamma-1)*e*rho  (Et = e*rho)
    if (i.eq.1) then
      f_ip1 = f_ip1 + gamma*(gamma-1.0_8)*Et(OPS_ACC4( 1, 0, 0))
      f_im1 = f_im1 + gamma*(gamma-1.0_8)*Et(OPS_ACC4(-1, 0, 0))
    else if (i.eq.2) then
      f_jp1 = f_jp1 + gamma*(gamma-1.0_8)*Et(OPS_ACC4( 0, 1, 0))
      f_jm1 = f_jm1 + gamma*(gamma-1.0_8)*Et(OPS_ACC4( 0,-1, 0))
    else
      f_kp1 = f_kp1 + gamma*(gamma-1.0_8)*Et(OPS_ACC4( 0, 0, 1))
      f_km1 = f_km1 + gamma*(gamma-1.0_8)*Et(OPS_ACC4( 0, 0,-1))
    endif

    ! Take derivatives
    deriv(1) = (f_ip1 - f_im1)/(x(OPS_ACC_MD1(1,1,0,0))-x(OPS_ACC_MD1(1,-1, 0, 0)))
    deriv(2) = (f_jp1 - f_jm1)/(x(OPS_ACC_MD1(2,0,1,0))-x(OPS_ACC_MD1(2, 0,-1, 0)))
    deriv(3) = (f_kp1 - f_km1)/(x(OPS_ACC_MD1(3,0,0,1))-x(OPS_ACC_MD1(3, 0, 0,-1)))

    ! Sum derivatives to get residual
    rhou_res(OPS_ACC_MD5(i,0,0,0)) = deriv(1) + deriv(2) + deriv(3)
  enddo

end subroutine
