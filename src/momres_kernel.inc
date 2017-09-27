subroutine momres_kernel( x, rho, rhou, rhou_res)

  implicit none

  ! Declare arguments 
  real (kind=8), INTENT(in), DIMENSION(3) :: x,rhou
  real (kind=8), INTENT(in), DIMENSION(1) :: rho
  real (kind=8), DIMENSION(3) :: rhou_res

  ! Local variables
  real (kind=8) f_ip1, f_im1, &
                f_jp1, f_jm1, &
                f_kp1, f_km1, &
                deriv(3)

  integer j

! loop through x,y,z momentum eqns and calculate the 3 terms to get gradients off
  do j = 1,3
! Build rho*u_i*u_j terms at stencil points
    f_ip1 = rhou(OPS_ACC_MD3(1,1,0,0)) *                     &
            rhou(OPS_ACC_MD3(j,1,0,0))/rho(OPS_ACC2(1,0,0))  !&
            !- Pij(1,j)  TODO: do properly, i+1 for stencil etc
    f_im1 = rhou(OPS_ACC_MD3(1,-1,0,0)) *                    &
            rhou(OPS_ACC_MD3(j,-1,0,0))/rho(OPS_ACC2(-1,0,0)) !&
            !- Pij(1,j)                                       
    f_jp1 = rhou(OPS_ACC_MD3(2,0,1,0)) *                     &
            rhou(OPS_ACC_MD3(j,0,1,0))/rho(OPS_ACC2(0,1,0))  !&
            !- Pij(2,j)
    f_jm1 = rhou(OPS_ACC_MD3(2,0,-1,0)) *                    &
            rhou(OPS_ACC_MD3(j,0,-1,0))/rho(OPS_ACC2(0,-1,0)) !&
            !- Pij(2,j)
    f_kp1 = rhou(OPS_ACC_MD3(3,0,0,1)) *                     &
            rhou(OPS_ACC_MD3(j,0,0,1))/rho(OPS_ACC2(0,0,1))  !&
            !- Pij(3,j)
    f_km1 = rhou(OPS_ACC_MD3(3,0,0,-1)) *                    &
            rhou(OPS_ACC_MD3(j,0,0,-1))/rho(OPS_ACC2(0,0,-1)) !&
            !- Pij(3,j)

    ! Take derivatives
    deriv(1) = (f_ip1 - f_im1)/(x(OPS_ACC_MD1(1,1,0,0))-x(OPS_ACC_MD1(1,-1, 0, 0)))
    deriv(2) = (f_jp1 - f_jm1)/(x(OPS_ACC_MD1(2,0,1,0))-x(OPS_ACC_MD1(2, 0,-1, 0)))
    deriv(3) = (f_kp1 - f_km1)/(x(OPS_ACC_MD1(3,0,0,1))-x(OPS_ACC_MD1(3, 0, 0,-1)))

    ! Sum derivatives to get residual
    rhou_res(OPS_ACC_MD4(j,0,0,0)) = deriv(1) !+ deriv(2) + deriv(3)
  enddo

end subroutine
