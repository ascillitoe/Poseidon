subroutine stress_kernel( x, rho, rhou, mul, tau)

  implicit none

  ! Declare arguments 
  real (kind=8), INTENT(in), DIMENSION(3) :: x,rhou
  real (kind=8), INTENT(in), DIMENSION(1) :: rho, mul
  real (kind=8), DIMENSION(6) :: tau

  ! Local variables
  real (kind=8) uip1, uim1, ujp1, ujm1, ukp1, ukm1, dudx(3,3), &
                div, S11, S22, S33, S12, S13, S23
  integer i,j

  do i = 1,3
    ! Calculate velocity derivatives
    uip1 = rhou(OPS_ACC_MD3(i, 1, 0, 0))/rho(OPS_ACC2( 1, 0, 0))
    uim1 = rhou(OPS_ACC_MD3(i,-1, 0, 0))/rho(OPS_ACC2(-1, 0, 0))
    ujp1 = rhou(OPS_ACC_MD3(i, 0, 1, 0))/rho(OPS_ACC2( 0, 1, 0))
    ujm1 = rhou(OPS_ACC_MD3(i, 0,-1, 0))/rho(OPS_ACC2( 0,-1, 0))
    ukp1 = rhou(OPS_ACC_MD3(i, 0, 0, 1))/rho(OPS_ACC2( 0, 0, 1))
    ukm1 = rhou(OPS_ACC_MD3(i, 0, 0,-1))/rho(OPS_ACC2( 0, 0,-1))

    ! i direction
    dudx(i,1) = (uip1-uim1) / &
                (   x(OPS_ACC_MD1(1,1,0,0)) -    x(OPS_ACC_MD1(1,-1,0,0)))
 
    ! j direction
    dudx(i,2) = (ujp1-ujm1) / &
                (   x(OPS_ACC_MD1(2,0,1,0)) -    x(OPS_ACC_MD1(2,0,-1,0)))
  
    ! k direction
    dudx(i,3) = (ukp1-ukm1) / &
                (   x(OPS_ACC_MD1(3,0,0,1)) -    x(OPS_ACC_MD1(3,0,0,-1)))
  enddo

  ! Calculate Sij
  div = dudx(1,1) + dudx(2,2) + dudx(3,3)
  S11 = dudx(1,1) - (div/3.0_8)
  S22 = dudx(2,2) - (div/3.0_8)
  S33 = dudx(3,3) - (div/3.0_8)
  S12 = 0.5_8*(dudx(1,2)+dudx(2,1))
  S13 = 0.5_8*(dudx(1,3)+dudx(3,1))
  S23 = 0.5_8*(dudx(2,3)+dudx(3,2))

  ! Calculate tauij
  tau(OPS_ACC_MD5(1,0,0,0)) = 2.0_8*mul(OPS_ACC4(0,0,0))*S11
  tau(OPS_ACC_MD5(1,0,0,0)) = 2.0_8*mul(OPS_ACC4(0,0,0))*S22
  tau(OPS_ACC_MD5(1,0,0,0)) = 2.0_8*mul(OPS_ACC4(0,0,0))*S33
  tau(OPS_ACC_MD5(1,0,0,0)) = 2.0_8*mul(OPS_ACC4(0,0,0))*S12
  tau(OPS_ACC_MD5(1,0,0,0)) = 2.0_8*mul(OPS_ACC4(0,0,0))*S13
  tau(OPS_ACC_MD5(1,0,0,0)) = 2.0_8*mul(OPS_ACC4(0,0,0))*S23

end subroutine
