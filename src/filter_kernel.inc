subroutine filter_kernel( rho, rhou, Et, rho_filt, rhou_filt, Et_filt)

  implicit none

  ! Declare arguments 
  real (kind=8), DIMENSION(3), INTENT(IN) :: rhou
  real (kind=8), DIMENSION(1), INTENT(IN) :: rho, Et
  real (kind=8), DIMENSION(3)             :: rhou_filt
  real (kind=8), DIMENSION(1)             :: rho_filt, Et_filt

  ! Average of 6 surrounding points
  rho_filt(OPS_ACC4(0,0,0)) = ( rho(OPS_ACC1(1,0,0)) + rho(OPS_ACC1(-1, 0, 0)) + &
                                rho(OPS_ACC1(0,1,0)) + rho(OPS_ACC1( 0,-1, 0)) + &
                                rho(OPS_ACC1(0,0,1)) + rho(OPS_ACC1( 0, 0,-1)) ) / 6.0_8

  do d = 1,3
    rhou_filt(OPS_ACC_MD5(d,0,0,0)) = ( rhou(OPS_ACC_MD2(d,1,0,0)) + rhou(OPS_ACC_MD2(d,-1, 0, 0)) + &
                                        rhou(OPS_ACC_MD2(d,0,1,0)) + rhou(OPS_ACC_MD2(d, 0,-1, 0)) + &
                                        rhou(OPS_ACC_MD2(d,0,0,1)) + rhou(OPS_ACC_MD2(d, 0, 0,-1)) ) / 6.0_8
  enddo

  Et_filt(OPS_ACC6(0,0,0)) = ( Et(OPS_ACC3(1,0,0)) + Et(OPS_ACC3(-1, 0, 0)) + &
                               Et(OPS_ACC3(0,1,0)) + Et(OPS_ACC3( 0,-1, 0)) + &
                               Et(OPS_ACC3(0,0,1)) + Et(OPS_ACC3( 0, 0,-1)) ) / 6.0_8

  ! Add little bit of average to nodal value
  rho(OPS_ACC1(0,0,0)) = (1.0_8-phi)*rho(OPS_ACC1(0,0,0)) + phi*rho_av
  do d = 1,3
    rhou(OPS_ACC_MD2(d,0,0,0)) = (1.0_8-phi)*rhou(OPS_ACC_MD2(d,0,0,0)) + phi*rhou_av(d)
  enddo
  Et(OPS_ACC3(0,0,0)) = (1.0_8-phi)*Et(OPS_ACC3(0,0,0)) + phi*Et_av

end subroutine
