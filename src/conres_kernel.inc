subroutine conres_kernel( x, rhou, rho_res, idx)

  implicit none

  ! Declare arguments 
  real (kind=8), INTENT(in), DIMENSION(3) :: x,rhou
  real (kind=8), DIMENSION(1) :: rho_res
  integer (kind=4), DIMENSION(3), INTENT(IN) :: idx

  ! Local variables
  real (kind=8) drhoudx(3)

  

  ! i direction
  drhoudx(1) = (rhou(OPS_ACC_MD2(1,1,0,0)) - rhou(OPS_ACC_MD2(1,-1,0,0))) / &
               (   x(OPS_ACC_MD1(1,1,0,0)) -    x(OPS_ACC_MD1(1,-1,0,0)))

  ! j direction
  drhoudx(2) = (rhou(OPS_ACC_MD2(2,0,1,0)) - rhou(OPS_ACC_MD2(2,0,-1,0))) / &
               (   x(OPS_ACC_MD1(2,0,1,0)) -    x(OPS_ACC_MD1(2,0,-1,0)))

  ! k direction
  drhoudx(3) = (rhou(OPS_ACC_MD2(3,0,0,1)) - rhou(OPS_ACC_MD2(3,0,0,-1))) / &
               (   x(OPS_ACC_MD1(3,0,0,1)) -    x(OPS_ACC_MD1(3,0,0,-1)))

  ! RHS
  rho_res(OPS_ACC3(0,0,0)) = drhoudx(1) !+ drhoudx(2) + drhoudx(3)

!if ( idx(2).eq.10 .and. idx(3).eq.10) then
!    write(*,'(3i3,3f17.9)') idx(1), idx(2), idx(3), rhou(OPS_ACC_MD2(1,1,0,0)), rhou(OPS_ACC_MD2(1,-1,0,0)), rho_res(OPS_ACC3(0,0,0))
!  endif

end subroutine
