subroutine profile_kernel( x, rho, rhou, Et, xprof, idx)

  implicit none

  ! Declare arguments 
  real (kind=8), DIMENSION(3), INTENT(IN) :: x, rhou
  real (kind=8), DIMENSION(1), INTENT(IN) :: rho, Et
  real (kind=8), DIMENSION(1) :: xprof
  integer, DIMENSION(3) :: idx

  xprof(OPS_ACC5(0,0,0)) = xprof(OPS_ACC5(0,0,0)) + rho(OPS_ACC2(0,0,0))

  if (idx(2).eq.3 .and. idx(3).eq.3) then
    write(11,'(2f17.9)') x(OPS_ACC_MD1(1,0,0,0)), rho(OPS_ACC2(0,0,0))
  endif


end subroutine
