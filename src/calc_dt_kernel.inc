subroutine calc_dt_kernel( x, rho, rhou, dtmin)

  implicit none

  real (kind=8), DIMENSION(3), INTENT(IN) :: x, rhou
  real (kind=8), DIMENSION(1), INTENT(IN) :: rho
  real (kind=8), DIMENSION(1) :: dtmin

  real (kind=8) uav(3), di, dj, dk, dti, dtj, dtk
  real (kind=8) soundspeed   ! TODO - calc soundspeed properly
  integer d

  soundspeed = 1.0_8

  ! Average velocities across stencil
  do d = 1,3
    uav(d)  = rhou(OPS_ACC_MD3(d,0,0,0))/rho(OPS_ACC2(0,0,0))
  enddo

  ! Grid spacings
  di = 0.5_8*(x(OPS_ACC_MD1(1,1,0,0))-x(OPS_ACC_MD1(1,-1, 0, 0)))
  dj = 0.5_8*(x(OPS_ACC_MD1(2,0,1,0))-x(OPS_ACC_MD1(2, 0,-1, 0)))
  dk = 0.5_8*(x(OPS_ACC_MD1(3,0,0,1))-x(OPS_ACC_MD1(3, 0, 0,-1)))

  ! Timestep for each stencil
  dti = cfl*di / (abs(uav(1))+soundspeed)  ! TODO - where to calc speed of sound
  dtj = cfl*dj / (abs(uav(2))+soundspeed)  ! TODO - where to calc speed of sound
  dtk = cfl*dk / (abs(uav(3))+soundspeed)  ! TODO - where to calc speed of sound

  ! Take min dt from i,j,k directions
  dtmin = min(dtmin,min(dti,dtj,dtk))

end subroutine
