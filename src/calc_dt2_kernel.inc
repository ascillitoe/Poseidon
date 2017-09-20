subroutine calc_dt2_kernel( x, rho, rhou, work)

  implicit none
  real (kind=8), DIMENSION(3), INTENT(IN) :: x, rhou
  real (kind=8), DIMENSION(1), INTENT(IN) :: rho
  real (kind=8), DIMENSION(1) :: work

  real (kind=8) uav(3),dti, dtj, dtk
  real (kind=8) soundspeed   ! TODO - calc soundspeed properly
  integer d

  soundspeed = 0.0_8

  ! Average velocities across stencil
  do d = 1,3
    uav(d) = 0.5_8*( (rhou(OPS_ACC_MD3(d,1,0,0))/rho(OPS_ACC2(1,0,0))) + &
                     (rhou(OPS_ACC_MD3(d,0,0,0))/rho(OPS_ACC2(0,0,0)))   )
  enddo

  ! Timestep for each stencil
  dti = cfl*(x(OPS_ACC_MD1(1,1,0,0))-x(OPS_ACC_MD1(1,0,0,0))) / &
        (abs(uav(1))+soundspeed)  ! TODO - where to calc speed of sound
  dtj = cfl*(x(OPS_ACC_MD1(2,0,1,0))-x(OPS_ACC_MD1(2,0,0,0))) / &
        (abs(uav(2))+soundspeed)  ! TODO - where to calc speed of sound
  dtk = cfl*(x(OPS_ACC_MD1(3,0,0,1))-x(OPS_ACC_MD1(3,0,0,0))) / &
        (abs(uav(3))+soundspeed)  ! TODO - where to calc speed of sound

  ! Take min dt from i,j,k directions
  work(OPS_ACC4(0,0,0)) = min(dti,dtj,dtk)

end subroutine
