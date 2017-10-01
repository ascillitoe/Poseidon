subroutine viscosity_kernel( rho, Et, mul )

  implicit none

  real (kind=8), DIMENSION(1), INTENT(IN) :: rho, Et
  real (kind=8), DIMENSION(1)             :: mul

  real (kind=8) Tref,S,muref,C1,p,T

  Tref  = 273.15_8
  S     = 110.4_8
  muref = 1.716d-5
  C1    = muref*(Tref+S)/Tref**(3.0_8/2.0_8)

  ! Static pressure and temp
  p = gamma*(gamma-1.0_8)*Et(OPS_ACC2(0,0,0))
  T = p/(R*rho(OPS_ACC1(0,0,0)))

  ! Laminar viscosity
  mul(OPS_ACC3(0,0,0)) = (C1*T**(3.0_8/2.0_8))/(T+S)
 

end subroutine
