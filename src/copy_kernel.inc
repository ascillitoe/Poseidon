subroutine copy_kernel( var, work, d)

  implicit none

  real (kind=8), DIMENSION(3), INTENT(IN) :: var
  real (kind=8), DIMENSION(1) :: work
  integer d

  work(OPS_ACC2(0,0,0)) = var(OPS_ACC_MD1(d,0,0,0))

end subroutine
