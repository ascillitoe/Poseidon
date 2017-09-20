subroutine min_dt_kernel( work, dtmin)

  implicit none
  real (kind=8), DIMENSION(1), INTENT(IN) :: work
  real (kind=8), DIMENSION(1) :: dtmin

  ! Take min dt from i,j,k directions
  dtmin = min(dtmin,work(OPS_ACC1(0,0,0)))

end subroutine
