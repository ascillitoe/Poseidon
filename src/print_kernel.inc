subroutine print_kernel( var, var2, idx)

  implicit none

  integer (kind=4), DIMENSION(3), INTENT(IN) :: idx
  real (kind=8), DIMENSION(1), INTENT(IN) :: var
  real (kind=8), DIMENSION(3), INTENT(IN) :: var2

  if ( idx(2).eq.3 .and. idx(3).eq.3) then
  write(*,'(3i3,2f17.9)') idx(1), idx(2), idx(3), var(OPS_ACC1(0,0,0)), var2(OPS_ACC_MD2(1,0,0,0))
  endif

end subroutine
