subroutine dvdx_kernel( x, v, dvdx)
  implicit none

  ! Declare arguments 
  real (kind=8), INTENT(in), DIMENSION(3) :: x, v
  real (kind=8), DIMENSION(9) :: dvdx

  ! Declare local variables
  integer f

  do f = 1,3
    ! i direction
    dvdx(OPS_ACC_MD3(f,0,0,0))   = (v(OPS_ACC_MD2(f,1,0,0)) - v(OPS_ACC_MD2(f,-1,0,0))) / &
                                   (x(OPS_ACC_MD1(1,1,0,0)) - x(OPS_ACC_MD1(1,-1,0,0)))

    ! j direction
    dvdx(OPS_ACC_MD3(f+3,0,0,0)) = (v(OPS_ACC_MD2(f,0,1,0)) - v(OPS_ACC_MD2(f,0,-1,0))) / &
                                   (x(OPS_ACC_MD1(2,0,1,0)) - x(OPS_ACC_MD1(2,0,-1,0)))

    ! k direction
    dvdx(OPS_ACC_MD3(f+6,0,0,0)) = (v(OPS_ACC_MD2(f,0,0,1)) - v(OPS_ACC_MD2(f,0,0,-1))) / &
                                   (x(OPS_ACC_MD1(3,0,0,1)) - x(OPS_ACC_MD1(3,0,0,-1)))
  enddo
end subroutine
