
!----------------
module atmosphere
!----------------

  use parameters

  implicit none

  real(kind=sp), dimension(:), allocatable :: o3prof,no2prof,airprof,aeroprof,atmosalt

contains

  subroutine create_layering(atmosalt)
    use global_variables, only : atmos_layer_thickness
    !TODO: Change this so that altitudes are the same as the input altitudes
    !      Maybe after the linear integration?
    implicit none
    real(kind=sp), dimension(:), intent(inout) :: atmosalt
    integer :: i
    atmos_layer_thickness = (ratm - req) / atmos_layers
    atmosalt(1) = 0.0_sp
    do i = 2,atmos_layers
      atmosalt(i) = atmosalt(i-1) + atmos_layer_thickness
    end do
  end subroutine create_layering



end module atmosphere
