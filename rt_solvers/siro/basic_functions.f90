!---------------------
module basic_functions
!---------------------

  ! these procedures don't require other modules
  ! (except parameters and global variables)

  use parameters
  use global_variables

contains

  !------------------------------------------------------
  function binarysearch(kumtn,nocells,random) result(ind)
  !------------------------------------------------------

    implicit none

    ! input/output
    integer, intent(in) :: nocells
    real(kind=sp), dimension(nocells),intent(in) :: kumtn
    real(kind=sp), intent(in) :: random

    ! local
    integer :: ind,k,first,last,apu

    if (random <= kumtn(1)) then
       k = 1
    else
       first = 1
       last = nocells
       apu = last-first
       do while(last-first >1 )
          k = (first+last)/2
          if (random <= kumtn(k)) then
             last = k
          else
             first = k
          end if
       end do
       k = last
    end if
    ind = k

  end function binarysearch

  function interpolate(x,y,N,x_interp) result(y_interp)
    implicit none
    integer, intent(in) :: N
    real(kind=sp), dimension(N), intent(in) :: x, y
    real(kind=sp), intent(in) :: x_interp
    real(kind=sp) :: y_interp
    integer :: ind1, ind2
    ind2 = binarysearch(x,N,x_interp)
    if (ind2 == 1) then
      y_interp = y(ind2)
    else
      ind1 = ind2 - 1
      y_interp = lin_interp(x(ind1),x(ind2),y(ind1),y(ind2),x_interp)
    end if
  end function interpolate

  function lin_interp(x1,x2,y1,y2,x) result(y)
    implicit none
    real(kind=sp), intent(in) :: x1, x2, y1, y2, x
    real(kind=sp) :: y
    y = y1 + (y2 - y1)/(x2 - x1) * (x - x1)/(x2 - x1)
  end function lin_interp

  !--------------------------------------------
  function phasef(index,costheta) result(phase)
  !--------------------------------------------

    implicit none

    ! input/output
    real(kind=sp), intent(in) :: costheta
    integer, intent(in) :: index
    integer :: ind
    real(kind=sp) :: phase

    if (index == 1) then                                       ! Rayleigh scattering phase function
       phase = (0.75_sp/(4.0_sp*pi))*(1.0_sp+costheta**2.0_sp) ! IQUV-kannassa

    elseif ((index == 2).or.(index == 3)) then                                   ! Mie scattering phase function
       !phase = (1.0_sp/(4.0_sp*pi))*(1.0_sp-g**2.0_sp)/ &
       !      ((1.0_sp+g**2.0_sp-2.0_sp*g*costheta)**1.5_sp)
       !Above is HG-phase function for aerosols, below is phasetable from file!
       ind = binarysearch(cosmie(1:mielength),mielength,costheta)
       !1/4pi on oikea!!!
       phase = 1.0_sp / (4.0_sp * pi) * phasetable1(ind,current_wl_ind)
       !write (*,*) costheta, phase
    else
       phase = 0.0_sp
    end if

  end function phasef

end module basic_functions
