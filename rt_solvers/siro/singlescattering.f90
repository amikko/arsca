!----------------------
module singlescattering
!----------------------

  use atmosphere
  use parameters
  use routines, only : get_nearest_altitude
  implicit none

contains

  !---------------------------------------------------
  subroutine singlelimb(single,o3cr,no2cr,aircr,aercr)
  !---------------------------------------------------

    use global_variables, only : satx,saty,satz,detx,dety,detz

    implicit none

    ! input output
    real(kind=sp), dimension(atmos_layers), intent(in) :: o3cr,no2cr
    real(kind=sp), intent(in) :: aircr,aercr
    real(kind=sp), intent(out) :: single

    ! local
    real(kind=sp) :: scatx,scaty,scatz,r,r0,abs1
    real(kind=sp) :: vakio1,vakio2,vakio3,vakio4
    logical :: inflag
    integer :: index, n

    vakio1 = step*100000.0_sp
    vakio2 = step*100000.0_sp
    vakio3 = step*aircr*100000.0_sp
    vakio4 = step*aercr*100000.0_sp

    single = 0.0_sp
    abs1 = 0.0_sp

    scatx = satx
    scaty = saty
    scatz = satz

    inflag = .true.

    r0 = ratm + 1.0_sp

    r = ratm

    n=1

    do while ((r < ratm) .or. inflag)

       !write(*,*) (r-req)

       ! r = altitude earth center (origo)

       scatx = scatx + step*detx
       scaty = scaty + step*dety
       scatz = scatz + step*detz

       r = sqrt(scatx*scatx+scaty*scaty+scatz*scatz)

       if (r < ratm) then

          if (r > r0) inflag = .false.

          index = get_nearest_altitude(r)

          if (index < 0) then
             index = 1
          elseif (index == 0 .or. index == 1) then
             index = 2
          elseif (index > atmos_layers) then
             index = atmos_layers
          end if

!          write(*,*) o3prof(index),o3cr(index)
!          write(*,*) airprof(index),aircr

          abs1 = abs1 + o3prof(index)*o3cr(index)*vakio1 + no2prof(index)*no2cr(index)*vakio2 + &
               airprof(index)*vakio3 + aeroprof(index)*vakio4



!          write(*,*) abs1
!          stop

          r0 = r

!          write(*,*) n!, abs1
          call leg(single,scatx,scaty,scatz,abs1,o3cr,no2cr,aircr,aercr)
          n=n+1

!          stop

       end if ! end if r > ratm

       !write(*,*) n

    end do ! end while r < ratm or inflag

    !write(*,*) single

    !stop

    !write(*,*) n

  end subroutine singlelimb


  !-------------------------------------------------------------------
  subroutine leg(single,scatx,scaty,scatz,abs1,o3cr,no2cr,aircr,aercr)
  !-------------------------------------------------------------------

    use global_variables, only : sunx,suny,sunz,detx,dety,detz
    use basic_functions, only : phasef

    implicit none

    ! input/output
    real(kind=sp), dimension(atmos_layers), intent(in) :: o3cr,no2cr
    real(kind=sp), intent(in) :: scatx,scaty,scatz,abs1,aircr,aercr
    real(kind=sp), intent(inout) :: single

    ! local
    real(kind=sp) :: r,x2,y2,z2,r2,singtn,abs2,cosi
    real(kind=sp) :: vakio1,vakio2,vakio3,vakio4
    integer :: index

    vakio1 = step*100000.0_sp
    vakio2 = step*100000.0_sp
    vakio3 = step*aircr*100000.0_sp
    vakio4 = step*aercr*100000.0_sp

    cosi = sunx*detx + suny*dety + sunz*detz

    abs2 = 0.0_sp

    r = sqrt(scatx*scatx+scaty*scaty+scatz*scatz)

    x2 = scatx + step*sunx
    y2 = scaty + step*suny
    z2 = scatz + step*sunz
    r2 = sqrt(x2*x2+y2*y2+z2*z2)

    do while ((r2 < ratm) .and. (r2 > req))

       index = get_nearest_altitude(r2)
       if (index < 0) then
          index = 1
       elseif (index == 0 .or. index == 1) then
          index = 2
       elseif (index > atmos_layers) then
          index = atmos_layers
       end if
       abs2 = abs2 + o3prof(index)*o3cr(index)*vakio1 + no2prof(index)*no2cr(index)*vakio2 + &
            airprof(index)*vakio3 + aeroprof(index)*vakio4
       x2 = x2 + step*sunx
       y2 = y2 + step*suny
       z2 = z2 + step*sunz
       r2 = sqrt(x2*x2+y2*y2+z2*z2)
    end do

    !write(*,*) r2,req

    if (r2 > req) then
       index = get_nearest_altitude(r)
       if (index < 0) then
          index = 1
       elseif (index == 0 .or. index == 1) then
          index = 2
       elseif (index > atmos_layers) then
          index = atmos_layers
       end if
       singtn = airprof(index)*vakio3*phasef(1,cosi) + aeroprof(index)*vakio4*phasef(2,cosi)
       singtn = singtn*exp(-abs1)*exp(-abs2)
       single = single + singtn
    end if

  end subroutine leg


end module singlescattering
