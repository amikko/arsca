
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

  !-----------------------------
  subroutine create_atmosphere()
  !-----------------------------
    use global_variables, only : atmos_layer_thickness

    implicit none
    integer :: n
    real(kind=sp) :: dummy

    open(unit=10,file='input/atmosphere_profiles/air_atmos_ALTIUS.dat',status='old')

    ! atmosphere altitude 0-100 km
    !write(*,*) atmos_layers
    do n = 1,atmos_layers
       !atmosalt(n) = 0.1_sp*real((n-1),kind=sp)
        read(10,*) atmosalt(n), dummy
        !write(*,*) n
    end do
    close(10)
    atmos_layer_thickness = atmosalt(2) - atmosalt(1)
    !write(*,*) atmosalt
    !stop

    !-------!
    ! Ozone !
    !-------!

    if (climatology == 0 ) then ! no ozone

       o3prof = 0.0_sp

    else if (climatology == 1) then ! osiris+lowtran (tropic)

       !open(unit=10,file='input/atmosphere_profiles/o3_test_clima.dat',status='old')
       !open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_ozone.dat',status='old')
       !open(unit=10,file='input/atmosphere_profiles/o3_clima_equator_jan_2008.dat',status='old')
       open(unit=10,file='input/atmosphere_profiles/o3_clima_equator_june_2004.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, o3prof(n), dummy, dummy, dummy, dummy
       end do
       close(10)
       o3prof(1) = 0.0_sp

    else if (climatology == 2) then ! osiris+lowtran (mid summer)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_ozone.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, o3prof(n), dummy, dummy, dummy
       end do
       close(10)
       o3prof(1) = 0.0_sp

    else if (climatology == 3) then ! osiris+lowtran (mid winter)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_ozone.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, o3prof(n), dummy, dummy
       end do
       close(10)
       o3prof(1) = 0.0_sp

    else if (climatology == 4) then ! osiris+lowtran (antarctica)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_ozone.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, dummy, o3prof(n), dummy
       end do
       close(10)
       o3prof(1) = 0.0_sp

    else if (climatology == 5) then ! osiris+lowtran (arctic)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_ozone.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, dummy, dummy, o3prof(n)
       end do
       close(10)
       o3prof(1) = 0.0_sp

    else if (climatology == 6) then ! ALTIUS ozone
      open(unit=10,file='input/atmosphere_profiles/o3_atmos_ALTIUS.dat',status='old')
      do n=1,atmos_layers
         read(10,*) dummy, o3prof(n)
      end do
      close(10)
      o3prof(1) = 0.0_sp
    end if

    !else if (climatology == -9999) then ! Erkki's function
    !  do n=1,1001
    !     z = atmosalt(n)
    !     if (z <= 50.0_sp) then
    !        o3prof(n) = 4.5e+11_sp*(z-10.0_sp)**3.0_sp*exp(-z*0.25_sp)
    !        if (z <= 12.0_sp) then
    !           o3prof(n) = 5.95e+11_sp*exp(-z/10.0_sp)
    !        end if
    !     else
    !        o3prof(n) = 1.23_sp*8.7e+10_sp*exp(-(z-50.0_sp)/4.34_sp)+8.0e+8_sp*exp(-(z-83.0_sp)**2.0_sp/25.0_sp)
    !     end if
    !  end do
    !  o3prof(1) = 0.0_sp

    !-----!
    ! NO2 !
    !-----!

    if (climatology == 0) then ! no NO2

       no2prof = 0.0_sp

    else if (climatology == 1) then ! osiris + lowtran (tropic)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_no2.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, no2prof(n), dummy, dummy, dummy, dummy
       end do
       close(10)
       no2prof(1) = 0.0_sp

    else if (climatology == 2) then ! osiris + lowtran (mid summer)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_no2.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, no2prof(n), dummy, dummy, dummy
       end do
       close(10)
       no2prof(1) = 0.0_sp

    else if (climatology == 3) then ! osiris + lowtran (mid winter)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_no2.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, no2prof(n), dummy, dummy
       end do
       close(10)
       no2prof(1) = 0.0_sp

    else if (climatology == 4) then ! osiris + lowtran (antarctica)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_no2.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, dummy, no2prof(n), dummy
       end do
       close(10)
       no2prof(1) = 0.0_sp

    else if (climatology == 5) then ! osiris + lowtran (arctic)

       open(unit=10,file='input/atmosphere_profiles/osiris_and_lowtran_no2.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, dummy, dummy, no2prof(n)
       end do
       close(10)
       no2prof(1) = 0.0_sp

    else if (climatology == 6) then ! ALTIUS no2 (actually aerosols)
      open(unit=10,file='input/atmosphere_profiles/aer_atmos_ALTIUS.dat',status='old')
      do n=1,atmos_layers
         read(10,*) dummy, no2prof(n)
      end do
      close(10)
      no2prof(1) = 0.0_sp
    end if


    !-------------!
    ! Neutral Air !
    !-------------!

    if (climatology == 0 .or. climatology == 1) then ! ecmwf + lowtran (equator)

       !open(unit=10,file='input/atmosphere_profiles/ecmwf_and_lowtran_air.dat',status='old')
       !open(unit=10,file='input/atmosphere_profiles/air_clima_equator_jan_2008.dat',status='old')
       open(unit=10,file='input/atmosphere_profiles/air_clima_equator_june_2004.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, airprof(n), dummy, dummy, dummy, dummy
       end do
       close(10)
       airprof(1) = 0.0_sp

    else if (climatology == 2) then ! ecmwf + lowtran (mid summer)

       open(unit=10,file='input/atmosphere_profiles/ecmwf_and_lowtran_air.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, airprof(n), dummy, dummy, dummy
       end do
       close(10)
       airprof(1) = 0.0_sp

    else if (climatology == 3) then ! ecmwf + lowtran (mid winter)

       !open(unit=10,file='input/atmosphere_profiles/ecmwf_and_lowtran_air.dat',status='old')
       open(unit=10,file='input/atmosphere_profiles/air_clima_mid_june_2008.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, airprof(n), dummy, dummy
       end do
       close(10)
       airprof(1) = 0.0_sp

    else if (climatology == 4) then ! ecmwf + lowtran (antarctica)

       open(unit=10,file='input/atmosphere_profiles/ecmwf_and_lowtran_air.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, dummy, airprof(n), dummy
       end do
       close(10)
       airprof(1) = 0.0_sp

    else if (climatology == 5) then ! ewmwf + lowtran (arctic)

       open(unit=10,file='input/atmosphere_profiles/ecmwf_and_lowtran_air.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, dummy, dummy, dummy, dummy, airprof(n)
       end do
       close(10)
       airprof(1) = 0.0_sp

    else if (climatology == 6) then ! ALTIUS air
      open(unit=10,file='input/atmosphere_profiles/air_atmos_ALTIUS.dat',status='old')
      do n=1,atmos_layers
         read(10,*) dummy, airprof(n)
      end do
      close(10)
      airprof(1) = 0.0_sp

    end if

    !----------!
    ! Aerosols !
    !----------!

    if (climatology == 0) then ! no aerosols

       aeroprof = 0.0_sp

    else if (climatology == 1) then
       do n=1,atmos_layers
          !if (atmosalt(n)<34.3_sp) then
          !   aeroprof(n) = -0.75_sp*atmosalt(n)+26.2_s ! equator 2008 june
          if (atmosalt(n)<33.0_sp) then
             aeroprof(n) = -0.5_sp*atmosalt(n)+17_sp ! equator 2004 june
          else
             aeroprof(n) = 0.0_sp
          end if
       end do


    else if (climatology == 2 .or. climatology == 3 .or. climatology == 4 .or. climatology == 5) then ! lowtran aerosol

       open(unit=10,file='input/atmosphere_profiles/lowtranaeroext525.dat',status='old')
       do n=1,atmos_layers
          read(10,*) dummy, aeroprof(n)
       end do
       close(10)
       aeroprof = aeroprof*1.0_sp/(1.0e+5_sp*3.1e-10_sp/0.525_sp) ! extinction -> number density
       aeroprof(1) = 0.0_sp

    else if (climatology == 6) then ! ALTIUS aerosols
      open(unit=10,file='input/atmosphere_profiles/aer_atmos_ALTIUS.dat',status='old')
      do n=1,atmos_layers
         read(10,*) dummy, aeroprof(n)
      end do
      close(10)
      aeroprof(1) = 0.0_sp
    end if

    !These are for debugging gas profiles
    !write(*,*) 'Profiles:'
    !write(*,*) 'O3:'
    !write(*,*) o3prof
    !write(*,*) 'NO2:'
    !write(*,*) no2prof
    !write(*,*) 'air:'
    !write(*,*) airprof
    !write(*,*) 'aer:'
    !write(*,*) aeroprof

  end subroutine create_atmosphere


end module atmosphere
