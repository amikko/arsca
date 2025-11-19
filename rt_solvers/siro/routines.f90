
!--------------
module routines
!--------------

  use parameters

  implicit none

  contains

    !---------------------------
    subroutine read_parameters()
    !---------------------------
    !
    ! Read namelist parameters (siro.nml)

      implicit none

      namelist /param/ noph,usepolar,userefrac,atmos_layers,maxnoord,maxnolay, &
           minweight,req,ratm,step,brdf_reflection,BRF_FILENAME,AER_FILENAME

      open(unit=11,file='siro_settings.nml',status='old')
      read(11,nml=param)
      close(11)

    end subroutine read_parameters


    !-----------------------
    subroutine get_lengths()
    !-----------------------

      implicit none

      integer :: dummy,erro
      real(kind=sp) :: dummy2

      ! check how long the miefile is
      open(unit=12,file=AER_FILENAME,status='old')
      erro = 0
      mielength = 0
      do while (erro==0)
         dummy2=0.0_sp
         read(12,*,iostat=erro) dummy2
         if (dummy2<1.0_sp) mielength = mielength+1
      end do
      close(12)
    end subroutine get_lengths


    !--------------------------
    subroutine read_mie_files()
    !--------------------------

      use global_variables
      use basic_functions, only : phasef

      implicit none

      integer :: i,k,wl_idx
      real(kind=sp) :: cumtablemax, P11_integral,P12_integral,P21_integral,P22_integral, sinmie
      ! Initialise global polarisation variables, inat = naturally polarised light
      inat(1) = 0.5_sp
      inat(2) = 0.5_sp
      inat(3) = 0.0_sp
      inat(4) = 0.0_sp

      ! read mie files
      ! note that the following isn't a loop, because of the intrinsic
      ! difficulty of interpreting strings
      call mie_helper(1,AER_FILENAME)

      wl_idx = current_wl_ind ! should be 1 always
      do i=1,mielength

         phasetable1(i,wl_idx) = (mie_table(i,wl_idx,1) + mie_table(i,wl_idx,6))

      end do
      deltaphasef1(1) = phasetable1(2,wl_idx)-phasetable1(1,wl_idx)
      !deltaphasef2(1) = phasetable2(2)-phasetable2(1)
      deltamy(1) = cosmie(2)-cosmie(1)
      write(*,*) mielength, wl_idx
      cumtable1(1,wl_idx) = 0.5_sp*(phasetable1(1,wl_idx)+phasetable1(2,wl_idx))*deltamy(1)
      !cumtable2(1) = 0.5_sp*(phasetable2(1)+phasetable2(2))*deltamy(1)

      sinmie = sin(0.5_sp*(acos(cosmie(2)) + acos(cosmie(1))))
      P11_integral = 0.5_sp*sinmie*(mie_table(1,wl_idx,1) + mie_table(2,wl_idx,1))*deltamy(1)
      P12_integral = 0.5_sp*sinmie*(mie_table(1,wl_idx,2) + mie_table(2,wl_idx,2))*deltamy(1)
      P21_integral = 0.5_sp*sinmie*(mie_table(1,wl_idx,5) + mie_table(2,wl_idx,5))*deltamy(1)
      P22_integral = 0.5_sp*sinmie*(mie_table(1,wl_idx,6) + mie_table(2,wl_idx,6))*deltamy(1)

      do i=2,mielength-1
         deltaphasef1(i) = phasetable1(i+1,wl_idx)-phasetable1(i,wl_idx)
         !deltaphasef2(i) = phasetable2(i+1)-phasetable2(i)
         deltamy(i) = cosmie(i+1)-cosmie(i)
         sinmie = sin(0.5_sp*(acos(cosmie(i+1)) + acos(cosmie(i))))
         !write (*,*) sinmie
         P11_integral = P11_integral + 0.5_sp*sinmie*(mie_table(i,wl_idx,1)+mie_table(i+1,wl_idx,1))*deltamy(i)
         P12_integral = P12_integral + 0.5_sp*sinmie*(mie_table(i,wl_idx,2)+mie_table(i+1,wl_idx,2))*deltamy(i)
         P21_integral = P21_integral + 0.5_sp*sinmie*(mie_table(i,wl_idx,5)+mie_table(i+1,wl_idx,5))*deltamy(i)
         P22_integral = P22_integral + 0.5_sp*sinmie*(mie_table(i,wl_idx,6)+mie_table(i+1,wl_idx,6))*deltamy(i)
         cumtable1(i,wl_idx) = cumtable1(i-1,wl_idx)+ 0.5_sp*(phasetable1(i,wl_idx)+phasetable1(i+1,wl_idx))*deltamy(i)
      end do
      cumtablemax = cumtable1(mielength-1,wl_idx)
      do i=1,mielength-1
        cumtable1(i,wl_idx) = cumtable1(i,wl_idx) / cumtablemax
      end do
      ! this here normalizes the polarzation matrices.

      !write(*,*) P11_integral


      !end do !wl_idx
    end subroutine read_mie_files

    subroutine mie_helper(k,filename)

      !this subroutine loads the 4x4 muller matrices for each of the wavelengths
      !Author: Antti Mikkonen
      !16.1.2019
      use global_variables
      implicit none
      integer :: i,k
      character(len=35) :: filename

      open(unit=15,file=filename,status='old',form='formatted')

      do i=1,mielength
        read(15,*) cosmie(i),mie_table(i,k,1),mie_table(i,k,2),mie_table(i,k,3), &
        mie_table(i,k,4),mie_table(i,k,5),mie_table(i,k,6),mie_table(i,k,7), &
        mie_table(i,k,8),mie_table(i,k,9),mie_table(i,k,10),mie_table(i,k,11), &
        mie_table(i,k,12),mie_table(i,k,13),mie_table(i,k,14),mie_table(i,k,15), &
        mie_table(i,k,16)
      end do
      close(15)

    end subroutine mie_helper

    !---------------------------------
    subroutine read_waves(wl,wlfactor)
    !---------------------------------

    ! Read wavelenghts

      implicit none

      ! input/output
      real(kind=sp),dimension(nosirowl),intent(out) :: wl,wlfactor

      ! local
      integer :: n

      open(unit=12,file=wlfile,status='old')
      do n=1,nosirowl
         read(12,*) wl(n)
      end do
      close(12)

      ! all the inputs are nanometres now
      ! wl = wl/10.0_sp !conversion from angstroms to nanometres

      ! calculate wlfactor
      do n=1,nosirowl
         wlfactor(n) = 1.0e-6_sp * (64.328_sp + 29498.1_sp / (146.0_sp-1e+6_sp/(wl(n)*wl(n))) &
              + 255.4_sp / (41.0_sp-1e+6_sp/(wl(n)*wl(n))) )
      end do

    end subroutine read_waves

    !-----------------------------------
    subroutine read_altitudes(altitudes)
    !-----------------------------------

      implicit none

      ! input/output
      real(kind=sp), dimension(nosiroalt), intent(out) :: altitudes

      ! local
      integer :: n

      open(unit=12,file=altfile,status='old')
      do n=1,nosiroalt
         read(12,*) altitudes(n)
      end do
      close(12)

    end subroutine read_altitudes

    !-------------------------
    subroutine sun_direction()
    !-------------------------
    ! Calculate direction of Sun at tangent point
    ! (from zenith and azimuth angles)

      use global_variables, only : sunx,suny,sunz

      implicit none

      sunx = cos(pi/180.0_sp*azimuth)*sin(pi/180.0_sp*zenith)
      suny = sin(pi/180.0_sp*azimuth)*sin(pi/180.0_sp*zenith)
      sunz = cos(pi/180.0_sp*zenith)

    end subroutine sun_direction


    !----------------------
    subroutine brdf_setup(brf_M,brf_zen_out,brf_cdf_zen_out,brf_cdf_azi_out)
    !----------------------
      implicit none

      real(kind=sp), dimension(:,:,:,:,:,:), intent(in) :: brf_M
      real(kind=sp), dimension(:), intent(in) :: brf_zen_out
      real(kind=sp), dimension(:,:,:), intent(out) :: brf_cdf_zen_out
      real(kind=sp), dimension(:,:,:,:), intent(out) :: brf_cdf_azi_out
      integer :: i,j,k,l
      real(kind=sp) :: sum
      real(kind=sp), dimension(len_brf_zen_in,len_brf_wavelengths,len_brf_zen_out) :: brf_cdf_zen_out_temp

      write (*,*) len_brf_zen_in, len_brf_wavelengths, len_brf_zen_out, len_brf_azi_out

      do i = 1,len_brf_zen_in ! the incoming zeniths
        do j = 1,len_brf_wavelengths ! the wavelengths
          do k = 1,len_brf_zen_out ! out_zeniths
            brf_cdf_azi_out(i,j,k,1) = brf_M(i,j,k,1,1,1) + brf_M(i,j,k,1,2,2)
            do l = 2,len_brf_azi_out ! out_azis
              brf_cdf_azi_out(i,j,k,l) = brf_cdf_azi_out(i,j,k,l-1) + brf_M(i,j,k,l,1,1) + brf_M(i,j,k,l,2,2)
            end do
            sum = brf_cdf_azi_out(i,j,k,len_brf_azi_out)

            brf_cdf_zen_out_temp(i,j,k) = sum
            do l = 1,len_brf_azi_out
              brf_cdf_azi_out(i,j,k,l) = brf_cdf_azi_out(i,j,k,l) / sum
            end do
            !now brf_cdf_azi_out is normalized
          end do
        end do
      end do

      do i = 1,len_brf_zen_in
        do j = 1,len_brf_wavelengths
          brf_cdf_zen_out(i,j,1) = brf_cdf_zen_out_temp(i,j,1)
          do k = 2,len_brf_zen_out
            brf_cdf_zen_out(i,j,k) = brf_cdf_zen_out(i,j,k-1) + brf_cdf_zen_out_temp(i,j,k)
          end do
          do k = 1,len_brf_zen_out
            !this equalizes the probabilities over the sphere. Maybe there's a
            !proper formula for this too!
            brf_cdf_zen_out(i,j,k) = brf_cdf_zen_out(i,j,k) * sin(brf_zen_out(k))
          end do

          sum = brf_cdf_zen_out(i,j,len_brf_zen_out)

          do k = 1,len_brf_zen_out
            brf_cdf_zen_out(i,j,k) = brf_cdf_zen_out(i,j,k) / sum
          end do
          !write (*,*) brf_cdf_zen_out(i,j,:)
          !now brf_cdf_zen_out is normalized
        end do
      end do
    end subroutine brdf_setup

    subroutine get_closest_value(val,arr,amt,idx)
      implicit none
      real(kind=sp), dimension(:), intent(in) :: arr
      real(kind=sp), intent(in) :: val
      integer, intent(in) :: amt
      integer, intent(out) :: idx

      integer :: temp_idx,i
      real(kind=sp) :: min_dist
      min_dist = 1e14_sp ! practically infinite

      idx = 0

      do i=1,amt

        if (abs(val - arr(i)) < min_dist) then
          min_dist = abs(val - arr(i))
          idx = i
        end if
      end do
      if (idx == 0) then
        write(*,*) 'Error in get_closest_value!'
        write(*,*) amt
        write(*,*) val
        write(*,*) min_dist
        write(*,*) arr
        stop
      end if

    end subroutine get_closest_value

    !-------------------------------------
    subroutine make_coordinates(tanheight)
    !-------------------------------------
    ! Calculate satellite coordinates according to
    ! satellite altitude and tangent point altitude

      use global_variables, only : satx,saty,satz,detx,dety,detz

      implicit none

      ! input
      real(kind=sp), intent(in) :: tanheight

      satz = req + tanheight
      saty = 0.0_sp
      satx = -1.0_sp * sqrt((req+satalt)**2.0_sp-satz**2.0_sp)

      detx = 1.0_sp
      dety = 0.0_sp
      detz = 0.0_sp

    end subroutine make_coordinates


    subroutine tables_general (kumtn,nocells,process,startweight,x,y,z,dx,dy,dz, &
      pathlayer,wl_abs_xsec,wl_sca_xsec,wlfact,abs_prof,sca_prof,transmitted_weight, &
      ground_los)
      use global_variables, only : satx,saty,satz,detx,dety,detz,current_wl
      use atmosphere
      use refra

      implicit none

      !I/O
      real(kind=sp), dimension(maxtable), intent(out) :: kumtn
      integer, intent(out) :: nocells
      integer, dimension(maxtable), intent(out) :: pathlayer
      real(kind=sp), dimension(:,:), intent(in) :: wl_abs_xsec, wl_sca_xsec
      real(kind=sp), dimension(:,:), intent(in) :: abs_prof, sca_prof
      real(kind=sp), dimension(maxtable,nosir), intent(out) :: process
      real(kind=sp), dimension(maxtable), intent(out) :: startweight
      real(kind=sp), dimension(maxtable), intent(out) :: x,y,z,dx,dy,dz
      real(kind=sp), intent(out) :: transmitted_weight
      logical, intent(out) :: ground_los
      real(kind=sp), intent(in) :: wlfact
      ! local
      real(kind=sp) :: r,tauabs,tausir,tauex,kumabs,kumsir,kumex,scalew
      real(kind=sp) :: x_test,y_test,z_test,dxnew,dynew,dznew
      integer :: j,nohalf,lay

      real(kind=sp) :: step_cm,vakio2,vakio3,vakio4,cm_in_km
      real(kind=sp) :: total_dist, step_length, sat_r, small_step_cm, small_step
      integer :: alt_idx, abs_idx, sca_idx, step_diminisher
      logical :: in_atmos

      ! This was added by Antti on 4.5.2023 since tables_general doesn't
      ! take into account the last small and incomplete step. Thus the contribution
      ! of bottom layers of the atmosphere aren't handled properly in Earth viewing
      ! geometries. This can cause up to 10% biases in polarization parameters
      ! in SWIR wavelengths.
      ! Since this is a precomputed table, this shouldn't cause large compuation
      ! overhead and this should ease the last step problem. However, this doesn't
      ! fix it exactly.

      step_diminisher = 4

      cm_in_km = 100000.0_sp
      step_cm = step * cm_in_km
      small_step_cm = step_cm / step_diminisher
      small_step = step / step_diminisher

      x_test = satx
      y_test = saty
      z_test = satz

      in_atmos = .false.
      ground_los = .false.
      total_dist = 0.0_sp !this makes sure that the search ends if we completely
      !miss the atmosphere
      step_length = sqrt(detx**2.0_sp + dety**2.0_sp + detz**2.0_sp)
      sat_r = sqrt(satx**2.0_sp + saty**2.0_sp + satz**2.0_sp)


      do while(.not.in_atmos)
        x_test = x_test + detx
        y_test = y_test + dety
        z_test = z_test + detz
        r = sqrt(x_test**2.0_sp+y_test**2.0_sp+z_test**2.0_sp)
        total_dist = total_dist + step_length
        if ((req < r) .and. (r < ratm)) then
          in_atmos = .true.
        end if
        if (total_dist > 2.0_sp * sat_r + 2.0_sp * ratm) then
          !the cutoff-distance is absurdly long, but those wasted cycles are
          !what you deserve for missing the whole atmosphere
          write (*,*) 'Check the input parameters! The line-of-sight beam doesn"t hit the atmosphere!'
          write (*,*) 'Halting Siro...'
          stop
        end if
      end do

      x(1) = x_test
      y(1) = y_test
      z(1) = z_test
      dx(1) = detx
      dy(1) = dety
      dz(1) = detz
      kumabs = 0.0_sp
      kumsir = 0.0_sp
      j = 1

      do while(in_atmos)

        r = sqrt(x(j)**2.0_sp+y(j)**2.0_sp+z(j)**2.0_sp)

        alt_idx = get_nearest_altitude(r)
        !write (*,*) r,alt_idx
        if (alt_idx > atmos_layers) then
          alt_idx = atmos_layers
        !elseif (index == 0 .or. index == 1) then
        !   index = 2
        elseif (alt_idx <= 0) then
          alt_idx = 1
        end if

        tauabs = 0.0_sp
        do abs_idx=1,noabs
          tauabs = tauabs + abs_prof(alt_idx,abs_idx) * wl_abs_xsec(alt_idx,abs_idx) * small_step_cm
        end do

        kumabs = kumabs + tauabs

        tausir = 0.0_sp
        do sca_idx=1,nosir
          tausir = tausir + sca_prof(alt_idx,sca_idx) * wl_sca_xsec(alt_idx,sca_idx) * small_step_cm
          !write(*,*) sca_prof(alt_idx,sca_idx), wl_sca_xsec(alt_idx,sca_idx), step_cm
          process(j,sca_idx) = tausir
        end do
        !stop
        kumsir = kumsir + tausir
        !write(*,*) kumsir
        !process handles the decisions on the process of scattering (after
        !deciding to scatter), so that pdf needs to be normalized
        process(j,1:nosir) = process(j,1:nosir) / process(j,nosir)

        tauex = tauabs + tausir

        kumex = kumabs + kumsir

        ! uniform sampling
        startweight(j) = tausir * exp(-kumex)
        !write(*,*) startweight(j)
        ! above line corresponds to the fact that tausir is the weight of the photon
        ! which scatters from that point. The cumulative transmittance at that point
        ! is exp(-kumex)
        ! TODO: is uniform samplin the best choice for this situation?
        ! The pdf really depends on which layers scatter the most
        kumtn(j) = real(j,kind=sp)
        !kumtn(j) = startweight(j)
        !kumtn(j) = 1.0_sp - exp(-kumex)

        lay = int(r-req) + 1
        !write(*,*) j, startweight(j), tausir, x(j), y(j), z(j)
        pathlayer(j) = lay

        call calcrefra(x(j),y(j),z(j),dx(j),dy(j),dz(j),dxnew,dynew,dznew,wlfact)
        !write(*,*)dxnew, dynew, dznew
        x(j+1) = x(j) + small_step*(dx(j)+dxnew)/2.0_sp
        y(j+1) = y(j) + small_step*(dy(j)+dynew)/2.0_sp
        z(j+1) = z(j) + small_step*(dz(j)+dznew)/2.0_sp

        !write(*,*) x(j+1),y(j+1),z(j+1)

        dx(j+1) = dx(j)
        dy(j+1) = dy(j)
        dz(j+1) = dz(j)

        j = j + 1

        r = sqrt(x(j)*x(j)+y(j)*y(j)+z(j)*z(j))
        !write(*,*) x(j),y(j),z(j), r
        if ((r > ratm) .or. (req > r)) then
          in_atmos = .false.
          if (req > r) then
            ground_los = .true.
          end if
        end if
      end do

      transmitted_weight = exp(-kumex)

      nocells = max(1,j-1)
      !write(*,*) nocells
      ! The following 10 rows have been disabled by Antti on 11.7.2019 due to nadir problems
      ! tables subroutine considers the first order scattering. If the ground is in los,
      ! then the ground contribution is computed separately in siro.f90
      !if (ground_los) then
      !  !if the ground is in the line of sight, then we sample that too in the
      !  !startweight and kumtn!
      !  nocells = nocells + 1
      !  startweight(j) = transmitted_weight
      !  !kumtn(j) = real(j,kind=sp) / (1.0_sp - transmitted_weight)
      !  kumtn(j) = 1.0_sp
      !  !kumtn(j) = startweight(j)
      !  pathlayer(j) = 1
      !end if
      scalew = kumtn(nocells)

      ! if lisatty 21.4.1999 - jos koko rata hyvin korkealla kumtn(nocells)=0
      if (scalew > 0.0_sp) then
      kumtn(1:nocells) = kumtn(1:nocells)/scalew
      else
      kumtn(nocells) = 1.0_sp
      end if
      !write (*,*) startweight(1:nocells)
      !write (*,*) scalew
      do j = 1,nocells
         startweight(j) = startweight(j)*scalew
      end do
      !do j = 1,nocells
      !  write(*,*) x(j), y(j), z(j)
      !end do
      !scalew = 1.0_sp - exp(-1.0_sp*(kumex))
      !write(*,*) startweight(1:nocells)
      !stop
      !stop
      if (current_wl > 761.04) then
        !write (*,*) transmitted_weight
        !write (*,*) kumtn
        !write (*,*) startweight
        !stop
      end if
    end subroutine tables_general

    subroutine direct_transmissivity(dx_end,dy_end,dz_end,source_angular_radius,direct_los)
      use global_variables, only : sunx,suny,sunz

      implicit none

      real(kind=sp), intent(in) :: dx_end, dy_end, dz_end, source_angular_radius
      logical, intent(out) :: direct_los

      real(kind=sp) :: inner_prod, angle


      inner_prod = dx_end * sunx + dy_end * suny + dz_end * sunz
      angle = acos(inner_prod) / pi * 180.0_sp
      !if beams point in the same direction, then radiance is directly transmitted

      if (angle < source_angular_radius) then
        direct_los = .true.
      else
        direct_los = .false.
      endif

    end subroutine direct_transmissivity

    !---------------------------------------------------------------------------------------------------------------------
    subroutine tables_limb (tanheight,kumtn,nocells,process,startweight,x,y,z,dx,dy,dz,pathlayer,o3cr,no2cr,aircr,aercr,wlfact)
    !---------------------------------------------------------------------------------------------------------------------

      use global_variables, only : satx,saty,satz,detx,dety,detz
      use atmosphere
      use refra

      implicit none

      ! input/output
      real(kind=sp), intent(in) :: tanheight
      real(kind=sp), dimension(atmos_layers),intent(in) :: o3cr,no2cr
      real(kind=sp), intent(in) :: aircr,aercr
      integer, intent(out) :: nocells
      real(kind=sp), dimension(maxtable), intent(out) :: kumtn
      real(kind=sp), dimension(maxtable,nosir), intent(out) :: process
      real(kind=sp), dimension(maxtable), intent(out) :: startweight
      real(kind=sp), dimension(maxtable), intent(out) :: x,y,z,dx,dy,dz
      integer, dimension(maxtable), intent(out) :: pathlayer
      real(kind=sp), intent(in) :: wlfact

      ! local
      real(kind=sp), dimension(maxtable) :: halfx,halfy,halfz,halfdx,halfdy,halfdz
      real(kind=sp) :: r,tauabs,tausir,tauex,kumabs,kumsir,kumex,scalew
      integer :: j,nohalf,lay

      real(kind=sp) :: vakio1,vakio2,vakio3,vakio4
      integer :: index

      vakio1 = step*100000.0_sp
      vakio2 = step*100000.0_sp
      vakio3 = step*aircr*100000.0_sp
      vakio4 = step*aercr*100000.0_sp

      ! Coordinates and propagation direction at tangent point
      halfx(1) = 0.0_sp
      halfy(1) = 0.0_sp
      halfz(1) = req + tanheight
      halfdx(1) = 1.0_sp
      halfdy(1) = 0.0_sp
      halfdz(1) = 0.0_sp
      r = halfz(1)

      ! Trace the second half of the LOS
      j=1
      do while (r <= ratm)
         call calcrefra(halfx(j),halfy(j),halfz(j),halfdx(j),halfdy(j),halfdz(j), &
              halfdx(j+1),halfdy(j+1),halfdz(j+1),wlfact)
         j = j+1
         halfx(j) = halfx(j-1) + step*(halfdx(j-1) + halfdx(j))/2.0_sp
         halfy(j) = halfy(j-1) + step*(halfdy(j-1) + halfdy(j))/2.0_sp
         halfz(j) = halfz(j-1) + step*(halfdz(j-1) + halfdz(j))/2.0_sp
         r = sqrt(halfx(j)**2.0_sp + halfy(j)**2.0_sp + halfz(j)**2.0_sp)
      end do

      nohalf = j - 1

      ! Fill first half
      do j=1,nohalf
         x(j) = -1.0_sp*halfx(nohalf+2-j)
         y(j) = -1.0_sp*halfy(nohalf+2-j)
         z(j) = halfz(nohalf+2-j)
         dx(j) = halfdx(nohalf+2-j)
         dy(j) = halfdy(nohalf+2-j)
         dz(j) = -1.0_sp*(halfdz(nohalf+2-j))
      end do

      ! Fill second half
      do j=nohalf+1,2*nohalf+1
         x(j) = halfx(j-nohalf)
         y(j) = halfy(j-nohalf)
         z(j) = halfz(j-nohalf)
         dx(j) = halfdx(j-nohalf)
         dy(j) = halfdy(j-nohalf)
         dz(j) = halfdz(j-nohalf)
      end do

      nocells = 2*nohalf + 1

      ! calculate probabilitites

      kumabs = 0.0_sp
      kumsir = 0.0_sp

      do j = 1,nocells

         r = sqrt(x(j)**2.0_sp+y(j)**2.0_sp+z(j)**2.0_sp)

         index = get_nearest_altitude(r)

         if (index > atmos_layers) then
            index = atmos_layers
         elseif (index == 0 .or. index == 1) then
            index = 2
         elseif (index < 0) then
            index = 1
         end if

         tauabs = o3prof(index)*o3cr(index)*vakio1 + no2prof(index)*no2cr(index)*vakio2

         kumabs = kumabs + tauabs

         tausir = airprof(index)*vakio3 + aeroprof(index)*vakio4

         process(j,1) = airprof(index)*vakio3
         process(j,2) = tausir
         !process(j,2) = aeroprof(index)*vakio4

         kumsir = kumsir + tausir

         process(j,1:nosir) = process(j,1:nosir) / process(j,nosir)



         tauex = tauabs + tausir

         kumex = kumabs + kumsir

         ! uniform sampling
         startweight(j) = tausir*exp(-kumex)

         kumtn(j) = real(j,kind=sp)

         lay = int(r-req) + 1

         pathlayer(j) = lay

      end do

      scalew = kumtn(nocells)

      do j = 1,nocells
         kumtn(j) = kumtn(j)/scalew
         startweight(j) = startweight(j)*scalew
      end do

    end subroutine tables_limb

    !----------------------------------------------------------------------
    subroutine first_scattering_point(kumtn,nocells,process,startweight, &
         x,y,z,dx,dy,dz,pathlayers,phox,phoy,phoz,dirx,diry,dirz,weight, &
         scatgas,photlayer,transmitted_weight,ground_los,ground_scattering)
    !----------------------------------------------------------------------

      use basic_functions, only : binarysearch

      implicit none

      ! input/output
      real(kind=sp), dimension(:,:), intent(in) :: process
      real(kind=sp), dimension(maxtable), intent(in) :: startweight,kumtn
      real(kind=sp), dimension(maxtable), intent(in) :: x,y,z,dx,dy,dz
      integer, dimension(maxtable), intent(in) ::  pathlayers
      integer, intent(in) :: nocells
      real(kind=sp), intent(in) :: transmitted_weight
      logical, intent(in) :: ground_los
      logical, intent(out) :: ground_scattering


      real(kind=sp), intent(out) :: weight
      real(kind=sp), intent(out) :: dirx,diry,dirz
      real(kind=sp), intent(out) :: phox,phoy,phoz
      integer, dimension(:), intent(out) :: photlayer
      integer, intent(out) :: scatgas

      ! local
      integer :: ind,j
      real(kind=sp) :: random

      !this following if is quite ugly, but it saves an unnecessary call of
      !random if ground_los = .false.

      !if (ground_los) then
      !  call random_number(random)
      !  !this is the similar sort of sampling as done in the siro main photon loop
      !  if (random < transmitted_weight) then
      !    ind = nocells
      !    ground_scattering = .true.
      !  else
      !    call random_number(random)
      !    ind = binarysearch(kumtn,nocells,random)
      !    ground_scattering = .false.
      !  end if
      !else
        call random_number(random)
        ind = binarysearch(kumtn,nocells,random)
        ground_scattering = .false.
      !end if

      !The next 7 rows have been disabled by Antti in 11.7.2019, because of nadir bugs
      !Same reasoning as in tables subroutine: the ground scattering is handled
      !call random_number(random)
      !ind = binarysearch(kumtn,nocells,random)
      !if (ground_los .and. ind == nocells) then
      !  ground_scattering = .true.
      !else
      !  ground_scattering = .false.
      !end if

      phox = x(ind)
      phoy = y(ind)
      phoz = z(ind)
      dirx = dx(ind)
      diry = dy(ind)
      dirz = dz(ind)

      !The next 4 rows have been disabled by Antti in 11.7.2019, because of nadir bugs
      !if (ground_scattering) then
      !  weight = transmitted_weight
      !else
        weight = startweight(ind)
      !end if


      do j=1,ind
         photlayer(pathlayers(j)) = photlayer(pathlayers(j))+1
      end do

      if (.not.ground_scattering) then
        call random_number(random)
        j=1
        do while(process(ind,j) < random)
           j=j+1
        end do
        scatgas = j
      else
        scatgas = 10 !the process index for ground scattering (doesn't do nothing (which means nothing))
      end if

    end subroutine first_scattering_point


    !------------------------------------------------------------------------------------------
    subroutine point_sun(phox,phoy,phoz,dirx,diry,dirz,noscat,weight, &
         detw,cummatrix,detvec,wl_abs_xsec, &
         wl_sca_xsec,wlfact,abs_prof,sca_prof)
    !------------------------------------------------------------------------------------------

      use atmosphere
      use global_variables, only : sunx,suny,sunz,inat
      use basic_functions, only : phasef
      use polarisation, only : forwardrota,phasematrix,maxvecmulti

      implicit none

      ! input/output
      real(kind=sp), dimension(maxnoord),intent(inout) :: detw
      real(kind=sp), intent(in) :: phox,phoy,phoz,dirx,diry,dirz
      real(kind=sp), intent(in) :: weight,wlfact
      integer, intent(in) :: noscat
      real(kind=sp), dimension(:,:), intent(in) :: wl_abs_xsec, wl_sca_xsec
      real(kind=sp), dimension(:,:), intent(in) :: abs_prof, sca_prof

      ! polarization input/output
      real(kind=sp), dimension(4,4), intent(inout) :: cummatrix
      real(kind=sp), dimension(4,maxnoord), intent(inout) :: detvec

      ! local
      real(kind=sp), save :: r,singtn,scatr,ksum,sirx,siry,sirz,absolut,cosi
      real(kind=sp), dimension(maxtable),save :: x,y,z
      real(kind=sp), allocatable, dimension(:), save :: phases,kex
      integer, allocatable, dimension(:),save :: pslayer
      integer,save :: j,k,lay,ord,noce,index,l,abs_idx,sca_idx
      logical,save :: groundhit,half
      real(kind=sp),save :: x1,y1,z1,rin,u,ptulo,limit,rinx,riny,rinz, &
           vakio, cm_in_km, photon_tol, lvec_norm

      ! for polarization..
      real(kind=sp),dimension(4,4),save :: rota,phasem
      real(kind=sp),dimension(4),save :: vec,vec2,sirv,lvec

      !TODO: Check if this and ground_sun could be refactored into something else
      !as they look quite similar...

      !$OMP THREADPRIVATE (r,singtn,scatr,ksum,sirx,siry,sirz,absolut,cosi)
      !$OMP THREADPRIVATE (x,y,z,phases,kex,pslayer,j,k,lay,ord,noce,index)
      !$OMP THREADPRIVATE (groundhit,half,x1,y1,z1,rin,u,ptulo,limit,rinx,riny,rinz)
      !$OMP THREADPRIVATE (vakio,sca_idx)
      if (.not.allocated(phases)) then
        allocate(phases(nosir))
        allocate(kex(nosir))
        allocate(pslayer(maxnolay))
      end if
      cm_in_km = 100000.0_sp
      vakio = step * cm_in_km

      ! scattering orders above maxnoord will be summed to order maxnoord
      if (noscat > maxnoord) then
         ord = maxnoord
      else
         ord = noscat
      end if

      scatr = sqrt(phox*phox+phoy*phoy+phoz*phoz)

      half = .false.

      pslayer = 0

      x1 = phox
      y1 = phoy
      z1 = phoz

      ! "iterate path subroutine" begins
      rinx = phox
      riny = phoy
      rinz = phoz

      call unit_vector(rinx,riny,rinz,rin)

      u = req / rin

      if ((1.0_sp-u**2.0_sp) < epsilon1) then
         limit = sin(refangle)
      else
         limit = sqrt(1.0_sp-u**2.0_sp)*cos(refangle) + u*sin(refangle)
      end if

      ptulo = -1.0_sp*(rinx*sunx+riny*suny+rinz*sunz)

      sirx = sunx
      siry = suny
      sirz = sunz

      if (ptulo > limit) then
         groundhit = .true.         !point definitely in shadow
      else

         if (half) then
            x(1) = x1 + 0.5_sp*step*sunx
            y(1) = y1 + 0.5_sp*step*suny
            z(1) = z1 + 0.5_sp*step*sunz
         else
            x(1) = x1 + step*sunx
            y(1) = y1 + step*suny
            z(1) = z1 + step*sunz
         end if

         r = sqrt(x(1)**2+y(1)**2+z(1)**2)

         noce = 1

         do while ((r < ratm) .and. (r > req))
            noce = noce + 1
            x(noce) = x(noce-1) + step*sunx
            y(noce) = y(noce-1) + step*suny
            z(noce) = z(noce-1) + step*sunz
            r = sqrt(x(noce)**2+y(noce)**2+z(noce)**2)
         end do

         noce = noce - 1

         if (r > req) then
            groundhit = .false.
         else
            groundhit = .true.
         end if

      end if


      if (.not.groundhit) then

         absolut = 0.0_sp

         ! follow the path through the atmosphere
         do j = 1,noce

            r = sqrt(x(j)**2 + y(j)**2 + z(j)**2)

            index = get_nearest_altitude(r)

            if (index < 0) then
               index = 1
            elseif (index == 0 .or. index == 1) then
               index = 2
            elseif (index > atmos_layers) then
               index = atmos_layers
            end if

            do abs_idx = 1,noabs
              absolut = absolut + abs_prof(index,abs_idx) * wl_abs_xsec(index,abs_idx) * vakio
            end do
            do sca_idx = 1,nosir
              absolut = absolut + sca_prof(index,sca_idx) * wl_sca_xsec(index,sca_idx) * vakio
            end do

            lay = int(r-req) + 1

            pslayer(lay) = pslayer(lay) + 1

         end do

         singtn = 0.0_sp

         cosi = sirx*dirx + siry*diry + sirz*dirz

         index = get_nearest_altitude(scatr)

         if (index == 0 .or. index == 1) then
            index = 2
         elseif (index < 0) then
            index = 1
         elseif (index>atmos_layers) then
            index = atmos_layers
         end if

         ksum = 0.0_sp
         do sca_idx = 1,nosir
           kex(sca_idx) = sca_prof(index,sca_idx) * wl_sca_xsec(index,sca_idx) * vakio
           !write(*,*) sca_prof(index,sca_idx), wl_sca_xsec(index,sca_idx)
           ksum = ksum + kex(sca_idx)
           !TODO: This will cause problems: phasef has an input of sca_idx,
           !but there's no guarantee that the scattering kernels are used in
           !the same order as in the input. They may also use same kernels,
           !but with different profiles and cross-sections
           phases(sca_idx) = phasef(sca_idx,cosi) * kex(sca_idx)
           !if (noscat == 1) then
          !   write (*,*) cosi,phasef(sca_idx,cosi)
           !end if
         end do


         do k=1,nosir
            singtn = singtn + phases(k)/ksum
         end do

         singtn = singtn*exp(-absolut)*weight

         !$OMP CRITICAL
         detw(ord) = detw(ord) + singtn

         ! update polarisation variables

         if (usepolar) then

            call forwardrota(dirx,diry,dirz,sirx,siry,sirz,rota)

            sirv = 0.0_sp

            do k=1,nosir
               !write(*,*) 'cosi: ', cosi
               call phasematrix(k,cosi,phasem)
               !write(*,*) 'phasem: ', phasem
               call maxvecmulti(phasem,inat,vec)
               !write(*,*) 'vec: ', vec

               do l=1,4
                  !write(*,*) 'k:', k
                  sirv(l) = sirv(l) + kex(k)/ksum*vec(l)
               end do

            end do

            !if (sirv(1) > 1) then
            !  write (*,*) 'new instance'
            !  write (*,*) singtn
            !  write (*,*) sirv
            !end if

            !write(*,*) 'rota: ', rota
            call maxvecmulti(rota,sirv,vec2)
            call maxvecmulti(cummatrix,vec2,lvec)


            ! Antti commented the following 10 rows 7.3.2023. This should be
            ! fixed in some other place than here! Probably related to
            ! These photon might pack too much a punch
            ! The phenomenon isn't too common, so we'll just scale them to the
            ! photon_tol
            !photon_tol = 10.0_sp
            !lvec_norm = sqrt(lvec(1) * lvec(1) + lvec(2) * lvec(2) + lvec(3) * lvec(3) &
            !    + lvec(4) * lvec(4))
            !if (lvec_norm > photon_tol) then
              !write (*,*) lvec
            !  do l=1,4
            !    lvec(l) = photon_tol / lvec_norm * lvec(l)
            !  end do
            !end if

            do l=1,4
               detvec(l,ord) = detvec(l,ord) + lvec(l)*exp(-absolut)*weight
            end do

         end if
        !$OMP END CRITICAL

      end if   !.not.groundhit

    end subroutine point_sun


    !-----------------------------------------------------------------------
    subroutine scattering_angle(dirx,diry,dirz,scatgas,kosini,cosphi,sinphi)
    !-----------------------------------------------------------------------

      use global_variables
      use basic_functions
      implicit none

      ! input/output
      real(kind=sp), intent(inout) :: dirx,diry,dirz
      real(kind=sp), intent(out) :: kosini,cosphi,sinphi
      integer, intent(in) :: scatgas

      ! local
      real(kind=sp) :: ran1,ran2,ran3,sini,x1,y1,z1,norm,x2,y2,z2,phi
      integer :: ind

      ! parameters for the H-G phase function
      !real(kind=sp) :: aaero,baero,caero


      !aaero = (1.0_sp+g*g) / (2.0_sp*g)
      !baero = -1.0_sp*(1.0_sp-g)**2 / (2.0_sp*g)
      !caero = -2.0_sp*g / (1.0_sp+g)


      if (scatgas == 1) then
         ! NOTE: This is a weirdly efficient way to sample Rayleigh scattering
         ! direction. This isn't an approximation but exact solution.
         call random_number(ran1)
         if (ran1 <= 0.75_sp) then
            call random_number(kosini)
         else
            call random_number(ran1)
            call random_number(ran2)
            call random_number(ran3)
            kosini = max(ran1,ran2,ran3)
         end if
         call random_number(ran1)
         if (ran1 <= 0.5_sp) kosini = -1.0_sp*kosini
      else
         call random_number(ran1)
         ind = binarysearch(cumtable1(1:mielength,current_wl_ind),mielength,ran1)
         kosini = cosmie(ind)

         !The above is for scattering using the mie tables
         !The below is for Henyey-Greenstein scattering of aerosols
         !kosini = 1.0_sp/(1.0_sp + caero*ran1)
         !kosini = aaero + baero*kosini*kosini

      end if

      !write(*,*) kosini

      sini = sqrt(1.0_sp - kosini*kosini)

      x1 = diry ! Miksi tämä tehdään?
      y1 = -dirx
      z1 = 0.0_sp
      norm = sqrt(x1*x1+y1*y1+z1*z1)

      if (norm < 0.001_sp) then ! Tämä voi olla biasin syynä!! T. Antti
         x1 = dirz
         y1 = 0.0_sp
         z1 = -dirx
         norm = sqrt(x1*x1+y1*y1+z1*z1)
      end if

      x1 = x1/norm
      y1 = y1/norm
      z1 = z1/norm

      ! simo: the follwing was changed in the other siro version??!?!

      !x2 = y1*dirz - z1*diry
      !y2 = z1*dirx - x1*dirz
      !z2 = x1*diry - y1*dirx

      x2 = -1.0_sp*y1*dirz + z1*diry
      y2 = -1.0_sp*z1*dirx + x1*dirz
      z2 = -1.0_sp*x1*diry + y1*dirx

      !at this point, the another angle is drawn. The phi angle is independent
      !of a theta-dependent scattering kernel. Phi is around the theta = 0-axis.
      call random_number(ran1)
      phi = 2.0_sp*pi*ran1
      cosphi = cos(phi)
      sinphi = sin(phi)

      dirx = kosini*dirx + sini*(x1*cosphi+x2*sinphi)
      diry = kosini*diry + sini*(y1*cosphi+y2*sinphi)
      dirz = kosini*dirz + sini*(z1*cosphi+z2*sinphi)

    end subroutine scattering_angle


    function get_nearest_altitude(r) result(index)
      use global_variables, only : atmos_layer_thickness
      implicit none
      integer :: index
      real(kind=sp), intent(in) :: r
      !TODO: This function needs to be generalized into a variable thickness
      !layers in the next iteration of Siro
      index = nint((r-req)/atmos_layer_thickness)
    end function get_nearest_altitude

    !----------------------------------------------------------------------------
    subroutine newtrace(phox,phoy,phoz, dirx,diry,dirz,x,y,z, &
     dx,dy,dz,kumtn,nocells,scalew,pathlayers,down,wlfact,wl_abs_xsec, &
     wl_sca_xsec,abs_prof,sca_prof)
    !----------------------------------------------------------------------------

      use atmosphere
      use refra

      implicit none

      ! input/output
      real(kind=sp), intent(in) :: phox,phoy,phoz,dirx,diry,dirz
      real(kind=sp), intent(out) :: scalew
      real(kind=sp), intent(in) :: wlfact
      real(kind=sp), dimension(:,:), intent(in) :: wl_abs_xsec, wl_sca_xsec
      real(kind=sp), dimension(:,:), intent(in) :: abs_prof, sca_prof
      real(kind=sp), dimension(maxtable), intent(out) :: x,y,z,dx,dy,dz
      real(kind=sp), dimension(maxtable), intent(out) :: kumtn
      integer, dimension(maxtable) :: pathlayers
      integer, intent(out) :: nocells
      logical, intent(out) :: down

      ! local
      integer :: j,lay,index,abs_idx,sca_idx
      real(kind=sp) :: r,tauex,vakio,dxnew,dynew,dznew,cm_in_km

      cm_in_km = 100000.0_sp
      vakio = step * cm_in_km

      down = .false.

      j = 1

      r = sqrt(phox**2+phoy**2+phoz**2)

      if ((r-req) < 1.0e-6_sp) then
         x(1) = phox + 0.5_sp*step*dirx
         y(1) = phoy + 0.5_sp*step*diry
         z(1) = phoz + 0.5_sp*step*dirz
      else
         x(1) = phox + step*dirx
         y(1) = phoy + step*diry
         z(1) = phoz + step*dirz
      end if

      dx(1) = dirx
      dy(1) = diry
      dz(1) = dirz

      r = sqrt(x(j)**2.0_sp+y(j)**2.0_sp+z(j)**2.0_sp)

      if (r < req) then

         x(1) = phox
         y(1) = phoy
         z(1) = phoz
         scalew = 0.0_sp   ! scalew=0 = > reflection
         kumtn(1) = 1.0_sp
         nocells = 1

      elseif (r >= ratm) then

         x(1) = phox
         y(1) = phoy
         z(1) = phoz
         scalew = 1.0_sp
         kumtn(1) = 1.0_sp
         nocells = 1

      else

         tauex = 0.0_sp

         do while ((r <= ratm).and.(r > req))
           !write(*,*) 'r'
           !write(*,*) r
            index = get_nearest_altitude(r)
            !write (*,*) 'index:'
            !write (*,*) index
            if (index == 0 .or. index == 1) then
               index = 2
            elseif (index < 0) then
               index = 1
            elseif (index>atmos_layers) then
               index = atmos_layers
            end if

            do abs_idx = 1,noabs
              tauex = tauex + abs_prof(index,abs_idx) * wl_abs_xsec(index,abs_idx) * vakio
            end do
            do sca_idx = 1,nosir
              tauex = tauex + sca_prof(index,sca_idx) * wl_sca_xsec(index,sca_idx) * vakio
            end do

            kumtn(j) = 1.0_sp - exp(-1.0_sp*(tauex))

            lay = int(r-req) + 1
            pathlayers(j) = lay

            ! Take a new step
            call calcrefra(x(j),y(j),z(j),dx(j),dy(j),dz(j),dxnew,dynew,dznew,wlfact)

            x(j+1) = x(j) + step*(dx(j)+dxnew)/2.0_sp
            y(j+1) = y(j) + step*(dy(j)+dynew)/2.0_sp
            z(j+1) = z(j) + step*(dz(j)+dznew)/2.0_sp

            dx(j+1) = dx(j)
            dy(j+1) = dy(j)
            dz(j+1) = dz(j)

            j = j + 1

            r = sqrt(x(j)*x(j)+y(j)*y(j)+z(j)*z(j))

         end do

         ! nocells=j-1
         nocells = max(1,j-1)

         scalew = kumtn(nocells)

         ! if lisatty 21.4.1999 - jos koko rata hyvin korkealla kumtn(nocells)=0
         if (scalew > 0.0_sp) then
            kumtn(1:nocells) = kumtn(1:nocells)/scalew
         else
            kumtn(nocells) = 1.0_sp
         end if

      end if

      if (r <= req) down = .true.


    end subroutine newtrace


    !-----------------------------------------------------------------------------
    subroutine scattering_point(kumtn,nocells,scalew,pathlayers,x,y,z,dx,dy,dz, &
     photlayer,weight,phox,phoy,phoz,dirx,diry,dirz,scatgas,wl_abs_xsec, &
     wl_sca_xsec,wlfact,abs_prof,sca_prof)
    !------------------------------------------------------------------------------

      use basic_functions, only : binarysearch
      use atmosphere

      implicit none

      ! input/output
      real(kind=sp), dimension(maxtable), intent(in) :: x,y,z,dx,dy,dz
      real(kind=sp), dimension(maxtable), intent(in) :: kumtn
      real(kind=sp), dimension(:,:), intent(in) :: wl_abs_xsec, wl_sca_xsec
      real(kind=sp), dimension(:,:), intent(in) :: abs_prof, sca_prof
      real(kind=sp), intent(out) :: phox,phoy,phoz
      real(kind=sp), intent(out) :: dirx,diry,dirz
      real(kind=sp), intent(inout) :: weight, wlfact

      real(kind=sp), intent(in) :: scalew
      integer, dimension(maxnolay), intent(inout) :: photlayer
      integer, dimension(maxtable), intent(in) :: pathlayers
      integer, intent(in) :: nocells
      integer, intent(out) :: scatgas

      ! local
      real(kind=sp) :: random,r,tausir,tauabs
      integer :: ind,i,index,abs_idx,sca_idx
      real(kind=sp), allocatable, dimension(:) :: process

      allocate(process(nosir))
      ! draw a point from kumtn
      call random_number(random)

      ind = binarysearch(kumtn,nocells,random)
      !write(*,*) 'ind:', ind
      r = sqrt(x(ind)**2.0_sp+y(ind)**2.0_sp+z(ind)**2.0_sp)
      !write(*,*) 'r = ', r
      index = get_nearest_altitude(r)

      if (index == 0 .or. index == 1) then
         index = 2
      elseif (index < 0) then
         index = 1
      elseif (index>atmos_layers) then
         index = atmos_layers
      end if

      tausir = 0.0_sp

      do sca_idx = 1,nosir
        tausir = tausir + sca_prof(index,sca_idx) * wl_sca_xsec(index,sca_idx)
        process(sca_idx) = tausir
      end do

      do i=1,nosir
         process(i) = process(i)/tausir
      end do

      call random_number(random)

      i = 1
      do while(process(i) < random)
         i = i+1
      end do

      scatgas = i
      !write(*,*) 'scatgas: ',scatgas
      ! update weight, photlayer and photon coordinates

      !TODO: integrate tauabs linearly through the atmosphere
      !in the future iterations, this should be done with quadraturic integrals
      tauabs = 0.0_sp
      do abs_idx = 1,noabs
        tauabs = tauabs + abs_prof(index,abs_idx) * wl_abs_xsec(index,abs_idx)
      end do

      weight = weight*scalew*tausir/(tausir+tauabs)

      do i=1,ind
         photlayer(pathlayers(i)) = photlayer(pathlayers(i)) + 1
      end do

      phox = x(ind)
      phoy = y(ind)
      phoz = z(ind)

      dirx = dx(ind)
      diry = dy(ind)
      dirz = dz(ind)

    end subroutine scattering_point


    !------------------------------------------------------------------------------
    subroutine ground_sun(phox,phoy,phoz,noscat,weight,detw, &
         cummatrix,detvec,abs_prof,sca_prof,wl_abs_xsec,wl_sca_xsec,diroldx,diroldy,diroldz, &
         brf_zen_in,brf_zen_out,brf_azi_out,brf_M,brf_wavelengths)
    !------------------------------------------------------------------------------
      ! Called when computing the contribution of radiation from the reflection
      ! off the surface.
      use atmosphere
      use global_variables, only : sunx,suny,sunz,inat,brdf_reflection
      use polarisation, only : maxvecmulti

      implicit none

      ! input/output
      real(kind=sp), dimension(maxnoord) :: detw
      real(kind=sp), intent(in) :: phox,phoy,phoz
      real(kind=sp), intent(in) :: diroldx,diroldy,diroldz
      real(kind=sp), intent(in) :: weight
      real(kind=sp), dimension(:,:), intent(in) :: wl_abs_xsec, wl_sca_xsec
      real(kind=sp), dimension(:,:), intent(in) :: abs_prof, sca_prof
      integer, intent(in) :: noscat
      real(kind=sp),dimension(:), intent(in) :: brf_zen_in, brf_zen_out, brf_azi_out, brf_wavelengths
      real(kind=sp),dimension(:,:,:,:,:,:), intent(in) :: brf_M
      ! polarisation
      real(kind=sp), dimension(4,4) :: cummatrix, phasem
      real(kind=sp), dimension(4,maxnoord) :: detvec

      ! local
      real(kind=sp), dimension(maxtable),save :: x,y,z
      real(kind=sp),save :: sirx,siry,sirz,limit,ptulo,rinx,riny,rinz
      real(kind=sp),save :: absolut,r,singtn,scatr,cosi,x1,y1,z1,rin,u
      real(kind=sp),save :: vakio, cm_in_km
      integer, allocatable, dimension(:),save :: pslayer
      logical, save :: half,groundhit
      integer, save :: ord,j,lay,noce,index,l,abs_idx,sca_idx
      ! polarisation
      real(kind=sp), dimension(4), save :: vec
      real(kind=sp) :: ref
      !$OMP THREADPRIVATE(x,y,z,sirx,siry,sirz,limit,ptulo,rinx,riny,rinz)
      !$OMP THREADPRIVATE(absolut,r,singtn,scatr,cosi,x1,y1,z1,rin,u,vakio)
      !$OMP THREADPRIVATE(pslayer,half,groundhit,ord,j,lay,noce,index)
      if (.not.allocated(pslayer)) then
        allocate(pslayer(maxnolay))
      end if

      cm_in_km = 100000.0_sp
      vakio = step * cm_in_km

      if (noscat > maxnoord) then
         ord = maxnoord
      else
         ord = noscat
      end if

      half = .true.

      pslayer = 0

      x1 = phox
      y1 = phoy
      z1 = phoz

      rinx = phox
      riny = phoy
      rinz = phoz

      !write(*,*) phox, phoy, phoz

      call unit_vector(rinx,riny,rinz,rin)

      u = req / rin
      !write (*,*) 'u: ', u
      ! This is probably about figuring out if a surface point is visible due to
      ! refraction in a case where it would be shadowed otherwise
      if ((1.0_sp-u**2.0_sp) < epsilon1) then
         limit = sin(refangle)
      else
         limit = sqrt(1.0_sp-u**2.0_sp)*cos(refangle) + u*sin(refangle)
      end if

      ptulo = -1.0_sp*(rinx*sunx+riny*suny+rinz*sunz)

      sirx = sunx
      siry = suny
      sirz = sunz

      if (ptulo > limit) then

         groundhit = .true. ! point definitely in shadow

      else

         if (half) then
            x(1) = x1 + 0.5_sp*step*sunx
            y(1) = y1 + 0.5_sp*step*suny
            z(1) = z1 + 0.5_sp*step*sunz
         else
            x(1) = x1 + step*sunx
            y(1) = y1 + step*suny
            z(1) = z1 + step*sunz
         end if

         r = sqrt(x(1)**2+y(1)**2+z(1)**2)

         noce = 1

         do while ((r < ratm) .and. (r > req))
            noce = noce + 1
            x(noce) = x(noce-1) + step*sunx
            y(noce) = y(noce-1) + step*suny
            z(noce) = z(noce-1) + step*sunz
            r = sqrt(x(noce)**2+y(noce)**2+z(noce)**2)
         end do

         noce = noce - 1

         if (r > req) then
            groundhit = .false.
         else
            groundhit = .true.
         end if

      end if

      if (.not.groundhit) then

         absolut = 0.0_sp

        ! follow the path through the atmosphere
         do j = 1,noce

            r = sqrt(x(j)**2.0_sp+y(j)**2.0_sp+z(j)**2.0_sp)

            index = get_nearest_altitude(r)

            if (index < 1) then
               index = 1
            elseif (index == 0 .or. index == 1) then
               index = 2
            elseif (index > atmos_layers) then
               index = atmos_layers
            end if

            do abs_idx = 1,noabs
              absolut = absolut + abs_prof(index,abs_idx) * wl_abs_xsec(index,abs_idx) * vakio
            end do
            do sca_idx = 1,nosir
              absolut = absolut + sca_prof(index,sca_idx) * wl_sca_xsec(index,sca_idx) * vakio
            end do

            lay = int(r-req) + 1

            if (lay > 0) then
               pslayer(lay) = pslayer(lay) + 1
            end if

         end do


         !compute the contribution of the reflection
         singtn = 0.0_sp
         scatr = sqrt(phox*phox+phoy*phoy+phoz*phoz)
         cosi = (sirx*phox+siry*phoy+sirz*phoz)/scatr
         if (brdf_reflection) then
           call brdf_get_reflectivity(ref,phox,phoy,phoz,sirx,siry,sirz,diroldx, &
              diroldy,diroldz,brf_zen_in,brf_zen_out,brf_azi_out,brf_M,phasem,brf_wavelengths)
           singtn = cosi/pi*ref
           !write (*,*) singtn
         else
           singtn = cosi/pi*albedo
         end if
         singtn = singtn*exp(-absolut)*weight

         !write(*,*) 'singtn: ', singtn

         ! update detw
         !$OMP CRITICAL
         detw(ord) = detw(ord) + singtn
         !write (*,*) 'ground_sun_detw: ', detw(ord)
         if (usepolar) then
            ! update polarisation parameters
            call maxvecmulti(cummatrix,inat,vec)
            do l=1,4
               detvec(l,ord) = detvec(l,ord)+vec(l)*singtn
            end do
         end if
           !$OMP END CRITICAL
      end if ! .not.groundhit

    end subroutine ground_sun


    !-----------------------------------------------------------------------
    subroutine ground_hit(x,y,z,nocells,pathlayers,phox,phoy,phoz,photlayer)
    !-----------------------------------------------------------------------
      ! sets up the layering properly after hitting the ground
      implicit none

      ! input / output
      real(kind=sp), dimension(maxtable), intent(in) :: x,y,z
      real(kind=sp), intent(out) :: phox,phoy,phoz
      integer, intent(in) :: nocells
      integer, dimension(maxtable), intent(in) :: pathlayers
      integer, dimension(maxnolay), intent(inout) :: photlayer

      ! local
      real(kind=sp) :: norm,deep,scale
      integer :: i

      phox = x(nocells)
      phoy = y(nocells)
      phoz = z(nocells)

      norm = sqrt(phox**2+phoy**2+phoz**2)
      deep = norm - req

      scale = req / norm + 1.0e-7_sp
      phox = phox*scale
      phoy = phoy*scale
      phoz = phoz*scale

      if (deep > (step/2.0_sp)) photlayer(1) = photlayer(1) + 1

      do i=1,nocells-1
         photlayer(pathlayers(i)) = photlayer(pathlayers(i)) + 1
      end do

    end subroutine ground_hit


    !--------------------------------------------------------------------
    subroutine reflection(phox,phoy,phoz,weight,dirx,diry,dirz,cummatrix, &
      brf_zen_in, brf_zen_out, brf_azi_out, brf_cdf_zen_out, brf_cdf_azi_out, brf_M, brf_wavelengths)
    !--------------------------------------------------------------------
      ! Called when a traced photon is scattered off the surface into a new
      ! direction.
      ! Models lambertian reflection with albedo for backward MC
      ! Modified for a BRDF

      use polarisation, only : phasematrix,maxmulti
      use global_variables, only : brdf_reflection
      implicit none

      ! input/output
      real(kind=sp), intent(inout) :: phox,phoy,phoz
      real(kind=sp), intent(inout) :: weight
      real(kind=sp), intent(inout) :: dirx,diry,dirz
      real(kind=sp), dimension(4,4), intent(inout) :: cummatrix

      real(kind=sp), dimension(:) :: brf_zen_in, brf_zen_out, brf_azi_out, brf_wavelengths
      real(kind=sp), dimension(:,:,:) :: brf_cdf_zen_out
      real(kind=sp), dimension(:,:,:,:) :: brf_cdf_azi_out
      real(kind=sp), dimension(:,:,:,:,:,:) :: brf_M

      !local
      integer :: i,j
      logical :: baddir
      real(kind=sp) :: ran1,ptulo,norm,rnorm,ref,diroldx,diroldy,diroldz
      ! for polarization..
      real(kind=sp), dimension(4,4) :: phasem,matrix

      rnorm = sqrt(phox*phox+phoy*phoy+phoz*phoz)

      ! let phox,y,z, be a unit vector for a while
      phox = phox/rnorm
      phoy = phoy/rnorm
      phoz = phoz/rnorm

      diroldx = dirx
      diroldy = diry
      diroldz = dirz

      if (brdf_reflection) then
        !write (*,*) dirx,diry,dirz
        call brdf_sample_direction(phox,phoy,phoz,dirx,diry,dirz, &
           brf_zen_in, brf_zen_out, brf_azi_out, brf_cdf_zen_out, brf_cdf_azi_out)
        !write (*,*) dirx,diry,dirz

      else !sample a lambertian direction for the photon
        baddir = .true.
        do while(baddir)

           call random_number(ran1)
           dirx = ran1 - 0.5_sp

           call random_number(ran1)
           diry = ran1 - 0.5_sp

           call random_number(ran1)
           dirz = ran1 - 0.5_sp

           norm = sqrt(dirx*dirx+diry*diry+dirz*dirz)

           ptulo = dirx*phox + diry*phoy + dirz*phoz

           if ((norm < 0.5_sp).and.(ptulo > 0.0_sp)) then
              baddir = .false.
           end if
        end do

      end if
      ptulo = dirx*phox + diry*phoy + dirz*phoz
      norm = sqrt(dirx*dirx+diry*diry+dirz*dirz)
      dirx = dirx/norm
      diry = diry/norm
      dirz = dirz/norm

      ! add a small number to req to ensure that photon lies above surface
      phox = phox*(req+0.001_sp)
      phoy = phoy*(req+0.001_sp)
      phoz = phoz*(req+0.001_sp)

      if (brdf_reflection) then
        !compute angles, get reflectance from the brdf, update weight
        call brdf_get_reflectivity(ref,phox,phoy,phoz,dirx,diry,dirz,diroldx, &
          diroldy,diroldz,brf_zen_in,brf_zen_out,brf_azi_out,brf_M,phasem,brf_wavelengths)
        weight = ref*weight*2.0_sp*ptulo/norm
      else
        weight = albedo*weight*2.0_sp*ptulo/norm
      end if

      if (usepolar) then
         ! update polarisation matrix
         ! BRDF phase matrix is obtained in the above brdf_get_reflectivity call
         if (.not.brdf_reflection) then
           call phasematrix(10,1.0_sp,phasem)
         endif
         call maxmulti(cummatrix,phasem,matrix)
         do i=1,4
            do j=1,4
               cummatrix(i,j) = matrix(i,j)
            end do
         end do
       end if

    end subroutine reflection


    !-----------------------------
    subroutine unit_vector(x,y,z,r)
    !-----------------------------

      ! Returns vector (x,y,z) scaled to unit length and the original length

      implicit none

      real(kind=sp), intent(inout) :: x,y,z
      real(kind=sp), intent(out) :: r

      r = sqrt(x**2.0_sp + y**2.0_sp + z**2.0_sp)

      if (r > 0.0_sp) then
         x = x/r
         y = y/r
         z = z/r
      else
         write(*,*) 'Error in Unitvector, r = ',r
         write(*,*) 'x = ',x
         write(*,*) 'y = ',y
         write(*,*) 'z = ',z
         stop
      end if

    end subroutine unit_vector

    subroutine brdf_sample_direction(phox,phoy,phoz,dirx,diry,dirz, &
       brf_zen_in, brf_zen_out, brf_azi_out, brf_cdf_zen_out, brf_cdf_azi_out)
      use global_variables, only : brf_wl_idx
      use basic_functions, only : binarysearch
      implicit none
      real(kind=sp), intent(in) :: phox,phoy,phoz
      real(kind=sp), intent(inout) :: dirx,diry,dirz
      real(kind=sp), dimension(:) :: brf_zen_in, brf_zen_out, brf_azi_out
      real(kind=sp), dimension(:,:,:) :: brf_cdf_zen_out
      real(kind=sp), dimension(:,:,:,:) :: brf_cdf_azi_out

      !locals
      real(kind=sp) :: rand_zen, rand_azi
      real(kind=sp) :: diroldx,diroldy,diroldz
      real(kind=sp) :: phonx, phony, phonz
      real(kind=sp) :: rnorm, dotp, zen_in
      real(kind=sp) :: zen_out, azi_out
      integer :: ind_zen, ind_azi, brf_zen_idx

      diroldx = dirx
      diroldy = diry
      diroldz = dirz

      rnorm = sqrt(phox*phox+phoy*phoy+phoz*phoz)
      phonx = phox / rnorm
      phony = phoy / rnorm
      phonz = phoz / rnorm

      dotp = (-diroldx) * phonx + (-diroldy) * phony + (-diroldz) * phonz

      zen_in = acos(dotp)
      call get_closest_value(zen_in,brf_zen_in,len_brf_zen_in,brf_zen_idx)

      call random_number(rand_zen)
      call random_number(rand_azi)

      ind_zen = binarysearch(brf_cdf_zen_out(brf_zen_idx,brf_wl_idx,:),len_brf_zen_out,rand_zen)
      ind_azi = binarysearch(brf_cdf_azi_out(brf_zen_idx,brf_wl_idx,ind_zen,:),len_brf_azi_out,rand_azi)
      zen_out = brf_zen_out(ind_zen)
      azi_out = brf_azi_out(ind_azi)

      !we got the scattering angles. now to figure out where they point.

      call angles_to_directions(phonx, phony, phonz, diroldx, diroldy, diroldz, zen_out, azi_out, dirx, diry, dirz)

    end subroutine brdf_sample_direction

    subroutine brdf_get_reflectivity(ref,phox,phoy,phoz,dirx,diry,dirz,diroldx,diroldy,&
      diroldz,brf_zen_in,brf_zen_out,brf_azi_out,brf_M,phasem,brf_wavelengths)
      use global_variables, only : brf_wl_idx, current_wl
      implicit none

      real(kind=sp), intent(in) :: phox,phoy,phoz
      real(kind=sp), intent(in) :: dirx,diry,dirz
      real(kind=sp), intent(in) :: diroldx,diroldy,diroldz
      real(kind=sp), intent(out) :: ref
      real(kind=sp), intent(out) :: phasem(4,4)
      real(kind=sp),dimension(:), intent(in) :: brf_zen_in, brf_zen_out, brf_azi_out, brf_wavelengths
      real(kind=sp),dimension(:,:,:,:,:,:), intent(in) :: brf_M

      !Ota vaan mullerista (1,1 + 2,2) / 2

      !locals
      real(kind=sp) :: phonx, phony, phonz, fwd_x, fwd_y, fwd_z, flat_x, flat_y, flat_z
      real(kind=sp) :: rnorm, dotp, zen_in, cpx, cpy, cpz, dotp2, azi_out, zen_out, temp
      real(kind=sp) :: dirnx, dirny, dirnz, diroldnx, diroldny, diroldnz, wl_dist, wl_scal
      integer :: zen_out_idx, azi_out_idx, zen_in_idx, i, j, another_brf_wl_idx
      logical :: flat_extrap

      rnorm = sqrt(phox*phox+phoy*phoy+phoz*phoz)
      phonx = phox / rnorm
      phony = phoy / rnorm
      phonz = phoz / rnorm

      rnorm = sqrt(dirx*dirx+diry*diry+dirz*dirz)
      dirnx = dirx / rnorm
      dirny = diry / rnorm
      dirnz = dirz / rnorm

      rnorm = sqrt(diroldx*diroldx+diroldy*diroldy+diroldz*diroldz)
      diroldnx = diroldx / rnorm
      diroldny = diroldy / rnorm
      diroldnz = diroldz / rnorm

      !zen_in_idx
      dotp = (-diroldnx) * phonx + (-diroldny) * phony + (-diroldnz) * phonz
      zen_in = acos(dotp)
      call get_closest_value(zen_in,brf_zen_in,len_brf_zen_in,zen_in_idx)

      !zen_out_idx
      dotp = dirnx * phonx + dirny * phony + dirnz * phonz
      !if (dotp > 1.0_sp) then
      !  dotp = 1.0_sp
      !end if
      zen_out = acos(dotp)
      call get_closest_value(zen_out,brf_zen_out,len_brf_zen_out,zen_out_idx)

      !azi_out_idx

      call flatten(phonx, phony, phonz, -diroldnx, -diroldny, -diroldnz, fwd_x, fwd_y, fwd_z)
      call flatten(phonx, phony, phonz, dirnx, dirny, dirnz, flat_x, flat_y, flat_z)
      call cross_product(phonx, phony, phonz, fwd_x, fwd_y, fwd_z, cpx, cpy, cpz)
      dotp = fwd_x * flat_x + fwd_y * flat_y + fwd_z * flat_z
      if (dotp < -1.0_sp) then
        ! This might happen in weird situations.
        dotp =  -1.0_sp
      end if
      dotp2 = cpx * flat_x + cpy * flat_y + cpz * flat_z
      if (dotp2 > 0) then
        azi_out = acos(dotp)
      else
        azi_out = 2 * pi - acos(dotp)
      end if
      !write (*,*) 'azi'
      !write (*,*) azi_out


      call get_closest_value(azi_out,brf_azi_out,len_brf_azi_out,azi_out_idx)
      ! Periaatteessa ota mullerista elementit (1,1 + 2,2) / 2, mutta nyt
      ! riittää tämä, koska brf_M ei sisällä polarisaatiotietoa

      !To intepolate linearly, we need the closest wl index and the one next to it
      !If it is outside, then we'll just extrapolate using a flat spectrum
      !write (*,*) brf_wavelengths(brf_wl_idx), current_wl
      if (brf_wavelengths(brf_wl_idx) < current_wl) then
        another_brf_wl_idx = brf_wl_idx + 1
      else
        another_brf_wl_idx = brf_wl_idx - 1
      end if
      !write (*,*) another_brf_wl_idx, brf_wl_idx
      !stop
      if (another_brf_wl_idx < 1 .or. another_brf_wl_idx > len_brf_wavelengths) then
        flat_extrap = .true.
      else
        flat_extrap = .false.
        wl_dist = abs(current_wl - brf_wavelengths(brf_wl_idx))
        wl_scal = wl_dist / abs(brf_wavelengths(brf_wl_idx) - brf_wavelengths(another_brf_wl_idx))
        !write (*,*) wl_scal
      end if

      if (flat_extrap) then
        ref = 0.5_sp * ( brf_M(zen_in_idx, brf_wl_idx, zen_out_idx, azi_out_idx, 1, 1) &
        + brf_M(zen_in_idx, brf_wl_idx, zen_out_idx, azi_out_idx, 2, 2) )
      else
        ref = 0.5_sp * (1 - wl_scal) * ( brf_M(zen_in_idx, brf_wl_idx, zen_out_idx, azi_out_idx, 1, 1) &
        + brf_M(zen_in_idx, brf_wl_idx, zen_out_idx, azi_out_idx, 2, 2) ) + &
        0.5_sp * wl_scal * ( brf_M(zen_in_idx, another_brf_wl_idx, zen_out_idx, azi_out_idx, 1, 1) &
        + brf_M(zen_in_idx, another_brf_wl_idx, zen_out_idx, azi_out_idx, 2, 2) )
      end if
      do i=1,4
        do j=1,4
          if (flat_extrap) then
            phasem(i,j) = brf_M(zen_in_idx, brf_wl_idx, zen_out_idx, azi_out_idx, i, j)
          else
            phasem(i,j) = (1 - wl_scal) * brf_M(zen_in_idx, brf_wl_idx, zen_out_idx, azi_out_idx, i, j) &
            + wl_scal * brf_M(zen_in_idx, brf_wl_idx, zen_out_idx, azi_out_idx, i, j)
          end if
        end do
      end do
      ! NOTE: We assume the surface to be non-polarizing!
      ! That means that in this case we'll just compute:
      temp = 0.5_sp * (phasem(1,1) + phasem(2,2))
      phasem(1,1) = 0.5_sp * temp
      phasem(1,2) = 0.5_sp * temp
      phasem(2,1) = 0.5_sp * temp
      phasem(2,2) = 0.5_sp * temp
      !write (*,*) phasem
      !stop

    end subroutine brdf_get_reflectivity

    subroutine angles_to_directions(surf_normx, surf_normy, surf_normz, inc_dirx, inc_diry, &
       inc_dirz, out_zenith, out_azimuth, out_dirx, out_diry, out_dirz)
      implicit none
      real(kind=sp), intent(in) :: surf_normx, surf_normy, surf_normz, inc_dirx, inc_diry, inc_dirz, out_zenith, out_azimuth
      real(kind=sp), intent(out) ::  out_dirx, out_diry, out_dirz

      !locals
      real(kind=sp) :: rot_x, rot_y, rot_z, fwd_x, fwd_y, fwd_z, dotp
      real(kind=sp) :: ax_x, ax_y, ax_z

      call flatten(surf_normx, surf_normy, surf_normz, inc_dirx, inc_diry, inc_dirz, fwd_x, fwd_y, fwd_z)

      call general_rotation(surf_normx, surf_normy, surf_normz, fwd_x, fwd_y, fwd_z, out_azimuth, rot_x, rot_y, rot_z)
      call cross_product(surf_normx, surf_normy, surf_normz, rot_x, rot_y, rot_z, ax_x, ax_y, ax_z)
      call general_rotation(ax_x, ax_y, ax_z, surf_normx, surf_normy, surf_normz, out_zenith, out_dirx, out_diry, out_dirz)

    end subroutine angles_to_directions

    subroutine cross_product(a1,a2,a3,b1,b2,b3,res1,res2,res3)
      implicit none
      real(kind=sp), intent(in) :: a1,a2,a3,b1,b2,b3
      real(kind=sp), intent(out) :: res1, res2, res3
      res1 = a2 * b3 - a3 * b2
      res2 = -(a1 * b3 - a3 * b1)
      res3 = a1 * b2 - a2 * b1
    end subroutine cross_product

    subroutine general_rotation(ax1,ax2,ax3,v1,v2,v3,theta,res1,res2,res3)
      ! rotates the vector v about ax theta much.
      ! needs normalized vectors!
      ! source:  Rodrigues Rotation Wikipedia page
      implicit none
      real(kind=sp), intent(in) :: ax1,ax2,ax3,v1,v2,v3,theta
      real(kind=sp), intent(out) :: res1, res2, res3

      !locals
      real(kind=sp) :: axv1, axv2, axv3, dotp
      call cross_product(ax1,ax2,ax3,v1,v2,v3,axv1, axv2, axv3)
      dotp = ax1 * v1 + ax2 * v2 + ax3 * v3
      res1 = v1 * cos(theta) + axv1 * sin(theta) + ax1 * dotp * (1 - cos(theta))
      res2 = v2 * cos(theta) + axv2 * sin(theta) + ax2 * dotp * (1 - cos(theta))
      res3 = v3 * cos(theta) + axv3 * sin(theta) + ax3 * dotp * (1 - cos(theta))
    end subroutine general_rotation

    subroutine flatten(normx, normy, normz, v1, v2, v3, res1, res2, res3)
      !flattens a vector v onto the surface with surface normal norm
      implicit none
      real(kind=sp), intent(in) :: normx, normy, normz, v1, v2, v3
      real(kind=sp), intent(out) :: res1, res2, res3

      !locals
      real(kind=sp) :: dotp, rnorm
      dotp = normx * v1 + normy * v2 + normz * v3

      ! if 1 + dotp ~= 0, then the vectors are parallel/antiparallel. In that case, the fwd
      ! computation yields a zero vector. We'll fix it with the following check:
      if (abs(1 - abs(dotp)) < epsilon1) then
        res1 = 0.0_sp
        res2 = 0.0_sp
        res3 = 1.0_sp
      else
        res1 = v1 - normx * dotp
        res2 = v2 - normy * dotp
        res3 = v3 - normz * dotp
      end if

      rnorm = sqrt(res1*res1+res2*res2+res3*res3)
      res1 = res1 / rnorm
      res2 = res2 / rnorm
      res3 = res3 / rnorm

    end subroutine flatten

    !----------------------------------------------------------------
    subroutine write_output(single,unscaw,unscavec,output_I,output_sI,output_sQ,output_sU,output_sV,wl,noalt,altitudes,outfilename)
    !----------------------------------------------------------------

      use parameters

      implicit none

      ! input/output
      integer, intent(in) :: noalt
      real(kind=sp),dimension(maxnoord,nosirowl,noalt), intent(in) :: output_I,output_sI,output_sQ,output_sU,output_sV
      real(kind=sp),dimension(nosirowl,noalt), intent(in) :: single
      real(kind=sp),dimension(nosirowl), intent(in) :: wl
      real(kind=sp),dimension(nosirowl,noalt), intent(in) :: unscaw
      real(kind=sp),dimension(4,nosirowl,noalt), intent(in) :: unscavec
      integer,dimension(noalt), intent(in) :: altitudes
      character (len = 128), intent(in) :: outfilename

      ! local
      integer :: n,m
      character(len=100) :: clima
      clima = "radiance"

      !write(filename,'(A4,I3.3,A4,I3.3,A4,I3.3)') '_alb', nint(albedo*100.0_sp), '_zen', nint(zenith), &
      !     '_azi', nint(azimuth)

      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.siro',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), unscaw(n,m), output_I(1,n,m), output_I(2,n,m), output_I(3,n,m), &
                 output_I(4,n,m),output_I(5,n,m), output_I(6,n,m)
                 !write(*,*) altitudes(m), wl(n), unscaw(n,m), output_I(1,n,m), output_I(2,n,m), output_I(3,n,m), &
                !      output_I(4,n,m),output_I(5,n,m), output_I(6,n,m)
         end do
      end do
      close(10)

      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.param',status='unknown')
      write(10,'(f4.1)') version
      write(10,'(i7)') noph
      write(10,'(f5.1)') zenith
      write(10,'(f5.1)') azimuth
      write(10,'(f5.2)') albedo
      write(10,'(f5.2)') g
      write(10,'(f6.1)') satalt
      write(10,'(i2)') climatology
      close(10)

      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.psiroI',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), unscavec(1,n,m), output_sI(1,n,m), output_sI(2,n,m), output_sI(3,n,m), &
                 output_sI(4,n,m),output_sI(5,n,m), output_sI(6,n,m)
         end do
      end do
      close(10)
      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.psiroQ',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), unscavec(2,n,m), output_sQ(1,n,m), output_sQ(2,n,m), output_sQ(3,n,m), &
                 output_sQ(4,n,m),output_sQ(5,n,m), output_sQ(6,n,m)
         end do
      end do
      close(10)
      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.psiroU',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), unscavec(3,n,m), output_sU(1,n,m), output_sU(2,n,m), output_sU(3,n,m), &
                 output_sU(4,n,m),output_sU(5,n,m), output_sU(6,n,m)
         end do
      end do
      close(10)
      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.psiroV',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), unscavec(4,n,m), output_sV(1,n,m), output_sV(2,n,m), output_sV(3,n,m), &
                 output_sV(4,n,m),output_sV(5,n,m), output_sV(6,n,m)
         end do
      end do
      close(10)


    end subroutine write_output

    !----------------------------------------------------------------
    subroutine write_errors(detw_rsd,I_rsd,Q_rsd,U_rsd,V_rsd,wl,noalt,altitudes,outfilename)
    !----------------------------------------------------------------

      use parameters

      implicit none

      ! input/output
      integer, intent(in) :: noalt
      real(kind=sp),dimension(nosirowl,noalt), intent(in) :: detw_rsd,I_rsd,Q_rsd,U_rsd,V_rsd
      real(kind=sp),dimension(nosirowl), intent(in) :: wl
      integer,dimension(noalt), intent(in) :: altitudes
      character (len = 128), intent(in) :: outfilename

      ! local
      integer :: n,m
      character(len=100) :: clima
      clima = "errors"

      !write(filename,'(A4,I3.3,A4,I3.3,A4,I3.3)') '_alb', nint(albedo*100.0_sp), '_zen', nint(zenith), &
      !     '_azi', nint(azimuth)

      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.siro',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), detw_rsd(n,m)
                 !write(*,*) altitudes(m), wl(n), unscaw(n,m), output_I(1,n,m), output_I(2,n,m), output_I(3,n,m), &
                !      output_I(4,n,m),output_I(5,n,m), output_I(6,n,m)
         end do
      end do
      close(10)

      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.param',status='unknown')
      write(10,'(f4.1)') version
      write(10,'(i7)') noph
      write(10,'(f5.1)') zenith
      write(10,'(f5.1)') azimuth
      write(10,'(f5.2)') albedo
      write(10,'(f5.2)') g
      write(10,'(f6.1)') satalt
      write(10,'(i2)') climatology
      close(10)

      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.psiroI',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), I_rsd(n,m)
         end do
      end do
      close(10)
      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.psiroQ',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), Q_rsd(n,m)
         end do
      end do
      close(10)
      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.psiroU',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), U_rsd(n,m)
         end do
      end do
      close(10)
      open(unit=10,file='output/'//trim(clima)//trim(outfilename)//'.psiroV',status='unknown')
      do m = 1,noalt
         do n = 1,nosirowl
            write(10,*) altitudes(m), wl(n), V_rsd(n,m)
         end do
      end do
      close(10)


    end subroutine write_errors


    !-------------------------------------
    function interp1(x,y,xin) result(yout)
    !-------------------------------------

      ! linear interpolation (from M. Laine)

      implicit none

      real(kind=8) :: x(:), y(:), xin
      real(kind=8) :: yout

      integer :: n, il, iu, im
      real(kind=8) :: x1, x2, y1, y2

      n = size(x,1)
      if (size(y,1) .ne. n) then
         stop 'wrong lengts in interp1'
      end if

      ! find location by bisection
      il=1
      iu=n
      do while (iu-il>1)
         im = floor(real(iu+il,kind=sp)/2.0_sp)
         if (xin>x(im)) then
            il=im
         else
            iu=im
         end if
      end do
      x1 = x(il)
      x2 = x(iu)
      y1 = y(il)
      y2 = y(iu)

      ! linear interpolation
      yout = y1 + (y2-y1)/(x2-x1)*(xin-x1)

    end function interp1


    !----------------------
    subroutine write_info()
    !----------------------

      implicit none

      write(*,*) ' '
      write(*,*) ' --------------------'
      write(*,*) ' Running Siro V2.2.1 '
      write(*,*) ' --------------------'
      write(*,*) ' '
      if (climatology==1) then
         write(*,'(2x,a34)') 'Climatology:        OSIRIS+Lowtran tropic'
      elseif (climatology==2) then
         write(*,'(2x,a47)') 'Climatology:        OSIRIS+Lowtran mid-latitude summer'
      elseif (climatology==3) then
         write(*,'(2x,a47)') 'Climatology:        OSIRIS+Lowtran mid-latitude winter'
      else if (climatology==4) then
         write(*,'(2x,a41)') 'Climatology:        OSIRIS+Lowtran antarctica'
      else if (climatology==5) then
         write(*,'(2x,a41)') 'Climatology:        OSIRIS+Lowtran arctic'
      end if

      write(*,'(2x,a7,12x,f5.1)') 'Zenith: ', zenith
      write(*,'(2x,a8,12x,f5.1)') 'Azimuth: ', azimuth
      write(*,'(2x,a7,11x,f5.1)') 'Albedo: ', albedo
      write(*,'(2x,a18,i8)') 'Number of photons: ', noph
      !if (usetemp) then
      !   write(*,'(2x,a39)') 'Temperature:        climatology profile'
      !else
      !   write(*,'(2x,a28)') 'Temperature:        constant'
      !end if
      if (crossec == 1) then
         write(*,'(2x,a43)') 'Cross sections:     GOMOS/OSIRIS resolution'
      elseif (crossec == 2) then
         write(*,'(2x,a48)') 'Cross sections:     GOMOS bright limb resolution'
      end if
      if (userefrac) then
         write(*,'(2x,a22)') 'Refraction:         on'
      else
         write(*,'(2x,a23)') 'Refraction:         off'
      end if

      write(*,*) ' '
      write(*,*) ' -------------------------------'

    end subroutine write_info


  end module routines
