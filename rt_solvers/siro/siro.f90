
!-----------
program siro
!-----------
!
! Siro 2.2.
!
! Fully spherical backward Monte Carlo simulator of the photon
! paths through the atmosphere
!
! Liisa Oikarinen (1965-2002), Simo Tukiainen, Erkki Kyrola, Antti Mikkonen
!
! (c) Finnish Meteorological Institute, 9.9.2008 - 22.5.2020

  use parameters
  use routines
  use atmosphere
  use singlescattering
  use polarisation, only : unitmatrix,polaris
  use global_variables, only : current_wl_ind, current_wl,brf_wl_idx ! these is ONLY for
  use netcdf_io
  !$ use omp_lib
  ! wavelength-dependent muller matrices t. Antti

  implicit none

  ! variables to be allocated dynamically
  real(kind=sp), dimension(:,:), allocatable :: singlescat          ! Analytic Solution Of Single Scattered Radiance
  real(kind=sp), dimension(:,:,:), allocatable :: detw,output_I,output_sI,output_sQ,output_sU,output_sV
  real(kind=sp), dimension(:,:,:,:), allocatable :: detvec
  real(kind=sp), dimension(:,:), allocatable :: detw_rsd, detw_rm, detw_rm_c
  real(kind=sp), dimension(:,:,:), allocatable :: detvec_rsd, detvec_rm, detvec_rm_c
  real(kind=sp), dimension(:,:), allocatable :: unscaw      !if unscattered (transmitted) radiation is measured
  real(kind=sp), dimension(:,:,:), allocatable :: unscavec  !polarized unscattered (transmitted) radiation
  real(kind=sp), dimension(:,:,:), allocatable :: out_unscavec  !polarized unscattered (transmitted) radiation
  real(kind=sp), dimension(:), allocatable :: wl
  real(kind=sp), dimension(:), allocatable :: wlfactor       ! wl-dependent part for refraction
  real(kind=sp), dimension(:,:), allocatable :: incident_stokes
  real(kind=sp), dimension(:), allocatable :: pres_prof, temp_prof
  real(kind=sp), dimension(:,:), allocatable :: abs_prof
  real(kind=sp), dimension(:,:), allocatable :: sca_prof
  !real(kind=sp), dimension(:,:,:), allocatable :: sca_xsec, abs_xsec !these are defined in global_variables.f90
  real(kind=sp), dimension(:,:), allocatable :: wl_sca_xsec, wl_abs_xsec
  real(kind=sp), dimension(:,:), allocatable :: solar_dir
  real(kind=sp), dimension(:,:), allocatable :: sat_position
  real(kind=sp), dimension(:,:), allocatable :: sat_direction
  integer, dimension(:), allocatable :: altitudes
  ! The dimensions of the BRF are 4,3,90,90,4,4
  !                               input_zenith_angles
  !                                 wavelengths
  !                                   output_zenith_angles
  !                                      output_azimuth_angles
  !                                         muller_rows
  !                                           muller_columns
  real(kind=sp), dimension(:,:,:,:,:,:), allocatable :: brf_M
  real(kind=sp), dimension(:), allocatable :: brf_zen_in
  real(kind=sp), dimension(:), allocatable :: brf_wavelengths
  real(kind=sp), dimension(:), allocatable :: brf_zen_out
  real(kind=sp), dimension(:), allocatable :: brf_azi_out
  real(kind=sp), dimension(:,:,:), allocatable :: brf_cdf_zen_out
  real(kind=sp), dimension(:,:,:,:), allocatable :: brf_cdf_azi_out

  ! other variables
  real(kind=sp), dimension(maxtable) :: loskumtn
  real(kind=sp), dimension(maxtable) :: kumtn                 ! cumulative scattering probability
  real(kind=sp), dimension(maxtable) :: losstartweight
  real(kind=sp), allocatable, dimension(:,:) :: losprocess
  real(kind=sp), dimension(maxtable) :: losx,losy,losz
  real(kind=sp), dimension(maxtable) :: losdx,losdy,losdz
  real(kind=sp), dimension(maxtable) :: x,y,z                 ! LOS coordinates
  real(kind=sp), dimension(maxtable) :: dx,dy,dz
  real(kind=sp) :: phox,phoy,phoz                             ! photon coordinates
  real(kind=sp) :: dirx,diry,dirz                             ! photon direction
  real(kind=sp) :: scalew,random
  real(kind=sp) :: weight                                     ! photon intensitivity
  real(kind=sp) :: tanheight                                  ! tangent point altitude
  real(kind=sp) :: wlfact
  real(kind=sp) :: source_angular_radius
  real(kind=sp) :: x0 !for error estimation
  real(kind=sp), dimension(4) :: xv
  ! for polarization
  real(kind=sp), dimension(4,4) :: cummatrix                  ! cumulative polarisation matrix for one photon
  real(kind=sp) :: diroldx,diroldy,diroldz
  real(kind=sp) :: costeta,cosalpha,sinalpha

  integer, dimension(maxtable) :: loslayers
  integer, dimension(maxtable) :: pathlayers
  integer :: losnocells,nocells,ph,n,m,scatgas,noscat,ind,ord,glz,i
  integer, allocatable, dimension(:) :: photlayer
  logical :: down
  logical :: ground_los !is ground in the LOS
  logical :: direct_los !is source in the LOS
  logical :: ground_scattering !if the photon's first scattering point is on ground
  real(kind=sp) :: transmitted_weight
  integer :: DEBUG_scattering_counter
  integer :: wl_progress_counter
  integer :: discarded_photons, noph_wl
  real(kind=sp), dimension(4,1) :: temp_detvec
  real(kind=sp), dimension(1) :: temp_detw

  !get_command_argument is 0-indexed, so the 0 is the name of the executable and
  !1 is the input file name
  call get_command_argument(1,INPUT_FILE_NAME)
  call get_command_argument(2,OUTPUT_FILE_NAME)

  !write(*,*) INPUT_FILE_NAME

  call read_parameters() ! read input parameters from siro.nml file

  call get_lengths() ! get mielength

  call nc_get_lengths() ! get nosirowl, nosiroalt, nosir and noabs

  !if (usepolar) call read_mie_files() ! read mie LUT if we use polarization
  !Nowadays, for the custom aerosol scattering kernel we need to read these too
  !for the unpolarized phase function
  call read_mie_files()

  ! allocate memory:
  allocate(singlescat(nosirowl,nosiroalt))
  allocate(unscaw(nosirowl,nosiroalt))
  allocate(detw(maxnoord,nosirowl,nosiroalt))
  if (MCERRORS) then
    allocate(detw_rsd(nosirowl,nosiroalt))
    allocate(detw_rm(nosirowl,nosiroalt))
    allocate(detw_rm_c(nosirowl,nosiroalt))
  end if
  allocate(output_I(maxnoord,nosirowl,nosiroalt))
  allocate(output_sI(maxnoord,nosirowl,nosiroalt))
  allocate(output_sQ(maxnoord,nosirowl,nosiroalt))
  allocate(output_sU(maxnoord,nosirowl,nosiroalt))
  allocate(output_sV(maxnoord,nosirowl,nosiroalt))
  allocate(wl(nosirowl))
  allocate(wlfactor(nosirowl))
  allocate(incident_stokes(nosirowl,4))
  allocate(altitudes(nosiroalt)) !This is 'altitudes' only as a legacy; it should be 'geometries'
  allocate(unscavec(4,nosirowl,nosiroalt))
  allocate(out_unscavec(4,nosirowl,nosiroalt))
  allocate(detvec(4,maxnoord,nosirowl,nosiroalt))
  allocate(detvec_rsd(4,nosirowl,nosiroalt))
  allocate(detvec_rm(4,nosirowl,nosiroalt))
  allocate(detvec_rm_c(4,nosirowl,nosiroalt))
  allocate(losprocess(maxtable,nosir))
  allocate(photlayer(maxnolay))

  allocate(sca_prof(atmos_layers,nosir))
  allocate(abs_prof(atmos_layers,noabs))
  allocate(pres_prof(atmos_layers))
  allocate(temp_prof(atmos_layers))
  allocate(atmosalt(atmos_layers))

  allocate(sca_xsec(atmos_layers,nosirowl,nosir))
  allocate(abs_xsec(atmos_layers,nosirowl,noabs))
  allocate(wl_sca_xsec(atmos_layers,nosir))
  allocate(wl_abs_xsec(atmos_layers,noabs))

  allocate(solar_dir(nosiroalt,3))
  allocate(sat_position(nosiroalt,3))
  allocate(sat_direction(nosiroalt,3))
  atmosalt=0.0_sp

  if(brdf_reflection) then
    call nc_get_brdf_lengths()
  end if
  allocate(brf_M(len_brf_zen_in,len_brf_wavelengths,len_brf_zen_out,len_brf_azi_out,4,4))
  allocate(brf_zen_in(len_brf_zen_in))
  allocate(brf_wavelengths(len_brf_wavelengths))
  allocate(brf_zen_out(len_brf_zen_out))
  allocate(brf_azi_out(len_brf_azi_out))
  allocate(brf_cdf_zen_out(len_brf_zen_in,len_brf_wavelengths,len_brf_zen_out))
  allocate(brf_cdf_azi_out(len_brf_zen_in,len_brf_wavelengths,len_brf_zen_out,len_brf_azi_out))

  !TODO: Tarkista: sämpläys oikein
  !TODO: Tarkista: get_reflectivity oikein
  !TODO: Tarkista: pure albedot vs brdf-jatke
  !TODO: HUOM! Pitikö brdf_inputissa korjata sinillä!
  ! TODO: tarkista että input luetaan oikein!
  ! TODO: Muokkaa: reflection, ground_sun
  ! TODO: Se sirontafunktio lambertilaiselle pinnalle?
  ! TODO: Selvitä, missä kaikkialla oletataan lambertilaisuus!!!
  ! TODO: uuden suunnan sämpläys cdf:ien avulla
  ! TODO: Polarisaatioon tarvitaan phasematrix, point_sun
  ! NOTE: Nykyisessä muodossaan maasta heijastuva säde ei siroa enää. Tämä voi olla ongelma.
  ! Se oisi ratkaistavissa first_scattering_point:ssa, mutta nyt se on vähän jännittävä
  !write(*,*) nosirowl, nosiroalt

  call nc_read_waves(wl,wlfactor,incident_stokes,source_angular_radius) ! read wavelenghts
  !write(*,*) 'ran nc_read_Waves'
  call create_layering(atmosalt)

  ! create atmosphere profiles: temperature, pressure, scatterers, absorbers and their cross-sections
  call nc_create_atmosphere(atmosalt, sca_prof, abs_prof, pres_prof, temp_prof, &
                  sca_xsec, abs_xsec)


  if (brdf_reflection) then
    call nc_get_brdf(brf_zen_in,brf_wavelengths,brf_zen_out,brf_azi_out,brf_M)
    call brdf_setup(brf_M,brf_zen_out,brf_cdf_zen_out,brf_cdf_azi_out)
    write(*,*) 'No polarization effects of the surface are taken into account due to unverified implementation'
    !This requires changes in polaris functions and routines where you sample from BRDF. 
    !write(*,*) 'Uncomment the polaris call around line number 420 to enable it.'
  end if

  call nc_read_albedo()
  !write(*,*) sca_xsec(:,1,1)
  !write(*,*) nosir
  !write(*,*) 'ran nc_create_atmosphere'
  call nc_get_geometries(solar_dir,sat_position,sat_direction)
  !write(*,*) 'ran_nc_get_geomentries'
  !write(*,*) 'got geometries', solar_dir, sat_position, sat_direction

  !TODO: Add information on all relevant subroutines which you've modified/created
  !TODO: Check out the libradtran netcdf/MOPSMAP phase matrix format,
  !TODO: implement that shit if necessary
  !TODO: otherwise, implement netcdf phase matrix/scattering kernel input

  !TODO: output file format should include specifics on which data it was ran on
  !such as checksums on input/scatterer data. No data is needed for this, except
  !the description field. Also the settings used in the solver and checksum of
  !the program itself and the git commit from which that solver was built from.

  !TODO: netcdf output file setup
  !this actually will happen later -v
  !TODO: photon paths to the output file

  !NOTE: As Siro requires a modernizing rewrite, then these TODOs mentioned above
  ! should be carried out as a part of it.

  call write_info()  ! display information about the simulation on the screen

  ! init some variables:
  output_I = 0.0_sp
  output_sI = 0.0_sp
  output_sQ = 0.0_sp
  output_sU = 0.0_sp
  output_sV = 0.0_sp
  singlescat = 0.0_sp
  detvec = 0.0_sp
  cummatrix = 0.0_sp
  unscaw = 0.0_sp
  unscavec = 0.0_sp

  do m = 1,nosiroalt ! through altitudes
     altitudes(m) = m
     sunx = solar_dir(m,1)
     suny = solar_dir(m,2)
     sunz = solar_dir(m,3)

     satx = sat_position(m,1)
     saty = sat_position(m,2)
     satz = sat_position(m,3)

     detx = sat_direction(m,1)
     dety = sat_direction(m,2)
     detz = sat_direction(m,3)

     ind = 1

     write(*,*) ' '
     write(*,'(2x,a6,1x,i2,1x,a1,1x,i2,2x,a15,f5.1,1x,a2,1x,a13,f5.1,a2)') 'Layer:', m, '/',nosiroalt,'Tangent height:', &
          tanheight, 'km', 'Temperature: ', 0.0_sp,' K'
     write(*,*) ' '

     wl_progress_counter = nosirowl / 10

     do n = 1,nosirowl ! through wavelengths

        current_wl = wl(n)
        if (brdf_reflection) then
          call get_closest_value(current_wl,brf_wavelengths,len_brf_wavelengths,brf_wl_idx)
        end if
        losnocells = 0
        loslayers = 0
        detw = 0.0_sp

        wl_sca_xsec = sca_xsec(:,n,:)
        wl_abs_xsec = abs_xsec(:,n,:)

        if (modulo(n,wl_progress_counter) == 0) then
          write(*,*) n,' out of ', nosirowl, ' wavelengths done'
        end if
        wlfact = wlfactor(n)

        call tables_general (loskumtn,losnocells,losprocess,losstartweight,losx,losy,losz, &
             losdx,losdy,losdz,loslayers,wl_abs_xsec,wl_sca_xsec,wlfact,abs_prof,sca_prof, &
             transmitted_weight,ground_los) ! discretisize LOS

        !if the LOS beam hits the light source
        call direct_transmissivity(losdx(losnocells),losdy(losnocells),losdz(losnocells), &
        source_angular_radius,direct_los)
        if (direct_los) then
          unscaw(n,m) = transmitted_weight
          if (usepolar) then
            do i=1,4
              unscavec(i,n,m) = transmitted_weight * incident_stokes(n,i)
            end do
          end if
        end if

        if (ground_los) then

          !The first scattering event is from the ground
          call ground_hit(losx,losy,losz,losnocells,loslayers,phox,phoy,phoz,photlayer)
          !write (*,*) detw(:,n,m)

          if (usepolar) call unitmatrix(cummatrix)

          call ground_sun(phox,phoy,phoz,1,transmitted_weight,temp_detw(:),cummatrix, &
            temp_detvec(:,:),abs_prof,sca_prof,wl_abs_xsec,wl_sca_xsec,detx,dety,detz, &
            brf_zen_in,brf_zen_out,brf_azi_out,brf_M,brf_wavelengths)
          unscaw(n,m) = temp_detw(1)
          if (usepolar) unscavec(:,n,m) = temp_detvec(:,1)

          temp_detw = 0.0_sp
          temp_detvec = 0.0_sp

        end if

        !FIRST scattering pointissa on handlattava se, että osuuko biimi maahan!
        !!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) &  FIRSTPRIVATE(loskumtn,losnocells,losprocess,losstartweight, &
        !!$OMP losx,losy,losz,losdx,losdy,losdz,loslayers,o3cr,no2cr,aircr,aercr,n,m,wlfact,noph) &
        !!$OMP SHARED (detw,detvec,totalphots)
        ! The above is the original one, but the one below is the new one. The new one
        ! uses a bit more memory per thread, but the old one didn't work for detvec
        !$OMP PARALLEL DEFAULT(FIRSTPRIVATE) &
        !$OMP SHARED (detw,detvec)
        !call random_seed()

        photlayer = 0
        pathlayers = 1

        !$OMP DO
        do ph = 1,noph ! photon loop

           photlayer = photlayer + 1
           if (usepolar) call unitmatrix(cummatrix)

           noscat = 1
           !if (.true.) then
            call first_scattering_point(loskumtn,losnocells,losprocess,losstartweight,losx,losy,losz,losdx, &
                     losdy,losdz,loslayers,phox,phoy,phoz,dirx,diry,dirz,weight,scatgas,photlayer, &
                     transmitted_weight,ground_los,ground_scattering)

                     !if some layers might not scatter at all, then we run into trouble if
                     !we simulate photons with 0 weight. Let's skip them altogether
                     !Added by Antti 28.6.2019

                ! The first scattering event is from the atmosphere
               ! calculate total extinction to the top of the atmopshere (towards sun)
               call point_sun(phox,phoy,phoz,dirx,diry,dirz,noscat,weight,detw(:,n,m), &
                    cummatrix,detvec(:,:,n,m),wl_abs_xsec, &
                    wl_sca_xsec,wlfact,abs_prof,sca_prof)

               diroldx = dirx
               diroldy = diry
               diroldz = dirz

               ! calcuate scattering angle
               call scattering_angle(dirx,diry,dirz,scatgas,costeta,cosalpha,sinalpha)

               if (usepolar) then ! calculate polarization

                  call polaris(dirx,diry,dirz,diroldx,diroldy,diroldz,costeta,cosalpha,sinalpha,scatgas,cummatrix)

               end if

                 !write(*,*) 0, weight

           do while (weight > minweight) ! follow the photon

              ! calculate new path
              call newtrace(phox,phoy,phoz,dirx,diry,dirz,x,y,z, &
                   dx,dy,dz,kumtn,nocells,scalew,pathlayers,down,wlfact,wl_abs_xsec, &
                   wl_sca_xsec,abs_prof,sca_prof)

              if (down) then ! photon is going towards ground

                 call random_number(random)

                 if (random < scalew) then ! scattering event happened before hitting ground

                    ! calculate phox,phoy,phoz and dirx,diry,dirz
                    call scattering_point(kumtn,nocells,1.0_sp,pathlayers,x,y,z,dx,dy,dz,photlayer,weight, &
                         phox,phoy,phoz,dirx,diry,dirz,scatgas,wl_abs_xsec, &
                         wl_sca_xsec,wlfact,abs_prof,sca_prof)


                    noscat = noscat + 1

                    call point_sun(phox,phoy,phoz,dirx,diry,dirz,noscat,weight,detw(:,n,m), &
                         cummatrix,detvec(:,:,n,m),wl_abs_xsec, &
                         wl_sca_xsec,wlfact,abs_prof,sca_prof)

                    diroldx = dirx
                    diroldy = diry
                    diroldz = dirz

                    call scattering_angle(dirx,diry,dirz,scatgas,costeta,cosalpha,sinalpha)

                    if (usepolar) then ! calculate polarization
                       call polaris(dirx,diry,dirz,diroldx,diroldy,diroldz,costeta,cosalpha,sinalpha,scatgas,cummatrix)
                    end if

                 else ! photon hits ground

                    call ground_hit(x,y,z,nocells,pathlayers,phox,phoy,phoz,photlayer)

                    noscat = noscat + 1

                    diroldx = dirx
                    diroldy = diry
                    diroldz = dirz

                    call ground_sun(phox,phoy,phoz,noscat,weight,detw(:,n,m),cummatrix,&
                    detvec(:,:,n,m),abs_prof,sca_prof,wl_abs_xsec,wl_sca_xsec,diroldx,diroldy,diroldz,&
                    brf_zen_in,brf_zen_out,brf_azi_out,brf_M,brf_wavelengths)

                    call reflection(phox,phoy,phoz,weight,dirx,diry,dirz,cummatrix, &
                    brf_zen_in, brf_zen_out, brf_azi_out, brf_cdf_zen_out, brf_cdf_azi_out, brf_M,brf_wavelengths)

                    if (usepolar.and.brdf_reflection) then
                       ! calculate polarization, but only if the surface is non-Lambertian!
                       ! NOTE by Antti on 21.4.2022: We'll just have a nonpolarizing surface now so we don't have to update the
                       ! cumulative polarization matrix.
                       !call polaris(dirx,diry,dirz,diroldx,diroldy,diroldz,costeta,cosalpha,sinalpha,10,cummatrix)
                    end if
                 end if

              else
                 ! calculate phox,phoy,phoz and dirx,diry,dirz
                 call scattering_point(kumtn,nocells,scalew,pathlayers,x,y,z,dx,dy,dz,photlayer,weight, &
                      phox,phoy,phoz,dirx,diry,dirz,scatgas,wl_abs_xsec, &
                      wl_sca_xsec,wlfact,abs_prof,sca_prof)

                 !if (weight < minweight) then
                 !   cycle
                 !end if

                 noscat = noscat + 1

                 call point_sun(phox,phoy,phoz,dirx,diry,dirz,noscat,weight,detw(:,n,m),&
                      cummatrix,detvec(:,:,n,m),wl_abs_xsec, &
                      wl_sca_xsec,wlfact,abs_prof,sca_prof)

                 diroldx = dirx
                 diroldy = diry
                 diroldz = dirz

                 call scattering_angle(dirx,diry,dirz,scatgas,costeta,cosalpha,sinalpha)

                 if (usepolar) then ! calculate polarization
                    call polaris(dirx,diry,dirz,diroldx,diroldy,diroldz,costeta,cosalpha,sinalpha,scatgas,cummatrix)
                 end if

              end if

           end do !photon order of scattering trace loop
           !DEBUG_scattering_counter = DEBUG_scattering_counter + noscat

           if(MCERRORS) then
             ! temporarily save current running means
             detw_rm_c(n,m) = detw_rm(n,m)
             do i = 1, 4
               detvec_rm_c(i,n,m) = detvec_rm(i,n,m)
             end do

             ! update running means (of contributions summed over scattering orders)
             x0 = 0.0_sp
             do ord = 1, maxnoord
               x0 = x0 + detw(ord,n,m)
             end do
             detw_rm(n,m) = detw_rm_c(n,m) + (x0 - detw_rm_c(n,m))/ph

             xv = 0.0_sp
             do ord = 1, maxnoord
               xv(1) = xv(1) + detvec(1,ord,n,m) + detvec(2,ord,n,m)
               xv(2) = xv(2) + detvec(1,ord,n,m) - detvec(2,ord,n,m)
               xv(3) = xv(3) + detvec(3,ord,n,m)
               xv(4) = xv(4) + detvec(4,ord,n,m)
             end do
             do i = 1, 4
               detvec_rm(i,n,m) = detvec_rm_c(i,n,m) + (xv(i) - detvec_rm_c(i,n,m))/ph
             end do

             ! update running standard deviations (of contributions summed over scattering orders)
             detw_rsd(n,m) = detw_rsd(n,m) + (x0 - detw_rm_c(n,m))*(x0 - detw_rm(n,m))

             do i = 1, 4
               detvec_rsd(i,n,m) = detvec_rsd(i,n,m) + (xv(i) - detvec_rm_c(i,n,m))*(xv(i) - detvec_rm(i,n,m))
             end do
           end if ! end MCERRORS
        end do ! end of photon loop

        !$OMP END DO
        !$OMP END PARALLEL

        do ord=1,maxnoord
            output_sI(ord,n,m) = (detvec(1,ord,n,m) + detvec(2,ord,n,m)) / noph
            output_sQ(ord,n,m) = (detvec(1,ord,n,m) - detvec(2,ord,n,m)) / noph
            output_sU(ord,n,m) = - detvec(3,ord,n,m) / noph
            output_sV(ord,n,m) = detvec(4,ord,n,m) / noph
        end do
        output_I(:,n,m) = detw(:,n,m) / noph

        out_unscavec(1,n,m) = unscavec(1,n,m) + unscavec(2,n,m)
        out_unscavec(2,n,m) = unscavec(1,n,m) - unscavec(2,n,m)
        out_unscavec(3,n,m) = - unscavec(3,n,m)
        out_unscavec(4,n,m) = unscavec(4,n,m)


     end do ! end of wavelength loop

  end do ! end of altitude loop
  write (*,*) 'Siro done simulating ', nosirowl, ' wavelenghts in ', nosiroalt, ' geometries.'

  call write_output(singlescat,unscaw,out_unscavec,output_I,output_sI,output_sQ, &
  output_sU,output_sV,wl,nosiroalt,altitudes,OUTPUT_FILE_NAME)

  if(MCERRORS) then
    call write_errors(SQRT(ABS(detw_rsd)) / (noph * noph), &
    SQRT(ABS(detvec_rsd(1,:,:))) / (noph * noph), &
    SQRT(ABS(detvec_rsd(2,:,:))) / (noph * noph), &
    SQRT(ABS(detvec_rsd(3,:,:))) / (noph * noph), &
    SQRT(ABS(detvec_rsd(4,:,:))) / (noph * noph), &
    wl,nosiroalt,altitudes,OUTPUT_FILE_NAME)
  end if

end program siro
