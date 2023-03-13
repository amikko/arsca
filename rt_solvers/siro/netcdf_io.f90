module netcdf_io

  !This module contains routines for data I/O via NetCDF4 files.
  use parameters
  use netcdf
  use basic_functions
  implicit none
  character (len = 128) :: INPUT_FILE_NAME
  character (len = 128) :: OUTPUT_FILE_NAME

  !group names
  character (len = *), parameter :: MEDI_GRP_NAME = "medium"
  character (len = *), parameter :: INST_GRP_NAME = "instrument"
  character (len = *), parameter :: SOUR_GRP_NAME = "source"
  character (len = *), parameter :: BOUN_GRP_NAME = "boundary"

  !dimension names
  character (len = *), parameter :: MPOS_DIM_NAME = "n_medium_position"
  character (len = *), parameter :: WL_DIM_NAME = "n_input_wavelength"
  character (len = *), parameter :: ABS_DIM_NAME = "n_absorber"
  character (len = *), parameter :: SCA_DIM_NAME = "n_scatterer"
  character (len = *), parameter :: GEOM_DIM_NAME = "n_los"
  character (len = *), parameter :: SRC_DIM_NAME = "n_source"
  character (len = *), parameter :: COORD_DIM_NAME = "n_coordinate"
  character (len = *), parameter :: STOKES_DIM_NAME = "n_stokes"

  !atmospheric profile variables
  character (len = *), parameter :: MPOS_VAR_NAME = "position"
  character (len = *), parameter :: ABS_VAR_NAME = "absorber"
  character (len = *), parameter :: SCA_VAR_NAME = "scatterer"
  character (len = *), parameter :: AXSEC_VAR_NAME = "absorbing_cross_section"
  character (len = *), parameter :: SXSEC_VAR_NAME = "scattering_cross_section"

  !instrument variables
  character (len = *), parameter :: INSTPOS_VAR_NAME = "position"
  character (len = *), parameter :: VIEWVEC_VAR_NAME = "view_vector"

  !light source variables
  character (len = *), parameter :: WL_VAR_NAME = "input_wavelength"
  character (len = *), parameter :: INCDIR_VAR_NAME = "incident_direction"
  character (len = *), parameter :: INCSTOKES_VAR_NAME = "incident_stokes"
  character (len = *), parameter :: SOURANGR_VAR_NAME = "source_angular_radius"

  !scattering variables
  character (len = *), parameter :: REFKERPAR_VAR_NAME = "reflection_kernel_parameter"

  integer, parameter :: NDIMS = 8 !parameter on vain käännösajan vakio; C define

  !The output file variables
  character (len = *), parameter :: RADIANCE_VAR_NAME = "radiance"
  character (len = *), parameter :: STATUS_VAR_NAME = "status"

  ! These are for BRDF input, t. Antti on 21.5.2020
  character (len = *), parameter :: BRF_VAR_ZEN_IN = "zen_inc"
  character (len = *), parameter :: BRF_VAR_WL = "wavelength"
  character (len = *), parameter :: BRF_VAR_ZEN_OUT = "zen_ref"
  character (len = *), parameter :: BRF_VAR_AZI_OUT = "azi_ref"
  character (len = *), parameter :: BRF_VAR_MULLER = "muller"
  character (len = *), parameter :: BRF_DIM_ZEN_IN = "zenith_inc"
  character (len = *), parameter :: BRF_DIM_WL = "wavelengths"
  character (len = *), parameter :: BRF_DIM_ZEN_OUT = "zenith_ref"
  character (len = *), parameter :: BRF_DIM_AZI_OUT = "azimuth_ref"



contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop 2
    end if
  end subroutine check

  subroutine nc_create_atmosphere(atmosalt, sca_prof, abs_prof, pres_prof, temp_prof, sca_xsec, abs_xsec)
    implicit none
    !the input/output
    real(kind=sp), dimension(:), intent(in) :: atmosalt
    real(kind=sp), dimension(:), intent(out) :: pres_prof, temp_prof
    real(kind=sp), dimension(:,:), intent(out) :: sca_prof, abs_prof
    real(kind=sp), dimension(:,:,:), intent(out) :: sca_xsec, abs_xsec
    !the first dimension is the atmosalts and the second one is the nosir/noabs

    !local
    integer :: ncid, prof_groupid, alti_varid, pres_varid, temp_varid, sca_varid, abs_varid, sxsec_varid, axsec_varid
    integer :: altsize_in
    real(kind=sp), dimension(:), allocatable :: alt_in, pres_in, temp_in
    real(kind=sp), dimension(:,:), allocatable :: sca_prof_in, abs_prof_in, alt_in_2d
    real(kind=sp), dimension(:,:,:), allocatable :: sca_xsec_in, abs_xsec_in
    integer :: dimids(NDIMS)
    integer :: i, j, k
    character (len = nf90_max_name) :: dimname

    !reserve the memory for input data
    call get_dimids(dimids)
    call check( nf90_open(INPUT_FILE_NAME, NF90_NOWRITE, ncid) )
    call check ( nf90_inquire_dimension(ncid, dimids(1), dimname, altsize_in) )
    allocate(alt_in_2d(altsize_in, 3))
    allocate(sca_prof_in(altsize_in, nosir))
    allocate(abs_prof_in(altsize_in, noabs))
    allocate(sca_xsec_in(altsize_in, nosirowl, nosir))
    allocate(abs_xsec_in(altsize_in, nosirowl, noabs))

    !obtain the varids of all the profiles
    call check( nf90_inq_ncid(ncid, MEDI_GRP_NAME, prof_groupid) )
    call check( nf90_inq_varid(prof_groupid, MPOS_VAR_NAME, alti_varid) )
    call check( nf90_inq_varid(prof_groupid, SCA_VAR_NAME, sca_varid) )
    call check( nf90_inq_varid(prof_groupid, ABS_VAR_NAME, abs_varid) )
    call check( nf90_inq_varid(prof_groupid, SXSEC_VAR_NAME, sxsec_varid) )
    call check( nf90_inq_varid(prof_groupid, AXSEC_VAR_NAME, axsec_varid) )

    !call check( nf90_get_var(prof_groupid, alti_varid, alt_in_2d) )
    call read_2d_data_from_variable(alt_in_2d, (/altsize_in, 3/), prof_groupid, alti_varid)
    call read_2d_data_from_variable(sca_prof_in, (/altsize_in, nosir/), prof_groupid, sca_varid)
    call read_2d_data_from_variable(abs_prof_in, (/altsize_in, noabs/), prof_groupid, abs_varid)
    call read_3d_data_from_variable(sca_xsec_in, (/altsize_in, nosirowl, nosir/), prof_groupid, sxsec_varid)
    call read_3d_data_from_variable(abs_xsec_in, (/altsize_in, nosirowl, noabs/), prof_groupid, axsec_varid)
    call check( nf90_close(ncid) )
    allocate(alt_in(altsize_in))
    alt_in(:) = alt_in_2d(:,1)


    do i = 1,altsize_in
      !The radius of Earth is taken into account elsewhere
      alt_in(i) = alt_in(i) - req
    end do

    !write(*,*) alt_in

    do i = 1,atmos_layers
      do j = 1,nosir
        sca_prof(i,j) = interpolate(alt_in,sca_prof_in(:,j),altsize_in,atmosalt(i))
        do k = 1,nosirowl
          sca_xsec(i,k,j) = interpolate(alt_in,sca_xsec_in(:,k,j),altsize_in,atmosalt(i))
        end do
      end do
      do j = 1,noabs
        abs_prof(i,j) = interpolate(alt_in,abs_prof_in(:,j),altsize_in,atmosalt(i))
        do k = 1,nosirowl
          abs_xsec(i,k,j) = interpolate(alt_in,abs_xsec_in(:,k,j),altsize_in,atmosalt(i))
        end do
      end do
    end do

  end subroutine nc_create_atmosphere

  subroutine nc_get_geometries(solar_dir, sat_position, sat_direction)
    implicit none
    !the input/output
    real(kind=sp), dimension(:,:), intent(out) :: solar_dir, sat_position, sat_direction

    !local
    real(kind=sp), dimension(1,3) :: temp_solar_dir
    integer :: ncid, sour_groupid, inst_groupid, satpos_varid, viewv_varid, incdir_varid
    integer :: i, j

    !reserve the memory for input data
    call check( nf90_open(INPUT_FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_ncid(ncid, INST_GRP_NAME, inst_groupid) )
    call check( nf90_inq_ncid(ncid, SOUR_GRP_NAME, sour_groupid) )

    !obtain the varids of all the profiles
    call check( nf90_inq_varid(inst_groupid, INSTPOS_VAR_NAME, satpos_varid) )
    call check( nf90_inq_varid(inst_groupid, VIEWVEC_VAR_NAME, viewv_varid) )
    call check( nf90_inq_varid(sour_groupid, INCDIR_VAR_NAME, incdir_varid) )

    call read_2d_data_from_variable(sat_position,(/nosiroalt,3/),inst_groupid,satpos_varid)
    call read_2d_data_from_variable(sat_direction,(/nosiroalt,3/),inst_groupid,viewv_varid)
    !call read_2d_data_from_variable(solar_dir,(/nosiroalt,3/),sour_groupid,incdir_varid)
    call read_2d_data_from_variable(temp_solar_dir,(/1,3/),sour_groupid,incdir_varid) !the first solar dir workaround
    write(*,*) 'Using only the first solar direction!'
    ! Antti of 21.09.2020: Why is this above situation here?
    call check( nf90_close(ncid) )

    !first solar dir workaround begins here
    do i=1,nosiroalt
      do j=1,3
        solar_dir(i,j) = temp_solar_dir(1,j)
      end do
    end do
    !first solar dir workaround ends here

    do i=1,nosiroalt
      do j=1,3
        solar_dir(i,j) = - solar_dir(i,j)
        !incident direction isn't solar dir exactly, but completely opposite!
      end do
    end do

  end subroutine nc_get_geometries

  subroutine nc_get_lengths()
    !This function obtains the various array sizes from the simulation
    !configuration NetCDF file.

    !Sizes obtained in this function: nosirowl, nosiroalt, nosir, noabs
    character (len = nf90_max_name) :: dimname
    integer :: ncid
    integer :: wl_datasize, geom_datasize, abs_datasize, sca_datasize
    integer :: dimids(NDIMS)
    call get_dimids(dimids)
    call check( nf90_open(INPUT_FILE_NAME, NF90_NOWRITE, ncid) )
    !write (*,*) dimids
    call check ( nf90_inquire_dimension(ncid, dimids(2), dimname, wl_datasize) )
    call check ( nf90_inquire_dimension(ncid, dimids(5), dimname, geom_datasize) )
    call check ( nf90_inquire_dimension(ncid, dimids(3), dimname, abs_datasize) )
    call check ( nf90_inquire_dimension(ncid, dimids(4), dimname, sca_datasize) )
    nosirowl = wl_datasize
    nosiroalt = geom_datasize !even though this isn't a tanget altitude anymore
    !the old variable name is still used for compatibility
    noabs = abs_datasize
    nosir = sca_datasize
    !write (*,*) wl_datasize, geom_datasize, abs_datasize, sca_datasize
    call check( nf90_close(ncid) )
  end subroutine nc_get_lengths

  subroutine nc_read_albedo()
    use parameters
    implicit none
    integer :: ncid, phase_groupid, refl_kernel_param_varid
    call check( nf90_open(INPUT_FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_ncid(ncid, BOUN_GRP_NAME, phase_groupid) )
    call check( nf90_inq_varid(phase_groupid, REFKERPAR_VAR_NAME, refl_kernel_param_varid) )
    call check( nf90_get_var(phase_groupid, refl_kernel_param_varid, albedo) )
    call check( nf90_close(ncid) )

  end subroutine nc_read_albedo

  !------------------------------------
  subroutine nc_read_waves(wl,wlfactor,incident_stokes,source_angular_radius)
  !------------------------------------

  ! Read wavelenghts

    implicit none

    ! input/output
    real(kind=sp),dimension(:),intent(out) :: wl,wlfactor
    real(kind=sp),dimension(:,:), intent(out) :: incident_stokes
    real(kind=sp), intent(out) :: source_angular_radius

    ! local
    integer :: n, sour_groupid, wl_varid, incstokes_varid, source_radi_varid, ncid

    call check( nf90_open(INPUT_FILE_NAME, NF90_NOWRITE, ncid) )
    !write(*,*) 'opened'
    call check( nf90_inq_ncid(ncid, SOUR_GRP_NAME, sour_groupid) )
    call check( nf90_inq_varid(sour_groupid, WL_VAR_NAME, wl_varid) )
    call check( nf90_inq_varid(sour_groupid, INCSTOKES_VAR_NAME, incstokes_varid) )
    call check( nf90_inq_varid(sour_groupid, SOURANGR_VAR_NAME, source_radi_varid) )
    !write(*,*) 'varid inq ok'
    call read_1d_data_from_variable(wl, nosirowl, sour_groupid, wl_varid)
    !write(*,*) wl
    call read_2d_data_from_variable(incident_stokes, (/nosirowl, 4/), sour_groupid, incstokes_varid)
    !write(*,*) 'incstokes ok'
    call check( nf90_get_var(sour_groupid, source_radi_varid, source_angular_radius) )
    !write(*,*) 'vars got'
    call check( nf90_close(ncid) )
    ! calculate wlfactor
    do n=1,nosirowl

       wlfactor(n) = 1.0e-6_sp * (64.328_sp + 29498.1_sp / (146.0_sp-1e+6_sp/(wl(n)*wl(n))) &
            + 255.4_sp / (41.0_sp-1e+6_sp/(wl(n)*wl(n))) )
    end do

  end subroutine nc_read_waves

  !subroutine save_radiance_data(unscaw,unscavec,I,sI,sQ,sU,sV)
    !TODO: netcdf output

    !Captain's log, stardate 13.06.19: Even though netcdf output would be very useful,
    !it is not sufficient in today's hectic world of satellite proposals.

  !  implicit none
    !integer, intent(in) !toivottavasti tämä ei ollut mikään liian tärkee; en muista yhtään
  !  real(kind=sp), intent(in), dimension(maxnoord,nosirowl,nosiroalt) :: I,sI,sQ,sU,sV
  !  integer :: ncid, radiance_varid
  !  call nf90_open(OUTPUT_FILE_NAME, NF90_WRITE, ncid)             ! open existing netCDF dataset

  !  call nf90_inq_varid(ncid, RADIANCE_VAR_NAME, radiance_varid)      ! get variable IDs

  !  call nf90_put_var()        ! provide new values for variables, if any
  !  call nf90_close(ncid)
  !end subroutine save_radiance_data

  !subroutine save_description_data()

  !end subroutine save_description_data

  subroutine get_dimids(dimids)
    implicit none
    integer, dimension(NDIMS), intent (inout) :: dimids
    integer :: ncid
    call check( nf90_open(INPUT_FILE_NAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_dimid(ncid, MPOS_DIM_NAME, dimids(1)) )
    call check( nf90_inq_dimid(ncid, WL_DIM_NAME, dimids(2)) )
    call check( nf90_inq_dimid(ncid, ABS_DIM_NAME, dimids(3)) )
    call check( nf90_inq_dimid(ncid, SCA_DIM_NAME, dimids(4)) )
    call check( nf90_inq_dimid(ncid, GEOM_DIM_NAME, dimids(5)) )
    call check( nf90_inq_dimid(ncid, COORD_DIM_NAME, dimids(6)) )
    call check( nf90_inq_dimid(ncid, STOKES_DIM_NAME, dimids(7)) )
    call check( nf90_inq_dimid(ncid, SRC_DIM_NAME, dimids(8)) )
    call check( nf90_close(ncid) )
  end subroutine get_dimids
  !The three following subroutines are not beautiful, efficient or at least
  !civilized. However, they were the only way which I could devise to circumvent
  !the weird behaviour of netcdf4 in fortran for inverting the matrix dimensions
  !in reading the file. You can notice this in starting the read at (j,i) and (k,j,i)
  !instead the correct way around. Be wary of these sort of phenomenon in the
  !future and may Erdős have mercy on our souls.
  !
  !                                 Sincerely,
  !                                        Antti
    subroutine read_1d_data_from_variable(data,data_size,ncid,varid)
      implicit none
      real(kind=sp),dimension(:),intent(out) :: data
      integer, intent(in) :: data_size
      integer, intent(in) :: ncid,varid
      integer :: i
      real(kind=sp),dimension(1) :: temp_data
      do i = 1,data_size
        call check( nf90_get_var(ncid, varid, temp_data, (/i/), (/1/)))
        data(i) = temp_data(1)
      end do
    end subroutine read_1d_data_from_variable

    subroutine read_2d_data_from_variable(data,data_size,ncid,varid)
      implicit none
      real(kind=sp),dimension(:,:),intent(out) :: data
      integer,dimension(2), intent(in) :: data_size
      integer, intent(in) :: ncid,varid
      integer :: i,j
      real(kind=sp),dimension(1,1) :: temp_data
      do i = 1,data_size(1)
        do j = 1,data_size(2)
          call check( nf90_get_var(ncid, varid, temp_data, (/j,i/), (/1,1/)))
          data(i,j) = temp_data(1,1)
        end do
      end do
    end subroutine read_2d_data_from_variable

    subroutine read_3d_data_from_variable(data,data_size,ncid,varid)
      implicit none
      real(kind=sp),dimension(:,:,:),intent(out) :: data
      integer,dimension(3), intent(in) :: data_size
      integer, intent(in) :: ncid,varid
      integer :: i,j,k
      real(kind=sp),dimension(1,1,1) :: temp_data
      do i = 1,data_size(1)
        do j = 1,data_size(2)
          do k = 1,data_size(3)
            call check( nf90_get_var(ncid, varid, temp_data, (/k,j,i/), (/1,1,1/)))
            data(i,j,k) = temp_data(1,1,1)
          end do
        end do
      end do
    end subroutine read_3d_data_from_variable

! BRDF reading from a file.
! Antti, 21.5.2020
! This is hardcoded and bad practice.

subroutine nc_get_brdf(brf_zen_in, brf_wavelengths, brf_zen_out, brf_azi_out, brf_M)
  implicit none
  real(kind=sp),dimension(:),intent(out) :: brf_zen_in
  real(kind=sp),dimension(:),intent(out) :: brf_wavelengths
  real(kind=sp),dimension(:),intent(out) :: brf_zen_out
  real(kind=sp),dimension(:),intent(out) :: brf_azi_out
  real(kind=sp),dimension(:,:,:,:,:,:),intent(out) :: brf_M
  integer :: i,j,k,l,m,n
  integer :: ncid, n_zen_in, n_wl, n_zen_out, n_azi_out
  integer :: varid_zen_in, varid_wl, varid_zen_out, varid_azi_out, varid_M
  real(kind=sp),dimension(1,1,1,1,1,1) :: temp_data

  ! these are also hardcoed in the beginning of siro.f90
  n_zen_in = len_brf_zen_in
  n_wl = len_brf_wavelengths
  n_zen_out = len_brf_zen_out
  n_azi_out = len_brf_azi_out
  call check( nf90_open(BRF_FILENAME, NF90_NOWRITE, ncid) )

  call check( nf90_inq_varid(ncid, BRF_VAR_ZEN_IN, varid_zen_in) )
  call check( nf90_inq_varid(ncid, BRF_VAR_WL, varid_wl) )
  call check( nf90_inq_varid(ncid, BRF_VAR_ZEN_OUT, varid_zen_out) )
  call check( nf90_inq_varid(ncid, BRF_VAR_AZI_OUT, varid_azi_out) )
  call check( nf90_inq_varid(ncid, BRF_VAR_MULLER, varid_M) )

  call read_1d_data_from_variable(brf_zen_in,n_zen_in,ncid,varid_zen_in)
  call read_1d_data_from_variable(brf_wavelengths,n_wl,ncid,varid_wl)
  call read_1d_data_from_variable(brf_zen_out,n_zen_out,ncid,varid_zen_out)
  call read_1d_data_from_variable(brf_azi_out,n_azi_out,ncid,varid_azi_out)
  ! then, the monstrosity...
  write (*,*) 'Reading the BRDF (this might take a while)...'
  do i = 1,n_zen_in
    write(*,*) 'i:',i
    do j = 1,n_wl
      !write(*,*) 'j:',j
      do k = 1,n_zen_out
        !write(*,*) 'k:',k
        do l = 1,n_azi_out
          !write(*,*) 'l:',l
          do m = 1,4
            !write(*,*) 'm:',m
            do n = 1,4
              !write(*,*) 'n:',n

              call check( nf90_get_var(ncid, varid_M, temp_data, (/n,m,l,k,j,i/), (/1,1,1,1,1,1/)))
              brf_M(i,j,k,l,m,n) = temp_data(1,1,1,1,1,1)
            end do
          end do
        end do
      end do
    end do
  end do
  call check( nf90_close(ncid) )
end subroutine nc_get_brdf

  subroutine nc_get_brdf_lengths()
    !This function obtains the various array sizes from BRDF file

    character (len = nf90_max_name) :: dimname
    integer :: ncid
    integer :: zenin_datasize, wl_datasize, zenout_datasize, aziout_datasize
    integer :: zenin_dimid, wl_dimid, zenout_dimid, aziout_dimid

    call check( nf90_open(BRF_FILENAME, NF90_NOWRITE, ncid) )
    call check( nf90_inq_dimid(ncid, BRF_DIM_ZEN_IN, zenin_dimid) )
    call check( nf90_inq_dimid(ncid, BRF_DIM_WL, wl_dimid) )
    call check( nf90_inq_dimid(ncid, BRF_DIM_ZEN_OUT, zenout_dimid) )
    call check( nf90_inq_dimid(ncid, BRF_DIM_AZI_OUT, aziout_dimid) )
    call check ( nf90_inquire_dimension(ncid, zenin_dimid, dimname, zenin_datasize) )
    call check ( nf90_inquire_dimension(ncid, wl_dimid, dimname, wl_datasize) )
    call check ( nf90_inquire_dimension(ncid, zenout_dimid, dimname, zenout_datasize) )
    call check ( nf90_inquire_dimension(ncid, aziout_dimid, dimname, aziout_datasize) )
    len_brf_zen_in = zenin_datasize
    len_brf_wavelengths = wl_datasize !even though this isn't a tanget altitude anymore
    !the old variable name is still used for compatibility
    len_brf_zen_out = zenout_datasize
    len_brf_azi_out = aziout_datasize
    !write (*,*) wl_datasize, geom_datasize, abs_datasize, sca_datasize
    call check( nf90_close(ncid) )
    !write(*,*) len_brf_zen_in, len_brf_wavelengths, len_brf_zen_out, len_brf_azi_out
    !stop
  end subroutine nc_get_brdf_lengths

end module netcdf_io
