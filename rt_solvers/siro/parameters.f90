
!----------------
module parameters
!----------------

  implicit none

  logical :: MCERRORS = .false.

  integer, parameter :: &
       sp = selected_real_kind(12,50)         ! accuracy of floats

  real(kind=sp), parameter :: &                ! siro version
       version = 3.0_sp

  ! do not change
  integer, parameter :: &
       !crossize       = 5201, & !original value !test values
       !o3crossize     = 8944, & !originally commented !test values
       !no2crossize    = 8953, & !originally commented
       !o3crossize     = 5201, & !original value
       !no2crossize    = 5201, & !original value
       crossize       = 11, & ! ALTIUS
       o3crossize     = 11, & ! ALTIUS
       no2crossize    = 11, & ! this is actually the aerosol absorption cross-section for ALTIUS
       maxtable       = 80000
       !mielength      = 210                     ! length of mie file
       !mielength      = 5400
       !mielength      = 251
       !mielength      = 51 ! This is for raysca validation


  ! do not change
  real(kind=sp), parameter :: &
       pi             = 3.141592653589793238_sp, &   ! number of phi
       epsilon1       = 1.0e-20_sp, &
       refangle       = 2.617993877991494e-2_sp   ! refraction angle
       !selvitä mihin refangle perustuu!


  ! Most of these parameters are now read from siro_settings.nml file
  integer :: &
       noph, &          !photon count for each wavelength
       atmos_layers, &  !internal atmosphere layer discretization
       maxnolay, &      !LOS discretization count
       maxnoord      !output file count

  integer, dimension(2) :: mielength !lengths of the aerosol files

  real(kind=sp) :: &
       minweight, &     !minimum weight of a photon before elimination
       req, &           !radius of the ground surface
       ratm, &          !radius of the top of the atmosphere
       step             !integration step for scattering events

  logical, save :: &
       userefrac, &
       usepolar, &
       brdf_reflection

  character(len=256) :: BRF_FILENAME
  character(len=256) :: AER_FILENAME
  character(len=256) :: AER2_FILENAME

  integer :: nosirowl, nosiroalt, nosir, noabs
  integer :: len_brf_zen_in, len_brf_zen_out, &
  len_brf_azi_out, len_brf_wavelengths

  real(kind=sp) :: albedo

end module parameters
