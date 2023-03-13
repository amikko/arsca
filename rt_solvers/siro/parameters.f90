
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
       maxtable       = 20000
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


  ! These parameters are now read from siro_settings.nml file
  integer :: &
       noph, &          !photon count for each wavelength
       atmos_layers, &  !internal atmosphere layer discretization
       maxnolay, &      !LOS discretization count
       maxnoord, &      !output file count
       mielength        !length of miefile

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

  ! mie files
  character(len=*), parameter :: miefile1 = 'input/miefiles/Continew500.dat'
  character(len=*), parameter :: miefile2 = 'input/miefiles/BgStrnew500.dat'
  character(len=*), parameter :: miefile1_altius = 'input/miefiles/aer_ALTIUS_ind01.dat'
  character(len=*), parameter :: miefile2_altius = 'input/miefiles/aer_ALTIUS_ind02.dat'
  character(len=*), parameter :: miefile3_altius = 'input/miefiles/aer_ALTIUS_ind03.dat'
  character(len=*), parameter :: miefile4_altius = 'input/miefiles/aer_ALTIUS_ind04.dat'
  character(len=*), parameter :: miefile5_altius = 'input/miefiles/aer_ALTIUS_ind05.dat'
  character(len=*), parameter :: miefile6_altius = 'input/miefiles/aer_ALTIUS_ind06.dat'
  character(len=*), parameter :: miefile7_altius = 'input/miefiles/aer_ALTIUS_ind07.dat'
  character(len=*), parameter :: miefile8_altius = 'input/miefiles/aer_ALTIUS_ind08.dat'
  character(len=*), parameter :: miefile9_altius = 'input/miefiles/aer_ALTIUS_ind09.dat'
  character(len=*), parameter :: miefile10_altius = 'input/miefiles/aer_ALTIUS_ind10.dat'
  character(len=*), parameter :: miefile11_altius = 'input/miefiles/aer_ALTIUS_ind11.dat'
  ! Hello. You're here probably to set up a wavelength-dependent scattering.
  ! Be sure to verify that you're loading the precisely same files as you're
  ! inputting. Hopefully you'll do better than me. If not, please increment this
  ! counter:
  ! Total work time wasted due to this file: 80 hours

  integer, parameter :: &
  climatology    = 6,& !This value is a development era artifact, remove this asap
  crossec = 4!This value is a development era artifact, remove this asap

  integer :: nosirowl, nosiroalt, nosir, noabs
  integer :: len_brf_zen_in, len_brf_zen_out, &
  len_brf_azi_out, len_brf_wavelengths

  real(kind=sp) :: albedo

  real(kind=sp), parameter :: &
  g = 0.1, & !This value is a development era artifact, remove this asap
  azimuth = 90.0, &!This value is a development era artifact, remove this asap
  zenith = 30.0, &!This value is a development era artifact, remove this asap
  satalt = 800.0

  ! wavelength and altitutde files
  character(len=256) :: altfile
  character(len=256) :: wlfile


end module parameters
