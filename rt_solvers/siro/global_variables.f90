
!----------------------
module global_variables
!----------------------

  use parameters

  implicit none

  real(kind=sp) :: satx = 0.0_sp, saty = 0.0_sp, satz = 0.0_sp
  real(kind=sp) :: sunx = 0.0_sp, suny = 0.0_sp, sunz = 0.0_sp
  real(kind=sp) :: detx = 0.0_sp, dety = 0.0_sp, detz = 0.0_sp

  real(kind=sp),dimension(o3crossize,11) :: o3cross
  real(kind=sp),dimension(no2crossize,11) :: no2cross
  real(kind=sp),dimension(crossize) :: aircross, aercross
  real(kind=sp),dimension(o3crossize) :: o3crosswl
  real(kind=sp),dimension(no2crossize) :: no2crosswl
  real(kind=sp),dimension(crossize) :: aircrosswl, aercrosswl

  real(kind=sp),dimension(:,:,:),allocatable :: sca_xsec
  real(kind=sp),dimension(:,:,:),allocatable :: abs_xsec

  ! added for polarization
  real(kind=sp), dimension(4), save :: inat
  real(kind=sp), dimension(maxtable), save :: cosmie, mie1p1,mie1p2,mie1p3,mie1p4,mie2p1,mie2p2,mie2p3,mie2p4
  real(kind=sp), dimension(maxtable,11), save :: cumtable1, cumtable2, phasetable1, phasetable2
  real(kind=sp), dimension(maxtable), save :: deltaphasef1, deltaphasef2, deltamy
  real(kind=sp), dimension(maxtable,1,16), save :: mie_table

  !                                 ^-- these 11's stand for each of the wavelengths

  real(kind=sp) :: atmos_layer_thickness
  integer, parameter :: current_wl_ind = 1
  integer :: brf_wl_idx
  real(kind=sp) :: current_wl
  logical :: direct_beam
  ! if brdf_reflection is 0, use lambertian. Otherwise use brdf

end module global_variables
