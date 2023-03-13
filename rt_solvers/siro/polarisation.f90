!------------------
module polarisation
!------------------

  use parameters
  use global_variables

  implicit none

  ! this is just some kind of epsilon value..
  real(kind=sp), parameter :: numacc = 5.0e-7_sp

  ! depolarisation of dry air, Young 1980
  ! parameter(depol=0.0279)
  ! for limb model comparisons depol=0
  real(kind=sp), parameter :: depol = 0.0_sp



contains

  !------------------------------
  subroutine maxmulti(r,s,matrix)
    !------------------------------

    implicit none

    ! input/output
    real(kind=sp), dimension(4,4), intent(in) :: r,s
    real(kind=sp), dimension(4,4), intent(out) :: matrix
    ! local
    integer :: i,j

    do i=1,4
       do j=1,4
          matrix(i,j) = r(i,1)*s(1,j)+r(i,2)*s(2,j)+r(i,3)*s(3,j)+r(i,4)*s(4,j)
       end do
    end do

  end subroutine maxmulti

  !-----------------------------------------------
  subroutine forwardrota(ox,oy,oz,nx,ny,nz,matrix)
    !-----------------------------------------------

    ! L.O. 11.12.1998
    ! Calculates rotation matrix for rotation of polarisation parameters from
    ! scatteringplane defined by o and n to the reference plane (x-y-plane)
    ! by clockwise rotation around -1*o
    ! Input:
    !   ox,y,z = old propagation direction (in backward terminology)
    !   nx,y,z = new propagation direction (in backward terminology)
    ! Output:
    !   matrix = rotation matrix

    implicit none

    ! input/output
    real(kind=sp), intent(in) :: ox,oy,oz,nx,ny,nz
    real(kind=sp), dimension(4,4), intent(out) :: matrix
    ! local
    real(kind=sp) :: ptulo,proj1,proj2,norma,cosfi,sinfi
    logical :: cond1,cond2,cond3


    ! special cases where rotation matrix=unity (rot. angle 0 or pi)
    ! o parallel to z-axis
    cond1 = abs(abs(oz)-1.0_sp) < numacc
    ! o and n or -n parallel
    cond2 = ((ox-nx)**2+(oy-ny)**2+(oz-nz)**2) < numacc
    cond3 = ((ox+nx)**2+(oy+ny)**2+(oz+nz)**2) < numacc
    if (cond1.or.cond2.or.cond3) then
       call unitmatrix(matrix)
    else
       ! calculate cosfi and sinfi
       ptulo = ox*nx+oy*ny+oz*nz
       proj1 = 1.0_sp/sqrt(1.0_sp-ptulo**2)
       proj2 = 1.0_sp/sqrt(ox**2+oy**2)
       norma = proj1*proj2
       cosfi = norma*(nx*oy-ny*ox)
       sinfi = norma*(-nz+oz*ptulo)
       call rotmatrix(cosfi,sinfi,matrix)
    end if

  end subroutine forwardrota

  !----------------------------
  subroutine unitmatrix(matrix)
  !----------------------------

    implicit none

    ! output
    real(kind=sp), dimension(4,4), intent(out) :: matrix
    ! local
    integer :: i,j

    do i=1,4
       do j=1,4
          if (i.eq.j) then
             matrix(i,j) = 1.0_sp
          else
             matrix(i,j) = 0.0_sp
          end if
       end do
    end do

  end subroutine unitmatrix

  !----------------------------
  subroutine IQUV2IIUV(matrix)
  !----------------------------

    implicit none

    ! output
    real(kind=sp), dimension(4,4), intent(out) :: matrix
    ! local
    integer :: i,j

    call unitmatrix(matrix)
    matrix(1,1) = 0.5_sp
    matrix(2,1) = 0.5_sp
    matrix(1,2) = 0.5_sp
    matrix(2,2) = -0.5_sp

  end subroutine IQUV2IIUV

  !---------------------------------------------
  subroutine phasematrix(process,costeta,matrix)
    !---------------------------------------------
    ! L.O. 8.12.1998
    ! Output: matrix=phase matrix of "process"
    ! Input:
    !  process = scattering process, 1=Rayleigh, 2-3=Aerosol, 10=reflection
    !  costeta = cosine of scattering (reflection) angle

    use global_variables, only : brdf_reflection
    implicit none

    integer, intent(in) :: process
    real(kind=sp) :: costeta,matrix(4,4),basis_change(4,4),tempmat(4,4)
    real(kind=sp) :: coeff,p1,p2,p3,p4,gamma
    real(kind=sp) :: p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16
    integer :: aero

    if (process == 1)  then
       ! Rayleigh scattering with anisotropy (gamma.ne.0),
       ! ks. Chandrasekhar p. 49

       ! depol is currently set to zero
       gamma = depol/(2.0_sp-depol)
       coeff = 3.0_sp/(8.0_sp*pi)*(1.0_sp/(1.0_sp+2.0_sp*gamma))

       matrix(1,1) = coeff*(costeta**2*(1.0_sp-gamma) + gamma)
       matrix(2,1) = gamma*coeff
       matrix(3,1) = 0.0_sp
       matrix(4,1) = 0.0_sp

       matrix(1,2) = gamma*coeff
       matrix(2,2) = coeff
       matrix(3,2) = 0.0_sp
       matrix(4,2) = 0.0_sp

       matrix(1,3) = 0.0_sp
       matrix(2,3) = 0.0_sp
       matrix(3,3) = coeff*costeta*(1.0_sp-gamma)
       matrix(4,3) = 0.0_sp

       matrix(1,4) = 0.0_sp
       matrix(2,4) = 0.0_sp
       matrix(3,4) = 0.0_sp
       matrix(4,4) = coeff*costeta*(1.0_sp-3.0_sp*gamma)

    elseif ((process == 2).or.(process == 3)) then
       ! Mie phase matrix
       !if (process == 2) then
      !    aero = 1
       !else
      !    aero = 2
       !end if

       !call readmietable(aero,costeta,p1,p2,p3,p4)

       !matrix(1,1)=p2
       !matrix(2,1)=0.0_sp
       !matrix(3,1)=0.0_sp
       !matrix(4,1)=0.0_sp

       !matrix(1,2)=0.0_sp
       !matrix(2,2)=p1
       !matrix(3,2)=0.0_sp
       !matrix(4,2)=0.0_sp

       !matrix(1,3)=0.0_sp
       !matrix(2,3)=0.0_sp
       !matrix(3,3)=p3
       !matrix(4,3)=-p4

       !matrix(1,4)=0.0_sp
       !matrix(2,4)=0.0_sp
       !matrix(3,4)=p4
       !matrix(4,4)=p3

       !Above is the original
       !Below is the ALTIUS wavelength-dependent scattering kernel

       call readmietable_altius(costeta,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, &
       p12,p13,p14,p15,p16)

       tempmat(1,1) = p1
       tempmat(1,2) = p2
       tempmat(1,3) = p3
       tempmat(1,4) = p4

       tempmat(2,1) = p5
       tempmat(2,2) = p6
       tempmat(2,3) = p7
       tempmat(2,4) = p8

       tempmat(3,1) = p9
       tempmat(3,2) = p10
       tempmat(3,3) = p11
       tempmat(3,4) = p12

       tempmat(4,1) = p13
       tempmat(4,2) = p14
       tempmat(4,3) = p15
       tempmat(4,4) = p16

       matrix = 1/(2*pi)*tempmat


    elseif (process == 10) then
       ! Reflection from ground

       matrix(1,1)=0.5_sp
       matrix(2,1)=0.5_sp
       matrix(3,1)=0.0_sp
       matrix(4,1)=0.0_sp

       matrix(1,2)=0.5_sp
       matrix(2,2)=0.5_sp
       matrix(3,2)=0.0_sp
       matrix(4,2)=0.0_sp

       matrix(1,3)=0.0_sp
       matrix(2,3)=0.0_sp
       matrix(3,3)=0.0_sp
       matrix(4,3)=0.0_sp

       matrix(1,4)=0.0_sp
       matrix(2,4)=0.0_sp
       matrix(3,4)=0.0_sp
       matrix(4,4)=0.0_sp
     else
       write(*,*) 'Error: scattering process out of range !!!'
    end if

  end subroutine phasematrix

  !------------------------------------------------
  subroutine readmietable(aero,costeta,p1,p2,p3,p4)
  !------------------------------------------------
    !L.O. 26.2.1999
    !Output: p1,p2,p3,p4 = elements of Mie phase matrix
    !                      (by linear interpolation from miefile)
    !Input: aero= aerosol "type" (1 or 2)
    !       costeta=cosine of scattering (reflection) angle

    use basic_functions, only : binarysearch

    implicit none

    integer, intent(in) :: aero
    real(kind=sp) :: costeta,p1,p2,p3,p4
    integer :: i1,i2
    real(kind=sp) :: d1,d2

    i1 = binarysearch(cosmie,mielength,costeta)
    i2 = i1+1

    d1 = (cosmie(i2)-costeta)/(cosmie(i2)-cosmie(i1))
    d2 = (costeta-cosmie(i1))/(cosmie(i2)-cosmie(i1))

    if (aero == 1) then
       p1 = mie1p1(i1)*d1+mie1p1(i2)*d2
       p2 = mie1p2(i1)*d1+mie1p2(i2)*d2
       p3 = mie1p3(i1)*d1+mie1p3(i2)*d2
       p4 = mie1p4(i1)*d1+mie1p4(i2)*d2
    elseif (aero == 2) then
       p1 = mie2p1(i1)*d1+mie2p1(i2)*d2
       p2 = mie2p2(i1)*d1+mie2p2(i2)*d2
       p3 = mie2p3(i1)*d1+mie2p3(i2)*d2
       p4 = mie2p4(i1)*d1+mie2p4(i2)*d2
    else
       write(*,*) 'error in Readmietable'
    end if


  end subroutine readmietable

  subroutine readmietable_altius(costeta,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11, &
    p12,p13,p14,p15,p16)
    !apologies for this abomination of a function, Fortran isn't my first
    !language, t. Antti
    use basic_functions, only : binarysearch

    implicit none

    real(kind=sp) :: costeta,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16
    integer :: i1,i2,wl_ind
    real(kind=sp) :: d1,d2

    wl_ind = current_wl_ind ! this a global variable, set in the siro main loop
    i1 = binarysearch(cosmie,mielength,costeta)
    i2 = i1+1
    d1 = (cosmie(i2)-costeta)/(cosmie(i2)-cosmie(i1))
    d2 = (costeta-cosmie(i1))/(cosmie(i2)-cosmie(i1))
    p1 = mie_table(i1,wl_ind,1) * d1 + mie_table(i2,wl_ind,1) * d2
    p2 = mie_table(i1,wl_ind,2) * d1 + mie_table(i2,wl_ind,2) * d2
    p3 = mie_table(i1,wl_ind,3) * d1 + mie_table(i2,wl_ind,3) * d2
    p4 = mie_table(i1,wl_ind,4) * d1 + mie_table(i2,wl_ind,4) * d2
    p5 = mie_table(i1,wl_ind,5) * d1 + mie_table(i2,wl_ind,5) * d2
    p6 = mie_table(i1,wl_ind,6) * d1 + mie_table(i2,wl_ind,6) * d2
    p7 = mie_table(i1,wl_ind,7) * d1 + mie_table(i2,wl_ind,7) * d2
    p8 = mie_table(i1,wl_ind,8) * d1 + mie_table(i2,wl_ind,8) * d2
    p9 = mie_table(i1,wl_ind,9) * d1 + mie_table(i2,wl_ind,9) * d2
    p10 = mie_table(i1,wl_ind,10) * d1 + mie_table(i2,wl_ind,10) * d2
    p11 = mie_table(i1,wl_ind,11) * d1 + mie_table(i2,wl_ind,11) * d2
    p12 = mie_table(i1,wl_ind,12) * d1 + mie_table(i2,wl_ind,12) * d2
    p13 = mie_table(i1,wl_ind,13) * d1 + mie_table(i2,wl_ind,13) * d2
    p14 = mie_table(i1,wl_ind,14) * d1 + mie_table(i2,wl_ind,14) * d2
    p15 = mie_table(i1,wl_ind,15) * d1 + mie_table(i2,wl_ind,15) * d2
    p16 = mie_table(i1,wl_ind,16) * d1 + mie_table(i2,wl_ind,16) * d2


  end subroutine readmietable_altius

  !---------------------------------------
  subroutine rotmatrix(cosfi,sinfi,matrix)
    !---------------------------------------

    ! L.O. 8.12.1998
    ! Output: matrix=rotation matrix
    ! Input:
    !    cosfi,sinfi = cos(rotation angle), sin(rotation angle)
    !   cos2fi,sin2fi = cos(2*rotation angle), sin(2*rotation angle)

    implicit none

    real(kind=sp) :: cosfi,sinfi,cos2fi,sin2fi,matrix(4,4)
    real(kind=sp) :: cosfisq,sinfisq

    cosfisq = cosfi**2
    sinfisq = sinfi**2
    cos2fi = 2.0_sp*cosfi**2-1.0_sp
    sin2fi = 2.0_sp*sinfi*cosfi

    matrix(1,1) = cosfisq
    matrix(2,1) = sinfisq
    matrix(3,1) = -1.0_sp*sin2fi
    matrix(4,1) = 0.0_sp

    matrix(1,2) = sinfisq
    matrix(2,2) = cosfisq
    matrix(3,2) = sin2fi
    matrix(4,2) = 0.0_sp

    matrix(1,3) = 0.5_sp*sin2fi
    matrix(2,3) = -0.5_sp*sin2fi
    matrix(3,3) = cos2fi
    matrix(4,3) = 0.0_sp

    matrix(1,4) = 0.0_sp
    matrix(2,4) = 0.0_sp
    matrix(3,4) = 0.0_sp
    matrix(4,4) = 1.0_sp

  end subroutine rotmatrix

  !-----------------------------------------
  subroutine maxvecmulti(matrix,vec,outcome)
  !-----------------------------------------

    ! L.O. 11.12.1998
    ! Input:  a matrix and vector vec
    ! Output: outcome=matrix*vec

    implicit none

    real(kind=sp), dimension(4,4),intent(in) :: matrix
    real(kind=sp), dimension(4),intent(in) :: vec
    real(kind=sp), dimension(4),intent(out) :: outcome
    integer :: i

    do i=1,4
       outcome(i) = matrix(i,1)*vec(1)+matrix(i,2)*vec(2)+ &
            matrix(i,3)*vec(3)+matrix(i,4)*vec(4)
    end do

  end subroutine maxvecmulti


  !--------------------------------------------------------------------------------
  subroutine polaris(nx,ny,nz,ox,oy,oz,costeta,cosalpha,sinalpha,process,cummatrix)
    !--------------------------------------------------------------------------------
    ! L.O. 17.12.1998.
    ! Updates polarisation parameters (cummatric and pweight) of a photon
    ! Input:  n,o = new and old propagation direction (backwards)
    !        teta = scattering angle
    !        alpha = rotation after scattering in forward terminology
    !        cummatrix = product of rotation and phase matrixes of all
    !        former collisions
    ! Output: cummatrix = cummatrix*L(alpha)*S(teta)*L(beta) (beta=rotation
    !        before coollision in forward terminology, L=rotation, S=phasem.)

    !        pweight=pweight*1/P(teta), P=phase function of current collision
    ! Version 8.11.99

    use basic_functions, only : phasef

    implicit none

    !input/output
    real(kind=sp), intent(in) :: nx,ny,nz,ox,oy,oz,costeta,cosalpha,sinalpha
    integer, intent(in) :: process
    real(kind=sp), intent(inout) :: cummatrix (4,4)
    ! local
    !real(kind=sp), parameter :: epsilon = 1.0e-20_sp ! check this?!
    ! Antti @ 2.9.2022: Checked the epsilon and it requires some tuning.
    ! In this case, it probably should be about 0.1 degree?
    real(kind=sp), parameter :: epsilon = 0.1 / 180.0_sp * pi
    logical :: cond1,cond2
    real(kind=sp), dimension(4,4) :: phasem,rot1,rot2,matrix,matrix1,oldcm
    real(kind=sp) :: cosbeta,sinbeta,nim,pweight,ph,norm
    integer :: i,j

    call phasematrix(process,costeta,phasem)
    do i=1,4
       do j=1,4
          oldcm(i,j)=cummatrix(i,j)
       end do
    end do
    ! exceptions
    ! n=+-o, rotation angles zero
    cond1 = ((ox-nx)**2+(oy-ny)**2+(oz-nz)**2) < epsilon**2 ! old = new; sirotaan etusuuntaan
    cond2 = ((ox+nx)**2+(oy+ny)**2+(oz+nz)**2) < epsilon**2 ! old = -new; sirotaan takasuuntaan


    if (cond1.or.cond2) then

       call Maxmulti(cummatrix,phasem,matrix)
    else
       if (abs(abs(oz)-1.0_sp) < epsilon) then
          ! o parallel to +-z-axis (beta=+-pi/2, rot matrix symm to shift of pi)
          cosbeta = 0.0_sp
          sinbeta = 1.0_sp
          ! n parallel to +z-axis
       else if (abs(nz-1.0_sp) < epsilon) then
          norm = sqrt(1.0_sp-oz**2) ! Mitä tässä tapahtuu t. Antti
          cosbeta = -1.0_sp*ox/norm
          sinbeta = oy/norm
          !     n parallel to -z-axis
       else if (abs(nz+1.0_sp) < epsilon) then
          norm = sqrt(1.0_sp-oz**2)
          cosbeta = -1.0_sp*ox/norm
          sinbeta = -1.0_sp*oy/norm
       else ! on jotain muutakin kuin zetaa oldissa ja newssa
          nim = ox*ny-oy*nx
          if (abs(nim) < epsilon) then !
             cosbeta = 0.0_sp
             sinbeta = 1.0_sp   !   vai pitaisiko olla -1 ?
          else
             cosbeta = -1.0_sp*cosalpha*sqrt(ox**2+oy**2)/sqrt(nx**2+ny**2)
             sinbeta = cosbeta*(oz-nz*costeta)/nim
          end if
       end if

       call rotmatrix(cosalpha,sinalpha,rot1)
       call rotmatrix(cosbeta,sinbeta,rot2)

       call maxmulti(phasem,rot2,matrix) ! M = PR_2

       call maxmulti(rot1,matrix,matrix1) ! M1 = R_1PR_2

       call maxmulti(cummatrix,matrix1,matrix)

    end if

    !write (*,*) cummatrix
    ! "if" below  just for safety, polaris should be called only for process<=nosir

    pweight = 0.0_sp ! simo added

    if (process <= nosir) then
       ph = phasef(process,costeta)

       if (abs(ph-0.0_sp) > epsilon) then
          pweight = 1.0_sp/ph
       end if
    else ! process == 10, i.e. surface reflection
         !
      pweight = 1.0_sp
    end if
    do i=1,4
       do j=1,4
          cummatrix(i,j)=matrix(i,j)*pweight
       end do
    end do
    !if (cummatrix(1,1) > 100) then
      ! write (*,*) 'abnormally high cummatrix detected'
      ! write (*,*) 'process'
      ! write (*,*) process
      ! write (*,*) 'conds'
      ! write (*,*) cond1
      ! write (*,*) cond2
      ! write (*,*) 'phasem'
      ! write (*,*) phasem
      ! write (*,*) 'old cummatrix'
      ! write (*,*) oldcm
      ! write (*,*) 'new cummatrix'
      ! write (*,*) cummatrix
      ! write (*,*) 'rot1'
      ! write (*,*) rot1
      ! write (*,*) 'rot2'
      ! write (*,*) rot2
      ! write (*,*) 'angles'
      ! write (*,*) cosalpha
      ! write (*,*) sinalpha
      ! write (*,*) cosbeta
      ! write (*,*) sinbeta
      ! write (*,*) 'direction'
      ! write (*,*) ox
      ! write (*,*) oy
      ! write (*,*) oz
      ! write (*,*) nx
      ! write (*,*) ny
      ! write (*,*) nz
      ! write (*,*) 'scaling'
      ! write (*,*) phasef(process,costeta)
      !stop
    !end if
    !write(*,*) cummatrix
    !write(*,*) '----------------'
  end subroutine polaris



end module polarisation
