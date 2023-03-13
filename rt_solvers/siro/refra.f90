
!-----------
module refra
!-----------

  use parameters
  
  implicit none

  ! refraction parameters
  real(kind=sp) :: dhscale = 6.5_sp
  
contains
  
  !---------------------------------------------------------------
  subroutine calcrefra(rx,ry,rz,dx,dy,dz,dnewx,dnewy,dnewz,wlfact)
  !---------------------------------------------------------------
    !
    ! Refraction calculation

    implicit none
    
    ! input/output
    real(kind=sp), intent(in) :: rx,ry,rz
    real(kind=sp), intent(in) :: dx,dy,dz
    real(kind=sp), intent(in) :: wlfact
    real(kind=sp), intent(out) :: dnewx,dnewy,dnewz
    ! local
    real(kind=sp) r1,r2,z1,z2,n1,n2,n,rx12,ry12,rz12,rx012,ry012, &
         rz012,pistetulo,pituus,cosa,sinap,cosap,apu,dnorm

    if (userefrac) then
       
       r1 = sqrt(rx**2.0_sp + ry**2.0_sp + rz**2.0_sp)
       
       r2 = sqrt((rx+dx*step)**2.0_sp+(ry+dy*step)**2.0_sp+(rz+dz*step)**2.0_sp)
       
       z1 = r1 - req
       z2 = r2 - req
       
       n1 = wlfact*rindex(z1)
       n2 = wlfact*rindex(z2)
       n = 1.0_sp + 0.0_sp - n2 + n2**2.0_sp - n2**3.0_sp + n2**4.0_sp + n1 - n1*n2 + n1*n2**2.0_sp - n1*n2**3.0_sp
       
       rx12 = rx+dx*step/2.0_sp
       ry12 = ry+dy*step/2.0_sp
       rz12 = rz+dz*step/2.0_sp
       
       pistetulo = rx12*dx+ry12*dy+rz12*dz
       pituus = sqrt(rx12*rx12+ry12*ry12+rz12*rz12)
       
       if (pistetulo > 0.0_sp) then
          rx012 = rx12/pituus
          ry012 = ry12/pituus
          rz012 = rz12/pituus
       else
          rx012 = -1.0_sp*rx12/pituus
          ry012 = -1.0_sp*ry12/pituus
          rz012 = -1.0_sp*rz12/pituus
       end  if
       
       cosa = abs(pistetulo) / pituus
       if (cosa > 1.0_sp) then
          cosa = 1.0_sp
       end if
       
       sinap = n*sqrt(1.0_sp-cosa**2.0_sp)
       
       if (abs(sinap) > 1.0_sp) then
          sinap = sign(1.0_sp,sinap)
       end if
       
       cosap = sqrt(1.0_sp-sinap**2.0_sp)
       
       if (abs(cosap) > 1.0_sp) then
          cosap = sign(1.0_sp,cosap)
       end if
       
       if (cosa < epsilon1) then
          dnewx = dx
          dnewy = dy
          dnewz = dz
       elseif ((1.0_sp-sinap) < epsilon1) then
          apu = 1.0_sp / sqrt(1.0_sp-cosa**2.0_sp)
          dnewx = apu*(dx-cosa*rx012)
          dnewy = apu*(dy-cosa*ry012)
          dnewz = apu*(dz-cosa*rz012)
       else
          dnewx = (cosap-n*cosa)*rx012 + n*dx
          dnewy = (cosap-n*cosa)*ry012 + n*dy
          dnewz = (cosap-n*cosa)*rz012 + n*dz
       end if
       
       dnorm = sqrt(dnewx**2.0_sp + dnewy**2.0_sp + dnewz**2.0_sp)
       dnewx = dnewx/dnorm
       dnewy = dnewy/dnorm
       dnewz = dnewz/dnorm
       
    else
       
       dnewx = dx
       dnewy = dy
       dnewz = dz
       
    end if

  end subroutine calcrefra
  
  !------------------------------
  function rindex(z) result(rind)
  !------------------------------
    
    implicit none
    
    real(kind=sp), intent(in) :: z
    real(kind=sp) :: rind
    
    if (z < 0.0_sp) then
       rind = 1.0_sp
    elseif (z < 11.0_sp) then
       rind = (1.0_sp+0.0_sp+z*0.062519_sp)*exp(-z/dhscale)
    else
       rind = 1.6877_sp*exp(-z/dhscale)
    end if
    
  end function rindex
  
  !---------------------------------------------
  subroutine tangent(rx,ry,rz,px,py,pz,tx,ty,tz)
  !---------------------------------------------
    
    implicit none
    
    real(kind=sp) :: rx,ry,rz,px,py,pz,tx,ty,tz,ptulo,alpha,beta,gamma,delta,d1,d2,d3
    
    ptulo = px*rx+py*ry+pz*rz
    
    if (abs(ptulo) < epsilon1) then
       tx = px
       ty = py
       tz = pz
    else if (abs(rz) < epsilon1) then
       if (abs(ry) < epsilon1) then
          if (abs(py).lt.epsilon1) then
             tx = 0.0_sp
             ty = 0.0_sp
             tz = 1.0_sp
          else
             tx = 0.0_sp
             ty = 1.0_sp/sqrt(1.0_sp+(pz/py)**2.0_sp)
             tz = ty*pz/py
          end if
       else
          alpha=px*ry-py*rx
          if (abs(alpha) < epsilon1) then
             tx = 0.0_sp
             ty = 0.0_sp
             tz = 1.0_sp
          else
             beta = pz*ry+rx**2*pz/ry
             gamma = beta/alpha
             delta = gamma**2.0_sp+1.0_sp+(rx/ry)**2.0_sp
             tx = sqrt(1.0_sp/delta)
             ty = -1.0_sp*tx*rx/ry
             tz = tx*gamma
          end if
       end if
    else
       d1 = py*rz-pz*ry
       d2 = pz*rx-px*rz
       d3 = px*ry-py*rx
       alpha = d2-ry/rz*d3
       if (abs(alpha) < epsilon1) then
          tx = 0.0_sp
          ty = 1.0_sp/sqrt(1.0_sp+(ry/rz)**2.0_sp)
          tz = -1.0_sp*ry/rz*ty
       else
          alpha = -1.0_sp*(d1-rx/rz*d3)/alpha
          beta = 1.0_sp+(rx/rz)**2.0_sp+(1.0_sp+(ry/rz)**2.0_sp)*alpha**2.0_sp
          beta = beta+2.0_sp*rx*ry/rz**2.0_sp*alpha
          tx = 1.0_sp/sqrt(beta)
          ty = alpha*tx
          tz = -1.0_sp*(rx/rz*tx+ry/rz*ty)
       end if
    end if
    
    ptulo = tx*px+ty*py+tz*pz
    
    if (ptulo < 0.0_sp) then
       tx = -1.0_sp*tx
       ty = -1.0_sp*ty
       tz = -1.0_sp*tz
    end if
    
  end subroutine tangent
  
end module refra
