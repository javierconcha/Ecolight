C     Last change:  LKS   7 Aug 2008    6:54 pm
      SUBROUTINE ztozeta
C 
C     core routine on file ztozeta.f
C
C     called by INISHRAD [MAIN->RADIANCE->INISHRAD->ZTOZETA]
C
C     calls routines BEAMCZ (in this file), QAG (in file quadrat.f)
C
C     This routine first computes the OPTICAL depths zeta 
c     (nondimensional) that correspond to the GEOMETRICAL depths z 
c     (in meters) where output is requested.  The equation
c
c        dzeta = c(z)dz
c
c     is integrated, where c(z) gives the beam attenuation (in 1/meter)
c     at GEOMETRIC depth z.
c
c     A look up table of z to zeta values is also computed for use in
c     routine sumphas, when integrating the Riccati equations (when
c     optical depths must be converted to geometric depths before
c     calling the abscat routine).
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /CMISC/ imisc(30),fmisc(30) 
c
c     lookup table for z to zeta grid
      Common /Cztozeta/ nzvals, zetavals(mxnzvals),zvals(mxnzvals)
      integer nzvals
c
c     temporary storage for the numerical integration routine qag:
      parameter (mxlimit=25,mxlenw=4*mxlimit)
      dimension iwork(mxlimit),workq(mxlenw)
c 
c     the function beamcz(z) returns c(z), for numerical integration
      EXTERNAL beamcz
c
      data ierrp/0/
      save ierrp
C 
      nz = imisc(4)
c
c     parameters for qag
      epsabs = 1.0e-4
      epsrel = 1.0e-5
      key = 1
      limit = mxlimit
      lenw = mxlenw
!      nzvals = mxnzvals 
C
!     check to see if we have enough z memory to store ALL wavelen depths
      nzAvail = mxnzvals - nz
      ndz = min(int(float(nzAvail)/z(nz)), 10)  !use 0.1m res or max res given avail pts
!!      print *, 'ndz: ',ndz, nzAvail, nz, z(nz), 1./float(ndz), 
!!     1 float(nzAvail)/z(nz)

c     Compute zeta values at the z levels where output is requested
      ii = 1
      zeta(ii) = beamcz(z(1))*z(1)
      DO iz=2,nz
c
cccccc  Canned routine call for integration of c(z)dz  cccccccccccccccccc
!     get estimate of dtau for layer iz
      call qag(beamcz,z(iz-1),z(iz),epsabs,epsrel,key,value,
     1  abserr,neval,ier,limit,lenw,last,iwork,workq)
c
!     loop thru to find zetas in dtau resolution
      z1 = z(iz-1)
      zvals(ii) = z1
!     determine how many steps you need to make to have dtau=0.1 resolution
      ndtau = int(ndz *  (z(iz) - z1)) + 1
      dz = (z(iz) - z1)  / float(ndtau) 
      Do idtau = 1, ndtau
        z2 = z1 + dz  
        if(idtau.eq.ndtau) z2 = z(iz)
        call qag(beamcz, z1,z2, epsabs,epsrel,key,dval,
     1    abserr,neval,ier,limit,lenw,last,iwork,workq)
        ii = ii + 1
        zvals(ii) = z2
        zetavals(ii) = zetavals(ii-1) + dval
        z1 = z2
!       report if there is an error, on first error only     
        if(ier.ne.0 .and. ierrp.eq.0) then
          write(10,109)
          ierrp=1
        endif
      enddo
      nzvals = ii
      zeta(iz) = zeta(iz-1) + value
      end do
!!      print *, 'nzvals: ',nzvals
C
      RETURN
  109 format(/5x,'Routine "QAG" had a difficult time integrating',
     1' dzeta/c in routine zetatoz.'/5x,
     2'This error message does not mean the results are necessarily',
     3' inaccurate, so'/5x,
     4'ingore the error message from "QAG" unless the',
     5' output is clearly wrong.')
      END 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function beamcz(depth)
c
c     This function returns c(z), given depth as the GEOMETRIC depth z 
c     in meters, for use in the numerical integration routine.  It is
c     ASSUMED that abscat is expecting geometric depth.
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /CMISC/ imisc(30),fmisc(30) 
c
c     temporary local storage:
      dimension acomp(mxcomp),bcomp(mxcomp)
c
c     get the absorption and scattering coefs at geometric depth z
c
      ncomp = imisc(6)
      wavelen = fmisc(13)
c
c     insert the call to the desired "ab" routine:
      call abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
c
c     Check to be sure that the total attenuation is NOT zero (trivial case)
      if((atotal + btotal).le.0.0) then
        write(10,*)
        write(10,*) 'Total attenuation is ZERO.  Light field does not ',
     &             'change with depth.'
        call HERR("ZTOZETA","trivial case")  !stop run
      endif
c
      beamcz = atotal + btotal
c
      return
      end
