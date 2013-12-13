C     Last change:  LKS  27 Feb 2008    2:10 pm
      SUBROUTINE zetatoz 
c 
C     core routine on file zetatoz.f
C
C     called by INISHRAD [MAIN->RADIANCE->INISHRAD->ZETAOZ]
C
C     calls routines beamc1 (in this file) and QAG (in file quadrat.f)
C 
C     This routine computes the GEOMETRICAL depths z (in meters) that 
c     correspond to the OPTICAL depths zeta (nondimensional) where 
c     output is requested.  The equation
c
c         dz = dzeta/c(zeta)
c
c     is integrated, where c(zeta) gives the beam attenuation in
c     1/meters at OPTICAL depth zeta.
c
c     (This routine is called only when ECOLIGHT is run in
c     monochromatic mode with iop = 1)
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
c     temporary storage for the numerical integration routine qag
      parameter (mxlimit=25,mxlenw=4*mxlimit)
      dimension iwork(mxlimit),work(mxlenw)
c
c     function beamc1(zeta) returns 1.0/c(zeta)
      EXTERNAL beamc1
c
      data ierrp/0/
      save ierrp
c
c     parameters for numerical integration 
      epsabs = 1.0e-4
      epsrel = 1.0e-5
      key = 1
      limit = mxlimit
      lenw = mxlenw
c
      nz = IMISC(4) 
c
      z(1) = zeta(1)*beamc1(zeta(1))
c     
      DO iz=2,nz
c
cccccc Canned Routine Call for numerical integration ccccccccccccccc
c
      call qag(beamc1,zeta(iz-1),zeta(iz),epsabs,epsrel,key,value,
     1  abserr,neval,ier,limit,lenw,last,iwork,work)

      if(ier.ne.0 .and. ierrp.eq.0) then
        write(10,109)
        ierrp=1
  109 format(/5x,'Routine "QAG" had a difficult time integrating',
     1' dzeta/c in routine zetatoz.'/5x,
     2'This error message does not mean the results are necessarily',
     3' inaccurate, so'/5x,
     4'ingore the error message from "QAG" unless the Hydrolight',
     5' output is clearly wrong.')
      endif
c
cccccc End canned routine call
c
      z(iz) = z(iz-1) + value
      end do
C 
      RETURN
      END 
c
c================================================================
      function beamc1(depth)
c
c     this function returns 1.0/c(zeta), given depth as the OPTICAL
c     depth zeta, for use in the integration routine.  It is ASSUMED
c     that the "ab" routine is expecting optical depth.
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /CMISC/ imisc(30),fmisc(30) 
c
c     temporary local storage:
      dimension acomp(mxcomp),bcomp(mxcomp)
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
        call HERR("ZETATOZ","trivial case")  !stop run
      endif
c
      beamc1 = 1.0/(atotal + btotal)
c
      return
      end


