C     Last change:  LKS   7 Aug 2008    6:49 pm
      SUBROUTINE RADIANCE
c
C     core routine on file radiance.f
c
c     called by MAIN
c
c     calls INISHRAD, BOTMBC, RICCATI, RADW, RADZETA, MATxMAT, VECxMAT
C
c     This routine solves for the radiance as a function of mu band and
c     depth for the current wavelength.  These calculations are
c     flowcharted in L&W Fig. 8.2 (as revised).  Various auxillary quantities
c     (such as water-leaving and surface-reflected radiances in the air)
c     are computed and saved for later use.
c 
      INCLUDE "DIMENS_XL.INC"
C
      common /Cbbopt/ ibbopt(mxcomp), bbfrac(mxcomp), 
     1                BfrefPL(mxcomp), Bf0PL(mxcomp), BfmPL(mxcomp)
      integer ibbopt
      real bbfrac, BfrefPL, Bf0PL, BfmPL
C 
      COMMON /CRTS/ Rzw(mxmu,mxmu,mxz),Twz(mxmu,mxmu,mxz),
     1              Rzb(mxmu,mxmu,mxz),Sptwz(mxmu,mxz),Smtbz(mxmu,mxz)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Crhotau/ rho(mxmu,mxmu),tau(mxmu,mxmu),
     1                 betatP(mxmu,mxmu,mxcomp),betatM(mxmu,mxmu,mxcomp)
      COMMON /CBOTBC/ Rmb(mxmu,mxmu)
      COMMON /Crtsurf/ raw(mxmu,mxmu),taw(mxmu,mxmu),rwa(mxmu,mxmu),
     1                 twa(mxmu,mxmu)
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)

c     radiance arrays are stored as follows:
c        RADsky(mu) = total incident sky radiance
c        RADwla(mu) = Lw(mu,z=a) = water-leaving radiance
c        RADrsa(mu) = Lreflsky(mu,z=a) = surface-reflected sky rad
c        RADupa(mu) = Lu(mu,z=a) = total upwelling rad in air
c           note: RADupa = RADwla + RADrsa, i.e., Lu = Lw + Lreflsky in air
c        RADupz(mu,z) = Lup(mu,z) = total upwelling rad in water
c        RADdnz(mu,z) = Ldn(mu,z) = total downwelling rad in water
      COMMON /Crad/ RADsky(mxmu),RADwla(mxmu),RADrsa(mxmu),RADupa(mxmu),
     1              RADupz(mxmu,mxz),RADdnz(mxmu,mxz)
      COMMON /Cpirrad/ npirad,izirad(mxz)
                          
      COMMON /CMISC/ imisc(30),FMISC(30) 
      COMMON /Cfrstcls/ iabscat,iqasky,iradian,iradanal,iradxcl
c
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
c
      dimension nuphas(mxcomp)
c
      DATA nurad/11/
c
c     temporary local storage:
      dimension RADupw(mxmu)
      integer errL, izerr 
      save
C
C     **********  INITIALIZATION  **********
C
C     Open the needed files on first call to radiance
C
      if (iradian.eq.1) then
c
c     Open files of phase functions:
*     IFF not using bbfrac(component) to select phase function
*		[ibbopt = -1 if no phase function needed (i.e., CDOM) ]
*		[ibbopt =  0 if phase function read from file ]
*		[ibbopt =  1,2,3 if phase function selected by bb/b ratio ]
*     [note: unit numbers are still 49+icomp (but some unit #s are skipped]
c
      ncomp = imisc(6)
      npfiles = iabs(ncomp)
      DO i=1,npfiles
         IF(ibbopt(i).ne.0) goto 10
         nuphas(i) = 49+i
         OPEN(nuphas(i),file=pfname(i),status='old')
 10   end do
c
c     the output file for radiances will be opened in subroutine radanal
c
      endif
c 
c     Initialize for this wavelength
c
      CALL INISHRAD(nuphas)
c
      CALL varzmax
c
      imisc(16) = nurad
C
      nmu = imisc(1)
c
C     **********  BEGIN COMPUTATIONS  **********
C 
c*****NOTE:  The model solves for the TOTAL radiances, not
c     for the direct and diffuse parts separately.  However,
c     for values in air (at z=a), direct and diffuse parts are computed
c     for later use in computing Rrs = Lw/Ed, etc.
C
C++++ COMPUTE Rmb = R(m,b) for the desired bottom boundary condition
C
      CALL BOTMBC
c
C++++ Integrate the Ricatti equations to get Rzw = R(zeta,w), Twz = T(w,zeta), etc
C
      CALL RICCATI
C
C++++ Compute the radiances L+(w) = Ld(w) and L-(w) = Lu(w) at zeta = w
c     (just below the air-water surface) using Eqns. (8.98) and (8.102)
C
      CALL RADw
C 
C++++ Compute the radiances Ld = L+(z) and Lu = L-(z) at all
c     interior depths, w .lt. zeta .le. m, using Eqns. (8.105) and (8.106)
C
      CALL RADzeta
c 
C++++ Compute the upward total radiances Lu(a) = L-(a) = RADupa just
c     above the sea surface, using Eq. (8.107).  These radiances
c     include both the water-leaving radiance and reflected sky
c     radiance.
c
c     first term:
      do i=1,nmu
         RADupw(i) = RADupz(i,1)
      enddo
      call VECxMAT(RADupw,twa,nmu,nmu,mxmu,RADwla)
c
c     second term:
      call VECxMAT(RADsky,raw,nmu,nmu,mxmu,RADrsa)
c
      DO i=1,nmu
         RADupa(i) = RADwla(i) + RADrsa(i)
      end do
C
C     **********  END OF RADIANCE COMPUTATIONS  **********

C     ***  Ensure all radiances are non-negative (numerical errors may
C          occur when strong sources are present with an inf deep bottom
c          in radiances near the horizon  ***
      nz = izirad(npirad)
      errL = 0 
      izerr= 1
      Do i=1,nmu
       If(RADwla(i).lt.0 .or. RADupa(i).lt.0) errL = 1
       If(RADwla(i).lt.0) RADwla(i)=0.0
       If(RADupa(i).lt.0) RADupa(i)=0.0
       Do iz=nz,1,-1
         If(RADupz(i,iz).lt.0 .or. RADdnz(i,iz).lt.0) then
           errL = 1
           izerr = iz
         Endif
         If(RADupz(i,iz).lt.0) RADupz(i,iz)=0.0
         If(RADdnz(i,iz).lt.0) RADdnz(i,iz)=0.0
       Enddo
      Enddo
      If(errL.gt.0) write(10,100) z(izerr)
c
c     Close all files.
      DO i=1,ncomp
         close(nuphas(i))
      end do
c

      return
  100 format(//,'WARNING:  some calculated radiances were negative.',
     1  /5x,'Results below ',f6.1,' m may not be correct. ')
      end
