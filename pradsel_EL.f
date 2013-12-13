C     Last change:  LKS   1 Jul 2008    4:40 am
      subroutine pradsel
c 
c     core routine on file pradsel.f
c 
c     called by RADANAL [MAIN->RADANAL->PRADSEL]
c 
c     This routine prints zenith and nadir radiances.
c     Radiance-irradiance ratios are also computed and printed.
c 
      INCLUDE "DIMENS_XL.INC"
c
c     radiance arrays are stored as follows:
c        RADsky(mu) = (old RAD0Pa) = L+(z=a;mu) = total incident sky radiance
c        RADwla(mu) = Lw(mu,z=a) = water-leaving radiance
c        RADrsa(mu) = (old RAD0Ma) = Lreflsky(mu,z=a) = surface-reflected sky rad
c        RADupa(mu) = (old RADMa) = Lu(mu,z=a) = total upwelling rad in air
c           note: RADupa = RADwla + RADrsa, i.e., Lu = Lw + Lreflsky in air
c        RADupz(mu,z) =(old RADMz) = L-(mu,z) = total upwelling rad in water = Lu(mu,z)
c        RADdnz(mu,z) = (old RADPz) = L+(mu,z) = total downwelling rad in water = Ld(mu,z)
      COMMON /Crad/ RADsky(mxmu),RADwla(mxmu),RADrsa(mxmu),RADupa(mxmu),
     1              RADupz(mxmu,mxz),RADdnz(mxmu,mxz)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cirrad/ Eou(0:mxz),Eod(0:mxz),Eu(0:mxz),Ed(0:mxz), 
     1                fMUu(0:mxz),fMUd(0:mxz),fMUtot(0:mxz),R(0:mxz),
     2                E2(0:mxz)
      COMMON /Cpirrad/ ipirad,izirad(mxz)
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)

      COMMON /CMISC/  IMISC(30),FMISC(30) 
c 
c     declare temp vars
      integer nmu, jwave, nzc
      real wavelen
c
      nmu = imisc(1)
      wavelen = fmisc(13)
c
c     printout down the the max depth where quantities were
c     computed by solving the RTE.  This depth was determined in varzmax.
      jwave = imisc(11)
      nzc = min(ipirad, indexz(jwave))
c
c     zeta = a (in the air, just above the surface)
c 
c     radup is the TOTAL upward radiance, including sky radiance
c     that is reflected upward by the water surface.
c     radw is the "water-leaving" radiance, which does not contain
c     the reflected sky radiance.
c
      radup = radupa(nmu)
      radw = RADwla(nmu)
      raddn = RADsky(nmu)

      Rrstot = radup/Ed(0)
      Q = Eu(0)/radup
      Rrs = radw/Ed(0)
c
      write(10,200) wavelen
      write(10,202) radup,raddn,Rrstot,Q,radw,Rrs
c 
c     depths w .le. zeta .le. m
c 
      do iiz=1,nzc
         iz = izirad(iiz)
         radup = RADupz(nmu,iz)
         raddn = RADdnz(nmu,iz)
         Rrs = radup/Ed(iz)
         Q = Eu(iz)/radup
         write(10,204)iz,zeta(iz),z(iz),radup,raddn,Rrs,Q
      end do
c 
      return
C 
  200 FORMAT(//2x,'Selected Radiances (units of W/m^2 sr nm) and',
     1' Radiance-Irradiance Ratios at ',f6.1,' nm'//,
     1'   iz   zeta     z       Lu(z)        Ld(z)        ',
     3'Lu/Ed      Q = Eu/Lu      Lw(z)     Rrs = Lw/Ed'/
     48x,'(m)',39x,
     5'(1/sr)',8x,'(sr)',8x,'(nadir)      (1/sr)'/)
  202 FORMAT(11X,'in air  ',1P,6E13.3/)
  204 FORMAT(I5,2F7.2,1P,4E13.3) 
      END 
