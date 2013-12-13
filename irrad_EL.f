C     Last change:  LKS   5 Apr 2008    5:38 pm
      subroutine irrad
c 
c     core routine on file irrad.f
c
c     called by RADANAL [MAIN->RADANAL->IRRAD]
c 
c     This routine computes various irradiances using the
c     azimuthally averaged radiances; mean cosines and the
c     irradiance reflectance are also computed and printed.

c     Irradiances are _computed_ at all zeta levels, for possible use in 
c     computing K-functions, etc., but _printout_ is only at selected 
c     depths, as specified in routine initial. 
c 
c     The zero element of arrays holds the values for zeta = a (in the air)
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
     1  fMUu(0:mxz),fMUd(0:mxz),fMUtot(0:mxz),R(0:mxz),E2(0:mxz)
      COMMON /Cpirrad/ ipirad,izirad(mxz)
      COMMON /Cmisc/ imisc(30),fmisc(30)
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)
c
c     declare temp vars
      integer nmu, jwave, nz
      real wavelen
c
      nmu = imisc(1)
      wavelen = fmisc(13)
c      nz = imisc(4)
c     compute the irradiances down the the max depth where they were
c     computed by solving the RTE.  This depth was determined in varzmax.
      jwave = imisc(11)
      nz = indexz(jwave) 
c
c     compute quantities in the air (at zeta = a) 
c 
      sum1 = 0. 
      sum2 = 0. 
      sum3 = 0.
      sum4 = 0.
      do i=1,nmu
         radM = RADupa(i)
         radP = RADsky(i)
         Omegamu = omega(i)
         sum1 = sum1 + radM*Omegamu
         sum2 = sum2 + radP*Omegamu
         sum3 = sum3 + radM*fmu(i)*Omegamu
         sum4 = sum4 + radP*fmu(i)*Omegamu
      end do
c 
      Eou(0) = sum1
      Eod(0) = sum2
      Eu(0) = sum3
      Ed(0) = sum4
c 
      Eo = Eou(0) + Eod(0)
c
      fMUu(0) = Eu(0)/Eou(0)
      fMUd(0) = Ed(0)/Eod(0) 
      fMUtot(0) = (Ed(0) - Eu(0))/Eo
      R(0) = Eu(0)/Ed(0) 
c 
      if(imisc(9).ge.0) then
c         *write in air values
		write(10,200) wavelen
		write(10,203) Eou(0),Eod(0),Eo,Eu(0),Ed(0),fMUu(0),
     1					fMUd(0),fMUtot(0),R(0)
	else
c	   ** minimal printout selected
c          *write in air values
		write(10,300) wavelen
		write(10,303) Eo,Eu(0),Ed(0),fMUu(0),
     1					fMUd(0),fMUtot(0),R(0)
	endif
c
c     Compute quantities within the water (w, ..., zeta, ..., m)
c
      do iz=1,nz

      sum1 = 0. 
      sum2 = 0. 
      sum3 = 0.
      sum4 = 0.
      sum5 = 0.
c 
c     compute irradiances
c 
      do i=1,nmu
         radM = RADupz(i,iz)
         radP = RADdnz(i,iz)
!      write(65,*) 'irrad: ',iz, i, radM, radP
         Omegamu = omega(i)
         sum1 = sum1 + radM*Omegamu
         sum2 = sum2 + radP*Omegamu
         sum3 = sum3 + radM*fmu(i)*Omegamu
         sum4 = sum4 + radP*fmu(i)*Omegamu
         sum5 = sum5 + (3.0*fmu(i)*fmu(i) - 1.0)*(radM + radP)*Omegamu
      end do
      Eou(iz) = sum1
      Eod(iz) = sum2
      Eu(iz) = sum3
      Ed(iz) = sum4
      E2(iz) = 0.5*sum5   !E2 is used for Raman calculations
c
      Eo = Eou(iz) + Eod(iz)
      fMUu(iz) = Eu(iz)/Eou(iz) 
      fMUd(iz) = Ed(iz)/Eod(iz) 
      fMUtot(iz) = (Ed(iz) - Eu(iz))/Eo
      R(iz) = Eu(iz)/Ed(iz) 
c 
c     check for printout
      if(imisc(9).ge.0) then
c         *write in water values
         iprint = 0
         do iiz=1,ipirad
           if(iz.eq.izirad(iiz)) iprint = 1
         end do
	   if(iprint.ne.0) write(10,202) iz,zeta(iz),z(iz),Eou(iz), 
     1	Eod(iz),Eo,Eu(iz),Ed(iz),fMUu(iz),fMUd(iz),fMUtot(iz),R(iz)
	else
c	   ** minimal printout selected
c         *write in water values
         iprint = 0
         do iiz=1,ipirad
           if(iz.eq.izirad(iiz)) iprint = 1
        end do
	   if(iprint.ne.0) write(10,302) iz,zeta(iz),z(iz), 
     1	Eo,Eu(iz),Ed(iz),fMUu(iz),fMUd(iz),fMUtot(iz),R(iz)
      endif
c
      end do
C 
      RETURN
C
  200 format(///2x,'Irradiances (units of W/m^2 nm), Mean Cosines',
     1' (Mubars), and Irradiance Reflectance at ',f6.1,' nm'//
     2'   iz   zeta   z(m)',8x,'Eou',12x,'Eod',13x,'Eo',
     3 13x,'Eu',13x,'Ed',7x,'MUBARu   MUBARd    MUBAR',6x,'R = Eu/Ed'/)
  202 FORMAT(I5,2F7.2,1P,5E15.4,0P,3F9.4,1P,E15.4)
  203 FORMAT(11X,'in air',2X,1P,5E15.4,0P,3F9.4,1P,E15.4/)
c
  300 format(//2x,'Irradiances (units of W/m^2 nm), Mean Cosines',
     1' (Mubars), and Irradiance Reflectance at ',f6.1,' nm'//
     2'   iz   zeta   z(m)',8x,'Eo',
     3 13x,'Eu',13x,'Ed',7x,'MUBARu   MUBARd    MUBAR',6x,'R = Eu/Ed'/)
  302 FORMAT(I5,2F7.2,1P,3E15.4,0P,3F9.4,1P,E15.4)
  303 FORMAT(11X,'in air',2X,1P,3E15.4,0P,3F9.4,1P,E15.4/)
C 
      END 
