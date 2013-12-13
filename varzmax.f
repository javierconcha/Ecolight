C     Last change:  LKS   7 Aug 2008    6:46 pm
      Subroutine varzmax
C     
C     core routine on file varzmax.f

c     new version of max depth algorithm using abar

c     called by routine RADIANCE if ioptflag .ne. 0
C
c     This routine determines the maximum depth to be used in solving the
c     RTE at each wavelength.
c
c     if ioptflag = 0, solve the RTE to the same depth at each wavelength,
c                      as in done in Hydrolight; zmax = max z of printout
c     if ioptflag = 1, dynamically determine the max depth so that PAR will
c                      be accurately computed throughout the euphotic zone.
c                      zmax will be estimated from a min PAR value
c     not yet implemented:
c     if ioptflag = 2, dynamically determine the max depth from 1/Kd at
c                      each wavelength.  This option is for use when only
c                      the water-leaving radiance is of interest.

      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Ciop/ acoef(mxz,0:mxcomp),bcoef(mxz,0:mxcomp),
     1		      atten(mxz),albedo(mxz), bbcoef(mxz,0:mxcomp)
      COMMON /CEospl/ nspl,zspl(mxz),Eospl(mxz,mxwave),
     1                E2spl(mxz,mxwave),Edspl(mxz,mxwave)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 

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
C
C     common blocks modified here
      COMMON /Cmisc/ imisc(30),fmisc(30) 
      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cpkfcn/ ipkfcn,izkfcn(mxz) 
C
c     Common blocks associated with the PAR testing 
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)

c     common blocks for communication with the Gregg & Carder model
c     (routine GCEd):
      parameter (nlt=701)
      common /cgcirr1/ iein,iblw,jday,rlon,rlat,the,hr,pres,am,rh,
     1                 wv,wsm,ws,vi,ro3
      common /cgcirr2/ Eddirgc(nlt),Eddifgc(nlt),Edtotgc(nlt)

      dimension expabar(mxz,mxwave)
      integer ioptflag
C
      Data kall/0/
C
      save
C
C     Initialize this routine on first call
      If(kall.eq.0) then

c     *** Constants for calculating photons from energy
!         planck = 6.626e-34
!         speed = 2.998e8
!         fnmtom = 1.0e-9
!         fmoltoumol = 1.0e6
!         avagadro = 6.023e23
c        factor converts energy to number of photons
cc         factor = fnmtom*fmoltoumol/(planck*speed*avagadro)
cc         nmu = imisc(1)
cc         nwave = imisc(7)
         nz = imisc(4)
         nwskip = imisc(26)
c
c        Save the user-requested max depth
         nzinit = imisc(4)   !this is the depth index of the initial max depth
         indexz(0) = nzinit

c        set parmin for debugging; make a user input later
c      iGC = 0
c      if(iGC .ne. 0) then
c         PARMIN = 1.0     ! the min PAR in (micromol photons)/(m^2 s)

c     call the Gregg and Carder irradiance model to estimate the clear-sky,
c     in-water quantum Ed for use in estimating PAR(0)

c     iein : = 0 to compute irradiance (W/m2/nm)
c            = 1 to compute quanta (microEinst/m2/sec)
c     iblw : = 0 to compute above-surface values
c            = 1 to compute below-surface values
cc      iein = 1
cc      iblw = 1
cc      call GCed
c     reset to above-water, irradiance for subsequent calculations
cc      iein = 0
cc      iblw = 0
c     Edtotgc now contains spectral quantum Ed in (micromol phot)/(m^2 s)
c     at 1 nm resolution from 350 to 800 nm.
c     Use this from 400 to 700 to estimate PAR(0).  Note that PAREd0
c     will underestimate the true PAR computed from Eo; this makes
c     fPAR too small, so we go too deep, which gives a safety factor
cc      PAREd0 = 0.0
cc      do i=51,351  ! 350 to 700 nm
cc         PAREd0 = PAREd0 + Edtotgc(i)
cc      enddo
cc      tunefact = 1.0  ! a tuning factor for later use
cc      FPAR = tunefact*PARmin/PAREd0

cc      write(10,100) PARmin,PAREd0,FPAR
c  100 format(/' The user-requested minimum PAR value is PARmin =   ',
c     1 f8.2,' (micromol photons)/(m^2 s)'/
c     2' The clear-sky Ed(400-700) estimated PAR(0) value is',f8.2,
c     3' (micromol photons)/(m^2 s)'/
c     4' FPAR = PARmin/PAR(0) = ',1p,e9.2)
c      endif

      FPAR = fmisc(25)
**    FPAR = 0.2
!      write(10,105) FPAR
  105 format(/' The user-requested Fo value for use in estimating the ,'
     1,'max computation depth zo is Fo =',f5.2)
c
      ioptflag=imisc(25)
      kall = 1
      Endif
C
C     All subsequent calls begin here....
      jwave = imisc(11)

c     Cannot optimize zmax if using Dynamic Depth!!
      If(imisc(27).eq.1 .and. imisc(12).eq.0 .and. imisc(8).gt.0) then
        indexz(jwave) = imisc(4)
        zopt(jwave) = z(imisc(4))
      Elseif(ioptflag .eq. 0) then
c        use the user-requested max depth at each wavelength
         indexz(jwave) = nzinit
         zopt(jwave) = z(nzinit)
         return
      ElseIf(ioptflag.eq.1) then

c     determine zmax so that PAR will be accurately computed down to the 
c     requested minimum PAR value, based on FPAR at each wavelength

c        compute the depth integrated abs coef factor at this wavelength
c        note that acoef(z,0) is the total abs coef at the current lambda
         
         expabar(1,jwave) = 1.0
            do iz=2,nz
               expabar(iz,jwave) = expabar(iz-1,jwave) *
     1         exp(-0.5*(acoef(iz,0)+acoef(iz-1,0))*(z(iz) - z(iz-1)))
c           write(10,*)' z, expabar = ',z(iz),expabar(iz,jwave)
            enddo

      if(jwave.eq.1) then

c     special calculation for first wavelength
c     find the depth where exp(-az) becomes less than FPAR

      do iz=1,nz
        if(fpar .ge. expabar(iz,jwave)) go to 200
      enddo
      izindx = nz
      go to 202
  200 izindx = iz
  202 continue
      indexz(jwave) = izindx
      zest = z(izindx)
c     for the first wavelength, accept the initial estimate
      zopt(jwave) = zest
      zestprev = zest

      else

c     calculations for wavelengths 2 to nwave use info from the
c     previous wavelength to correct the abar estimate
c
c     find the actual fPAR depth on the previous wavelength
c     (Eo was extrapolated to z(nz) in radanal)
      Eo0 =exp(Eospl(1,jwave-nwskip))
      do iz=1,nz
         Eo = exp(Eospl(iz,jwave-nwskip))
         if(fpar .ge. (Eo/Eo0))
     1      go to 210
      enddo
      zfpar(jwave-nwskip) = z(nz)
      go to 211
  210 Eo1 = exp(Eospl(iz-1,jwave-nwskip))
      zfpar(jwave-nwskip) = z(iz-1) + (z(iz) - z(iz-1))*
     1 log(fpar*Eo0/Eo1)/
     2 log(Eo/Eo1)
  211 continue
c
c     get the abar-estimated max depth at the current wavelength
      do iz=1,nz
         if(fpar .ge. expabar(iz,jwave)) go to 220
      enddo
      indexz(jwave) = nz
      go to 221
  220 indexz(jwave) = iz
  221 continue
      zest = z(indexz(jwave))
!      ztemp = zest
c
c     use the actual and estimated fpar depths at the previous
c     wavelength to correct the present estimate

      zfinal = zest*zfpar(jwave-nwskip)/zestprev

c     save the initial estimate for use next time
      zestprev = zest
c
c     find the index of the final estimated max depth (this is
c     the next printout depth below the estimated depth)
      do iz=1,nz
         if(zfinal .le. z(iz)) go to 230
      end do
      indexz(jwave) = nz
      go to 231
  230 indexz(jwave) = iz
  231 continue
c     set zopt to the next printout depth below zfinal
      zopt(jwave) = z(indexz(jwave))

      endif
      endif

!!      write(10,240)jwave,zfpar(jwave-nwskip),ztemp,zopt(jwave),
!!     1             indexz(jwave)
  240 format(' jwave = ',i3,'  zfpar(prev) =',f7.3,'   zest = ',f7.3,
     1'   zopt = ',f7.3,'  indexz = ',i3)

      Return
      End