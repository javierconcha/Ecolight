C     Last change:  LKS  24 Nov 2008    8:45 am
      SUBROUTINE QASKY
c 
c     core routine on file qasky.f
c
c     called by INISHRAD [MAIN->RADIANCE->INISHRAD->QASKY]
c
c     calls sunang (included in this file) and the sky radiance and 
c           sky irradiance models 
c 
c     This routine computes the quad-averaged incident sky radiances
c
c     The particular model used to compute the sky (relative or
c     absolute) radiance is inserted via a call to SKYRAD.
c
c     The particular model used to compute the sky irradiance is
c     inserted via a call to SKYIRRAD.
c
c     The sky radiances returned from the skyrad model are always
c     re-scaled to give agreement with the direct and diffuse
c     irradiances returned by the skyirrad model.
c
c     NOTE: for consistency, sky and solar angles are always in 
c     radians for calls to sky models.
c
      INCLUDE "DIMENS_XL.INC"
c
C     The following common block array was written on 1/19/2000 
C	to store and printout the above-water direct and diffuse Ed's 
C     at each wavelength (block shared with excel.f only)
      common /CErik/ EdifOut(mxwave), EdirOut(mxwave)

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
      common /Csky/  iskyflag,skyspecs(mxnsky)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      common /cgcirr1/ iein,iblw,jday,rlon,rlat,the,GMThr,pres,am,rh,
     1                 wv,wsm,ws,vi,ro3
      common /cmisc/ imisc(30),fmisc(30)
c
      data kall/0/
      save
c
c     temporary local storage:
      dimension thetab(mxmu),skyback(mxmu)
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      if(kall.eq.0) then
c     preliminary calculations on the first call.

c     Note that the unnormalized PATTERN of the background sky
c     radiance is the same at each wavelength (for the default sky
c     radiance model, hcnrad), so the quad averaging calculations
c     need to be done only once.  Only the IRradiances depend on
c     wavelength.

      nmu = imisc(1)
      pi = fmisc(1)
      degrad = fmisc(2)
      radeg = fmisc(3)
c
c---------------------------------------------------------------------
c
c     Preliminary calculations for sky radiance and irradiance models.
c
c     Various parameters needed as input by the sky radiance and
c     irradiance models are contained in the skyspecs array, or in the
c     common blocks used by the particular routines.
c
      if(iskyflag .eq. 3) then
c
c     The run uses julian day, lat, long, and GMT to compute the
c     solar zenith angle.
c     NOTE: these values will override the defaults (set in SETDFLTS)
c     for computation of ozone concentration in the skyirrad model
c     gcirrad.
c
      jday = ifix(skyspecs(6))
      rlat = skyspecs(7)
      rlon = skyspecs(8)
      GMThr = skyspecs(9)
c
c     compute the corresponding solar zenith angle
c
      call sunang(jday,GMThr,radeg,rlon,rlat,suntheta,suna)
c
      skyspecs(1) = suntheta
c     NOTE:  the computed sun azimuth angle, suna, is referenced to
c     true north.  The sun azimuth angle in Hydrolight, sunphi, is
c     referenced to the downwind direction.  Therefore, the default 
c     value of sunphi is not replaced by suna.
c
      write(10,503) jday,rlat,rlon,GMThr
      write(10,504) suntheta,suna
c
c     check for bad input
c    
         if(suntheta.gt. 90.0) then
         call HERR("QASKY","computed solar zenith angle > 90 deg")  !stop run
         endif
c
      endif
c
c--------------------------------------------------------------------
c
c     Compute the quad-averaged, unnormalized background sky radiances
c
c*****USER INPUT*********
c     Set the sub-quad partitioning for quad averaging of the sky 
c     radiance distribution; nsubmu and nsubphi = 1 and 36 correspond
c     to 10 x 10 deg sky resolution, which is good
c     enough for most purposes.

      nsubmu = 1
      nsubphi = 36
C*****END USER INPUT*********

      dphi = 2.0*pi/float(nsubphi)
c
c     thetas is known either from input or from the call to sunang
      suntheta = degrad*skyspecs(1)
      sunphi = degrad*skyspecs(2)
c     all angles are now in radians
c
      do iu=1,nmu
         dmu = deltmu(iu)/float(nsubmu)
         if(iu.eq.1) then
           umumin = 0.5*dmu
         else
           umumin = bndmu(iu-1) + 0.5*dmu
         endif
c
c     integrate the sky radiance over quad (mu band) Qu = Q(iu)
c     The phi integral is over 0 to 2pi
      sum = 0.0
      do j=1,nsubphi
         skyphi = float(j-1)*dphi
      do i=1,nsubmu
         skymu = umumin + float(i-1)*dmu
c
c     Obtain the diffuse (background) sky radiance for the current
c     sky (theta, phi) and for the given solar (theta, phi).
c     Unnormalized or relative radiance values will be set to the proper
c     magnitudes below. 
c
c---------------------------------------------------------------------
c     Insert the call to the desired "skyrad" routine, of the form:
c     call skyrad(suntheta,sunphi,skytheta,skyphi, skyrad)
c
c     NOTE:  the "skyrad" routine can be any routine that
c     returns a sky radiance value in direction (skytheta,skyphi).
c
c     skyrad routines always expect angles in radians
      skytheta = acos(skymu)
c
      call skyrad(suntheta,sunphi,skytheta,skyphi,skyrad0)
c
      sum = sum + skyrad0*dmu
      end do  ! i
      end do  ! j
c
      skyback(iu) = sum*dphi/omega(iu)
      end do  ! iu

c     -----------------------------------------------------------
c
c     Compute the downwelling diffuse irradiance associated with
c     the quad-averaged, UNnormalized diffuse sky radiance
c
      Ednorm = 0.0
      do i=1,nmu
         Ednorm = Ednorm + skyback(i)*fmu(i)*omega(i)
      end do

c     determine which quad contains the sun:
c 
c     convert the boundary mu values to degrees
      do i=1,nmu
         thetab(i) = radeg*acos(bndmu(i))
      end do
c
c     determine the mu index of the quad containing the sun
c
      thetas = suntheta*radeg
      imus = nmu
      do i=1,nmu-1
      if(thetas.lt.thetab(i) .and. thetas.ge.thetab(i+1)) imus = i + 1
      end do
      if(thetas.gt.thetab(1)) imus = 1
c 
  200 continue
c 
c        get the exact theta value at quad center for printout
         if(imus.eq.1) then
            theq = 0.5*(90.0 + radeg*acos(bndmu(imus)))
         elseif(imus.eq.nmu) then
            theq = 0.
         else
            theq = 0.5*radeg*(acos(bndmu(imus-1)) + acos(bndmu(imus)))
         endif
 
      write(10,510) imus,theq
c
c     Insert the call to the desired "skyirrad" routine:
      call skyirrad(suntheta,sunphi, Eddif,Eddir)

      kall = 1
      return     !first kall just initializes and generated msgs
      endif

c     end of preliminary calculations on the first call to qasky
c---------------------------------------------------------------------
c
c     subsequent calls start here and use the same background sky
c     UNnormalized radiance pattern
c
c---------------------------------------------------------------------
c     Insert the call to the desired "skyirrad" routine:
      call skyirrad(suntheta,sunphi, Eddif,Eddir)
c
c--------------------------------------------------------------------
c
c     Scale the quad-averaged UNnormalized diffuse sky radiances so
c     that they yield the diffuse Ed value at this wavelength
c
      factor = 0.0
      if(Ednorm.ne.0.0) factor = Eddif/Ednorm
      do i=1,nmu
         RADsky(i) = skyback(i)*factor
      end do
c
c     Add the solar direct beam quad-averaged radiance to the
c     appropriate quad
c
      RADsky(imus) = RADsky(imus) + Eddir/(fmu(imus)*omega(imus))
c
!      write(10,*) '**Lsky:  ',(RADsky(i),i=1,nmu)
c     RADsky(mu) now contains the absolute quad-averaged sky
c     radiance distribution, where mu gives the VIEWING
c     direction.
c
c     As a check, compute the total Ed from the final quad-averaged
c     sky radiances, and see that it equals the total Ed value from
c     the sky irradiance model.
c
      Edfinal = 0.0
      do i=1,nmu
         Edfinal = Edfinal + radsky(i)*fmu(i)*omega(i)
      end do
c
ccc	*** Store Diffuse and Direct sky Ed's at this wavelength
         jwave = imisc(11)	!get index to current wavelength
	EdifOut(jwave) = Eddif		!store diffuse Ed (in air) into array
	EdirOut(jwave) = Eddir		!store direct Ed (in air) into array
ccc
c     Print out normalization info on first call iff idbug>0 (added printout requested)
!      if(kall.eq.0.and.imisc(9).gt.0) write(10,520) Edfinal,Eddif+Eddir
c
      return
c 
c     formats 
c 
  503 format(/5x,"The sun's zenith angle is computed from"/
     110x,'day of year = ',i4,' (1 is Jan 1)'/
     210x,'latitude    = ',f9.4,' degrees (N is positive)'/
     310x,'longitude   = ',f9.4,' degrees (E is positive)'/
     410x,'GMT         = ',f7.2,' hours')
  504 format(/5x,'The computed solar zenith angle is',f6.2,' degrees'//
     15x,'The computed solar azimuth angle is',f7.2,' degrees (relative 
     2to true north)'//
     35x,'The computed solar azimuth angle is not used by EcoLight, ',
     4/5x,'which uses an azimuthally averaged solar radiance.'/)
  510 format(/5x,'The sun is placed in quad Q(r) = Q(',i2,
     1') centered at theta = ',f6.3/)
  520 format(/5x'Check on the sky radiance initialization:'/
     1'      total Ed from integrating the sky radiance =',1p,e11.4/
     2'      total Ed from the sky irradiance model     =',e11.4)
      end 
c

