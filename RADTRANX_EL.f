C     Last change:  LKS  13 Nov 2008   11:03 am
      subroutine RADTRANX(suntheta,sunphi, Eddif,Eddir)
c
c     called by qasky
c
c     calls skyparam in incfiles.for to get atmos params
c
c     This is a driver routine for the Gregg and Carder clear-sky 
c     irradiance model (RADTRAN) and associated routines (reference:
c     Gregg and c     Carder, Limnol. Oceanogr., 35(8), 1657-1675, 1990).
c     The original model has been extended from 700 to 800 nm for use in
c     Hydrolight 5.0.
c
c     This routine calls the G&C model to get the direct and diffuse
c     spectral irradiances from 300 to 1000 nm at 1 nm resolution.
c     These irradiances are then band averaged if necessary, so that 
c     Eddif and Eddir are averages over the current wavelength band.
c
c     The G&C clear-sky irradiances are ajusted by the model of Kasten 
c     and Czeplak if the cloud fraction is greater than 0.25 (reference:
c     Kasten and Czeplak, Solar Energy, 24, 177-189, 1980)
c
c     Conversion from time and location to solar angle has already been
c     done in routine qasky, if necessary
c 
c*****NOTE:  routine lidata opens and reads file "..\data\gcirrad.txt"
c*****NOTE:  routine lidata_new opens and reads file "..\data\RADTRANX_dbase.txt"
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Csky/ iskyflag,skyspecs(mxnsky)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      common /cmisc/ imisc(30),fmisc(30) 
c     common blocks for communication with the Gregg & Carder model
c     (RADTRANX: routine GCEd):
      parameter (nlt=701)
      common /cgcirr1/ iein,iblw,jday,rlon,rlat,the,GMThr,pres,am,rh,
     1                 wv,wsm,ws,vi,ro3
      real pres,am,rh,wv,vi,wsm,ro3
      common /cgcirr2/ Eddirgc(nlt),Eddifgc(nlt),Edtotgc(nlt)
      COMMON/gcifname/gcifile
      character*32 gcifile
c
      data kall/0/
      save kall
c
      nwave = imisc(7)
c
      the = skyspecs(1)
      cloud  = skyspecs(10)

      if(kall.eq.0) then
          write(10,100)	!cloud  (HCNRAD already printed out the cloud fraction)
          call GCEd        !get msgs
          kall = 1
          return
      Endif
c     -------------------------------
c
c     Obtain the clear-sky spectral irradiances (350-700 nm at 1 nm 
c     resolution) incident onto the sea surface, using the Gregg and 
c     Carder model (sub GCEd):
c
      call GCEd
c
c     obtain the spectral irradiances at the current wavelength, or
c     average them over the current wavelength band
c
      if(nwave.eq.0) then
c
c        monochromatic run:  use the exact wavelength
c
         wavenm = fmisc(13)
c         ilam = ifix(wavenm + 0.5) - 349   ! for Gregg & Carder only
         ilam = ifix(wavenm + 0.5) - 299   ! for Gregg & Carder only
         Eddif = Eddifgc(ilam)
         Eddir = Eddirgc(ilam)
         Edtot = Edtotgc(ilam)
         write(10,532) wavenm,Eddif,Eddir,Edtot
c
      else
c
c     average the irradiances over the current wavelength band
c
         jwave = imisc(11)
         wave1 = waveb(jwave)
         wave2 = waveb(jwave+1)
c         ilam1 = ifix(wave1 + 0.5) - 349
c         ilam2 = ifix(wave2 + 0.5) - 349
         ilam1 = ifix(wave1 + 0.5) - 299
         ilam2 = ifix(wave2 + 0.5) - 299
         count = 0.0
         sumdif = 0.0
         sumdir = 0.0
         sumtot = 0.0
         do ilam=ilam1,ilam2
            count = count + 1.0
            sumdif = sumdif + Eddifgc(ilam)
            sumdir = sumdir + Eddirgc(ilam)
            sumtot = sumtot + Edtotgc(ilam)
         end do
         Eddif = sumdif/count
         Eddir = sumdir/count
         Edtot = sumtot/count
         write(10,533) wave1,wave2,Eddif,Eddir,Edtot
c
      endif
c
c     ---------------------------------
c
c     Adjust the clear-sky irradiances by the model of Kasten and Czeplak
c     if the cloud fraction is greater than 0.25.
c
      if(cloud .gt. 0.25) then
         Edtot = Edtot*(1.0 - 0.75*cloud**3.4)
         Eddif = Edtot*(0.3 + 0.7*cloud**2)
         Eddir = Edtot - Eddif
         write(10,534) Eddif,Eddir,Edtot
      endif
c
      return
c
c     formats
c
  100 format(5x,'Diffuse and Direct Sky irradiances are computed ',
     6'using the "RADTRANX" (RADTRAN eXtended for 300-1000 nm) model'/
     25x,'The RADTRAN clear-sky irradiances are adjusted for cloudines',
     3's by the model of Kasten and Czeplak')
c      using'//5x,
c     4'cloud = cloud fraction (0 for clear to 1 for heavy overcast)
c     5 =',f5.2)
  532 format(/2x'RADTRANX clear-sky spectral irradiances at',
     1f6.1,' nm:'/10x,'Ed(diffuse) =',1pe11.3,5x,'Ed(direct) =',
     2e11.3,5x,'Ed(total) =',e11.3/)
  533 format(/2x'RADTRANX clear-sky spectral irradiances averaged over',
     1' the band from',f6.1,' nm to',f6.1,' nm:'/
     110x,'Ed(diffuse) =',1p,e11.3,5x,'Ed(direct) =',e11.3,5x,
     2'Ed(total) =',e11.3)
  534 format(2x'Spectral irradiances after adjusting for cloud cover:
     1'/10x,'Ed(diffuse) =',1p,e11.3,5x,'Ed(direct) =',e11.3,5x,
     2'Ed(total) =',e11.3/)
c
      end
c
c**********************************************************************
c
      subroutine GCEd
c
c     user "skyirrad" model on file RADTRANX.f
c
c     This routine (along with the associated subroutines and data file
c     RADTRANX_dbase.txt) computes the direct and diffuse spectral irradiances
c     incident onto the sea surface, using the model of Gregg and
c     Carder, Limnol & Oceanogr, 35(8), pp 1657-1675, 1990.
c
c     This code is a minor rewrite of the program "radtran", which was
c     kindly made available by Ken Carder.  Changes made in the rewrite
c     are as follows:
c
c    1.  main program "radtran" changed to subroutine "GCEd"
c     2.  input parameters are supplied by subroutine via
c        common block /cgcirr1/, rather than read from unit 5
c     3.  output irradiances are returned to gcirrad via common block
c         /cgcirr2/, rather than being written to a file
c     4.  double precision pi & rad values changed to single precision
c     5.  routine sunang is now on file qasky
c     6.  the data file gcirrad.txt has been extended to 800 nm, with
c         corresponding changes in the dimension statements (nlt = 452 now)
c
c     Input parameters needed by the Gregg & Carder model are as follows:
c
c     iein : = 0 to compute irradiance (W/m2/nm)
c            = 1 to compute quanta (uEin/m2/sec/nm)
c     iblw : = 0 to compute above-surface values
c            = 1 to compute below-surface values
c     jday : Julian day (year-day)
c     rlon,rlat : longitude (W negative), latitude (N positive)
c         Note: jday,rlon,rlat are used to compute the climatological
c         ozone value if ro3 = -99, and "the" if the = -99
c     the  : solar zenith angle; enter -99 to calculate from time and
c              earth position
c     hr   : time in hours GMT (= Eastern time plus 5) (not used if "the"
c            is given expicitly)
c     pres : pressure in inches mercury; enter 29.92 for standard
c     am   : air mass type (1-10); enter 1 for default (marine aerosols)
c     rh   : relative humidity; enter 80 for default
c     wv   : precipitable water (water vapor) in cm; enter 1.5 for default
c     wsm,ws : mean 24-hour wind speed, current wind speed (m/s); enter 
c          4.0 and 6.0 for defaults
c     vi   : visibility in km; enter 15 for default
c     ro3  : ozone in Dobson Units, or enter -99 for climatological ozone
c              values (computed from the day and Earth position)
c
      parameter (nlt=701)
      common /cgcirr1/ iein,iblw,jday,rlon,rlat,the,GMThr,pres,am,rh,
     1                 wv,wsm,ws,vi,ro3
      common /cgcirr2/ Edir(nlt),Edif(nlt),Ed(nlt)
      common /cmisc/ imisc(30),fmisc(30)
c
      integer lam(nlt)
      real Fobar(nlt),Fo(nlt),oza(nlt),ag(nlt),aw(nlt)
      real external get_o3_climat   !get_TOMS_zm
c
      data kall/0/
      save
c

c***** Initialize on first Kall
      If(kall.eq.0) then
          pi = fmisc(1)
          pi2 = 2.0*pi
          rad = fmisc(3)
c
c         Use user-specified total ozone value...
          if (ro3 .gt.0.0) then
           to3 = ro3
          else
c           ..or Climatology if ro3 is negative (flagged as -99)    
            write(10,504) jday,rlat,rlon
            to3 = get_o3_climat(jday, rlat, rlon)
          endif
c
c	Compute ozone scale height in cm
        sco3 = to3*1.0E-3
        if (the .eq. -99.0)then
         call sunang(jday,GMThr,rad,rlon,rlat,sunz,suna)
         theta = sunz
         the = sunz
        else
         theta = the
        endif
c
        write(10,600)'Parameter values used in the RADTRANX model:'
        write(10,601)'Solar zenith angle = ',theta,' degrees'
        write(10,601)'Pressure           = ',pres,' in. mercury'
      If(jday.lt.0) then
        write(10,604)'Day of year        =  ',
     1               'not set (mean earth-sun distance used)'
      Else
        write(10,603)'Day of year        = ',jday
      Endif
        write(10,601)'Air mass type      = ',am
        write(10,601)'Relative Humidity  = ',rh,' %'
        write(10,601)'Precipitable water = ',wv,' cm'
        write(10,601)'24-hr wind speed   = ',wsm,' m/sec'
        write(10,601)'Current wind speed = ',ws,' m/sec'
        write(10,601)'Visibility         = ',Vi,' km'
        write(10,601)'Total ozone        = ',to3,' Dobson units'
c
      if (iblw .eq. 0)then
      write(10,600)'Irradiances are computed just above the sea surface'
      else
      write(10,600)'Irradiance Just Below Surface'
      endif

        kall = 1
      endif
c
  504 format(/5x,'The ozone concentration is computed from ',
     4	'climatology using'/
     1	10x,'Day of Year = ',i4/
     2	10x,'latitude    = ',f9.4,' degrees (N is positive)'/
     3	10x,'longitude   = ',f9.4,' degrees (E is positive)')
 600  format(/5x, a)
 601  format(10x, a, f6.1, 1x,a)
 603  format(10x, a, i6, 1x,a)
 604  format(10x, 2a)
c
c***** later calls begin here

c  Read in extraterrestrial light (new format for RADTRAN
c  extended 300 to 1000 nm)
c
      call lidata_new(lam,Fobar,oza,ag,aw)
c
!  Correct for Earth-Sun distance
c     modify date if it is last day of a leap year
      if(jday.gt.365) jday=365.0
c
      If(jday.lt.0) then !use mean earth-sun distance values
        do nl = 1,nlt
          Fo(nl) = Fobar(nl)
        enddo
      Else               !correct for earth-sun distance based on jday
        do nl = 1,nlt
          Fo(nl) = Fobar(nl)*
     *     (1.0+1.67E-2*cos(pi2*(jday-3)/365.0))**2
        enddo
      Endif
c
c  Compute irradiances
c
      call atmodd(iblw,rad,lam,theta,oza,ag,aw,sco3,pres,wv,rh,am,
     * wsm,ws,Vi,Fo,Edir,Edif,Ed)
      if (iein .eq. 0)then   !radiant energy flux
       sumirr = 0.0
       do nl = 1,nlt
       sumirr = sumirr + Ed(nl)
       enddo
       write(10,*)
       if(kall.eq.0) write(10,602)'Total irradiance Ed(300-1000 nm) = ',
     1               sumirr,' W/m2'
 602  format(5x, a, f8.3, 1x,a)
      else    !quanta
       h = 6.6256E-34
       c = 2.998E8
       hc = 1.0/(h*c)
       sumirr = 0.0
       do nl = 1,nlt
        qlam = lam(nl)*1.0E-9
        Edir(nl) = Edir(nl)*qlam*hc/6.023E17
        Edif(nl) = Edif(nl)*qlam*hc/6.023E17
        Ed(nl) = Ed(nl)*qlam*hc/6.023E17
        sumirr = sumirr + Ed(nl)
       enddo
       write(10,602)' Total irradiance = ',sumirr,' microEinst/m2/sec'
      endif
c
      return
      end
c
c******************************************************************************
      subroutine lidata_new(lam,Fobar,oza,ag,aw)
c
c  Opens and reads light data for the radiative transfer model
c  of Gregg and Carder 1990.  This is the data of Table 1 in the L&O
c  paper, with values added to extend the model to 800 nm.c
c
c  New format for RADTRAN-X database (300 to 1000 nm) CDM May 2005
c
      parameter(nlt=701)
      character*50 title
      integer lam(nlt)
      real oza(nlt),ag(nlt),aw(nlt)
      real Hobar(nlt),Fobar(nlt)
c
c     **Common block shared only with setdflts.f (& listed in main.f)
c     **Contains the filename '..\data\RADTRANX_dbase.txt' written
c     **in a format agreeable with the system (UNIX/PC)
      COMMON/gcifname/gcifile
      character*32 gcifile
c
c  Light Data
c
      open(8,file=gcifile,status='old', form='formatted', err=666)
c  skip 10 header records
      do i=1,10
         read(8,10)title
c         write(10,fmt='(a30)') title
      enddo


      do nl = 1,nlt
       read(8,20)lam(nl),Hobar(nl),oza(nl),ag(nl),aw(nl)
c       write(10,fmt='(2i6,4f12.4)') nl,lam(nl),Hobar(nl),
c     1             oza(nl),ag(nl),aw(nl)
      enddo
      
      close(8)

      do nl = 1,nlt
       Fobar(nl) = Hobar(nl)*10.0   !convert W/cm2/um to W/m2/nm
      enddo
10    format(a50)
20    format(i4,1x,f9.5,f11.4,2f9.4)
c
      return
  666 call nofile(10, 'lidata_new', gcifile )   !err opening file
      end
c*******************************************************************************
      subroutine atmodd(iblw,rad,lam,theta,oza,ag,aw,sco3,p,wv,
     * rh,am,wsm,ws,vis,Fo,Edir,Edif,Ed)
c
c  Model for atmospheric transmittance of solar irradiance through
c  a maritime atmosphere.  Computes direct and diffuse separately.
c  Includes water vapor and oxygen absorption.
c
      parameter(nlt=701)

      common /cmisc/ imisc(30),fmisc(30) 

      integer lam(nlt)
      real Fo(nlt),oza(nlt),ag(nlt),Edir(nlt),Edif(nlt),Ed(nlt)
      real aw(nlt)
c      double precision rad
c
      data kall/0/
      save kall
c
      p0 = 29.92
	iprtout = imisc(9)
c
c  Compute atmospheric path lengths (air mass); pressure-corrected
      cosunz = cos(theta/rad)
c  Modified March, 1994 according to Kasten and Young 1989.
c      rex = -1.253
c      rtmp = (93.885-theta)**rex
c      rm = 1.0/(cosunz+0.15*rtmp)
      rex = -1.6364
      rtmp = (96.07995-theta)**rex
      rm = 1.0/(cosunz+0.50572*rtmp)
      rmp = p/p0*rm
      otmp = (cosunz*cosunz+44.0/6370.0)**0.5
      rmo = (1.0+22.0/6370.0)/otmp
c
c  Obtain aerosol parameters; simplified Navy aerosol model
      call navaer(rh,am,wsm,ws,vis,beta,alpha,wa,asymp)
c
      if(kall.eq.0) then
c     printout on first kall only and only if iprtout>0
        If(iprtout.gt.0) then
          write(10,600)'Aerosol parameters from the Navy aerosol model:'
          write(10,601)'alpha                    = ',alpha
          write(10,601)'beta                     = ',beta
          write(10,601)'asymmetry parameter      = ',asymp
          write(10,601)'single scattering albedo = ',wa
        Endif
      endif
 600	format(/5x, a)
 601  format(10x, a, f8.4, 1x,a)
c
      eta = -alpha
c   Forward scattering probability
      alg = alog(1.0-asymp)
      afs = alg*(1.459+alg*(.1595+alg*.4129))
      bfs = alg*(.0783+alg*(-.3824-alg*.5874))
      Fa = 1.0 - 0.5*exp((afs+bfs*cosunz)*cosunz)
c
c  Surface reflectance
      if (iblw .eq. 1)then
       call sfcrfl(rad,theta,ws,rod,ros)
      else
       rod = 0.0
       ros = 0.0
      endif
c
c  Compute spectral irradiance
      do nl = 1,nlt
c    Rayleigh, by Bird's method
       rlam = lam(nl)*1.0E-3
       tr = 1.0/(115.6406*rlam**4 - 1.335*rlam**2)
       rtra = exp(-tr*rmp)   !transmittance
c    Ozone
       to = oza(nl)*sco3   !optical thickness
       otra = exp(-to*rmo)   !transmittance
c   Aerosols
       ta = beta*rlam**eta
       if (lam(nl) .eq. 550)then
       if(kall.eq.0) write(10,602)
     1    'Aerosol optical thickness at 550 nm = ',ta
       endif
 602  format(/5x, a, f8.3, 1x,a)
       atra = exp(-ta*rm)
       taa = exp(-(1.0-wa)*ta*rm)
       tas = exp(-wa*ta*rm)
c   Oxygen/gases
       gtmp = (1.0 + 118.3*ag(nl)*rmp)**0.45
       gtmp2 = -1.41*ag(nl)*rmp
       gtra = exp(gtmp2/gtmp)
c   Water Vapor
       wtmp = (1.0+20.07*aw(nl)*wv*rm)**0.45
       wtmp2 = -0.2385*aw(nl)*wv*rm
       wtra = exp(wtmp2/wtmp)
c
c  Direct irradiance
       Edir(nl) = Fo(nl)*cosunz*rtra*otra*atra*gtra*wtra*(1.0-rod)
c
c   Diffuse irradiance
       dray = Fo(nl)*cosunz*gtra*wtra*otra*taa*0.5*
     *      (1.0-rtra**.95)
       daer = Fo(nl)*cosunz*gtra*wtra*otra*rtra**1.5*taa*Fa*
     *      (1.0-tas)
c
c  Total diffuse
       Edif(nl) = (dray + daer)*(1.0-ros)
c
       Ed(nl) = Edir(nl) + Edif(nl)
      enddo
c
      kall = 1
      return
      end
c**********************************************************************************
      subroutine navaer(rh,am,wsm,ws,vis,beta,alpha,wa,asymp)
c
c  Computes aerosol parameters according to a simplified version
c  of the Navy marine aerosol model.
c
      real a(3),ro(3),dndr(3),r(3)
      data ro/0.03,0.24,2.0/
      data r/0.1,1.0,10.0/
c
      rlam = 0.55
c
c  Relative humidity factor
      if (rh .ge. 100.0)rh = 99.9
      rnum = 2.0 - rh/100.0
      rden = 6.0*(1.0-rh/100.0)
      frh = (rnum/rden)**0.333
c
c  Size distribution amplitude components
      a(1) = 2000.0*am*am
      a(2) = 5.866*(wsm-2.2)
      if (a(2) .lt. 0.5)a(2) = 0.5
      a(3) = 0.01527*(ws-2.2)*0.05        !from Hughes 1987
      if (a(3) .lt. 1.4E-5)a(3) = 1.4E-5
c
c  Compute size distribution at three selected radii according to
c  Navy method
      do n = 1,3
       dndr(n) = 0.0
       do i = 1,3
        rden = frh*ro(i)
        arg = alog(r(n)/rden)*alog(r(n)/rden)
        rval = a(i)*exp(-arg)/frh
        dndr(n) = dndr(n) + rval
       enddo
      enddo
c
c  Least squares approximation
      sumx = 0.0
      sumy = 0.0
      sumxy = 0.0
      sumx2 = 0.0
      do n = 1,3
       sumx = sumx + alog10(r(n))
       sumy = sumy + alog10(dndr(n))
       sumxy = sumxy + alog10(r(n))*alog10(dndr(n))
       sumx2 = sumx2 + alog10(r(n))*alog10(r(n))
      enddo
      gama = sumxy/sumx2
      rlogc = sumy/3.0 - gama*sumx/3.0
      alpha = -(gama+3.0)
c
c  Compute beta
      cext = 3.91/vis
      beta = cext*rlam**alpha
c
c  Compute asymmetry parameter -- a function of alpha
      if (alpha .gt. 1.2)then
       asymp = 0.65
      else if (alpha .lt. 0.0)then
       asymp = 0.82
      else
       asymp = -0.14167*alpha + 0.82
      endif
c
c  Single scattering albedo at 550; function of RH
      wa = (-0.0032*am + 0.972)*exp(3.06E-4*RH)
c
      return
      end
c
c *****************************************************************************
      subroutine sfcrfl(rad,theta,ws,rod,ros)
c
c  Computes surface reflectance for direct (rod) and diffuse (ros)
c  components separately, as a function of theta, wind speed or
c  stress.
c
c      double precision rad
c
      rn = 1.341      !index of refraction of pure seawater
      roair = 1.2E3   !density of air g/m3
c
c  Foam and diffuse reflectance
      if (ws .gt. 4.0)then
       if (ws .le. 7.0)then
        cn = 6.2E-4 + 1.56E-3/ws
        rof = roair*cn*2.2E-5*ws*ws - 4.0E-4
       else
        cn = 0.49E-3 + 0.065E-3*ws
        rof = (roair*cn*4.5E-5 - 4.0E-5)*ws*ws
       endif
       rosps = 0.057
      else
       rof = 0.0
       rosps = 0.066
      endif
c
c  Direct
c   Fresnel reflectance for theta < 40, ws < 2 m/s
      if (theta .lt. 40.0 .or. ws .lt. 2.0)then
       if (theta .eq. 0.0)then
        rospd = 0.0211
       else
        rtheta = theta/rad
        sintr = sin(rtheta)/rn
        rthetar = asin(sintr)
        rmin = rtheta - rthetar
        rpls = rtheta + rthetar
        sinp = (sin(rmin)*sin(rmin))/(sin(rpls)*sin(rpls))
        tanp = (tan(rmin)*tan(rmin))/(tan(rpls)*tan(rpls))
        rospd = 0.5*(sinp + tanp)
       endif
      else
c   Empirical fit otherwise
       a = 0.0253
       b = -7.14E-4*ws + 0.0618
       rospd = a*exp(b*(theta-40.0))
      endif
c
c  Reflectance totals
      rod = rospd + rof
      ros = rosps + rof
c
      return
      end
