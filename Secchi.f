C     Last change:  LKS   7 Nov 2008    7:55 pm
C     MJM Version, 2008 October 13, 1355 EDT
C     MJM My modifications include robustness and slight rearrangement. 
C     MJM The strategy is now to calculate RWP.46b and RWP.38
C     MJM at each output depth. We re-write RWP.55  to be 
C     MJM zSD(alpha+K) = Gamma and realize that alpha and K each have
C     MJM a zSD in their denominator. Therefore the product zSD(alpha+K) has
C     MJM no explicit zSD. So, we calculate the LHS at each grid point (z instead of zsd). By definition,
C     MJM at zero depth it must be zero. It increases until at some grid point it is 
C     MJM greater than 0. This first grid point where it is greater than 0 and the 
C     MJM previous one now bracket zSD, where the equality is satisfied. I now abandon 
C     MJM all finesse and use bisection to find the root. 
C     MJM NOTE THAT: it may still fail if there is a lot of bioluminescence or no radiance from above the water surface.
C     MJM I have clearly assumed certain behavior for alpha and fkz. Give the canonical "secchi depth" problem, 
C     MJM these are probably OK. However, there should probably be a safety valve somewhere to cease after some 
C     MJM number of iterations so HL/EL can exit.
      subroutine Secchi
c
c     core routine on file Secchi.f
c
c     written by CDM, 10 Nov 04
c
c     This routine estimates the Secchi depth z_SD after
c     a Hydrolight run has been completed.
c     It is intended for use only with multi-wavelength runs
c     that span the visible spectrum from 400 to 700 nm (at least).
c
c     Calculations are based on R. W. Presendorfer, Secchi Disk Science,
c     Limnol. Oceanogr. 31(5), 909-926, 1986 (ref as RWP below)
c
      INCLUDE "DIMENS_XL.INC"
c
c     Cxcl holds arrays accumulated in routine storexcl
      COMMON /Cxcl/  nzxcl,zxcl(mxz),izxcl(mxz),
     1               axcl(mxz,mxwave,0:mxcomp),
     2               bxcl(mxz,mxwave,0:mxcomp),
     3               bbxcl(mxz,mxwave,0:mxcomp),
     4               Edxcl(0:mxz,mxwave),Euxcl(0:mxz,mxwave),
     5               Eodxcl(0:mxz,mxwave),Eouxcl(0:mxz,mxwave),
     5               Rxcl(0:mxz,mxwave),
     6               Raduxcl(0:mxz,mxwave),Radwxcl(mxwave)
c
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave)

      COMMON /CMISC/ IMISC(30),FMISC(30)

c     The CIE 1931 standard photopic luminosity function ybar(lambda)
c     Ref: L&W Table 2.1
      integer, parameter :: nCIE=40
      real, parameter :: fKM = 683.0
      real :: CIEwave(nCIE),CIEybar(nCIE)
      data CIEwave /380.0,390.0, 400.0, 410.0, 420.0, 430.0, 440.0,
     1  450.0, 460.0, 470.0, 480.0, 490.0, 500.0, 510.0, 520.0, 530.0,
     2  540.0, 550.0, 560.0, 570.0, 580.0, 590.0, 600.0, 610.0, 620.0,
     3  630.0, 640.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, 710.0,
     4  720.0, 730.0, 740.0, 750.0, 760.0, 777.0/

      data CIEybar /  0.0000,0.0001,0.0004,0.0012,0.0040,0.0116,0.0230,
     1  0.0380,0.0600,0.0910,0.1390,0.2080,0.3230,0.5030,0.7100,0.8620,
     2  0.9540,0.9950,0.9950,0.9520,0.8700,0.7570,0.6310,0.5030,0.3810,
     3  0.2650,0.1750,0.1070,0.0610,0.0320,0.0170,0.0082,0.0041,0.0021,
     4  0.0011,0.0005,0.0003,0.0001,0.0001,0.0000/

c     CIE ybar values for use with current wavelengths
      real :: ybarwave(mxwave), wl
c     PHOTOPIC beamc(z), Ed(z) and Kd(z)
      real :: alphaz(mxz),Edvz(mxz)
c MJM variables
      real, dimension(mxz)::alsumz,fkdz,zfunc
      integer::izlast,is_greater
      real:: shallow, deep, mid, v_s, v_d, v_m, a_s, a_d, a_m
      real:: f_s, f_d, f_m

c     The Secchi coupling constant Gamma of RWP Eq. (55).
c     A value of 8.0 corresponds to a disk of R = 0.85 photopic reflectance
c     seen against a background ocean with Rinf = 0.02, assuming that the
c     disk become invisible when the apparent contrast is about 0.015,
c     and when there is no surface effect (Script T = 1.0 in RWP).
c     See RWP Table 2 for Gamma values for other conditions.
      real, parameter :: Gamma = 8.0

!MJM      print*,'In secchi'
      nwave = imisc(7)

      i1 = 1
      izlast=0
      is_greater=0
      do iw=1,nwave
         wl = wave(iw)
         if(wl.le.CIEwave(1) .or. wl.gt.CIEwave(nCIE)) then
            ybarwave(iw) = 0.0
         else
            ybarwave(iw) = yinterp(i1, nCIE, wl, CIEwave, CIEybar)
         endif
      enddo
c
c     compute alpha(z), the photopic beam attenuation at the user
c     printout depths, using RWP Eq. (21) with the upwelling radiance Lu
c     (i.e., we are looking straight down)
c     compute the illuminance Edvz(z) (E in RWP) by RWP (15), using
c     the irradiance Ed(z)
c      
      do iz=1,nzxcl
         radsum = 0.0
         alsum = 0.0
         Edsum = 0.0
         do iwave=1,nwave
c       the denominator of RWP (21): (MJM: why not use a function integrator?)
           radsum = radsum + raduxcl(iz,iwave)*ybarwave(iwave)*
     1           (waveb(iwave+1)-waveb(iwave))
c       the numerator of RWP (21); c = a + b
           alsum = alsum + raduxcl(iz,iwave)*(axcl(iz,iwave,0) +
     1  bxcl(iz,iwave,0))*ybarwave(iwave)*(waveb(iwave+1)-waveb(iwave))
c       the illuminance by RWP (15)
           Edsum = Edsum + Edxcl(iz,iwave)*ybarwave(iwave)*
     1           (waveb(iwave+1)-waveb(iwave))
         end do
      
         alphaz(iz) = alsum/radsum  ! RWP 21
         Edvz(iz) = fKM*Edsum       ! RWP 15
         
!MJM     Now, accumulate the sums of a modified RWP 55         
         if (iz.eq.1) then ! depth = 0
            alsumz(iz)=0.d0 
            fkdz(iz)  =0.d0
         else 
! note that we really should call qtrap (Numerical recipes) or some 
! equivalent function integrator instead of just using this trapezoidal sum...
            alsumz(iz)=alsumz(iz-1)+ 0.5*(alphaz(iz) + alphaz(iz-1)) *       &
     1         (zxcl(iz)-zxcl(iz-1))
            fkdz(iz)=-1.0*alog(Edvz(iz)/Edvz(1))
         endif
!MJM     by definition, at the shallowest (0) depth, alsumz(1)+fkdz(1)=0.
!MJM     at some point it will exceed gamma - at that point we stop, and the
!MJM     secchi depth is between the first index where alsumz(iz)+fkdz(iz) is 
!MJM     > 0 and the previous one. 
         izlast=iz
         if (alsumz(iz)+fkdz(iz).ge.gamma) then
!MJM Set is_greater=1 iff we transition from < to > 0
            is_greater=1 
            exit ! loop
         endif 
         
      end do
!diag      write(10,*)'alphaz'
!diag      write(10,'(i5,3e20.8)')(iz,alphaz(iz),alsumz(iz),fkdz(iz),         &
!diag     1           iz=1,izlast)
!diag      write(10,*)'finished secchi part 1'
c

      if (is_greater.eq.0) then
!         write(10,*)"The current estimate of zSD is greater than "//     &
!     1       "the lowest output depth of ",zxcl(nzxcl)
         write(10,105) zxcl(nzxcl)
         return
      endif

c     get alpha(z_SD) and Kd(z_SD) from the spline interpolations
!MJM      alsd =  yinterp(i1, nzxcl, zSD, zxcl, alphaz)
!MJM      Edvsd = yinterp(i1, nzxcl, zSD, zxcl, Edvz)
!MJM  Now, since we are assuming a linear interpolation, then we know the
!MJM  know the form; however, the assumption is linear interp
!MJM  in Ed, not in fkdz. So, since we know the root is required
!MJM  to be in the known interval, throw away all finesse and 
!MJM  bisect:
      shallow= zxcl(izlast-1)
      deep   = zxcl(izlast)
      a_s = alsumz(izlast-1)
      a_d = alsumz(izlast)
      f_s = fkdz(izlast-1)
      f_d = fkdz(izlast)
      v_s = a_s + f_s - gamma ! < 0
      v_d = a_d + f_d - gamma ! > 0
      i1=izlast-1
!     epsilon=0.0001         !lks     relative eps
      epsilon=0.002         !lks     absolute eps to give 0.01m resolution
      do ;
         mid    =(shallow+deep)*0.5                   !middle depth for bisection
         alsd =  yinterp(i1, nzxcl, mid, zxcl, alphaz) ! interpolate alphaz
         alsumzz=alsumz(izlast-1)+0.5*(alphaz(izlast-1) + alsd)*                 &
     1             (mid-zxcl(izlast-1))    ! integrate from ilast-1 to mid
         Edvsd = yinterp(i1, nzxcl, mid, zxcl, Edvz) ! interpolate edvz
         fkzz = -1.0*alog(Edvsd/Edvz(1))             ! value at mid
         v_m = alsumzz + fkzz - gamma                ! "merit" function
!         write(10,*)'shallow, mid, deep = ',shallow, mid, deep
!         write(10,*)'v_s, v_m, v_d = ', v_s, v_m, v_d
!         write(10,*)shallow, gamma/(a_s/shallow + f_s/shallow)
!         write(10,*)deep, gamma/(a_d/deep + f_d/deep)
!         write(10,*)''
         if (v_m .le.0) then ! pull in shallow bound
            shallow=mid
            a_s = alsumzz
            f_s=fkzz
            v_s=v_m
         else                ! pull in deep  bound
            deep=mid
            a_d=alsumzz
            f_d=fkzz
            v_d=v_m
         endif
         val=Gamma/(alsumzz/mid + fkzz/mid) ! maybe test on val instead of distance between hi/low ?
!         if ((deep-shallow)/mid.le.epsilon) exit !loop, converged
         if (abs(deep-shallow).le.epsilon) exit !loop, converged
         
      enddo
      zSD=mid
!      write(10,*)'deep   = ',deep
!      write(10,*)'mid    = ',mid
!      write(10,*)'shalow = ',shallow
      val=Gamma/(alsumzz/mid + fkzz/mid)
!      write(10,'(a,2f9.4)')'Zsd,val = ',Zsd, val
      write(10,303)alsumzz/mid,fkzz/mid
      write(10,304)Zsd
!      write(10,*)'val=',val



!MJM      print*,'Done with Secchi'

  105 format(/5x,' WARNING: the current estimate for z_SD is greater tha
     1n the maximum output depth for the run: z_SD cannot be further rev
     2ised',/,
!     35x,' z_SD(estimate) =',f8.2,'    max output depth =',f8.2/
     35x,'    max output depth =',f8.2/
     45x,' For a better estimate, rerun HydroLight to a deeper depth')
  201 format(5x,' The Secchi depth estimated from the surface alpha and'
     1,' K is z_SD =',f8.2,' m')
 200  format(//' Secchi depth calculations:',//,5x,
     1' The near-surface photopic alpha and K are alpha =',f8.4,
     2' 1/m and K =',f8.4,' 1/m')
303   format(/,5x,
     1' The final depth-averaged photopic alpha and K are alpha =',f8.4,
     2' 1/m and K =',f8.4,' 1/m')
  304 format(5x,
     1' The final estimate for the Secchi depth is z_SD =',
     2f8.2,' m')

c
      return
      end


