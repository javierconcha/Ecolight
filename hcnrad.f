C     Last change:  LKS  25 Apr 2008    8:21 pm
      subroutine hcnrad(suntheta,sunphi,skytheta,skyphi, radnorm)
c
c     user "skyrad" routine on file hcnrad.f
c
c     called by qasky
c
c     "hcnrad" = Harrison and Coombes Normalized RADiances
c
c     This is the default "skyrad" routine; see the "skyrad.txt" file
c
c     This routine computes the normalized sky radiance by the semiempirical
c     model of Harrison and Coombes, Solar Energy, 41(4), pp 387-392, 1988.
c
c     CALLING ARGUMENTS:
C     INPUT
c	suntheta,sunphi = the (theta,phi) position of the sun
c	                 
c	skymu,skyphi = the (theta,phi) direction where the sky radiance 
c                    is to be computed
c
c     all angles are in radians, i.e.
c          0.0 .le. suntheta,skytheta .le. pi/2
c          0.0 .le. sunphi,skyphi .le. 2pi
c
c     OUTPUT:
c	radnorm = the normalized sky radiance, with units of 1/steradian, seen
c	  	  in direction (skytheta,skyphi).  This value must be multiplied
c		  by Ed(sky,lambda)/pi [units of W/(m**2 nm)] in order to
c		  convert it into an absolute spectral radiance; see Eq. (6)
c 		  of Harrison and Coombes.
c
c     psi is the angle between (suntheta,sunphi) and (skytheta,skyphi)
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Csky/ iskyflag,skyspecs(mxnsky)
c
      data kall/0/
      save kall
c
c     cloud is the cloud fraction (0.0 for a clear sky to 1.0 for a
c	a fully overcast sky)
c
      cloud = skyspecs(10)
      if(kall.eq.0) then
c     printout on first call
         write(10,100) cloud
c     data quality check on first call
      if(cloud.lt.0) then
        cloud = 0
        write(10,*)
     1    'WARNING: CLOUD fraction less than zero!  Set to zero.'
        skyspecs(10) = cloud
      elseif(cloud.gt.1) then
        cloud = 1
        write(10,*)
     1    'WARNING: CLOUD fraction greater than one!  Set to one.'
        skyspecs(10) = cloud
      endif
c
c     This model cannot be evaluated for the sun exactly at the horizon.
c     Move the sun slightly if user selected a sun angle of 90 degrees
         IF(suntheta.ge.2.*atan(1.)) THEN
           suntheta=(89./90.)*suntheta
           skyspecs(1) = 45./atan(1.)*suntheta
           write(10,'(/5x,2a)')'The Harrison and Coombes sky radiance',
     &           ' model cannot be evaluated at sunset.'
           write(10,'(5x,2a,F4.1,a)')'The solar zenith angle has been',
     &           ' adjusted to a position of ', 45./atan(1.)*suntheta,
     &           ' degrees.'
         endif
         kall = 1
      endif
c
      sunmu = cos(suntheta)
      skymu = cos(skytheta)
      cospsi = sunmu*skymu + sqrt(1.0 - sunmu*sunmu) * 
     1         sqrt(1.0 - skymu*skymu) * cos(sunphi - skyphi)
      psi = acos(cospsi)
c
c
c     compute the normalized radiance by Eqs (3)-(5) of Harrison and Coombes
c
      radclear = (1.63 + 53.7*exp(-5.49*psi) +
     1           2.04*cospsi*cospsi*sunmu) * (1.0 - exp(-0.19/skymu)) *
     2           (1.0 - exp(-0.53/sunmu))
      radcloud = 0.45 + 0.12*acos(sunmu) + 0.43*skymu +
     1           0.72*exp(-1.88*psi)
c
      radnorm = cloud*radcloud + (1.0 - cloud)*radclear
c
      return
c
  100 format(/5x,'Normalized sky radiances are computed using sky',
     1' model "HCNRAD" (Harrison and Coombes Normalized RADiances)'/
     28x,'The cloud fraction is cloud =',f5.1,' (0 for a clear sky to ',
     3'1 for a heavy overcast)')
      end
