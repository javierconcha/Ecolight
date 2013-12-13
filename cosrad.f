C     Last change:  LKS  21 Jan 2008   10:22 am
      SUBROUTINE cosrad(suntheta,sunphi,skytheta,skyphi, skyrad)
c 
c     user "skyrad" routine on file cosrad.f 
c
c*****NOTE:  This "skyrad" model is useful for studying simple
c     idealized situations such as a black sky, a uniform sky, or a
c     heavy overcast.  In general, a more realistic sky model should
c     be used.
c
c     This routine returns a diffuse (background) sky radiance given by 
c
c        Lsky(theta) = Lo*(1 + C*cos(skytheta))      (Eq. 1)
c
c     Note that C = 0 gives a uniform background sky
c               C = 2 gives a cardioidal sky
c               C = 1.25 is a good approximation to a heavy overcast
c                   (see Brunger, A. P. and F. C. Hooper, Solar Energy
c                    51(1), 53-64, 1993 for discussion)
c
c     Integration of Lsky*cos(skytheta) gives the plane irradiance 
c     Ed(diffuse):
c
c        Ed(diffuse) = 2*pi*Lo*(1/2 + C/3)
c
c     which sets the value of Lo, given Ed(diffuse) = Rsky*Ed(total)
c     from the user input.  
c
c     CALLING ARGUMENTS: 
c     INPUT:
c     suntheta and sunphi give the sun's position (in radians, relative to 
c          phi = 0.0 in the downwind direction and thetas = 0.0 being
c          the zenith direction)
c     skytheta and skyphi give the sky direction (in radians, relative to 
c          phi = 0.0 in the downwind direction and thetas = 0.0 being
c          the zenith direction) where the radiance is desired
c     NOTE: suntheta, sunphi, and skyphi are not used in this routine
c
c     OUTPUT:
c     skyrad is the sky radiance in direction (skytheta, skyphi)
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Csky/ iskyflag,skyspecs(mxnsky)
      common /cmisc/ imisc(30),fmisc(30)
c
      save rad0, C, kall 
c
      data kall/0/
c
c     other input for the cosine sky radiance distribution:
c
c     C is the cosine weight factor
c
c     rsky is the ratio of diffuse sky to total (sky + sun) plane 
c          irradiance:  rsky = Ed(diffuse)/(Ed(diffuse) + Ed(direct))
c          rsky = 0.0 for a black sky (sun only), rsky = 1.0 for
c          no sun visible 
c
c     Edtotal is the total (sky + sun) plane irradiance on the
c          water surface from above 
c 
c     initialization on the first call:
c
      if(kall.eq.0) then
c
         C   = skyspecs(3)
         rsky   = skyspecs(4)
         Edtotl = skyspecs(5)
c
         pi = fmisc(1) 
c 
         write(10,100) rsky,C,Edtotl 
         rad0 = rsky*Edtotl/(2.0*pi*(0.5 + C/3.0)) 
         kall = 1
      endif
c
c     subsequent calls start here:
c
      skyrad = rad0*(1.0 + C*cos(skytheta))
c 
      return
c 
c     format
c 
  100 format(/5x,'Normalized sky radiances are computed using sky',
     1' model "COSRAD"'/8x,'(COSine RADiance model)'//,
     28x,'COSRAD parameters:'/
     310x,'Rsky = ratio of sky to total plane irradiance =',0p,f6.3/
     410x,'C = cosine parameter =',f6.3/
     510x,'Edtotl = total (sun + sky) Ed in W/(m^2 nm) =',f6.3)
      end
