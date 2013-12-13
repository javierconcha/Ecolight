C     Last change:  LKS  21 Jan 2008   10:22 am
      SUBROUTINE cosirrad(suntheta,sunphi, Eddif,Eddir)
c 
c     user "skyirrad" routine on file cosirrad.f 
c
c     called by routine qasky
c
c*****NOTE:  This simple "skyirrad" model is intended for use with the 
c     simple "skyrad" model on file "cosrad.f" for single-wavelength
c     runs.  (This routine returns the same values of Eddif and Eddir 
c     at all wavelengths, and therefore should not be used for multiple 
c     wavelength runs.)
c
c     CALLING ARGUMENTS: 
c     INPUT:
c     suntheta and sunphi give the sun's position (in radians, relative 
c          to phi = 0.0 in the downwind direction and thetas = 0.0 being
c          the zenith direction)
c     NOTE: suntheta and sunphi are not used in this routine
c
c     OUTPUT:
c
c     Eddif is the diffuse (background sky) downwelling plane irradiance
c         incident onto the sea surface
c     Eddir is the direct (solar beam) downwelling plane irradiance
c         incident onto the sea surface
c
      INCLUDE "DIMENS_XL.INC"
c
      common /Csky/  iskyflag,skyspecs(mxnsky)
c
      save kall
c
      data kall/0/
c
c     user input:
c     rsky is the ratio of diffuse sky to total (sky + sun) plane 
c          irradiance:  rsky = Ed(diffuse)/(Ed(diffuse) + Ed(direct))
c          rsky = 0.0 for a black sky (sun only), rsky = 1.0 for
c          no sun visible
c     Edtotl is the total (sky + sun) plane irradiance incident onto
c          the sea surface 
      rsky   = skyspecs(4)
      Edtotl = skyspecs(5)
c
      Eddif = rsky*Edtotl
      Eddir = (1.0 - rsky)*Edtotl
c
      if(kall.eq.0) then
c     printout on the first call:
c
         write(10,100) Eddif,Eddir,Edtotl,rsky 
c
         kall = 1
      endif
c
      return
c 
c     format
c 
  100 format(/5x'Sky irradiances are computed using sky irradiance ',
     9'model "COSIRRAD"',
     1/8x,'(COSine IRRADiances; for use with the COSRAD model)'//10x,
     2'diffuse (sky) spectral plane irradiance Ed(diffuse) =  '
     3,  1p,e10.3,' W/(m^2 nm)'/ 
     410x,'direct (sun) spectral plane irradiance Ed(direct) =    ',
     5  1p,e10.3,' W/(m^2 nm)'/ 
     610x,'total (sun + sky) spectral plane irradiance Ed(total) =',
     71p,e10.3,' W/(m^2 nm)'/ 
     810x,'ratio of sky to total plane irradiance, Rsky =',0p,f6.3)
      end
