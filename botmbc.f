C     Last change:  LKS   1 Nov 2007    1:01 pm
      subroutine botmbc
c 
c     core routine on file botmbc.f
c 
c     This routine specifies the bottom boundary condition via the 
c     quad averaged rmb = r(m,b) array for the selected type of bottom
c     boundary, for the present wavelength.
c
c     if ibotm = 0 (infinitely deep water).  Set up and solve the eigenvalue
c                   problem of Eq. (9.50) and then use (9.76) to get
c                   r(z,infinity)
c
c     if ibotm .ge. 1 (finite depth at zeta = zeta(nz)), The bottom is
c                   assumed to be Lambertian.  (This is not a requirement
c                   of EcoLight, which can handle any BRDF in principle.
c                   A Lambertian bottom is used only to speed up the
c                   calculations.)
c
c          if ibotm = 1, the same reflectance R is used at each wavelength
c          if ibotm = 2, the wavelength-dependent reflectance R is obtained
c                   from Hydrolight standard format data file Rbottomdatafile
c
      INCLUDE "DIMENS_XL.INC"
c 
      COMMON /CBOTBC/ Rmb(mxmu,mxmu)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      common /CMISC/ imisc(30),fmisc(30)
 
      save
c 
      nmu = imisc(1)
      ibotm = imisc(12)
      pi = fmisc(1)
      wavel = fmisc(13)
c
      if(ibotm.eq.0) then
c
c        Infinitely deep water:
c        set up and solve the eigenvalue problem (9.50) and obtain
c        Rinfinity = r(m,infinity) by (9.76).
c        infbotm sets rmb = r(m,inf)
c 
         call infbotm
c
      elseif(ibotm.ge.1) then 
c
c        Finite depth water.  The quad-averaged BRRF is given
c        by Eq. 8.109.

c        Get the bottom irradiance reflectance
         Rbot = rbottom(ibotm,wavel)
c
c        define the r(m,b) array
         do ir=1,nmu
         value = (Rbot/pi)*fmu(ir)*Omega(ir)
           do iu=1,nmu
             rmb(ir,iu) = value
           end do
         end do
c
      endif

      return
      end