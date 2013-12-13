C     Last change:  CDM Jan 2012 to reformat output
      subroutine writerad
c 
c     core routine on file writerad.f
c
c     called from radanal.f only if iOptRad = 1
c
c     This routine writes a file (Lroot.txt) with the full radiance 
c     distribution on a format designed for ease of use in plotting
c     with IDL (this file is read by IDL routine read_E_Lfile.pro).

c     a depth of -1.0 labels values in air, just above the sea surface
c
c     Changes made Aug 2013 for v 5.2:
c     total in-air radiances are also decomposed to azimuthal averages of
c     L_sky, L_w, and L_sr

c     called by RADANAL [MAIN->RADANAL->WRITERAD]
c 
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Crad/ RADsky(mxmu),RADwla(mxmu),RADrsa(mxmu),RADupa(mxmu),
     1              RADupz(mxmu,mxz),RADdnz(mxmu,mxz)
c     radiance arrays are stored as follows:
c        RADsky(mu) = total incident sky radiance
c        RADwla(mu) = Lw(mu,z=a) = water-leaving radiance
c        RADrsa(mu) = Lreflsky(mu,z=a) = surface-reflected sky rad
c        RADupa(mu) = Lu(mu,z=a) = total upwelling rad in air
c           note: RADupa = RADwla + RADrsa, i.e., Lu = Lw + Lreflsky in air
c        RADupz(mu,z) = Lup(mu,z) = total upwelling rad in water
c        RADdnz(mu,z) = Ldn(mu,z) = total downwelling rad in water
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      COMMON /Cmisc/ imisc(30),fmisc(30)
      COMMON /Ctitle/ ititle
      Character ititle*120
      common /Cradfile/ Lrootname!, Erootname
      character Lrootname*120 !, Erootname*120
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)
c
      COMMON /Cpirrad/ ipirad,izirad(mxz)    !print out on irrad grid

      data kall/0/

c     on the first call, compute the directions of the quad
c     centers and write the file headers
c 
      dimension thedeg(mxmu)
      save
c
      nmu = imisc(1)
      jwave = imisc(11)
      nwave = imisc(7) 
      nuLfile = imisc(15)
      radeg = fmisc(3)

c
      if(kall.eq.0) then
c        get theta in degrees
         thedeg(1) = 0.5*(90.0 + radeg*acos(bndmu(1))) 
         do i=2,nmu-1
            thedeg(i) = 0.5*radeg*(acos(bndmu(i-1)) + acos(bndmu(i)))
         end do
         thedeg(nmu) = 0.

c        open the output file and write headers
      open(nuLfile,file=Lrootname, status = 'unknown')
c
      write(nuLfile,fmt='(2a)') 'ECOLIGHT Run Title: ',
     1                          ititle(1:lenstr(ititle))
      write(nuLfile,fmt='(a)') 'NOTE: This file is formatted for reading
     1 by IDL routine read_E_Lfile.pro v5.2'
      write(nuLfile,fmt='("Parameters: nz, ntheta, nlambda = ",
     1  3i5)')  ipirad+1,2*nmu,nwave
      write(nuLfile,150)

         kall = 1
      endif
c-----End of initialization on first call.
c
c     Subsequent calls start here:

c     write radiances in the air (at z = a = -1.0)
      depth = -1.0
      wavel = fmisc(13)
c
c     downwelling radiances (theta = 0 to 90)
      do iu=nmu,1,-1
         rad =RADsky(iu)
         write(nuLfile,152) depth,thedeg(iu),wavel,rad,rad,0.0,0.0
      enddo

c     upwelling radiances (theta = 90 to 180)
      do iu=1,nmu
          write(nuLfile,152) depth,180.0-thedeg(iu),wavel,RADupa(iu),
     * 0.0,RADwla(iu),RADrsa(iu)
      enddo

c     write radiances in the water (at z = user requested values)
c 
      do iiz=1,ipirad
         iz = izirad(iiz)
         depth = zgeo(iz)
c     downwelling radiances (theta = 0 to 90)
!poles      do iu=nmu,nmu
      do iu=nmu,1,-1
            rad = RADdnz(iu,iz)
            write(nuLfile,152) depth,thedeg(iu),wavel,rad
         enddo
c     upwelling radiances (theta = 90 to 180)
!poles      do iu=nmu,nmu
      do iu=1,nmu
            rad = RADupz(iu,iz)
            write(nuLfile,152) depth,180.0-thedeg(iu),wavel,rad
         enddo
      enddo     ! end izz loop
      return
c 
  150 format('Azimuthally averaged radiances L(z,theta,lambda) in W m^-2
     * sr^-1 nm^-1'/
     *'theta and phi are the directions of photon travel:'/
     *'theta = 0 to 90 is downwelling; theta = 90 to 180 is upwelling; (
     *theta = 180 is nadir-viewing; theta = 0 is zenith-viewing)'/
     *'depth = -1.0 labels values in air (just above the sea surface)'/
     *'L_sky is incident sky radiance; theta = 0 to 90 deg only; in air
     *only'/
     *'L_w is water-leaving radiance; theta = 90 to 180 deg only; in air
     * only'/
     *'L_sr is surface-reflected radiance; theta = 90 to 180 deg only; i
     *n air only'/
     *'   depth   theta  lambda  total radiance    L_sky         L_w
     *        L_sr'/
     *'    [m]    [deg]   [nm]   [W/(m^2 sr nm)]' )
  152 format(f8.2,2f8.1,4es15.5e3) !allows for 1.23456E-100
      END
