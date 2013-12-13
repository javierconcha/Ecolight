C     Last change:  LKS   7 Aug 2008    6:49 pm
      SUBROUTINE PNTGRID
c     
c     core routine on file pntgrid.f
c
c     called by INITIAL [MAIN->INITIAL->PNTGRID]
c
c     This routine prints out the (theta,depth,wavelength) grid values
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
!      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
      character*10 blank, line(8)
      data blank/'          '/
c
      nmu = imisc(1)
      nz = imisc(4)
      iop = imisc(5)
      nwave = imisc(7)
      radeg = fmisc(3)
c
c     write headers
c
      write(10,90)
c
      imax = max0(nmu,nz,nwave)
      do 100 i=1,imax
c
c     blank the line
c
         do j=1,8
            line(j) = blank
         end do
c
c     insert values, where they exist, into the blank line
c
         if(i.le.nmu) then
            write(line(1),fmt='(i10)') i
c
c           get exact theta values at cell centers
            if(i.eq.1) then
               theta = 0.5*(90.0 + radeg*acos(bndmu(1)))
            elseif(i.eq.nmu) then
               theta = 0.
            else
               theta = 0.5*radeg*(acos(bndmu(i-1)) + acos(bndmu(i)))
            endif
c
            write(line(2),fmt='(f10.2)') theta
            write(line(3),fmt='(f10.4)') fmu(i)
         endif

         if(i.le.nz) then
            write(line(4),fmt='(i10)') i
c           zeta contains the depths as read in, optical or geometric
            if(iop.eq.1) write(line(5),fmt='(f10.3)') zeta(i)
            if(iop.eq.0) write(line(6),fmt='(f10.3)') zeta(i)
         endif

         if(i.le.nwave) then
            write(line(7),fmt='(i10)') i
            write(line(8),fmt='(f10.2)') wave(i)
         endif
c
c     print the line
c
      write(10,fmt='(8a10)') line
 100  continue
      return
c
 90   format(/2x,'The AZIMUTHALLY AVERAGED radiance distribution is ',
     1'computed on the following (Theta, Depth, Wavelength) grid:'//
     2 9x,'I',5x,'Theta',6x,'Mu',11x,'K',6x,
     3 'zeta',7x,'z',11x,'L',4x,'Lambda'/
     4 15x,'(deg)',36x,'(m)',16x,'(nm)'/)
      end
