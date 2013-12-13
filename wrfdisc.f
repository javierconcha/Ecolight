C     Last change:  LKS   9 Dec 2011   10:44 am
      subroutine wrfdisc
c
c     core routine on file wrfdisc.f
c
c     called by MAIN
c
c     calls eval1D (in file dataintp.f) and p2aray
c
c     This routine discretizes the wavelength redistribution functions for
c     chlorophyll and CDOM fluorescence, and for Raman scattering.
c
c     The double integrations (CM 23 June 94) are evaluated using a wavelength
c     resolution of deltaw (= 1.0 nm) within each band
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave)
      COMMON /CsourceE/ SrcM(mxmu),SrcP(mxmu)
      Common /Csource0/ ibiolum,ichlfl,icdomfl,iraman, RamanEXP
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
      nwave = imisc(7)
!      idbug = imisc(9)
c
      idbug = 0
      if(idbug.ne.0) then
      write(10,500)
      do i=1,nwave
      write(10,502) i,waveb(i),waveb(i+1),waveb(i+1)-waveb(i)
      end do
      endif
c
c     Discretize the continuous wrf's for the desired inelastic 
c     scatter components.  Evaluate only the wrf's that are needed in
c     the run.
c
      deltaw = 1.0
      do j=2,nwave
         wj = waveb(j) + 0.5*deltaw
         nj = ifix((waveb(j+1) - waveb(j))/deltaw)
         factor = deltaw*deltaw/(waveb(j+1) - waveb(j))
         do i=1,j-1
            wi = waveb(i) + 0.5*deltaw
            ni = ifix((waveb(i+1) - waveb(i))/deltaw)
c
      if(ichlfl.ne.0) then
c
c     chlorophyll fluorescence wrf:
      sumchl = 0.0
      do jj=1,nj
         wavej = wj + float(jj-1)*deltaw
         do ii=1,ni
            wavei = wi + float(ii-1)*deltaw
            sumchl = sumchl + wrfchl(wavei,wavej)
         end do
      end do

      fijchl(i,j) = factor*sumchl
      endif
c
      if(icdomfl.ne.0) then
c
c     CDOM fluorescence wrf:
      sumcdom = 0.0
      do jj=1,nj
         wavej = wj + float(jj-1)*deltaw
         do ii=1,ni
            wavei = wi + float(ii-1)*deltaw
            sumcdom = sumcdom + wrfcdom(wavei,wavej)
         end do
      end do
      fijcdom(i,j) = factor*sumcdom
      endif
c
      if(iraman.ne.0) then
c
c     Raman scattering wrf:
      sumraman = 0.0
      do jj=1,nj
         wavej = wj + float(jj-1)*deltaw
         do ii=1,ni
            wavei = wi + float(ii-1)*deltaw
            sumraman = sumraman + wrframan(wavei,wavej)
         end do
      end do
      fijraman(i,j) = factor*sumraman
      endif
c
      end do   ! i loop
      end do   ! j
c
!      idbug=1
      if(idbug.ne.0) then
      if(ichlfl.ne.0) call p2aray(fijchl,nwave,nwave,mxwave,2,
     1  ' fijchl as discretized')
      if(icdomfl.ne.0) call p2aray(fijcdom,nwave,nwave,mxwave,2,
     1  ' fijcdom as discretized')
      if(iraman.ne.0) call p2aray(fijraman,nwave,nwave,mxwave,2,
     1  ' fijraman as discretized')
      endif
c   
      return
c
c     formats
c
  500 format(///2x,'Wavelength redistribution functions for inelastic',
     1' scatter are discretized over the following wavelength bands:'//
     2'      band     from       to     bandwidth'/
     315x,'(nm)',6x,'(nm)',6x,'(nm)'/)
  502 format(6x,i3,f11.1,f10.1,f9.1)
      end

