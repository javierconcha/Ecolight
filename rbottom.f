C     Last change:  LKS   7 Aug 2008    7:01 pm
      function rbottom(ibotm,wavel)
c
c     user routine on file Rbottom.f
c
C     called by BRDFbotm  [MAIN->RADAMPS->BOTMBC->qaBRRF->BRDFbotm->RBOTTOM]
c
c     This function returns the irradiance reflectance R = Eu/Ed for the
c     given bottom type and wavelength.  This routine is for use with
c     Lambertian bottom boundaries, as defined by the BRDFbotm on BRDFLamb.f 
c
c     if ibotm = 1 use the value of rflbot read from the input for all
c                  wavelengths
c     if ibotm = 2 read the bottom reflectance data from file Rbottomdatafile
c
c     This routine is designed for reading in bottom reflectance data,
c     which are assumed to be on the Hydrolight "Bottom-reflectance-data 
c     standard format" on file "Rbottomdatafile", where "Rbottomdatafile" 
c     is the path and filename as specified when running the front end 
c     program.  (Users can redefine that format by changing the statments 
c     in this routine that read the data file.)
c
c     When this routine is called for the first time with ibotm = 2, 
c     the data file is read and processed as follows:
c
c     1) the wavelength and corresponding reflectance value R is read for 
c        each wavelength where R was measured
c
c     2) these discrete-wavelength R values are linearly interpolated
c
c     On all subsequent calls, the spline is used to define R at any
c     wavelength that lies between the first and last wavelengths in the data
c     files.  Wavelengths less (greater) than the first (last) data 
c     wavelength are given the R value for the first (last) wavelength
c     measured. 
c
c     INPUT:
c        ibotm: flag for the bottom type
c        wavel: the wavelength in nm at which the reflectance is desired
c
c     OUTPUT:
c        rbottom:  the bottom reflectance R (nondimensional)
c
      INCLUDE "DIMENS_XL.INC"
c
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      common /Cmisc/ imisc(30),fmisc(30)
c
!     local variables
      integer kall, itype, nheadr, nh2p, ngrid
      real, allocatable,dimension(:) :: Rgrid, waveG
      real dum(1)      
      data kall/0/, waveold/-1./, wRef/-1/,itype/0/,nheadr/10/,nh2p/2/
c
c     Note on dynamic storage:  MUST save values between subroutine calls:
      save kall, ngrid, Rgrid
c
c-----  Initialization is needed on the first call only if the reflectance 
c       comes from a data file:
c
      if(kall.eq.0) then
c
        if(ibotm.eq.1) then
c         use the same reflectance at all wavelenths
          rbottom = fmisc(14)
          write(10,101) rbottom
        elseif(ibotm.gt.1) then
          write(10,100) datafiles(6)
c         open and read the data file
          ngrid = imisc(7)
          allocate(Rgrid(ngrid), waveG(ngrid))
          waveG(:ngrid) = wave(:ngrid)
          call getonGrid1(itype, nheadr,nh2p, ngrid,datafiles(6), 
     1                    wRef,aRef,waveG,Rgrid, dum)
        endif
        kall = 1
        return    !return after initializing
      endif
c-----  End of initialization on first call----------------------------
c
c-----  Subsequent calls start here:
      if(ibotm.eq.1) then
c        use the same reflectance at all wavelenths
         rbottom = fmisc(14)
         RETURN
      else
c       use the values gridded from data
        jwave = imisc(11)
        call checkiWav(jwave, wavel)
        rbottom = Rgrid(jwave)
        if(wavel .ne. waveold) then
          write(10,1030) wavel,rbottom
          waveold = wavel
        endif
        fmisc(14) = rbottom
        return
      endif
c
c     Formats:
  100 format(/5x,
     1'The bottom reflectance R is obtained from measured values read',
     2' from a file named:'/5x,a,/2x,
     4'(Linear interpolation is used to define R at all wavelengths)')
  101 format(/5x,
     1'The bottom reflectance R is contant at all wavelengths with ',
     2'the value: ',/10x, 'R = ',f6.3)
 1030 format(/5x,'The bottom reflectance at wavelength',f6.1,' nm is',
     1' R =', f6.3)
      end
