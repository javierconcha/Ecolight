C     Last change:  LKS   7 Aug 2008    7:28 pm
      Subroutine pureH2O(wavel, a,b,c)
c
c     core routine on file pureH2O.f
c
c     called by "ab" routine
c
c     This routine returns the absorption (a), scattering (b), and beam
c     attenuation (c) coefficients for "pure water", for use in IOP
c     models.  
c
c     On the first call, values of a and b at discrete wavelengths 
c     are read from file "pureH2O.txt".  "pureH2O.txt" can be any file
c     of "pure water" values.  The choice of file "SBH2Oab.txt" 
c     containing "Smith and Baker" data or file "PFH2Oab.txt" containing
c     "Pope and Fry" data is made in subroutine "setdflts.f".
c     The discrete wavelength a and b values are linearly interpolated, 
c     which are used in subsequent calls to define a, b, and c = a+b 
c     at any wavelength from 200 to 800 nm.
c
c     INPUT:
c        wavel: the wavelength in nm, where a, b, and c are requested
c
c     OUTPUT:
c        a:  the pure water absorption coef in 1/m
c        b:  the pure sea water scattering coef in 1/m
c        c:  the attenuation coef = a + b
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cmisc/ imisc(30),fmisc(30) 
c
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      integer kall, itype, nheadr, nh2p
      real, allocatable,dimension(:) :: agrid, bgrid
c
      data kall/0/, wRef/-1/, itype/2/, nheadr/10/, nh2p/2/
c     Note on dynamic storage:  MUST save values between subroutine calls:
      save kall, agrid, bgrid
c
c-----  Initialization on the first call ------------------------------
c
      if(kall.eq.0) then
        write(10,100) datafiles(0)
c
c       open and read the data file
        ngrid = imisc(7)
        allocate(agrid(ngrid), bgrid(ngrid))
        call getonGrid1(itype, nheadr,nh2p, ngrid,datafiles(0), 
     1                wRef,aRef,wave,agrid, bgrid)
        kall = 1
      endif
c-----  End of initialization on first call----------------------------
c
c-----  Subsequent calls start here:
      jwave = imisc(11)     !index to current wavelen
      call checkiWav(jwave, wavel)
      a = agrid(jwave)
      b = bgrid(jwave)
      c = a + b
      return
c
c     Formats:
  100 format(/5x,
     1'Pure water a and b values are obtained from the file:'/
     28x,a,/8x,
     3'(Linear interpolation is used to define values at all',
     4' wavelengths)')
      end

