C     Last change:  LKS   7 Aug 2008    7:09 pm
      function bstar(icomp, wavenm)
c
c     user routine on file bstar.f   
c
c     called by ABCASE1, and ABCASE2
c
c     calls LENSTR, and interpolation routines (in dataintrp.f).
c
c     This function defines the specific scattering coefficient,
c     b*(wavelength) in mg/m^3, for use by the "ab" routines.
c     If the power or linear models are used, bstar actually returns the
c     depth-independent, wavelenght-dependent part of b(waven,z) which may
c     not strictly be the specific scattering coefficient (i.e., if Xn is NOT 1,
c     bstar will not return the true "b*")
c
c     The b* specifications are:
c         ibstropt = -1:  Component is non-scattering (e.g., CDOM) b* = 0
c                     0:  user-supplied data file read to get b* values
c                     1:  Power Law used (including Gordon-Morel model)
c                     2:  Linear model (including Gould, et al. model)
c                     3:  Constant value, indep of wavelength or conc
c                     4:  Power Law used to calculate c (b=c-a)  !new in HE5
c
c             bstarRef =  reference wavelength for model
c             bstar0   =  bstar at reference wavelength
c             PLm      =  m for Power law
c             GAMm     =  m of GAM, the slope of the linear model
c             GAMi     =  i of GAM, the offset value in the linear model
c             Xn       =  power coefficient for the concentration (GAMn or PLn)
c
c     For "datafile" or PSM options, on the first call, this routine
c     reads in the data file which contains the b* data, and 
c     puts the data on the run wavel grid.  Subsequent calls return
c     b* at the given wavelength.
c
c     For the Power Law and Linear model options (CHL), the routines return b*
c
c
      INCLUDE "DIMENS_XL.INC"

c     Common blocks containing b* specifications for all components
      Common /Cbstar/ ibstropt, bstarRef, bstar0, CompN, Plm,
     1                GAMm, GAMi
      Common /Cbstar2/ bstarfile
!
      common/cmisc/ imisc(30),fmisc(30)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 

      Integer ibstropt(mxcomp)
      Real bstarRef(mxcomp), bstar0(mxcomp), CompN(mxcomp)
      Real PLm(mxcomp), GAMm(mxcomp),GAMi(mxcomp)
      Character*120 bstarfile(mxcomp+2)
c
      Character*120 datadir, digitdir, spreadir,
     1              phasedir, surfdir,bottdir, Pdir
      COMMON /Cdirnames/ datadir,digitdir,spreadir,
     1                   phasedir,surfdir,bottdir, Pdir

      real wavenm

!     local
      integer kall(mxcomp), nheadr, nh2p
      real bstargrid(mxwave,mxcomp+2), b0(mxcomp), zgrid(1)
      real, allocatable,dimension(:) :: bgrid
!
      data kall/mxcomp*0/,waveold/0/ 
      data itype/0/, nheadr/10/, nh2p/0/
!     itype=0 means abscissa is wavelen; nheadr is # header lines in file 
!     and nh2p is the number of header lines to be printed out
c
c     Note on dynamic storage:  MUST save data values between subroutine 
c     calls:
      save
c
      if(kall(icomp).eq.0) then
c     initialization on first call
c
        If(ibstropt(icomp) .eq. 0) then       !data file
c-------- INITIALIZE DATA FILE OPTION ----------
          write(10,100) icomp, bstarfile(icomp)
          ngrid = imisc(7)
          allocate(bgrid(ngrid))
          call getonGrid1(itype, nheadr,nh2p,ngrid,bstarfile(icomp), 
     1                    bstarRef(icomp),b0(icomp),wave,bgrid, zgrid)
          bstargrid(:ngrid, icomp) = bgrid(:ngrid)
          deallocate(bgrid)
!
        Elseif(ibstropt(icomp).eq. 1) then
c-------- INITIALIZE Power Law model OPTION ----------
          write(10,101) icomp,bstar0(icomp),CompN(icomp),
     1                  bstarRef(icomp),PLm(icomp)
c
        Elseif(ibstropt(icomp).eq. 2) then
c-------- INITIALIZE Linear model OPTION ----------
          write(10,102) icomp,bstar0(icomp),CompN(icomp),
     1                  bstarRef(icomp),GAMm(icomp), GAMi(icomp)
c
        Elseif(ibstropt(icomp).eq. 3) then
c-------- INITIALIZE CONSTANT OPTION ----------
          write(10,103) icomp,bstar0(icomp)
c
        Elseif(ibstropt(icomp).eq. 4) then
c-------- INITIALIZE Power Law model OPTION for computing C ----------
          write(10,104) icomp,bstar0(icomp),CompN(icomp),
     1                  bstarRef(icomp),PLm(icomp)
c
        Elseif(ibstropt(icomp).gt. 3 .or. ibstropt(icomp).lt.-1) then
          write(10,*) 'bstar: component #',icomp
          write(10,*) '       Invalid option', ibstropt(icomp)
          bstar = 1.0
          return
        Endif
        kall(icomp) = 1
      endif
c-----End of initialization on first call.
c
c     Subsequent calls start here:
       IF (waveold.ne.wavenm) THEN
            waveold = wavenm
       END IF
c
      If(ibstropt(icomp) .eq. 0) then       !data file will be read
c-------- EVALUATE DATA FILE OPTION ----------
        If(abs(wavenm-bstarRef(icomp)).le.1e-4) then
          bstar = b0(icomp)
        Else
          jwave = imisc(11)     !index to current wavelen
          call checkiWav(jwave, wavenm)
          bstar = bstargrid(jwave,icomp)
        Endif
!
      Elseif(ibstropt(icomp).eq.1 .or. ibstropt(icomp).eq.4) then
c-------- EVALUATE Power Law model OPTION ----------
          bstar = bstar0(icomp)*(bstarRef(icomp)/wavenm)**PLm(icomp)
c
      Elseif(ibstropt(icomp).eq.2) then
c-------- EVALUATE Linear model OPTION ----------
          bstar = bstar0(icomp)*(GAMm(icomp)*wavenm + GAMi(icomp))/
     1    (GAMm(icomp)*bstarRef(icomp) + GAMi(icomp))       
c
      Elseif(ibstropt(icomp).eq.3) then
          bstar = bstar0(icomp)
c-------- Return Dummy values for some components  ----------
      Elseif(ibstropt(icomp).eq.-1) then
          bstar= 0.0
      ELSE
          bstar= 1.0
      ENDIF

      return
 100  FORMAT(/5x,'Specific scattering for component',
     1 i2, ' read in from file: ',a)
 101  FORMAT(/5x,'Scattering for component',
     1 i2, ' calculated using a Power Law model',
     2 /10x, 'b = ',f6.4,' X^(',f6.4,') * (',f4.0,'/wavenm)^(',f6.4,')')
 102  FORMAT(/5x,'Scattering for component',
     1 i2, ' calculated using a Linear model:',
     2 /10x, 'b = ',f6.4,' X^(',f6.4,
     3	') * (m*wavenm + i)/(m*wave0 + i), where',
     4 /10x, 'wave0 = ',f4.0, 2x,'m = ',f6.4,2x,' and i = ',f6.4)
 103  FORMAT(/5x,'Scattering for component',
     1 i2, ' assumed constant:',
     2 /10x, 'b = ',f6.4)
 104  FORMAT(/5x,'Scattering for component',
     1 i2, ' calculated from b = c - a, using a Power Law model',
     2 /10x, 'c = ',f6.4,' X^(',f6.4,') * (',f4.0,'/wavenm)^(',f6.4,')')
      end


