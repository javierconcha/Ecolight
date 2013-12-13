C     Last change:  LKS  31 Aug 2013    6:41 pm
      SUBROUTINE INITIAL
c
c     user routine on file initial.f
c
c     This subroutine 
c
c        1) sets default values for various parameters
c        2) reads the ascii file of runtime input created
c           by the front end program or by a text editor(file 
c           Iroot.txt, which is assumed to be associated with 
c           Fortran Unit 5)
c        3) reads the ascii file of surface information for the
c           desired wind speed (file surfwind._U_, which was
c           created by a surfcode run)
c        4) selects the type and amount of printout for file
c           Proot.txt)
c        5) prints the initial information to file Proot.txt (unit 10)
c
c     called by MAIN
c
c     calls SETDFLTS, PNTGRID, LENSTR 
c
      INCLUDE "DIMENS_XL.INC"
c
c     temporary local storage:
      real Parmin, Parmax, PhiChl, Raman0, RamanXS
      real S, T, refr
      integer lenpfdir, lenbdir, lendata
      character*120 Prootname
      character*120 txtdatafile, CDOMdatafile, HSdatafile
      character rootname*32, sl*2, bsl*2
      Character ac9datafile*120,chlzdatafile*120,Rbottomdatafile*120
c     *****Variables for Check date and time *****
      character*10 tdate, zone,time, txtdate
      character*5 txttime 
      integer iv(8)

      DIMENSION skydata(mxnsky)
c
      common /CindexChl/ indexchl
      integer indexchl
c
      common /Cconc/  itype(mxcomp)
      common /Cbbopt/ ibbopt(mxcomp), bbfrac(mxcomp), 
     1                BfrefPL(mxcomp), Bf0PL(mxcomp), BfmPL(mxcomp)
      Common /Castar/ iastropt, astarRef, astar0, asgamma
      Common /Castar2/ astarfile
      Common /Cbstar/ ibstropt, bstarRef, bstar0, CompN, Plm,
     1                GAMm, GAMi
      Common /Cbstar2/ bstarfile
      Common /Cabac9n/ nac9files
c
c     Specific ab variables
      integer ibbopt, itype
      real bbfrac, BfrefPL, Bf0PL, BfmPL
      Integer iastropt(mxcomp+2), ibstropt(mxcomp)
      Real astarRef(mxcomp+2), astar0(mxcomp+2)
      Real asgamma(mxcomp+2)
      Real bstarRef(mxcomp), bstar0(mxcomp)
      Real CompN(mxcomp), PLm(mxcomp)
      Real GAMm(mxcomp),GAMi(mxcomp)
      Real coef1, coef2, coef3
      Character*120 astarfile(mxcomp+2), bstarfile(mxcomp+2)
      Character*120 ac9Filterfile, H2Odatafile

      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      character Lrootname*120
      common /Cradfile/ Lrootname
c     for use in PAR:
      character Erootname*120
      common /Efilename/ Erootname
C
C     Common blocks used or defined here
      COMMON /Csky/ iskyflag,skyspecs(mxnsky)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      Common /Csource0/ ibiolum,ichlfl,icdomfl,iraman, RamanEXP
      Common /CSourceBdata/  S0datafile    !read in initial, used in s0bdata
        character*120 S0datafile
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
c     common block for communication with the Gregg & Carder sky
c     irradiance model (routine gcirrad):
      common /cgcirr1/ iein,iblw,jday,rlon,rlat,the,GMThr,pres,am,rh,
     1                 wv,wsm,ws,vi,ro3
      real pres,am,rh,wv,vi,wsm,ro3

      COMMON /Ctitle/ ititle
      character ititle*120
      Character*120 datadir, digitdir, spreadir,
     1              phasedir, surfdir,bottdir, Pdir
      COMMON /Cdirnames/ datadir,digitdir,spreadir,
     1                   phasedir,surfdir,bottdir, Pdir
c
c     Blocks new in H5
      COMMON /CmodelOpts/ iparam(0:5)
      COMMON /CConstConc/ compconc(mxcomp)
      INTEGER iparam, iDynZ
      REAL compconc
c     
c        All the variables in /FRSTCLS/ are set = 1 when calling 
c        from the driver.  They allow the routine of that name
c        to initialize itself on the first call, but not on sub-
c        sequent calls.  This is necessary when the programs are
c        called for many wavelengths.
      COMMON /Cfrstcls/ iabscat,iqasky,iradamps,iradanal, iradxcl
c 
c     Set first-call parameters:
      data areset/0.0/, breset/0.0/
c
      iabscat = 1
      iqasky = 1
      iradamps = 1
      iradanal = 1
      iradxcl  = 1
c
c     Set various defaults (including the paths for various
c     directories):
      CALL SETDFLTS
c
      lendatadir = lenstr(datadir)
      lendigitdir = lenstr(digitdir)
      lenspreadir = lenstr(spreadir)
      lenPdir = lenstr(Pdir)
c
cccccc  BEGIN READING RUN-TIME INPUT  cccccccccccccccccccccccccccccccc
c
c     See Appendix A of the Users' Guide for additional documentation
c
c     Set the system slash or backslash to a string variable to aid in 
c     generating filenames
      sl='/ '
      sl=sl(1:1)
      bsl='\ '
      bsl=bsl(1:1)
C

c
c.....RECORD 1:  Run params and title, up to 120 characters
c
      read(5,*)iparam(0),Parmin,Parmax,PhiChl,Raman0,RamanXS,iDynZ,
     1         RamanEXP     
c     iparam(0) read by run.exe to det if need to compile (H5)
c       if iparam(0) = 0  <-- use STANDARD pre-compiled EXE
c       if iparam(0) = 1  <-- use USER routines in re-compiled EXE
c
c.....RECORD 2:  The run title
      read(5,fmt='(a)') ititle 
c     new in HE5 -- append run start time
c     ***Call internal system routine to get date and time
c     Optional time stamp routine for Lahey Fortran only
      call DATE_AND_TIME(tdate,time,zone,iv)
c     date: ccyymmdd  time: hhmmss.ss
c     iv:  1-ccyy, 2-mm, 3-dd, 4-UTC, 5-hh, 6-min, 7-sec, 8-msec
      txtdate = tdate(5:6)//char(47)//tdate(7:8)//char(47)//tdate(1:4)
      txttime = time(1:2)// ':' // time(3:4) 
c
      ititle = trim(ititle) // ' (' // txtdate // ' ' // txttime // ')'
c
c-----  Set wavelength range for PAR calculations  ---------------
      imisc(27) = iDynZ      
c-----  RAMAN SCATTERING COEFFICIENT  ----------------------------
c     Set the value of the Raman scattering coefficient for
c     water at some reference wavelength; these are the values 
c     of a_oR and lambda'_o in L&W Eq. 5.89.  (used in routine
c     shatram on file shat.f)
      fmisc(21) = RamanXS
      fmisc(22) = Raman0
!     RamanEXP is stored in CSource0
c-----  CHLOROPHYLL FLUORESCENCE QUANTUM EFFICIENCY  ----------
c     Set the value of the Chlorophyll fluorescence quantum
c     efficiency; this is the value of PhiC in L&W Eq. 5.102.
c     (used in routine wrfchl on file wrfdisc)
      fmisc(23) = PhiChl
c-----  CDOM FLUORESCENCE FUNCTION  ------------------------
c     The parameters defining the CDOM etaY function of L&W Eq. 5.101
c     are set in routine wrfcdom on file wrfdisc.f

c
c.....RECORD 3:  The root name for generating file names
c
      read(5,fmt='(a)') rootname 
c
c     generate the output file names from the default directory
c     names and the run-time file root name
      lenroot = lenstr(rootname)
c
c     the file for full digital output:
      Drootname = digitdir(1:lendigitdir) // 'D' // 
     1  rootname(1:lenroot) // '.txt'
c
c     the file for spreadsheet output on the single-wavelength format:
      Srootname = spreadir(1:lenspreadir) // 'S' //
     1  rootname(1:lenroot) // '.txt'
c
c     the file for spreadsheet output on the multiple-wavelength format:
      Mrootname = spreadir(1:lenspreadir) // 'M' //
     1  rootname(1:lenroot) // '.txt'
c
c     the file for general run printout:
      Prootname = Pdir(1:lenPdir) // 'P' //
     1  rootname(1:lenroot) // '.txt'

c     the file for the full radiance distribution:
      Lrootname = digitdir(1:lendigitdir) // 'L' //
     1  rootname(1:lenroot) // '.txt'
c
!c     the file for scalar irradiance and PAR output:
!      Erootname = Pdir(1:lenPdir) // 'E' //
!     1  rootname(1:lenroot) // '.txt'

c     Open unit 10 as printout file as begin writing output
      open(unit=10,file=Prootname, status='unknown')
c
c     write out Copyright info as header (routine at bottom of file)
      call cpright(10)
      write(10,fmt='(//"  RUN TITLE:  ",a)') ititle

c     Optional time stamp routine for Lahey Fortran only
      iout = 10
      call date_stamp(iout,0, 0, 0)  ! 0 indicated "start" of run
c
c.....RECORD 4a: Oputput options

c-----  PRINTOUT FOR DEBUGGING
c     Set idbug = -1, 0, 1, or 2 for increasing amounts of output for
c     debugging.  This is best done in specific routines where additional
c     output is desired:
c		-1:			minimal output produced
c		 0:			standard output produced
c		 2:			extensive output produced
c
c----- OUTPUT FOR SPREADSHEET ANALYSIS 
c     Set flags for writing output to files for Excel spreadsheet 
c     postprocessing
c
c     Flag for single-wavelength-format spreadsheet data for
c     multiwavelength runs:
c         iwrtSS1 = 0 if no single-wavelength-format spreadsheet file
c                     (file Srootname) is to be written
c         iwrtSS1 = 1 if a single-wavelength-format spreadsheet file
c                     is to be written for the first wavelength only
c         iwrtSS1 = 2 if the single-wavelength-format spreadsheet file
c                     is to be appended with a set of single-wavelength-
c                     format data for each wavelength in a multi-
c                     wavelength run
c
c     Flag for multi-wavelength-format spreadsheet data: 
c         iwrtSSM = 0 if no multi-wavelength-format spreadsheet file
c                     (file Mrootname) is to be written
c         iwrtSSM = 1 if a multi-wavelength-format spreadsheet file
c                     is to be written
c
c     Flag for radiance output data file (only used with HL, not EL): 
c         iOptRad = 0 if no radiance output data file is to be written
c         iOptRad = 1 if a radiance output data file is to be written
c
c     NwSkip signals whether to solve alternating wavelengths or not
c            =1 for every wavelength; = 2 for every other wavelength
c
c----- OUTPUT FOR RADIANCE ANALYSIS
c     Set flags for writing radiance in tabular form to a special output
c     file in the printout directory
      read(5,*) idbug, idigital, iwrtSS1, iwrtSSM, iOptRad, nwskip

c.....RECORD 4b: Model options
c     Set run params which control how incfiles_std handles run
c      iparam(1) = IOP model
c      iparam(2) = sky model
c      iparam(3) = irradat (iff=1, get incident irrad from datafile)
c      iparam(4) = chl component index
c      iparam(5) = CDOM component index
      read(5,*) (iparam(i),i=1,5)    !specifies models

c.....RECORD 5a: Number of components in the IOP model
c     ncomp is the number of scattering and absorbing components used in 
c     the "ab" routine.  The same number of previously discretized phase 
c     functions will be read from files 
c          pfname(i),i=1,ncomp

      read(5,*) ncomp, nconc
c
c     check for array dimensions being larger than allowed
      if(max(ncomp,nconc).gt.mxcomp) then
         write(6,fmt='(" ncomp =",i3," gt mxcomp =",i3)') ncomp,mxcomp
         write(10,fmt='(" ncomp =",i3," gt mxcomp =",i3)') ncomp,mxcomp
         call HERR("initial","increase model component LIMIT in UI")  !stop run
      endif
c
c.....RECORD 4c:  Component Concentrations
c     Ncomp is the number of components in the chosen "ab" model,
c     Nconc is the number of "ab" model components (each with its own dpf) PLUS
c     additional absorbing components that may be needed for fluorescence
c     calculations (if CDOM and/or CHL  are not part of the "ab" model being used)
c
      read(5,*) (compconc(i),i=1,nconc)  ! Constant component concentrations
c

      if(nwskip.le.1) nwskip=1  !make sure value is appropriate

      imisc(9)  = idbug			!sets how much output is put in proot file
      imisc(20) = iwrtSS1			!sets whether Sroot file is written
      imisc(21) = iwrtSSM			!sets whether Mroot file is written
      imisc(22) = idigital	              !sets whether Droot file is written
      imisc(24) = iOptRad                !sets whether Lroot file is written
c
c.....RECORD 5: IOP Specification
c
c.....RECORD 5b:  Specific Absorption parameters
c     For each component, a pair of  (a*,b*) data is read and corresponding 
c     pairs of (a*,b*) data file names are read.
c
c     Concentration and absorption specs are read for every nconc, to allow
c     for cases when fluoresence is included for a non-component 
c     (such as chl-fluorescence with the abacbb model)
c
c     The component concentration is specified by: 
c        itype = 0:  Constant component concentration profile
c                1:  User-subroutine called to give conc(z)
c                2:  Data file will be read to give conc(z)
c
c     The a* specifications are:
c         iastropt = -1:  Subroutine called for aCDOM (a* routine not called)
c                    -2:  aCDOM solved as a percentage of aw+achl in "ab" routine (a* not called)
c			 0:  user-supplied data file read to get a* values
c                     1:  Pope & Fry absorption model used (pure water only)
c                     2:  Smith and Baker absorption model used (pure water only)
c                     3:  Prieur-Sathyenranath-Morel model used (Chlorophyll only)
c                     4:  Exponential model used (CDOM only)
c                     5:  Prieur-Sathyendranath-Morel model used (CDOM only)
c
c         The other a* params are only used when iastropt=4 (exp model)
c         astarRef: reference wavelength for exp model
c         astar0:   astar at the reference wavelength
c         asgamma:  exponential decay constant
c
      Do i=1,nconc
          read(5,*) itype(i),
     1              iastropt(i), astarRef(i), astar0(i), asgamma(i)
      Enddo
c
c.....RECORD 5c: Specific Absorption filenames
c     Now read astar data file names to be used iff iastropt = 0
      Do ii=1, nconc
          read(5,'(a)') astarfile(ii)
c
          lendata = lenstr(astarfile(ii))
          iabspath = 0
          do i=1,lendata
              if(astarfile(ii)(i:i).eq.sl(1:1) .or.
     1            astarfile(ii)(i:i).eq.bsl(1:1)) iabspath = 1
          end do
          if(iabspath.eq.0) then
              astarfile(ii) = datadir(1:lendatadir) //
     1                     astarfile(ii)(1:lendata)
          endif
      Enddo
c
c.....RECORD 5d: Specific Scattering parameters
c     The b* specifications are:
c         ibstropt =  0:  user-supplied data file read to get b* values
c                     1:  Power Law used (Chlorophyll only)
c                     2:  Linear model (Chlorophyll only)
c                     3:  Constant value, indep of wavelength or conc
c                     4:  Power Law used to calculate c (b=c-a)  !new in HE5
c         Iff ibstropt = 1, 2, or 4 up to four more parameters are used (otherwise, -999)
c             bstarRef = reference wavelength for model
c             bstar0   = bstar at reference wavelength
c         If ibstropt = 1 or 4, 
c                          coef1 = m for Power law
c                          coef2 = n of PL, the power of the compoenent concentration, CompN
c         If ibstropt = 2, coef1 = m of GAM, the slope of the linear model, and
c                          coef2 = i of GAM, the offset value in the linear model
c                          coef3 = n of GAM, the power of the component concentration, CompN
c     CompN is the power of the Concentration in the model for b
c
      Do i=1,ncomp
          read(5,*) ibstropt(i), bstarRef(i), bstar0(i),
     1              coef1, coef2, coef3

          If(ibstropt(i).eq.1 .or. ibstropt(i).eq.4) then
              PLm(i)=coef1
              compN(i)=coef2
          Elseif(ibstropt(i).eq.2) then
              GAMm(i) = coef1
              GAMi(i) = coef2
              compN(i) = coef3
          Else
              compN(i) = 1.0
          Endif
      Enddo
c
c     Now read bstar data file names to be used iff ibstropt = 0
c
c
c.....RECORD 5e: Specific Scattering filenames
      Do ii=1, ncomp
          read(5,'(a)') bstarfile(ii)
c
          lendata = lenstr(bstarfile(ii))
          iabspath = 0
          do i=1,lendata
              if(bstarfile(ii)(i:i).eq.sl(1:1) .or. 
     1            bstarfile(ii)(i:i).eq.bsl(1:1)) iabspath = 1
          end do
          if(iabspath.eq.0) then
              bstarfile(ii) = datadir(1:lendatadir) // 
     1                     bstarfile(ii)(1:lendata)
          endif
      Enddo

c
c.....RECORD 5f:  phase function specifiers and filenames
c     read the component phase function specifiers
c
c     For each component these values are read
c     ibbopt= 0:  Phase function file explicitly specified (use pfname)
c             1:  Select phase function to match specified bb/b ratio
c             2:  bb/b specified by Hydroscat data (only for use w/ ac-9 model)
c             3:  Bf = bb/b given by Power Law
c
c     bbfrac: bb/b ratio used iff ibbopt=1
c
c     Bf*PL:      Power Law coefs for bb/b, wave_ref, Bf_0, and m, respectively
c     
      do i=1,ncomp
          read(5,*) ibbopt(i), bbfrac(i),
     1              BfrefPL(i), Bf0PL(i), BfmPL(i)
      enddo

c     Now read the names of the files containing the discretized phase 
c     functions corresponding to the IOP components used in the "ab" 
c     routine.  After reading each file name, concatinate it with the 
c     data directory name to create the full name for use in routine 
c     radamps, where these files are opened and read.
c
c     Note:  file name used if and only if ibbopt=0, otherwise dummy variable
      do i=1,ncomp
         read(5,fmt='(a)') pfname(i)
         lenpfname = lenstr(pfname(i))
         lenpfdir = lenstr(phasedir)
         pfname(i) = phasedir(1:lenpfdir) // 
     1               pfname(i)(1:lenpfname)
      end do
c
c.....RECORD GROUP 6:  Wavelength information
c
c     nwave = the number of wavelength BANDS at which the model is being run.
c
c        if nwave = 0, the next record gives the EXACT WAVELENGTH in nm, and values
c           of a and b that can be used to override the values given by the
c           "ab" routine.  For use of the a and b values as given by "ab", read in
c           negative values for areset and breset.  The sky spectral radiance at
c           the exact wavelength will be used (at 1 nm resolution).
c        if nwave .ge. 1, the next record gives the nwave+1 WAVELENGTH BAND
c           BOUNDARIES (in nm) for which the model is to be run.  a and b values as
c           returned by the "ab" routine at the band center will be used. 
c           The band-averaged sky radiance will be used.
c
      read(5,*) nwave
c
c     check for array dimensions being larger than allowed
      if(nwave.gt.mxwave) then
         write(6,fmt='(" nwave =",i3," gt mxwave =",i3)') nwave,mxwave
         write(10,fmt='(" nwave =",i3," gt mxwave =",i3)') nwave,mxwave
         call HERR("initial","increase wavelength LIMIT in UI")  !stop run
      endif
c
      if(nwave.eq.0) then
c     set up a 1 nm wide wavelength band
         read(5,*) wave(1),areset,breset
         nwave = 1
         waveb(1) = wave(1) - 0.5
         waveb(2) = wave(1) + 0.5
      else
         areset = -1.0
         breset = -1.0
         read(5,*) (waveb(iw),iw=1,nwave+1)
c        compute the nominal wavelengths (band center wavelengths) for this run
         do iw = 1,nwave
            wave(iw) = 0.5*(waveb(iw) + waveb(iw+1))
         end do
      endif
c
c.....RECORD 7:  Internal source and inelastic scatter flags
c
c     ibiolum is a flag for the inclusion/omission of bioluminescence:
c        if ibiolum = 0, there is no bioluminescence present
c        if ibiolum = 1, the run includes bioluminescence read from data file
c        if ibiolum = 2, the run includes bioluminescence from call to function
c     ichlfl is a flag for the inclusion/omission of chlorophyll fluorescence:
c        if ichlfl = 0, the is no chlorophyll fluorescence present
c        if ichlfl = 1, chlorophyll fluorescence is present
c     icdomfl is a flag for the inclusion/omission of CDOM fluorescence:
c        if icdomfl = 0, the is no CDOM fluorescence present
c        if icdomfl = 1, CDOM fluorescence is present
c     iraman is a flag for the inclusion/omission of Raman scattering:
c        if iraman = 0, the is no Raman scattering present
c        if iraman = 1, Raman scattering is present
c
c     ... and read the indexes to components containing CHL for routine SHAT
      read(5,*) ibiolum,ichlfl,icdomfl,iraman, indexchl
c           
c*****Note:  if isource = 0, the source-free RTE is solved (which runs
c     faster).  If bioluminescence or inelastic scatter are present, then
c     the appropriate source terms and equations are be included in the
c     solution of the full RTE.
      isource = ibiolum+ichlfl+icdomfl+iraman
c*****Note:  a single-wavelength run can include bioluminescence, but if
c     fluorescence or Raman scatter are included, then the run MUST 
c     include all contributing wavelengths less that the largest 
c     requested wavelength.
c
c.....RECORD 8:  Sky information
c
c     Parameters for computing the sky radiance distribution.
c
c     The format of this record depends on the sky model selected:
c
c     iskyflag = 1 if the "idealized" sky models are being used
c              = 2 for "realistic" sky models, with solar zenith 
c                  angle being specified
c              = 3 for "realistic" sky models, with time and location 
c                  being specified
c
c     Depending on iskyflag, various quantities are read into skydata,
c     for storage in skyspecs (for use in the skyrad and skyirrad 
c     models) as follows:
c
c     skyspecs(1) = suntheta = solar zenith angle (degrees)
c             (2) = sunphi = solar azimuthal angle (degrees, measured 
c                   counterclockwise from phi = 0 in the downwind 
c                   direction; sunphi = 0 by default)
c             (3) = C = cosine radiance parameter
c             (4) = Rsky = Ed(diffuse)/Ed(total) ratio
c             (5) = Edtotal = total (sky + sun) incident spectral 
c                   irradiance (W/(m^2 nm)) at the single wavelength
c             (6) = jday = Julian day (1.0 = Jan 1)
c             (7) = rlat = latitude in degrees (+ is north, - is south)
c             (8) = rlon = longitude in degrees (+ is east, - is west)
c             (9) = GMThr = GMT in hours (Pacific Std time + 8 hours)
c                   with minutes as a fraction, e.g. 21.5 for 21 hours 
c                   30 minutes
c     skyspecs(10) = cloud = cloud fraction (0.0 for a clear sky to 1.0 
c                    for a heavy overcast)
c
c	   pres     	sea-level atmos pressure (inches Hg)
c	   am		aerosol model type (1-10)
c	   rh		relative humidity (%)
c	   wv		water vapor content (cm)
c	   vi		average visibility (km)
c	   wsm		24-hour average windspeed (m/s)
c	   ro3		total ozone (dobson)
c

c*****NOTE:  other parameters required for the sky radiance computations
c     are set to default values in SETDFLTS, but could be read in here 
c     if desired 
c
c     "nsky" is the number of skydata values to be read in this record.
c
      read(5,*) iskyflag,nsky,(skydata(ii),ii=1,nsky)  
      read(5,*) fjday,rlat,rlon, pres,am,rh,wv,vi,wsm,ro3      !new in H5
c
c     store values used to set mean earth-sun distance and ozone 
c     (not used if iskyflag = 1)
      skyspecs(6) = fjday
      skyspecs(7) = rlat
      skyspecs(8) = rlon
      jday = ifix(fjday)   !temp fix

c     store data according to iskyflag:
c
      if(iskyflag .eq. 1) then
         suntheta = skydata(1)
         sunphi = skydata(2)
         C = skydata(3)
         Rsky = skydata(4)
         Edtotal = skydata(5)
c        store input in skyspecs:
         skyspecs(1) = suntheta
         skyspecs(2) = sunphi
         skyspecs(3) = C
         skyspecs(4) = Rsky
         skyspecs(5) = Edtotal
c
      elseif(iskyflag .eq. 2) then
         suntheta = skydata(1)
         sunphi = skydata(2)
         cloud = skydata(3)
c        store input in skyspecs:
         skyspecs(1) = suntheta
         skyspecs(2) = sunphi
         skyspecs(10) = cloud
c
      elseif(iskyflag .eq. 3) then
         GMThr = skydata(1)
         cloud = skydata(2)
         sunphi = skydata(3)
c        store input in skyspecs:
         skyspecs(2) = sunphi
         skyspecs(9) = GMThr
         skyspecs(10) = cloud
c
      else
         print*,' invalid sky model option:  iskyflag =',iskyflag
         print*,' skydata =',(skydata(ii),ii=1,nsky)
         call HERR("initial","invalid sky model")  !stop run
      endif
c
c.....RECORD 9:  windspeed and index of refraction
c
      read(5,*) windspd, refr, T, S  !new in HE5:  index of refr, temp, and salinity
c     if refr<0, T and S will be used to calculate index of refraction
c
c     pass the wind speed to the Gregg&Carder sky irradiance model
      wsm = windspd             ! mean 24-hour wind
      ws = windspd              ! current wind speed
c
c.....RECORD 10:  Bottom boundary specifications
c
c     ibotm = 0    for an infinitely deep, homogeneous layer of water
c                  below depth m = zeta(nz)
c     ibotm .ge. 1 for an opaque reflecting bottom at depth zeta(nz).  The bottom
c                  has a BRDF as defined in routine BRDFbotm.f.  For the default 
c                  Lambertian BRDF, the irradiance reflectance R is determined by 
c                  the bottom type and wavelength, as follows: 
c
c         if ibotm = 1 use the value of R = rflbot read from input record 5 
c                      for all wavelengths
c         if ibotm = 2 read the bottom reflectance from a HYDROLIGHT
c                      standard-format file of bottom reflectance data
c
c     rflbot....the bottom reflectance, (always read, but used only if
c               ibotm = 1):  0.0 .le. rflbot .le. 1.0
c  
      read(5,*) ibotm,rflbot
c
c.....RECORD 11:  Output-depth information
c 
C     iop......0 if zetanom AS READ contains GEOMETRIC depths in meters 
c              1 if zetanom AS READ contains OPTICAL depths 
c     nznom....the number of depths (zeta levels) where output is desired
c     zetanom(1) = w = 0.0, ...,zetanom(nznom) = m ...the depths where output is desired 
c
c*****NOTE:  On the assumption that the user will want to compute
c     K-functions, each user-requested output depth is used to define
c     a nearby depth.  Each such pair of depths is then used to
c     estimate the depth derivatives needed for computing K functions.
c     (Accurate estimates cannot be obtained unless the depths
c     are closely spaced.).  Thus the actual HydroLight run solves
c     the RTE at closely spaced PAIRS of depths, e.g. 0.00, 0.01, 
c     1.00, 1.01,...,10.00, 10.01, if the user requested output at
c     0.00, 1.00, ..., 10.00.  The default printout is set below to 
c     give printout only at the user-requested depths (iz = 1, 3, 5,...;
c     e.g. 0.00, 1.00,..., 10.00 in this example).  The K-functions
c     estimated for a give pair of depths are associated with the
c     average depths of the pair, e.g. 0.005, 1.005, ... in this 
c     example.
c
c     zetanom is the user-requested output depths
c
c     min(nznom,mxz) added to read ommand to protect array dimens
      READ(5,*,err=67) iop,nz,(zeta(iz),iz=1,min(nz,mxz))
c
c     check to see if the depths are monotonically increasing
  67  ierr = 0
      do iz=2,nz
         if(zeta(iz) .le. zeta(iz-1)) ierr = 1
      end do
      if(ierr .ne. 0) then
        write(10,fmt="(//' ERROR:',/' depths are not monotonically ',
     1'increasing')")
        write(10,fmt="(/' depths = ',(10f10.4))") (zeta(iz),iz=1,nz)
        close(10)
        call HERR("initial","invalid depths")  !stop run
      endif
c
c     NOTE:  the run MUST use geometric depths if internal sources are
c     present
      if(isource.ne.0 .and. iop.eq.1) then
         write(10,1040)
         call HERR("initial",
     1        "cannot use optic depth when inelastic sources present")  !stop run
      endif
c
c
c.....Read the names of the "HydroLight standard format" user data files,
c      which may be needed by various routines (if a data file will 
c      not be needed in a particular run, the front end writes a
c      "dummy" name as a placeholder).  The files are
c    
c      ac9datafile = the file of ac-9 data, which is read by
c                    routine "abac9"
c      chlzdatafile = the file of Chl(z) data, which is read by
c                    routine "chlzdata"
c      Rbottomdatafile = the file of R(wavelength) bottom reflectance
c                        data, which is read by routine "rbottom"
c
C
C******NOTE:  About the array datafiles in common block Cfilenames
C               filenames are stored with the following indices:
C        (0): abPureWater file
C        (1): "normal"   ac9 data file
C        (2): "filtered" ac9 data file (optional; used to specify CDOM in abac9)
C        (3): Hydroscat data file
C        (4): chl(z) data file
C        (5): CDOM data file
C        (6): bottom reflectance data file
C        (7 to ncomp+6):  [optional] data files for each component
C        (ncomp+7):  User-supplied irradiance data file
C
C     The filenames stored are set in initial
c
c
c.....RECORD 12:  read the name of the pure water data file
c
c     Either case uses the Smith and Baker scattering coefficients
c     for "pure sea water."  See the header records in files
c     ..\data\SBH20ab.txt and ..\data\PFH20ab.txt for references.
c    
      read(5,fmt='(a)') H2Odatafile       !ac-9 datafile
c
c     check the string H2Odatafile for a '/' or a '\'.  If either is found,
c     the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the data directory
c
      lendata = lenstr(H2Odatafile)
      iabspath = 0
C
      do i=1,lendata
         if(H2Odatafile(i:i).eq.sl(1:1) .or. 
     1      H2Odatafile(i:i).eq.bsl(1:1)) iabspath = 1
      end do
c
      if(iabspath.eq.0) then
         H2Odatafile = datadir(1:lendatadir) // 
     1                 H2Odatafile(1:lendata)
      endif
      datafiles(0) = H2Odatafile
C
c

c.....RECORD 13:  read the name of the ac-9 data files
c
      read(5,*) nAC9files     !the number of ac-9 files used in this run (1 or 2)
c    
      read(5,fmt='(a)') ac9datafile       !ac-9 datafile
c     check the string ac9datafile for a '/' or a '\'.  If either is found,
c     the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the data directory
c
      lendata = lenstr(ac9datafile)
      iabspath = 0
C
      do i=1,lendata
         if(ac9datafile(i:i).eq.sl(1:1) .or. 
     1      ac9datafile(i:i).eq.bsl(1:1)) iabspath = 1
      end do
c
      if(iabspath.eq.0) then
         ac9datafile = datadir(1:lendatadir) // 
     1                 ac9datafile(1:lendata)
      endif
      datafiles(1) = ac9datafile
c
c    
      read(5,fmt='(a)') ac9Filterfile     !ac-9 filtered datafile
c     check the string ac9datafile for a '/' or a '\'.  If either is found,
c     the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the data directory
c
      lendata = lenstr(ac9Filterfile)
      iabspath = 0
C
      do i=1,lendata
         if(ac9Filterfile(i:i).eq.sl(1:1) .or. 
     1      ac9Filterfile(i:i).eq.bsl(1:1)) iabspath = 1
      end do
c
      if(iabspath.eq.0) then
         ac9Filterfile = datadir(1:lendatadir) // 
     1                 ac9Filterfile(1:lendata)
      endif
      datafiles(2) = ac9Filterfile
c
      read(5,fmt='(a)') HSdatafile        !hydroscat datafile
c     check the string HSdatafile for a '/' or a '\'.  If either is found,
c     the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the data directory
c
      lendata = lenstr(HSdatafile)
      iabspath = 0
C
      do i=1,lendata
         if(HSdatafile(i:i).eq.sl(1:1) .or. 
     1      HSdatafile(i:i).eq.bsl(1:1)) iabspath = 1
      end do
c
      if(iabspath.eq.0) then
         HSdatafile = datadir(1:lendatadir) // 
     1                 HSdatafile(1:lendata)
      endif
      datafiles(3) = HSdatafile
c

c.....RECORD 11:  read the name of the Chl(z) and CDOM data files
c    
      read(5,fmt='(a)') chlzdatafile
c     check the string chlzdatafile for a '/' or a '\'.  If either is found,
c     the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the data directory
c
      lendata = lenstr(chlzdatafile)
      iabspath = 0
      do i=1,lendata
         if(chlzdatafile(i:i).eq.sl(1:1) .or. 
     1      chlzdatafile(i:i).eq.bsl(1:1))  iabspath = 1
      end do
c
      if(iabspath.eq.0) then
         chlzdatafile = datadir(1:lendatadir) // 
     1                  chlzdatafile(1:lendata)
      endif
      datafiles(4) = chlzdatafile
c    
      read(5,fmt='(a)') CDOMdatafile
c     check the string chlzdatafile for a '/' or a '\'.  If either is found,
c     the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the data directory
c
      lendata = lenstr(CDOMdatafile)
      iabspath = 0
      do i=1,lendata
         if(CDOMdatafile(i:i).eq.sl(1:1) .or. 
     1      CDOMdatafile(i:i).eq.bsl(1:1))  iabspath = 1
      end do
c
      if(iabspath.eq.0) then
         CDOMdatafile = datadir(1:lendatadir) // 
     1                  CDOMdatafile(1:lendata)
      endif
      datafiles(5) = CDOMdatafile
c

c.....RECORD 12:  read the name of the Rbottom(wavelength) data file
      read(5,fmt='(a)') rbottomdatafile
c
c     check the string rbottomdatafile for a '/' or a '\'.  If either is 
c     found, the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the "..\data\reflbotm" directory
c
      lenrdata = lenstr(rbottomdatafile)
      iabspath = 0
      do i=1,lenrdata
         if(rbottomdatafile(i:i).eq.sl(1:1) .or. 
     1      rbottomdatafile(i:i).eq.bsl(1:1)) iabspath = 1
      end do
c
      if(iabspath.eq.0) then
         lenbdir=lenstr(bottdir)
         rbottomdatafile = bottdir(1:lenbdir) // 
     1                     rbottomdatafile(1:lenrdata)
      endif
      datafiles(6) = rbottomdatafile
c
c.....RECORD 13:  read the names of the component data files
c     (might be needed by the various models
      Do j = 7, 6+ncomp
        read(5,fmt='(a)',end=32) txtdatafile
c
c     check the string rbottomdatafile for a '/' or a '\'.  If either is 
c     found, the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the "..\data\" directory
c
        lenrdata = lenstr(txtdatafile)
        iabspath = 0
        do i=1,lenrdata
          if(txtdatafile(i:i).eq.sl(1:1) .or. 
     1       txtdatafile(i:i).eq.bsl(1:1)) iabspath = 1
        end do
c
        if(iabspath.eq.0) then
          lenbdir=lenstr(datadir)
          txtdatafile = datadir(1:lendatadir) // 
     1                      txtdatafile(1:lenrdata)
        endif
        datafiles(j) = txtdatafile
      Enddo
c
C     NEW in version 4.2.1:  read in optional irradiance data file
      read(5,fmt='(a)',end=32) txtdatafile
c
c     check the string for a '/' or a '\'.  If either is 
c     found, the name is assumed to be an absolute path name.  Otherwise,
c     the name is assumed to be a file in the "..\data\" directory
c
      lenrdata = lenstr(txtdatafile)
      iabspath = 0
      do i=1,lenrdata
        if(txtdatafile(i:i).eq.sl(1:1) .or.
     1     txtdatafile(i:i).eq.bsl(1:1)) iabspath = 1
      end do
c
      if(iabspath.eq.0) then
        lenbdir=lenstr(datadir)
        txtdatafile = datadir(1:lendatadir) //
     1                    txtdatafile(1:lenrdata)
      endif
      datafiles(7 + ncomp) = txtdatafile
c
c     read in the bioluminescent source file    !new in HE5
      read(5,fmt='(a)',end=32) S0datafile

      goto 999

 32   print *, 'WARNING:  only ',j,' component datafiles were read.'
      write(10,*)'WARNING:  only ',j,' component datafiles were read.'
c


cccccc  END READING RUN-TIME INPUT  cccccccccccccccccccccccccccccccc
c
  999 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Now store the run-time input data read above, and transfer the
c     needed values of imisc and fmisc from the surfwind file.
c
c     Arrays imisc and fmisc hold integer and floating point numbers
c     for used by various routines, as follows:
c
c     imisc(1) = nmu = the number of theta bands in theta = 0 to pi/2
c     imisc(2) = nphi = the number of phi bands in phi = 0 to 2*pi
c     imisc(3) = nL = the number of L terms in Fourier decompositions
c     imisc(4) = nz = the number of depths where output is computed,
c                including closely spaced pairs for K functions
c     imisc(5) = iop = flag for optical (1) or geometric (0) depth
c     imisc(6) = ncomp = number of components in IOP model in "ab" routine
c     imisc(7) = nwave = number of wavelength bands
c     imisc(8) = isource = flag for presence of internal sources or 
c                inelastic scattering (1 = present, 0 = not present)
c     imisc(9) = idbug = flag for additional printout for debugging
c                (most such code has been removed in version 4.0)
c     imisc(10) = nhat = dimension for spectral amplitude arrays
c     imisc(11) = jwave = index of current wavelength band
c     imisc(12) = ibotm = flag for bottom type
c     imisc(13) = isweep = flag for sweep type in Riccati integration
c     imisc(14) = L = current Fourier L value
c     imisc(15) = nuLfile = Fortran unit number for Lroot.txt printout
c     imisc(16) = nurad = Fortran unit number for file D_name_.txt
c                 containing the digital output for the run
c     imisc(17) = numrays = number of initial rays started in surfwind
c                 ray-tracing code 
c     imisc(18) = nuscr1 = Fortran unit number for scratch file 1
c     imisc(19) = nuscr2 = Fortran unit number for scratch file 2
c     imisc(20) = iwrtss1 = flag for writing single-wavelength
c                 spreadsheet data (files written iff = 1)
c     imisc(21) = iwrtss2 = flag for writing multi-wavelength
c                 spreadsheet data (files written iff = 1)
c     imisc(22) = Digital output file option (file written iff = 1)
c     imisc(23) = nconc = the number of component concentrations needed
c			(this is only relavent if CDOM or fluorescence is present)
c     imisc(24) = iOptRad = flag to write radiancer Lroot file
c                 (files written iff = 1)
c     imisc(25) = iZmax = flag for specifying how zmax is determined at each wavelength:
c                 = 0 to use the max printout depth at all wavelengths
c                 = 1 to dynamically determine zmax from a min PAR value (as in E2.1)
c                 = 2 to dynamically determine zmax from 1/Kd, for rem. sens. studies only
c     imisc(26) = NwSkip = number of wavebands skipped; used to sel alternating
c                 wavebands IFF no inelastic scattering is included
c     imisc(27) = iDynZ : flags if Dynamic depth allowed when sources 
c                 and inf bottom opt are selected (improves calc accuracy)
c     imisc(28) to imisc(30) = unused
c
c     fmisc(1) = pi = 3.141952654...
c     fmisc(2) = degrad = 0.0174533... converts degrees to radians
c     fmisc(3) = radeg = 57.2958... converts radians to degrees
c     fmisc(4) = refr = 1.34 = index of refraction of water
c     fmisc(5) = rad48 = 0.842439... = asin(1/refr) = critical angle for
c                internal reflection
c     fmisc(6) = relerrs = relative error for ODE integration
c     fmisc(7) = abserrs = absolute error for ODE integration
c     fmisc(8) = areset = absorption coef or flag for routine abconst
c     fmisc(9) = breset = scattering coef or flag for routine abconst
c     fmisc(10) = deltazK = the depth increment for use in generating
c                 closely spaced pairs of depths for K function
c                 computations
c     fmisc(11) = acoef(nz) = the absorption coef at the deepest depth
c     fmisc(12) = bcoef(nz) = the scattering coef at the deepest depth
c     fmisc(13) = wavel = current nominal wavelength in nm
c     fmisc(14) = rflbot = the bottom reflectance
c     fmisc(15) = wndspeed = the wind speed in m/s
c     fmisc(16) = delta = parameter in L&W Eq. 4.39
c     fmisc(17) = epsilon = parameter in Eq. 4.39
c     fmisc(18) = znow = current depth in meters
c     fmisc(19) = atotal = total absorption coef the the current depth
c                 wavelength
c     fmisc(20) = btotal = total scattering coef at the current depth
c                 and wavelength
c     fmisc(21) = a0R = Raman scattering coefficient
c     fmisc(22) = wave0 = Raman scattering reference wavelength
c     fmisc(23) = PhiChl = chlorophyll fluorescence quantum efficiency
c     fmisc(24) = delbbTol = (no longer used)
c     fmisc(25) = F0 (used only by Ecolight code)
c     fmisc(26) = PARmin = min wavelength used for PAR calculations
c     fmisc(27) = PARmax = max wavelength used for PAR calculations
c     fmisc(28) = Salinity (ppt or PSU) 
c     fmisc(29) = Temperature (degC)
c     fmisc(30) = index of refractrion. If > 0, use same refr at all wavel,
c                 if < 0, calc as a func of wavel, S, T
c
c     Misc file names are stored as follows:
c
c     datafiles(0) = pureH20abc = the file containing the pure-water
c                    a, b, and c values to be used
c     datafiles(1) = abac9data = the file containing the user's ac-9 data
c     datafiles(2) = abac9Filtered = the file containing the user's FILTERED ac-9 data
c     datafiles(3) = user-supplied Hydroscat data file (for abac9 model)
c     datafiles(4) = chlzdatafile = the file containing the user's 
c                    chlorophyll vs depth profile
c     datafiles(5) = CDOMzdatafile = the file containing the user's aCDOM profile
c     datafiles(6) = Rbottomdatafile = the file containing the user's
c                    bottom reflectance vs wavelength data
c     datafiles(7) to datafiles(ncomp+6) = data files corresponding to each component
c
c     Set EL parameters
      iZmax = 0         !solve to max output depth at each wavelength
      F0 = 1.0          !solve all the way down

c     store values from runtime input
      imisc(4) = nz
      imisc(5) = iop
      imisc(6) = ncomp
      imisc(7) = nwave
      imisc(8) = isource
      imisc(11) = 1       !init to first wavelength
      imisc(12) = ibotm
      imisc(15) = 17
      imisc(23) = nconc
      imisc(25) = iZmax
      imisc(26) = nwskip
c
      fmisc(8) = areset
      fmisc(9) = breset
      fmisc(13) = wave(1)
      fmisc(14) = rflbot
      fmisc(15) = windspd
      fmisc(25) = F0
      fmisc(26) = PARmin
      fmisc(27) = PARmax
c
      fmisc(28) = S
      fmisc(29) = T
c
      fmisc(30) = refr
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Open the file of surface information for the chosen windspeed
c     and read in grid information from the header records of that
c     file. 
      wavelen = wave(1)
      call loadsurface(windspd, wavelen)
c
cccccc  PRINTOUT INITIALIZATIONS cccccccccccccccccccccccccccccccccccc

c-----initialize printout depths (EL/HL dependent)
      CALL initprintout(idbug)
      nz = imisc(4)  !value of nz updated in HL's initpntout
      do iz=1,nz     !used to be done in first call to inishamp
        zgeo(iz) = zeta(iz)
      end do
c
c-----Printout the (theta,phi,depth,wavelength) grid
      CALL PNTGRID
c
c-----Trigger IOP model first-call printout
      write(10,200) "IOP MODEL"
      depth = 0.0 
      wavelen =  wave(1)
      ncomp = imisc(6) 
      CALL abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
      write(10,201) "IOP MODEL"
c
c
c-----Printout the Ineleastic scattering and fluorescence
      write(10,200) "INELASTIC SCATTER AND INTERNAL SOURCE"
c      Print info to Prot.txt about whether inelastic scatter is included
       if(ibiolum.gt.0) then
         write(10,100)
         S0 = s0bioSub(0.0,0.0)
       endif 
c     Double check that indexchl is valid of chl-fluorescence is included
      if(ichlfl.ne.0 .and. (indexchl.eq.0 .or. indexchl.gt.nconc)) then
        ichlfl = 0
        write(10,101)
      elseif(ichlfl.eq.1) then
        write(10,102) PhiChl
c       make sure data file messages have been triggered
        call achlz(depth, wavelen, achl)        
      else
         write(10,103) 
      endif
      if(icdomfl.eq.1) then
         write(10,110)
c       make sure data file messages have been triggered
        call acdomsub(depth, wavelen, achl)
      else
         write(10,111)
      endif
       if(iraman.eq.1) then
         write(10,120) RamanXS, Raman0, RamanEXP
       else
         write(10,121) 
       endif
      write(10,201) "INELASTIC SCATTER AND INTERNAL SOURCE"


c-----Printout the ATMOSPHERIC MODEL VALUES
      write(10,200) "ATMOSPHERIC MODEL"
      CALL qasky
      write(10,201) "ATMOSPHERIC MODEL"

c
c-----Printout the SURFACE BC
      write(10,200) "AIR-WATER SURFACE"
      write(10,300) windspd
      If(fmisc(30).gt.0) then      !constant refr
          write(10,301) fmisc(30)
      Else                        !variable refr
          write(10,302) T, S
      Endif
      write(10,201) "AIR-WATER SURFACE"

c-----Printout the Seafloor BC
      write(10,200) "BOTTOM BOUNDARY"
      CALL PNTBC
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c.....add dynamic depth calc depths IFF iDynZ is selected
c     (this improves calc accuracy when inelastic sources and inf bottom 
c      are selected)  i(27)=iDynZ; i(12)=ibotm; i(8) = isource
      If(imisc(27).eq.1 .and. imisc(12).eq.0 .and. imisc(8).gt.0) then
        call getDynZ(0)
        if(imisc(27).eq.1) write(10,400) !DynZ may turn off DynZ if no avail nz
      Endif
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write(10,201) "BOTTOM BOUNDARY"
c

c-----printout the data file range warnings (gathered during above init)
      nuout = 10
      CALL RangeWarn(nuout)  

      return
c
  100 format(/5x,'A bioluminescent source is included in this run')
  101 format(/5x, 'WARNING: Chlorophyll fluorescence could not',
     1       ' be included because CHL is not specified!',
     2       /8x,'To include chl-fluorescence, rerun frontend.')
  102 format(/5x,'Chlorophyll fluorescence is included in this run',
     2       /8x,'Chl fluorescent quantum efficiency set to ',f6.3)
  103 format(/5x,'Chlorophyll fluorescence is not included in this run')
  110 format(/5x,'CDOM fluorescence is included in this run')
  111 format(/5x,'CDOM fluorescence is not included in this run')
  120 format(/5x,'Raman scattering is included in this run',
     1       /8x,'The Raman scattering cross section is ',1pe8.2,
     2' 1/m at reference wavelength ',0pf6.1, ' nm',
     3       /8x,'The Raman scattering cross section depends on the EXCI
     4TATION wavelength as (reference wavelength/excitation wavelength)^
     5',f4.2)
  121 format(/5x,'Raman scattering is not included in this run')
c
  200 format(//2x,"***** BEGIN ",a," SPECIFICATIONS *****")
  201 format(/2x,"***** END ",a," SPECIFICATIONS *****")
c
  300 format(/5x,'The air-water surface is for a wind speed of',
     1  f5.1,' m/s')
  301 format (/5x,'Index of refraction is constant for all wavelengths'
     1 ' with the value: ',f5.3)
  302 format (/5x,'Index of refraction will be calculated as a ',
     1 'function of wavelength, ',/8x,' for the given seawater temp ',
     2 f4.1, ' degC and salinity ',f4.1,' PSU')
  400 format(/5x,'The actual depth at which the infinitely deep bottom ',
     1 'boundary condition will be applied '/5x,'will be determined ',
     2 'dynamically from the IOPs to insure a good solution at the ',
     3 /5x,'requested output depths')
c
 1004 format(a)
 1005 format(//2x,'RUN TITLE:  ',a)
 1040 format(//'ERROR:'/' OPTICAL depth being used with internal',
     1' sources:  MUST use GEOMETRIC depth')
c
      end
c

