C     Last change:  LKS   6 Sep 2013   10:46 pm
      subroutine abacbb(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
c
c     Example user "ab" routine on file abacbb.f
c
c     called by INISHAMP, BEAMc1, BEAMcz, and RHOTAU
c
c     calls yinterp, pureh2o
c
c     New in v5:  calls function CheckRange to check that interpolation
c                 values are in bounds of the data.  If not, it returns
c                 the min/max in the data set (function in this file)
c
c     New in v5.2: bb data files can include or not include pure water bb values
c
c     INPUT:
c        depth: the depth in meters where the absorption and scattering
c               coefficients are requested
c        wavelen: the wavelength in nanometers where the absorption and
c                 scattering coefficients are requested
c        ncomp: the number of components expected
c
c     OUTPUT:
c        acomp:  acomp(i),i=1 to ncomp, is the absorption coefficient
c                of the ith component
c        bcomp:  bcomp(i),i=1 to ncomp, is the scattering coefficient
c                of the ith component
c        atotal: the total absorption coefficient (sum of the component
c                values)
c        btotal: the total scattering coefficient (sum of the component
c                values)
c
c
c     This routine will read and process up to TWO ac-9 data files (one
c     filtered, one not), and an optional HydroScat data file containing
c     bb data. If the filtered ac-9 data is included, ncomp=3 (CDOM added)
c
c     This routine models a and b with 2 independent components:
c          component 1 is pure water
c          component 2 is "everything else" (particles, CDOM, etc)
c
c     For component 1, the pure-water a and b values are retrieved
c           from calls to the routine pureH2O() 
c     
c     For component 2, the a and c values are read from a 
c           "ac9datafile" input data file 
c
c     This routine opens and reads the file "ac9datafile", where 
c     "ac9datafile" is the actual path and filename as specified when 
c     running the front end program.  If HydroScat data is included
c     the file "HSdatafile" will also be read.
c
c     This routine is designed for reading in WETLabs ac-9 data (or other
c     data on the same format, such as from an ac-s or from separate a
c     and c instruments after merging their data).
c
c     This routine reads a file of absorption (a) and beam
c     attenuation (c) data on the HydroLight "ac-9 standard format."
c     These values are assumed to have had the "pure-water" values 
c     subtracted out in the processing of the ac-9 data.
c
c     If HydroScat-format backscatter data are included for specifying bb,
c     the bb data file can include pure water bb values (ibbopt = +2)
c     or omit the pure water bb values (ibbopt = -2).  Pure water
c     values will either be subtracted or not from the data as
c     read to get the particle bb data for use in determining the
c     particle phase function.
c
c     When this routine is called for the first time, the data files are
c     read and processed as follows:
c     1) a and c values are read for the discrete wavelengths, and at
c        each discrete depth for component 2; bb may also be read
c     2) b = c - a is computed for each component for each depth and
c        wavelength
c
c     On all subsequent calls, linear interp is used to define a and b
c     (and bb, if needed) at any depth and wavelength.  
c
c**** NOTE: If the calling depth and wavelength are outside the domain of
c     the original data, the nearest measured values will be used (rather
c     than exiting with an error message), but this is clearly not going
c     to give a proper description of the water body in most cases.
c
      INCLUDE "dimens_XL.inc"
c
      Common /Csource0/ ibiolum,ichlfl,icdomfl,iraman, ramanEXP
      common /Cbbopt/ ibbopt(mxcomp), bbfrac(mxcomp), 
     1                BfrefPL(mxcomp), Bf0PL(mxcomp), BfmPL(mxcomp)
      integer ibbopt
      real bbfrac, BfrefPL, Bf0PL, BfmPL
      Common /Cabac9n/ nac9files
c       flag for calculating bb shared with inishamp routine (calc bb iff flag=1)
      Common /CbbCalc/ ibbCalc
c     ibbopt= 0:  Phase function file explicitly specified (use pfname)
c             1:  Select phase function to match specified bb/b ratio
c             2:  bb/b specified by bb data file (bb data includes pure water bb)
c            -2:  bb/b specified by bb data file (bb data does NOT include pure water bb)
c    ibbCalc= 0:  do not calc bb/b (z to zeta table is being constructed)
c             1:  calculate bb/b if needed (RTE is being solved)
c
      dimension acomp(*),bcomp(*)
c
      common/cmisc/ imisc(30),fmisc(30)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
c
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      Character*120 ac9datafile, ac9Filtered, HSdatafile
c
c     arrays for reading ac-9 type data
      real, allocatable ::  zdata(:),wavdat(:),ap(:,:),bp(:,:)
      real, allocatable ::  zdataF(:),wavdatF(:),apF(:,:)
      real, allocatable ::  zdataHS(:),wavdatHS(:),bbp(:,:)
      Integer nzdata(3), nwavdat(3)
c     arrays for storing data put onto run wavelength grid
      real, allocatable ::  apg(:,:), apFg(:,:), bpg(:,:), bbpg(:,:)
c     arrays for temporary processing of data 
      real, allocatable :: yDat(:), zDat(:)
      real, allocatable :: xgrid(:), ygrid(:), zgrid(:) !xgrid persists as wave()
      real, allocatable :: bbpure(:)  !pure water bb to be removed from HS data
c
      integer il(3)   !last index for starting point of linear interp
c
      data kall/0/, numcomp/2/, waveold/0.0/
      data iclose/0/, nheadr/10/, nh2p/2/
      data il/1,1,1/
c
      integer getndat, getnwavac
      external getndat, getnwavac

c     Note on dynamic storage:  MUST save values between subroutine calls:
      save
c
c-----  Begin Initialization on first call  ---------------------------
c
c     On the first call to abac9:
c       1) read the file of "HydroLight standard format" ac9 data 
c       2) process it to get a and b at the discrete depths and wavelengths 
c          of the original data
c       3) fit a 2-D polynomial to the discrete data.  This polynomial
c          will be used in the subsequent calls to get a and b at any
c          depth and wavelength by interpolation.
c
c     Note:  pure water a and b values will be obtained from routine
c            pureH2O
c
      nconc = imisc(23)
      nwave = imisc(7)
      jwave = imisc(11)
c
      if(kall.eq.0) then
c
c     initialization on first call
        ac9datafile = datafiles(1)
        ac9Filtered = datafiles(2)
        HSdatafile=datafiles(3)
        iop = imisc(5)

c         Test whether run is compatible with AC9 model
        if(numcomp.ne.ncomp) then
c           stop the run if the requested number of components is
c           not compatible with abac9, which is a two-component IOP
c           model
           write(10,102) ncomp, numcomp
           call HERR("abacbb","# components is not compatible")  !stop run
        endif
          if(iop.ne.0) then
c           stop the run if it is using optical depth; abac9 expects
c           geometric depth as input
            write(10,104) iop
           call HERR("abacbb","depth type is not compatible")  !stop run
          endif
c
c         Write description for this run
          write(10,100) ac9datafile(1:lenstr(ac9datafile))
          IF(nac9files.eq.2) write(10,201)
          IF(nconc.eq.4) write(10,202)
          IF(nconc.gt.ncomp) write(10,203)
c         trigger pure water messages
          call pureH2O(wavelen, awater,bwater,cwater)
c
c-------- ac9 data -------
c     open and read the ac-9 data file
      nudata = 42
      nzdata(1) = getndat(iclose,nudata,nheadr,nh2p,ac9datafile) - 1
      nwavdat(1) = getnwavac(nudata,nheadr) 
      allocate(zdata(nzdata(1)), wavdat(nwavdat(1)))
      allocate(ap(nzdata(1), nwavdat(1)), bp(nzdata(1), nwavdat(1)))
      call readac9(nudata,nzdata(1),nwavdat(1),zdata,wavdat,ap, bp)
c     generate warnings if data range doesn't cover run range
      CALL RangeCheck(1, zdata(1),zdata(nzdata(1)), 
     1                wavdat(1), wavdat(nwavdat(1)), ac9datafile )
!     put data on wavelength grid
      allocate(apg(nzdata(1), nwave), bpg(nzdata(1), nwave))
      allocate(ydat(nwavdat(1)), zdat(nwavdat(1)) )
      allocate(xgrid(nwave), ygrid(nwave), zgrid(nwave)) 
      xgrid = wave(1:nwave)
      itype = 2  !(fit 2 columns of data)
      xRef = 0.0
      Do i=1,nzdata(1)
        yDat = ap(i,:)
        zDat = bp(i,:)
        Call GridDat(itype,nwave,nwavdat(1),xref,yref,
     1               xgrid, ygrid,zgrid,
     1               wavdat,yDat,zDat)
        apg(i,:) = ygrid   !gridded ap values
        bpg(i,:) = zgrid   !gridded bp values
      Enddo
      deallocate(ap,bp, wavdat, ydat, zdat)
C
c-------- Filtered ac9 data -------
c     If the filtered data file is present, fit its data the same way
      IF(nac9files.eq.2) then
        if(nconc.lt.3) nconc=3            !double check that CDOM will be included in printout
        write(10,101) ac9Filtered(1:lenstr(ac9Filtered))
c       open and read the ac-9 filtered data file
        nudata = 42
        nzdata(2) = getndat(iclose,nudata,nheadr,nh2p,ac9Filtered) - 1
        nwavdat(2) = getnwavac(nudata,nheadr) 
        allocate(zdataF(nzdata(2)), wavdatF(nwavdat(2)))
        allocate(apF(nzdata(2), nwavdat(2)))
        call readac9F(nudata,nzdata(2),nwavdat(2),zdataF,wavdatF,apF)
c       generate warnings if data range doesn't cover run range
        CALL RangeCheck(1, zdataF(1),zdataF(nzdata(2)), 
     1                  wavdatF(1), wavdatF(nwavdat(2)),ac9Filtered )
!       put data on wavelength grid
        allocate(apFg(nzdata(2), nwave))
        allocate(ydat(nwavdat(2)), zdat(nwavdat(2)) )
        itype = 1  !(fit 1 column of data)
        Do i=1,nzdata(2)
          yDat = apF(i,:)
          Call GridDat(itype,nwave,nwavdat(2),xref,yref,
     1                 xgrid, ygrid,zgrid,
     1                 wavdatF,yDat,zDat)
          apFg(i,:) = ygrid   !gridded ap values
        Enddo
        deallocate(apF, wavdatF, ydat, zdat)
      ENDIF
c
c--------- DPF printout 
      call PNTDPF(1)
      call PNTDPF(2)
c
c-------- hydroscat data -------
c     initialize the hydroscatt data if and only if requested
      if(abs(ibbopt(2)).eq.2) then
c       open and read the Hydroscat data file
        nudata = 42
        nzdata(3) = getndat(iclose,nudata,nheadr,nh2p,HSdatafile) - 1
        nwavdat(3) = getnwavac(nudata,nheadr) 
        allocate(zdataHS(nzdata(3)), wavdatHS(nwavdat(3)))
        allocate(bbp(nzdata(3), nwavdat(3)))
        call readac9F(nudata,nzdata(3),nwavdat(3),zdataHS,wavdatHS,bbp)
c       generate warnings if data range doesn't cover run range
        CALL RangeCheck(1, zdataHS(1),zdataHS(nzdata(3)), 
     1                  wavdatHS(1), wavdatHS(nwavdat(3)),HSdatafile )
!       put data on wavelength grid
        allocate(bbpg(nzdata(3), nwave))
        allocate(ydat(nwavdat(3)), zdat(nwavdat(3)) )
        allocate(bbpure(nwave))

        if(ibbopt(2).gt.0) then         !subtract pure water
            Do ij=1,nwave
              call pureH2O(wave(ij),awater,bwater,cwater)
              bbpure(ij) = 0.5*bwater
            Enddo
        endif

        itype = 1  !(fit 1 column of data)
        Do i=1,nzdata(3)
          yDat = bbp(i,:)
          Call GridDat(itype,nwave,nwavdat(3),xref,yref,
     1                 xgrid, ygrid,zgrid,
     1                 wavdatHS,yDat,zDat)
         if(ibbopt(2).gt.0) then
           bbpg(i,:) = ygrid - bbpure  !subtract pure water from total bb data
         else
           bbpg(i,:) = ygrid ! data as read is particle bb
          endif
        Enddo
c
        deallocate(bbp, ydat, zdat, ygrid, zgrid, bbpure)
c
c       Call fit routines to get bbfrac for initial depth and wavelength
c       and load appropriate DPF file
        bpart =  yinterp(il(1), nzdata(1), depth,zdata,bpg(1,jwave) )
        bbhscat = yinterp(il(3),nzdata(3),depth,zdataHS,bbpg(1,jwave))
      Endif                  !end of Hydroscatt specs
c--------------------------------
      kall = 1
      endif                  !end of initialization
c-----  End of initialization on first call----------------------------
c
c-----  Subsequent calls start here:
c
c-------- pure water data -------
c     Pure water values.  If this call is at a new wavelength, call
c     "pureH2O" to get the pure water a and b values
      if(wavelen.ne.waveold) then				
         call pureH2O(wavelen, awater,bwater,cwater)
         waveold = wavelen
      endif    !end block for new wavelength
c
c-------- ac9 data -------
c     Get particle values from linear interp of gridded data
      apart =  yinterp(il(1), nzdata(1), depth, zdata, apg(1,jwave) )
      bpart =  yinterp(il(1), nzdata(1), depth, zdata, bpg(1,jwave) )

c-------- Filtered ac9 data -------
c     If the filtered data file is present, use it to get component 3
c     and to correct the absorption of component 2
      IF(nac9files.eq.2) then
        aCDOM= yinterp(il(2),nzdata(2), depth, zdataF, apFg(1,jwave)) !rev2
      ENDIF
c--------------------------------
c
      acomp(1) = awater
      bcomp(1) = bwater
      acomp(2) = apart
      bcomp(2) = bpart
      atotal = acomp(1) + acomp(2)
      btotal = bcomp(1) + bcomp(2)
c
c     Add aCDOM to printout iff available  
c        (for printout only)
      IF(nac9files.eq.2) then
c     from filtered data contribution iff it is included 
        acomp(3) = aCDOM
        bcomp(3) = 0.0
      ElseIf(icdomfl.gt.0) then
c     from fluorescence calculations (iff cdom fluorescence turned on)
        call acdomsub(depth,wavelen, aCDOM)
        acomp(3) = aCDOM
        bcomp(3) = 0.0
      
      ENDIF

c     add Chlorophyll data contribution iff it is included 
c        (for printout only)
      IF(nconc.eq.4) then
         call achlz(depth, wavelen, achl)
         acomp(4) = achl
         bcomp(4) = 0.0
      ENDIF
c
c     If Hydroscat option is selected...
C     AND not if we are just building up the z to zeta table
      if(abs(ibbopt(2)).eq.2 .and. ibbCalc.ge.1) then
        bbhscat = yinterp(il(3),nzdata(3),depth,zdataHS,bbpg(1,jwave))
        icomp = 2
c       calculate the local bb/b value for the particulates
c       **Check for b = 0
        If(bcomp(icomp).gt.0) then
c         ** calc new bb fraction
          temp = (bbhscat) / bcomp(icomp)   !water was removed from data if needed
c         **Make sure the value is positive
          IF( temp .gt. 0 ) then
            IF(ibbCalc.eq.2) then
              call selpfbb(temp,icomp)
            ENDIF
            bbfrac(icomp) = temp
          ENDIF
        Else
          bcomp(2) = 0
        Endif
      Endif  !end HydroScat block
c
 20   return
 
c     Formats:
  100 format(/5x,'The IOP routine "abacbb" is being used:'//5x,
     1'Absorption and beam attenuation coefficients will be read from ',
     2'a standard-format (unfiltered) ac data file named:'/
     38x,a,/5x,
     7'(2D interpolation in depth and wavelength is used to define',
     8' a and b at all depths and wavelengths)'//,
!     45x,'and will be used in a 2-component IOP model where'//
     58x,'component 1 is pure water'/
     68x,'component 2 is "everything else"' )
  101 format(//5x,'CDOM absorption will be obtained from a FILTERED ',
     2'ac data file named:'/ 8x,a,/)
  102 format(//'Error in sub abacbb:  ncomp =',i2,
     1' must equal',i2,' for this routine')
  104 format(//'Error in sub abacbb:  called with optical depth: ',
     1       'iop = ',i2) 
  110 format(//'Error in sub abacbb:  input data file contains',i5,
     1' depths; spline fitting requires at least 3 depths')
  112 format(//'Error in sub abacbb:  input data file contains',i5,
     1' wavelengths; spline fitting requires at least 3 wavelengths')
  201 format(/5x,'CDOM absorption values will be reported',
     1' as Component 3',/8x,'(only used for fluorescence calculations)')
  202 format(/5x,'Chlorophyll absorption values will be reported ',
     1'as Component 4'/,8x,'(only used for fluorescence calculations)')
  203 format(/5x,'Note:  Total absorption will be given by the sum of',
     1' components 1 and 2.',/5x,'Output for additional components are',
     2' provided solely for informational purposes',/5x,'to indicate',
     3' what values are used in the fluorescence calculations.'/)
      end
c
c     **************************************************
      subroutine readac9(nudata,nzdata,mxwavdat,
     1                   zdata,wavdat,ap, bp)
c     Reads TWO 2D data sets
      integer nudata, nzdata, mxwavdat
      real zdata(nzdata),wavdat(mxwavdat)
      real ap(nzdata,mxwavdat), bp(nzdata,mxwavdat)
c
      integer irec, nwavdat
      character*80 header, msg
      data numhead/10/
c
      rewind(nudata)
      do irec=1,numhead
         read(nudata,800) header
      end do
      read(nudata,*) nwavdat, (wavdat(i),i=1,nwavdat)
c     read the FILTERED absorption data records ONLY.  
      Do irec = 1, nzdata
        read(nudata,*,end=820) zdata(irec),(ap(irec,i),i=1,nwavdat),
     1                         (bp(irec,i),i=1,nwavdat)
c       convert ap and cp values to ap and bp = cp - ap
        bp(irec,:) = bp(irec,:) - ap(irec,:)
c       Data must be monotonically increasing
        If(irec.gt.1) then
          If(zdata(irec).le.zdata(irec-1)) then
            msg="Data abscissas not monotonically increasing." //
     1          "  Clean up your data and rerun."
            call HERR('readac9', msg)
          Endif
        Endif
c       *Check goodness of data (must be =>0)
        Do i=1,nwavdat
!          write(10,801) 'data: ',irec,i,zdata(irec),wavdat(i),
!     1                         ap(irec,i),bp(irec,i)
  801 format(2x,'***',a,2i6,2f8.2,1p,2E12.3)
          call checkDat2x(zdata(irec), wavdat(i), ap(irec,i))
          call checkDat2x(zdata(irec), wavdat(i), bp(irec,i))
        Enddo
      Enddo
c     entire data file has been read
  820 close(nudata)
      return 
  800 format(a)
  803 format(8x,a)
      end
c     **************************************************
      subroutine readac9F(nudata,nzdata,mxwavdat,
     1                   zdata,wavdat,ap)
c     Reads ONE 2D data set
      integer nudata, nzdata, mxwavdat 
      real zdata(nzdata),wavdat(mxwavdat)
      real ap(nzdata,mxwavdat)
c
      integer irec, nwavdat, numhead
      character*80 header, msg
      data numhead/10/
c
      rewind(nudata)
      do irec=1,numhead
         read(nudata,800) header
      end do
      read(nudata,*) nwavdat, (wavdat(i),i=1,nwavdat)
c     read the FILTERED absorption data records ONLY.  
      Do irec = 1, nzdata
        read(nudata,*,end=820) zdata(irec),(ap(irec,i),i=1,nwavdat)
c       Data must be monotonically increasing
        If(irec.gt.1) then
          If(zdata(irec).le.zdata(irec-1)) then
            msg="Data abscissas not monotonically increasing." //
     1          "  Clean up your data and rerun."
            call HERR('readac9F', msg)
          Endif
        Endif
c       *Check goodness of data (must be =>0)
        Do i=1,nwavdat
          call checkDat2x(zdata(irec), wavdat(i), ap(irec,i))
        Enddo
      Enddo
c     entire data file has been read
  820 close(nudata)
      return 
  800 format(a)
  803 format(8x,a)
      end

