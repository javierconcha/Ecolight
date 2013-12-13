C     Last change:  LKS  14 Nov 2009    7:32 pm
      subroutine abnewcase1(z,wavenm,ncomp, acomp,bcomp,atotal,btotal)
c
c     user routine on file abnewcase1.f   
c     Please see HE5 Tech Documentation for details of this IOP model
c
c     called by INISHAMP, BEAMc1, BEAMcz, and RHOTAU
c
c     calls pureh2o, astar, bstar, and abmodel
c
c     This routine returns the absorption and scattering coefficients
c     a and b (in 1/meter), given the wavelength (350 to 700 nm)
C     and GEOMETRIC depth z (in m), for use in simulation of case 1
c     waters using the new Morel model.
c
c     In general, the concentration of each component of a general N
c     component model will be determined via a single inclusion of the
c     an abmodel routine call.  The input to this module are the two
c     input variables: depth, and wavelen.  The output will be of the
c     the array compconc(icomp).  The names of these three variables
c     (depth, wavelen, and compconc(mxcomp)) and the output compconc 
c     array must be dimensioned.
c
c     Calls to the chlorophyll routine can still be made via a call to
c     the chlz routine.  Similarly, acdom can be solved by calling
c     acdom.  This is how the inelastic scattering routiens in shat.f 
c     will get the concentrations if needed.  However, within the "ab"
c     routine, users are strongly encouraged to use a single reference to  
c     abmodel at each depth/wavelenght pair -- it is safer (less prone
c     to programming errors by the user) and more consistent.
c
c     A component's specific absorption and scattering is calculated via
c     models provided below.
c
c     Descriptions of the components of this model:
c
c     For component 1, routine pureH2O is called to get the pure-water
c     absorption and scattering values.
c  
c     For components 2 and 3, a user-supplied subroutine (inserted via the 
c     call to abmodel) is called to obtain the chlorophyll 
c     concentration chlconc (in mg/m^3) at geometric depth z. 
c     Component 2 contains the particulate absorption and small particle scattering;
c     Component 3 contains the large particle scattering
c
c     For Component 4, co-varying absorption due to CDOM is calculated (non-scattering)
c
      INCLUDE "DIMENS_XL.INC"

!      COMMON /CConstConc/ compconc(mxcomp)
      common /Cconc/  itype(mxcomp)  
      integer itype
      real compconc(mxcomp), acd, temp

      dimension acomp(mxcomp),bcomp(mxcomp)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)             !EL version of common block
      common/cmisc/ imisc(30),fmisc(30)
      COMMON /Cpirrad/ npirad,izirad(mxz)
c
      real nu, Acapp, Ecapp, Acapp440, Ecapp440 
      real alpha, bpSmall, bpLarge
      real fmu, bndmu, omega, deltmu, zgeo, zeta   !type CGRID elements
c
      data kall/0/, numcomp/3/, waveold/0.0/
c
c     Note on dynamic storage:  MUST save data values between subroutine 
c     calls:
      save
!*********************************************************************
c
      if(kall.eq.0) then
c*****Initialization on the first call:
         write(10,100)
c        reset number of concentrations for printout so that CDOM can be
c        printed separately
         nconc = 4
         imisc(23) = nconc
c        check for a valid call to this routine
         iop = imisc(5)
         if(numcomp.ne.ncomp) then
            write(10,102) ncomp, numcomp
            call HERR("abNewCase1","# components is incompatable")  !stop run
         endif
         if(iop.ne.0) then
            write(10,104) iop
            call HERR("abNewCase1","depth type is incompatable")  !stop run
         endif
c
c        call the desired Chl(z) routine first time to trigger messages
         depth = 0.
         wavelen = wavenm
c       
c#### PURE WATER
         call pureH2O(440.0, aw440,bwater,cwater)
         aw440 = awater    !lks 070122
c--------- DPF printout
         CALL PNTDPF(1)
c
c#### CHLOROPHYLL
c     print out the chl-spec header info iff constant option was selected
c     (if function or datafile was specified the appropriate routine will
c     print out the header info)
      icomp = 2
!      IF(itype(icomp).eq.0) then		!constant with depth
!        write(10,201)  compconc(icomp)
!      ELSEIF (itype(icomp) .EQ. 2) then  !read datafile
!        call chlzdata(icomp, depth, compconc(icomp))
!      ELSE
        call achlz(depth, wavelen, temp)       !calling to spawn msgs
!      EndIf
c--------- DPF printout
      CALL PNTDPF(2)
c
      call abmodel(depth, compconc)    !print any stray msgs
c
c#### Printout CHL profile
         if(imisc(9).ge.0) then
           write(10,106)
           do iz=1,npirad
             depth = zgeo(izirad(iz))
             call abmodel(depth, compconc)
             write(10,107) depth,compconc(2)
           end do
         endif
         kall = 1
      endif
c
c-----End of initialization on first call.
c
c     Subsequent calls start here:
c
c     If this is the first call at this wavelength, call routine
c     pureH2O to get the pure water values and call getAEnewcase1 to get
c     ap coeff
c
      wavelen = wavenm
      if(wavelen.ne.waveold) then
c        get pure water IOPs
         call pureH2O(wavelen, awater,bwater,cwater)
c        get chlorophyll coeff
         call getAEnewcase1(Acapp, Ecapp, Acapp440, Ecapp440)
c        get CDOM coeff
         Xcdom = exp(-0.014*(wavelen-440.))
c        update waveold 
         waveold = wavenm
      end if
c     
c     Get the chlorophyll concentration at geometric depth z
c
c     insert a call to the desired Chl(z) routine:
      depth = z
c
c     The call to abmodel will return the array compconc(ncomp), the
c     concentration of all components based on selections made in the 
c     Hydrolight Frontend
      call abmodel(depth, compconc)
      chlconc = compconc(2)
c
c     get scattering chl-dep coefficients
!     calculate Junge slope
      If(chlconc .le. 0 .or. chlconc.gt.2.0) then  
        nu = 0.0  !if chl>2 or chl<=0, nu remains 0
      ElseIf(chlconc .lt. 0.02) then   !for 0<chl<0.02, use Mie theory value
        nu = -1.0 
      ElseIf(chlconc .le. 2.0) then  !for 0.02<chl<2, use model [nu (-1,0)]
        nu = 0.5 * (log10(chlconc) - 0.3)
      Else  !all cases should be defined above; just a catch-all
        nu = 0.0 
      EndIf
!     define fraction of small particles
      if(chlconc .gt. 0) then
        alpha = 0.855*(0.5 - 0.25*log10(chlconc)) 
        alpha = max(0.0, min(alpha,1.0))  !alpha must be contained in (0,1)  
      else
        alpha = 0 !not defined for chl=0
      endif
cc
c     particulate a by the Morel model.
      ap = Acapp*chlconc**(Ecapp)
c     CDOM absorption 
c      get chl-dep absorption coefficient at 440
      ap440 = Acapp440*chlconc**Ecapp440
c      In Bricaud, acd.ne.0 when chl=0; we changed so that simplifies to 
c      pure water IFF chl=0; only matters for very small chl values!
!     acd = 0.2 * (aw440 + ap440) * Xcdom    ! as written in Bricaud(98)
      acd = 0.2 * (        ap440) * Xcdom    ! as written in Tech Notes (a=aw if chl=0)
c
c     particulate b by the Morel and Mantosena model
      if(chlconc .gt. 0) then
        bp =  0.416* chlconc**(0.766) * (wavelen/550.)**nu
      else
        bp = 0 !not defined for chl=0
      endif
      bpSmall = bp * alpha 
      bpLarge = bp * (1 - alpha)
c
      acomp(1) = awater
      acomp(2) = ap 
      acomp(3) = 0.0 
      acomp(4) = acd
      bcomp(1) = bwater
      bcomp(2) = bpSmall
      bcomp(3) = bpLarge
      atotal = acomp(1) + acomp(2)            + acomp(4)
      btotal = bcomp(1) + bcomp(2) + bcomp(3)
c
      return
c
 100  format(/5x,'The IOP routine "abNewCase1" is being used:'//2x,
     1'Absorption and scattering coefficients for case 1 water ',
     1'(Morel et al. bio-optical models) are used.'//
     28x,'component 1 is pure water'/
     38x,'component 2 is small particles'/
     38x,'component 3 is large particles'/
     48x,'component 4 is CDOM')
 102  format(//'Error in sub abcase1:  ncomp =',i2,
     1' must equal',i2,' for this version')
 104  format(/'Error in sub abcase1:  called with optical depth: ',
     1'iop =',i2)
 106  format(/5x,'The chlorophyll values at the requested output',
     1' depths are:'//8x,'depth      Chl'/10x,'(m)    (mg/m^3)'/)
 107  format(3x,2f10.3)     
 201  format(/5x,'The chlorophyll concentration is ',
     1'constant with depth with a value of ',F8.2, ' (mg/m^3)')
      end
!_______________________________________________________________________

      Subroutine getAEnewcase1(Acapp, Ecapp, Acapp440, Ecapp440)
c     returns Acapp and Ecapp for the current wavelen and at 440nm
c
      INCLUDE "DIMENS_XL.INC"

c     IO vars
      real Acapp, Ecapp, Acapp440, Ecapp440
      integer jwave
c
      character header*1, filename*120
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
      integer nwav
      real, allocatable,save,dimension(:) :: Apcoef,Epcoef
      real, allocatable,dimension(:) :: waveDat,ApDat,EpDat
c
c     local vars
      real, save :: A440, E440
      integer i, j, iwav, i1, i2
      real x
c
      integer, save :: kall
      data kall/0/
c     Note on dynamic storage:  MUST save data values between subroutine 
c     calls:
!*********************************************************************
c
      nwave = imisc(7)
      if(kall.eq.0) then
c*****load IOP coefficients on the first call:
        filename = trim(datafiles(6+3))   !stored as 3rd comp
        write(10,101) trim(filename)
        open(34, file=filename,status = 'old', err=200)
        Do i=1,10
          read(34,*) header
        Enddo
        nwav = 0
        Do 
          read(34,*, end=10) x
          if (x.le.0) goto 10
          nwav = nwav+1
        Enddo
  10    allocate(waveDat(nwav), ApDat(nwav), EpDat(nwav))  !dim to #data
        allocate(Apcoef(nwave), Epcoef(nwave))             !dim to #wavel in run
!
        rewind(34)
        Do i=1,10
          read(34,*) header
        Enddo
        Do i=1,nwav    !load all data
           read(34,fmt=*, end=20) waveDat(i),ApDat(i),EpDat(i)
        Enddo          
!       interp data to wavel grid
  20    i1 = 1   !initialize data bound
        Do i=1,nwave + 1
          if(i.le.nwave) then
            wavelen = wave(i)
          else
            i1 = 1
            wavelen = 440.0
          endif
c         find place in wavelength array
          Do j=i1,nwav-1
            if(waveDat(j).ge.wavelen) goto 30
          Enddo
          j = nwav
c         find factor for interpolation
 30       iwav = j
          if(iwav.eq.1) then
            i1=1
            i2=2
            x=0.0
          elseif(iwav.ge.nwav) then
            i1=nwav-1
            i2=nwav
            x=1.0
          else
            i1 = iwav-1
            i2 = iwav
            x = (wavelen-waveDat(i1))/(waveDat(i2)-waveDat(i1))
         endif
c        interpolate in DB
         If(i.le.nwave) then
           ApCoef(i) = (1-x) * ApDat(i1)  +  x * ApDat(i2)
           EpCoef(i) = (1-x) * EpDat(i1)  +  x * EpDat(i2)
         Else     !store 440nm result
           A440 = (1-x) * ApDat(i1)  +  x * ApDat(i2)
           E440 = (1-x) * EpDat(i1)  +  x * EpDat(i2)
         Endif
        Enddo
        deallocate(waveDat, ApDat, EpDat)
        close(34)
        kall = 1      
      endif
c
c*****return A and E coeff on all calls
      jwave = imisc(11)
      Acapp = ApCoef(jwave)
      Ecapp = EpCoef(jwave)
      Acapp440 = A440
      Ecapp440 = E440
!      write(87,'(i3,5f12.5)')jwave,wave(jwave),
!     1                       Acapp,Ecapp,Acapp440,Ecapp440

      return
  200 call nofile(10,"abNewCase1",filename)

  101 format(/8x, 'chlorophyll AE filename: ',a/)

      end subroutine

