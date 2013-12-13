C     Last change:  LKS  16 Nov 2010    9:47 am
      subroutine abcase1(z,wavenm,ncomp, acomp,bcomp,atotal,btotal)
c
c     user routine on file abcase1.f   
c
c     called by INISHAMP, BEAMc1, BEAMcz, and RHOTAU
c
c     calls pureh2o, astar, bstar, and abmodel
c
c     This routine returns the absorption and scattering coefficients
c     a and b (in 1/meter), given the wavelength (350 to 700 nm)
C     and GEOMETRIC depth z (in m), for use in simulation of case 1
c     waters.
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
c     a call to the astar and bstar routines.  Both astar and bstar are
c     functions which can be called in the form (ex for chlorophyll):
c     	astarchl = astar(icomp, wavenm)
c     The Frontend specified options/params for astar/bstar are passed
c     via a common block from initial (where they are read) directly to
c     the function. 
c
c     In GENERAL the a* specifications are:
c         iastropt =  0:  user-supplied data file read to get a* values
c                     1:  Pope & Fry absorption model used 
c                             (pure water only; astar not called)
c                     2:  Smith and Baker absorption model used 
c                             (pure water only; astar not called)
c                     3:  Prieur-Sathyenranath-Morel model used (Chlorophyll only)
c                     4:  Exponential model used (CDOM only)
c     The b* specifications are:
c         ibstropt =  0:  user-supplied data file read to get b* values
c                     1:  Power Law used (Chlorophyll only)
c                     2:  Linear model (Chlorophyll only)
c                     3:  Gordon-Morel model used (Chlorophyll only; bstar not called)
c
c     This Case1 model ALWAYS uses the bio-optical models of Morel (Prog.
c     Oceanogr. 26, 263, 1991) and of Gordon and Morel (Lec. Notes
c     on Coastal and Estuarine Studies, vol 4, Springer, 114p, 1983)
c     to convert the chl concentration into particulate a and b values.
c     [i.e., for Chlorophyll iastropt is ALWAYS 3 and ibstropt is ALWAYS 3] 
c
c     Descriptions of the components of this model:
c
c     For component 1, routine pureH2O is called to get the pure-water
c     absorption and scattering values.
c  
c     For component 2, a user-supplied subroutine (inserted via the 
c     call to abmodel) is called to obtain the chlorophyll 
c     concentration chlconc (in mg/m^3) at geometric depth z.
c
c*****NOTE:  Component 1 is pure water, therefore phase function 1 MUST
c            be the pure water phase function (file pureh2o.dpf).
c
c            Component 2 is particles/pigments, therefore phase
c            function 2 MUST be a particle phase function (e.g., 
c            avgpart.dpf or ff10540.dpf) or be specified by a bb/b.
c
      INCLUDE "DIMENS_XL.INC"

!      COMMON /CConstConc/ compconc(mxcomp)
      common /Cconc/  itype(mxcomp)  
      integer itype
      real compconc(mxcomp), acd, temp

      dimension acomp(mxcomp),bcomp(mxcomp)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      common/cmisc/ imisc(30),fmisc(30)
      COMMON /Cpirrad/ npirad,izirad(mxz)
c
c
      data kall/0/, numcomp/2/, waveold/0.0/
c
c     Note on dynamic storage:  MUST save data values between subroutine 
c     calls:
      save
c
      if(kall.eq.0) then
c*****Initialization on the first call:
         write(10,100)
c
c        load and init pure water routine
         call pureH2O(wavenm, awater,bwater,cwater)
c--------- DPF printout
         CALL PNTDPF(1)
c
c        reset number of concentrations for printout so that CDOM can be
c        printed separately
         nconc = 3
         imisc(23) = nconc
c        check for a valid call to this routine
         iop = imisc(5)
         if(numcomp.ne.ncomp) then
            write(10,102) ncomp, numcomp
         call HERR("abcase1","# components is not compatible")  !stop run
         endif
         if(iop.ne.0) then
            write(10,104) iop
         call HERR("abcase1","depth type is not compatible")  !stop run
         endif
cc
c         insert a call to the desired Chl(z) routine:
          depth = 0.
          wavelen = wavenm
c
c#### CHLOROPHYLL
c     print out the chl-spec header info iff constant option was selected
c     (if function or datafile was specified the appropriate routine will
c     print out the header info)
      icomp = 2
      call achlz(depth, wavelen, temp)       !calling to spawn msgs
c
c     Get the chlorophyll-specific absorption coefficient
c     at 440 nm for later normalization (the Morel model for achl
c     requires that achl* be normalized to 1 at 440 nm)
c     ASTAR has two arguments:  component# and wavelength
      ac440 = astar(icomp, 440.0)
      If(abs(ac440).lt.1.e-7) then
        write(10,*) 'Astar for component 2 ',
     1              'is zero at reference wavelength 440 nm.'
        write(10,'(5x,a)') '(used to normalize absorption)',
     1              'for component 2 is zero.'
        write(10,*) 'Absorption for component 2 will be ',
     1              'set to zero for all depths and wavelengths!'
        write(10,*)
      Endif
c     
      bcomp(icomp) = bstar(icomp, wavelen)    !calling to trigger mgs
c--------- DPF printout
      CALL PNTDPF(2)
c
      call abmodel(depth, compconc)    !print any stray msgs
c
c#### Printout profiles of components
      if(imisc(9).ge.0) then
        write(10,106)
        do iz=1,npirad
          depth = zgeo(izirad(iz))
          call abmodel(depth, compconc)
          write(10,107) depth,compconc(2)
        end do
      endif
c 
      kall = 1
      endif
c
c-----End of initialization on first call.
c
c     Subsequent calls start here:
c
c     If this is the first call at this wavelength, call routine
c     pureH2O to get the pure water values and call astarchl to get
c     ac (ac is normalized to ac(440nm))
c
      wavelen = wavenm
      if(wavelen.ne.waveold) then
      
         call pureH2O(wavelen, awater,bwater,cwater)

c         Make sure avoid division by zero error
          If(abs(ac440).ge.1.e-7) then
          ac = astar(2, wavelen) / ac440
          Else
              ac = 0.0
          Endif
 
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
c     Now get the "particulate" a by the Morel model.
c     Note that the Morel model gives the TOTAL absorption (including the
c     water contribution); define ap = atotal - awater
c
cc      ap = (awater + 0.06*ac*chlconc**(0.65))*(1.0 +
cc     1      0.2*exp(-0.014*(wavelen - 440.))) - awater
      ap = 0.06*ac*chlconc**(0.65)
c
c     CDOM absorption modified from the PS-Morel model
c     aCDOM = 0.2 * ap(440nm) * exp(-0.014*(wavelen - 440.))
c     (ac is by defined ot be 1 by P&S)
      acd = 0.012*chlconc**(0.65)*exp(-0.014*(wavelen - 440.))
c
c     particulate b by the Gordon-Morel model
      bp = (550.0/wavelen)*0.30*chlconc**(0.62)
c
      acomp(1) = awater
      acomp(2) = ap
      acomp(3) = acd
      bcomp(1) = bwater
      bcomp(2) = bp
      atotal = acomp(1) + acomp(2) + acomp(3)
      btotal = bcomp(1) + bcomp(2)
c
      return
c
 100  format(/5x,'The IOP routine "abcase1" is being used:'//5x,
     1'Absorption and scattering coefficients for case 1 water ',
     1'(Gordon and Morel bio-optical models) are used.'//
     28x,'component 1 is pure water'/
     38x,'component 2 is pigmented particles'/
     38x,'component 3 is CDOM that covaries with chlorophyll')
 102  format(//'Error in sub abcase1:  ncomp =',i2,
     1' must equal',i2,' for this version')
 104  format(/'Error in sub abcase1:  called with optical depth: ',
     1'iop =',i2)
 106  format(/5x,'Chlorophyll concentrations at the requested output',
     1' depths'//8x,'depth      Chl'/8x,' (m)    (mg/m^3)'/)
 107  format(3x,2f10.3)     
 201  format(/5x,'The chlorophyll concentration is ',
     1'constant with depth with a value of ',F8.2, ' (mg/m^3)')
      end


