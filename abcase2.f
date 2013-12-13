C     Last change:  LKS  15 Nov 2010   10:31 pm
      subroutine abcase2(z,wavenm,ncomp, acomp,bcomp,atotal,btotal)
c
c     Example user routine on file abcase2
c
c     called by INISHAMP, BEAMc1, BEAMcz, and RHOTAU via ABSCAT
c
c     calls pureh2o, astar, bstar, and abmodel

c     This routine returns the absorption and scattering coefficients
c     a and b (in 1/meter), given the wavelength (350 to 800 nm)
C     and GEOMETRIC depth z (in m), for use in simulation of case 1 or
c     case 2 waters.
c
c     The water is modeled with 4 independent components:
c
c          component 1 is pure water
c          component 2 is chlorophyll-bearing particles (and covarying
c                      yellow matter)
c          component 3 is CDOM (yellow matter) that does NOT covary
c                      with the chlorophyll concentration; component
c                      3 is assumed to be nonscattering
c          component 4 is mineral particles
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
c     The a* specifications are:
c         iastropt = -1:  Subroutine called for aCDOM (a* routine not called)
c                    -2:  aCDOM solved as a percentage of aw+achl in "ab" routine (a* not called)
c                     0:  user-supplied data file read to get a* values
c                     1:  Pope & Fry absorption model used (pure water only)
c                     2:  Smith and Baker absorption model used (pure water only)
c                     3:  Prieur-Sathyenranath-Morel model used (Chlorophyll only)
c                     4:  Exponential model used (CDOM only)
c                     5:  Prieur-Sathyendranath-Morel model used (CDOM only)
c     The b* specifications are:
c         ibstropt =  0:  user-supplied data file read to get b* values
c                     1:  Power Law used (Chlorophyll only)
c                     2:  Linear model (Chlorophyll only)
c                     3:  Gordon-Morel model used (Chlorophyll only)
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
c     For component 3, abmodel is called to obtain the absorption 
c     by CDOM at depth z and wavelength wavenm.  The CDOM is assumed to
c     be non-scattering.
c
c     For component 4, abmodel is called to define the
c     wavelength and depth dependence of the mineral a and b
c
c*****NOTE:  Component 1 is pure water, therefore phase function 1 MUST
c            be the pure water phase function (file pureh2o.dpf).
c
c            Component 2 is particles/pigments, therefore phase
c            function 2 MUST be a particle phase function (e.g., 
c            avgpart.dpf or ff10540.dpf) or be specified by a bb/b.
c
c            Component 3 is nonscattering, so phase function 3 is 
c            irrevelant (its values are weighted by b3 = 0).  
c
c            Component 4 is mineral particles, so phase function 4
c            MUST be a particle phase function (e.g., ff11540.dpf,
c            which is for an index of refraction of 1.15, like
c            that of many mineral particles) or be specified by bb/b.
c
      INCLUDE "DIMENS_XL.INC"
c
      real compconc(mxcomp)  !local array 
      real astarcomp(mxcomp), bstarcomp(mxcomp)
      real astar1(mxcomp) !, bstar1(mxcomp)
      real temp, tempacdom

      dimension acomp(mxcomp),bcomp(mxcomp)
c     only zgeo is used (to printout concent profiles) in the following block
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      common/cmisc/ imisc(30),fmisc(30)
      COMMON /Cpirrad/ npirad,izirad(mxz)
c
c     Specific ab variables: only iastropt and ibstropt used here!
      Common /Castar/ iastropt, astarRef, astar0, asgamma
      Common /Cbstar/ ibstropt, bstarRef, bstar0, CompN, Plm,
     1                GAMm, GAMi
      Integer iastropt(mxcomp+2), ibstropt(mxcomp)
      Real astarRef(mxcomp+2), astar0(mxcomp+2)
      Real asgamma(mxcomp+2)
      Real bstarRef(mxcomp), bstar0(mxcomp), CompN(mxcomp)
      Real PLm(mxcomp), GAMm(mxcomp),GAMi(mxcomp)

c	itype indicates the method used to specify the component concentrations
c     = 0 if constant; =1 if subroutine, =2 if datafile
      common /Cconc/  itype(mxcomp)  
        integer itype
c
      data kall/0/, numcomp/4/, waveold/0.0/
c
c     Note on dynamic storage:  MUST save data values between subroutine 
c     calls:
      save

      if(kall.eq.0) then
c
c*****Initialization on the first call:
c
        write(10,100)
c       check for a valid call to this routine
         iop = imisc(5)
         if(numcomp.ne.ncomp) then
            write(10,102) ncomp, numcomp
         call HERR("abcase2","# components not compatible")  !stop run
         endif
         if(iop.ne.0) then
            write(10,104) iop
         call HERR("abcase2","depth type not compatible")  !stop run
         endif
c
c     Set up IOP models
c        call first time to trigger messages
         depth = 0.0
         wavelen = wavenm
c
c        load and init pure water routine
         call pureH2O(wavenm, awater,bwater,cwater)
c--------- DPF printout
         CALL PNTDPF(1)
c
c#### CHLOROPHYLL
c     print out the chl-spec header info iff constant option was selected
c     (if function or datafile was specified the appropriate routine will
c     print out the header info)
      icomp = 2
        call achlz(depth, wavelen, temp)       !calling to spawn msgs
      astarcomp(icomp) = astar(icomp, wavelen)     !calling to spawn msgs
      bstarcomp(icomp) = bstar(icomp, wavelen)     !calling to spawn msgs
c--------- DPF printout
      CALL PNTDPF(2)
c#### CDOM
c     print out the CDOM-spec header info iff constant option was selected
c     (if function or datafile was specified the appropriate routine will
c     print out the header info)
      icomp = 3
      call acdomsub(depth, wavelen, temp)           !calling to spawn msgs
      astarcomp(icomp) = astar(icomp, wavelen)      !calling to spawn msgs
c--------- DPF printout
      CALL PNTDPF(3)

c     #### loads compconc; delay call to control printout msgs
      call abmodel(depth, compconc)    !print any stray msgs

c#### MINERALS
c
c     print out the mineral-spec header info iff constant option was selected
c     (if function or datafile was specified the appropriate routine will
c     print out the header info)
      icomp = 4
      IF(itype(icomp).eq.0) then		!constant with depth
        write(10,401)  compconc(icomp)
      ELSEIF (itype(icomp) .EQ. 2) then  !read datafile
        call minzdata(icomp, depth, compconc(icomp))
      EndIf
      astarcomp(icomp) = astar(icomp, wavelen)
      bstarcomp(icomp) = bstar(icomp, wavelen)
c--------- DPF printout
      CALL PNTDPF(4)
c
c#### Printout profiles of components
c
      if(imisc(9).ge.0) then
        temp = wavenm   !handle if astarRef_CDOM==0, use first wavelength
        if(iastropt(3).gt.0) temp = astarRef(3)+0.5
        write(10,106) temp

        do iz=1,npirad
          depth = zgeo(izirad(iz))
          call abmodel(depth, compconc)
          call acdomsub(depth, temp, tempacdom)  !need aCDOM at refWave
          write(10,107)depth,(compconc(i),i=2,ncomp)
        end do
      endif
c
c        *** You may need to normalize astars and bstars
c            here depending on what models you are using!
      Do ii=2,ncomp
        If(iastropt(ii).eq.3) then	!use Prieur-Sathyendranath-Morel model
          astar1(ii) = astar(ii, astarRef(ii))
          If(abs(astar1(ii)).lt.1.e-7) then
            write(10,'(1x,a,i2,a,f4.0,a)') 'Astar for component ',
     1            ii,'is zero at reference wavelength ',
     2           astarRef(ii),' nm.'
            write(10,'(5x,a)') '(used to normalize absorption)',
     1            'for component ',ii,' is zero.'
            write(10,*) 'Absorption for component ',ii,' will be ',
     1            'set to zero for all depths and wavelengths!'
            write(10,*)
              Endif
          Endif
c         ** no model currently uses a normalized bstar
c        	bstar1(ii) = bstar(ii, bstarRef(ii))
         Enddo
c
      kall = 1
      endif
c
c-----End of initialization on first call.
c
c     Subsequent calls start here:
c
c       Save the routine input into dummy arrays
      depth = z
      wavelen = wavenm
c
c****** Note:  For multiwavelength ab models, it can be efficient to check
c       to see if the call is at the same wavelength as the previous call
c       or at a new one, and then perform any calculations that are needed
c       only once for a given wavelength.  This can be done via statements
c       of the form
c
      if(wavenm.ne.waveold) then
c     	 new-wavelength calculations here....
c
c       call routine pureH2O to get the pure water values
        call pureH2O(wavelen, awater,bwater,cwater)
c
c       call astar and bstar routines for each component (other than pure water)
c       using the form:  Xstar(icomp, wavenm)
        Do ii=2,ncomp
           astarcomp(ii) = astar(ii, wavelen)
           bstarcomp(ii) = bstar(ii, wavelen)
        Enddo

        waveold = wavenm
      endif
c
c**** Then perform the remainder of the a and b calculations for the current
c     depth and wavelength....
c
c     A call to abmodel will provide the concentration info for each
c     component:  compconc(1:ncomp)
      call abmodel(depth, compconc)
c
c     Store the pure water values into component 1
         acomp(1) = awater
         bcomp(1) = bwater
c
c     Call the astar and bstar routines to define the specific
c     scattering and absorption coefficients
c     Then use the concentrations to get the total a and b for
c     each component
c
      do ii=2, 4
c
c     **** if concentration is zero (or less than), set a and b to zero
       If(compconc(ii).le.0) then
              acomp(ii) = 0.0
              bcomp(ii) = 0.0
       Else
c     ----- calculate absorption ----
        If(iastropt(ii).eq.3) then	!use Prieur-Sathyendranath-Morel model
          If(abs(astar1(ii)).ge.1e-7) then   !make sure astar1 is nonzero
c           Get the Chl absorption, by the Morel model (L&W eq 3.27).
            acomp(ii) = 0.06 * astarcomp(ii)/astar1(ii)
     1                  * compconc(ii)**(0.65)
          Else
            acomp(ii) = 0.0
          Endif
        Else	!for all other options, use astarchl
            acomp(ii) = astarcomp(ii) * compconc(ii)
        Endif
c     ----- calculate scattering ----
        If(ibstropt(ii).le.0) then    !data file or nonscattering
c         Use the general form, b = conc * bstar
          bcomp(ii) = compconc(ii) * bstarcomp(ii)
        Else
c         Use the general form, b = conc**N * bstar
          bcomp(ii) = compconc(ii)**compN(ii) * bstarcomp(ii)
c         correct if c is returned as b from PowerLaw
          if(ibstropt(ii).eq.4) then 
            bcomp(ii) = bcomp(ii) - acomp(ii)
            if(bcomp(ii).lt.0) then 
               bcomp(ii) = 0.0
               write(10,500) ii, z 
            endif         
          endif         
        Endif

       Endif

      Enddo
c
c     Make sure b of CDOM is zero-ed out (b* should return as 0, so the next
c     line of code isn't strictly necessary)
      bcomp(3) = 0.0
c
c     Calculate Total a and b coefficients at this depth and wavelength
      atotal = 0.0
      btotal = 0.0
      do i=1,ncomp
         atotal = atotal + acomp(i)
         btotal = btotal + bcomp(i)
      enddo
      return
c
 100  format(/5x,'The IOP routine "abCASE2" is being used:'//5x,
     1'Absorption and scattering coefficients for case 1 or 2 water',
     2' are used.'//
     38x,'component 1 is pure water'/
     48x,'component 2 is pigmented particles (Chl)'/
     58x,'component 3 is CDOM (yellow matter)'/
     68x,'component 4 is mineral particles')
 102  format(//'Error in sub abCASE2:  ncomp =',i2,
     1' must equal',i2,' for this version')
 104  format(//'Error in sub abCASE2:  called with optical depth: ',
     1'iop =',i2)
 106  format(/5x,'Component concentration values at the requested',
     1' output depths'//,
     28x,'depth       Chl     aCDOM(',f5.1,')  Mineral'/,
     38x,' (m)      (mg/m^3)    (1/m)     (g/m^3)')
 107  format(3x,4f11.3)
 201  format(/5x,'The chlorophyll concentration is ',
     1'constant with depth with a value of ',F8.2, ' (mg/m^3)')
 301  format(/5x,'The CDOM absorption at each wavelength is ',
     1'constant with depth with a value of ',F8.2, ' (1/m)' )
 302  format(/5x,'CDOM absorption is modeled using'/
     18x,'aCDOM(lambda0) =',f7.3,' (1/m), independent of depth'/
     28x,'SCDOM          =',f7.4,' (1/nm)'/
     38x,'lambda0        =',f6.1,' (nm)') 
 303  format(/5x,'CDOM absorption will covary with chlorophyll ',
     1       'according to a Case 1 model')
 304  format(/5x,'A subroutine will be called to provide aCDOM as',
     1       ' a function of depth and wavelength')
 401  format(/5x,'The mineral concentration is ',
     1'constant with depth with a value of ',F8.2,' (g/m^3)' /)
 500  format(/5x, '** WARNING:  scattering for comp #', i2,
     1      ' at depth ', f6.2,/
     1      ' was calculated to be negative using c = b-a, ',
     3      'and b has been set to zero!'/)  
      end


