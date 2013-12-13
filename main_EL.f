      program main
c
c     core program on file main.f
c
c     This is the main program for running Ecolight Version 5.2
c
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     +                                                                   +
C     +    ECOLIGHT 5.2  Copyright (c) 2013 by Curtis D. Mobley           +
C     +                                                                   +
c     +   ECOLIGHT IS EXPERIMENTAL AND IS LICENSED "AS IS" WITHOUT        +
C     +   REPRESENTATION OF WARRANTY OF ANY KIND, EITHER EXPRESS OR       +
C     +   IMPLIED.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE     +
C     +   OF ECOLIGHT IS WITH THE USER.  ECOLIGHT IS NOT FAULT TOLERANT.  +
C     +                                                                   +
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c     This computer program, named ECOLIGHT and consisting of various main
c     programs and subroutines hereafter referred to collectively as 
c     "ECOLIGHT", is licensed to the User on a non-exclusive,
c     non-transferable basis in accordance with the License Agreement on 
c     file HE52\Documents\HE5LicenseAgreement.pdf.
c
c==========================================================================
c
c     The file "DIMENS.INC" contains parameter statements that give the
c     maximum array dimensions for the fundamental quantities such as
c     polar angle, depth, and wavelength.  All other dimensions can be
c     derived from these. DIMENS_XL.INC is the link to 
c     ..\COMMON\DIMENS_XL.INC allowing easier adaption to different
c     system file path conventions
c
      INCLUDE "DIMENS_XL.INC"
c
      INCLUDE "MAIN_COMMON_EL.INC"  !contains declaration of all COMMONS
c
c	Set the printout file unit # to be 10
      iout = 10
c
c     write out Copyright info as header (routine at bottom of file)
      call cpright(6)
c
      write(6,*) "ECOLIGHT IS RUNNING, PLEASE WAIT..."
      write(6,*) "    (To stop the run before normal termination,"
      write(6,*) "    kill this command window.  Output will be lost)"

c     Initialize the run:  read the runtime input and set various defaults
c
      CALL INITIAL

c     **** ECOLIGHT CODE  ****
c     Store user-requested max depth index to reset each wavelength
      nzinit = imisc(4)
c     ****
c     Discretize the wavelength redistribution functions if inelastic
c     scattering is included in the run
c
      if(ichlfl.ne.0 .or. icdomfl.ne.0 .or.iraman.ne.0) call wrfdisc
c
c     Report time of initialization
c     Optional time stamp routine for Lahey Fortran only
      call date_stamp(iout,2, 0, 0)  ! calculate initialization runtime
c
c     Loop over all wavelengths to solve the RTE:
      nwave = imisc(7)
      nwskip = imisc(26) ! =1 for every wavelength; = 2 for every other wavelength
      if (nwskip .ne. 1) write(10,'(/2x,a,i2,/)')
     1           'Alternating wavelengths solved, nwskip = ', nwskip
      do jwave=1,nwave, nwskip
           imisc(11) = jwave
           fmisc(13) = wave(jwave)
c
           write(6,203) jwave, nwave, wave(jwave)
           write(10,100) jwave,waveb(jwave),waveb(jwave+1),wave(jwave)
c
c          solve for the radiances at the current wavelength
c
c          **** next line is specific to ECOLIGHT (variable depth options)
           imisc(4) = nzinit
c
           call radiance
c
c          analyze the results at the current wavelength
c
           call radanal
c
c          Report time of each wavelength
c          Optional time stamp routine for Lahey Fortran only
           call date_stamp(iout,3, jwave, nwave)  ! calculate incremental runtime
c
      enddo  !jwave
c
c     Run last wavelength if not already run  !new in v5
c     (only happens if nwskip=2 and number of wavelengths in run is even)
      if(imisc(11).lt. nwave) then
           jwave = nwave
           write(6,203) jwave, nwave, wave(jwave)
           imisc(11) = jwave
           fmisc(13) = wave(jwave)
           write(10,100) jwave,waveb(jwave),waveb(jwave+1),wave(jwave)
c          solve for the radiances at the current wavelength
c          **** next line is specific to ECOLIGHT (variable depth options)
           imisc(4) = nzinit
           call radiance
c          analyze the results at the current wavelength
           call radanal
c          Report time of each wavelength
           call date_stamp(iout,3, jwave, nwave)  ! calculate incremental runtime
      endif
c
c     if we Alternated wavebands, use linear interpolation to fill in
c     the missing values
c
      if(nwskip.gt.1) call FillWave
c
c     If the run covers at least the 400-700 nm range, then compute
c     PAR, CIE chromaticity coordinates for various quantities, and
c     the Secchi depth.
         call par         !call for all cases to calc and store values
      if(waveb(1) .le. 400.0 .and. waveb(nwave+1) .ge. 700.0) then
c         call CIExyY
         call Color
         call Secchi
      endif  ! end PAR and other visible-band calculations
c
c     Write the multi-wavelength file for spreadsheet postprocessing,
c     if requested
      iwrtssM = imisc(21)
      if(iwrtssM .ne. 0) call wrtxcelM
c
c     Optional time stamp routine for Lahey Fortran only
      call date_stamp(iout,1, nwave, nwave)  ! 1 indicated "end" of run
c
      !close all file connections
      call CLOSEALL

c     end run
      print *, 'The EcoLight run is complete.'
      print *, 'You may now close this window.'

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  100 format(//2x,'* * * * * Output for wavelength band',i3,
     1' (',f5.1,' to ',f5.1,' nm; nominal wavelength =',f6.1,
     2' nm)  * * * * *')
 203  format(2x,'Now beginning wavelength ',i3, ' of ',i3,
     1      ' at nominal wavelength ',f6.1, ' nm')
      end
c

