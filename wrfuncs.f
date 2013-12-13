C     Last change:  LKS   5 Apr 2008    5:40 pm
c-------------------------------------------------------------------------
c
      function wrfchl(wave1,wave2)
c
c     This function evaluates the chlorophyll-fluorescence wavelength 
c     redistribution function, fC(wavelength1 to wavelength2), as defined
c     by Eq. (5.102), for use by routine wavedisc.
c
      COMMON /Cmisc/ imisc(30),fmisc(30)
      data kall/0/
      save
c
      if(kall.eq.0) then

c     chl fluorsence quantum efficiency:
      phichl = fmisc(23)
c
c     emission function parameters:
      sigmac = 10.6
      wavec0 = 685.0
      pi = fmisc(1)
      factor1 = 1.0/(sigmac*sqrt(2.0*pi))
      factor2 = 0.5/(sigmac*sigmac)
      kall = 1
      endif
c
c     excitation function, Eq. (5.103):
c
      if(wave1.lt.370.0 .or. wave1.gt.690.0) then
         gchl = 0.0
      else 
         gchl = 1.0
      endif
c
c     emission function, Eq. (5.104):
c
      hchl = factor1*exp(-factor2*(wave2 - wavec0)**2)
c
      wrfchl = phichl*gchl*hchl*wave1/wave2
c
      return
      end
c------------------------------------------------------------------------
c 
      function wrfcdom(wave1,wave2)
c
c     This function evaluates the CDOM-fluorescence wavelength 
c     redistribution function, fY(wavelength1 to wavelength2), as defined
c     by Eqs. (5.99) and (5.101), for use by routine wavedisc.
c
c*****NOTE:  The default parameter values for the CDOM "eta function" of
c     Eq. (5.101) correspond to those of Fig. 5.11 in Light and Water,
c     which is the same as Fig. 4a (sample FA7) of S. K. Hawes, "Quantum
c     fluorescence efficiencies of marine humic and fulvic acids," M.S.
c     Thesis, Marine Sci. Dept., Univ. of South Florida, 92 pages, 1992.
c     As shown in Hawes, there is considerable variability in these
c     parameter values for samples of CDOM extracted from Gulf of Mexico
c     waters, and almost no data are available from other natural
c     waters.
!     mod 9/29/04 by LKS to extend CDOM to work in UV range
c
      parameter(ndata=13)
      dimension wave(ndata),A0FA7(ndata)
c
c     data from Table 3 of Hawes (1992):
      data wave/310.,330.,350.,370.,390.,410.,430.,450.,470.,490.,
     1          510.,530.,550./
      data A0FA7/5.18e-5,6.34e-5,8.00e-5,9.89e-5,9.39e-5,10.48e-5,
     1       12.59e-5,13.48e-5,13.61e-5,9.27e-5,4.00e-5,1.00e-5,0.0/
      data A1/0.470/, A2/0.407/, B1/8.077e-4/,B2/-4.57e-4/
!      data wave1old/0.0/
c
      save
c
c     evaluate etaY for the given wavelengths
c
c     Check to see if we are at a new wavelength 
c     (cycle thru wavelenghts, strictly increasing as we go)
!      if(wave1.gt.wave1old) then  
c      get A0 for the current excitation wavelength
       if(wave1.ge.wave(1) .and. wave1.le.wave(ndata)) then        !only incl from 310 to 550nm
          A0 = yinterp(i1, ndata, wave1, wave, A0FA7)
          temp = (1.0/wave2 - A1/wave1 - B1)/(0.6*(A2/wave1 + B2))
          wrfcdom = A0*exp(-temp*temp)*wave1/wave2
       else
          A0 = 0.0
          wrfcdom = 0.0
       endif 
!       wave1old = wave1
!      endif
c      temp code below while debugging
          temp = (1.0/wave2 - A1/wave1 - B1)/(0.6*(A2/wave1 + B2))
          wrfcdom = A0*exp(-temp*temp)*wave1/wave2
!
      return
  88  format(2f8.1,1p,4E12.3)
      end
c------------------------------------------------------------------------
c
      function wrframan(wave1,wave2)
c
c     This function evaluates the Raman-scatter wavelength redistribution
c     function, fR(wavelength1 to wavelength2), as defined by Eqs. (5.92)
c     and (5.94) AS CORRECTED in the notes of 09 Nov 94.
c     Recoded 11 Dec 2011 and 16 Jan 2012 by CDM
c
      common /Cmisc/ imisc(30),fmisc(30)
c
      dimension A(4),fkappa(4),dkappa(4)
c
c     parameter values from Light and Water Table 5.3 (Walrafen's data)
c
      data A/0.41, 0.39, 0.10, 0.10/
      data fkappa/3250., 3425., 3530., 3625./
      data dkappa/210., 175., 140., 140./
      data kall/0/
c
      save
c
      if(kall.eq.0) then
      fl2 = 4.0*alog(2.0)
      factor = sqrt(fmisc(1)/fl2)*(A(1) + A(2) + A(3) + A(4))
      factor = 1.0/factor
      do j = 1,4
      dkappa(j) = 1.0/dkappa(j)
      end do
      kall = 1
      endif
c
c     compute kappa'' from wave1 and wave2
c
      fkappp = (1.0e7)*(1.0/wave1 - 1.0/wave2)
c
c     evaluate Eq. (5.92, corrected)
c
      fR = 0.0
      do j=1,4
      fR = fR+A(j)*dkappa(j)*exp(-fl2*((fkappp-fkappa(j))*dkappa(j))**2)
      end do
c
      wrframan = (1.0e7)*factor*fR/(wave2*wave2)

c     line added by CDM, 16 Jan 2012, to zero out very small wrframan
c     values, to prevent numerical problems far from the emission bands
      if(wrframan .lt. 1.0e-10) wrframan = 0.0

      close(63)
c
      return
      end


