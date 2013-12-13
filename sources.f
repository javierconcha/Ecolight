C     Last change:  LKS  31 Aug 2013    6:34 pm
      subroutine sumsrc(zetanow,znow)
c
c     core routine on file sources.f
c
c     This routine computes the total effective internal source term
c     for the current optical depth and wavelength.  The total is the
c     sum of contributions by bioluminescence, chlorophyll
c     fluorescence, CDOM fluorescence, and Raman scatter.

c     Depending on which source terms are included in the run,
c     the contribution of each is added to the total source term.
c
c     sumsrc is called by drtdzs

c     sumsrc calls:
c        achlz    to get absorption by chlorophyll (via a user-selected
c                 routine as determined in the GUI)
c        acdomsub to get absorption by CDOM (via a user-selected
c                 routine as determined in the GUI)
c        s0biosub to the the bioluminesce source (via a user supplied
c                 routine, e.g., s0biolum, as selected in the GUI)
c
c        eval1D   to perform interpolation of irradiances (eval1d is in
c                 file dataintp.f)
c
      INCLUDE "DIMENS_XL.INC"
      COMMON /CEospl/ nspl,zspl(mxz),Eospl(mxz,mxwave),
     1                E2spl(mxz,mxwave),Edspl(mxz,mxwave)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      COMMON /CsourceE/ SrcM(mxmu),SrcP(mxmu)
      Common /Csource0/ ibiolum,ichlfl,icdomfl,iraman, ramanEXP
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
c     temporary local storage:
      dimension sum1(mxmu),factor1(mxmu),Eozi(mxwave)
      data i1/1/
      save i1 !index of last depth Eo data used
c
      nmu = imisc(1)
      iop = imisc(5)
      jwave = imisc(11)
      pi = fmisc(1)
      wave0 = fmisc(22)

c     zero the source arrays srcm = S-t(mu) and srcp = S+t(mu).  These are the
c     azimuthally averaged source terms that appear in the Riccati equations.
c
      do i=1,nmu
        srcm(i) = 0.0
        srcp(i) = 0.0
      end do

c     get factor1(i) = 1/(4*pi*beam c*mu(i)) at the current wavelength
c     (c was just computed in the call to rhotau; the 4*pi is for isotropic
c     emission)
c
      atotal = fmisc(19)
      btotal = fmisc(20)
      do i=1,nmu
         factor1(i) = 1.0/(4.0*pi*(atotal + btotal)*fmu(i))
      enddo
c
c     allow for routines that may want optical depth as input
c
      if(iop.eq.1) then
         depth = zetanow
      else
         depth = znow
      endif
c
c     compute Eo at all previous wavelengths, at the current depth,
c     if any inelastic terms are present.
c     The stored linear spline information (from routine  radanal)
c     is used to interpolate between the requested output depths
c
      if((jwave.eq.1) .or. (ichlfl.eq.0 .and. icdomfl.eq.0 .and.
     1                      iraman.eq.0)) goto 999
      do iwave=1,jwave-1
        Eozi(iwave) =  yinterp(i1, nspl, depth, zspl, Eospl(1,iwave))
        Eozi(iwave) = exp(Eozi(iwave))
      enddo

c
c*****CHLOROPHYLL FLUORESCENCE*******************************************

      if(ichlfl.ne.0) then

c     This section computes the effective internal source term
c     for chlorophyll fluorescence.
c     (The formulation of CM 24 June 94, p 13, is used.)
c
c     loop over all wavelengths shorter than the current wavelength and
c     compute the chl fluorescence contribution from each
c
      sum0 = 0.0
      do iwave=1,jwave-1
c
c     Now compute the contribution from wavelength iwave, using the
c     chlorophyll-specific absorption, the chl concentration at this 
c     depth, and the previously discretized chlorophyll wavelength
c     redistribution function 
c
      wavei = wave(iwave)
c
c     The generic call to the Chl absorption routine (the exact
c     routine called is specified via the GUI in Incfiles.for):

      call achlz(depth, wavei, achl)

      sum0 = sum0 + Eozi(iwave)*achl*fijchl(iwave,jwave)
c
      enddo
c
         do i=1,nmu
            Srcm(i) = Srcm(i) + sum0*factor1(i)
         end do
      i = 1
      endif   ! end of if(ichlfl.ne.0) test
c
c*****CDOM FLUORESCENCE**************************************************

      if(icdomfl.ne.0) then
c
c     This section computes the effective internal source term
c     for CDOM fluorescence.
c     (The formulation of CM 24 June 94, p 13, is used.)

c     loop over all wavelengths shorter than the current wavelength and
c     compute the CDOM fluorescence contribution from each
c
      sum0 = 0.0
      do iwave=1,jwave-1
c
c     Now compute the contribution from wavelength iwave, using the
c     parameterized CDOM absorption and the previously discretized
c     CDOM wavelength redistribution function
 
      wavei = wave(iwave)
c
c     The generic call to the CDOM absorption routine (the exact
c     routine called is specified via the GUI in Incfiles.for):

      call acdomsub(depth,wavei, abscdom)
      
      sum0 = sum0 + Eozi(iwave)*abscdom*fijcdom(iwave,jwave)

      end do

c     add the CDOM contribution to the source terms
c
         do i=1,nmu
            Srcm(i) = Srcm(i) + sum0*factor1(i)
         end do
      i = 1
      endif   ! end of if(icdomfl.ne.0) test
  88  format(2x,a,2f6.1,1p,3E12.3)
c
c*****RAMAN SCATTER******************************************************

      if(iraman.ne.0) then
c
c     This section computes the internal source term for Raman scattering.
c
c*****NOTE:  The present treatment of Raman scattering uses an 
c     azimuthally averaged source function that is equivalent to
c     the source function described in Appendix A of Mobley, et al.
c     1993, Comparison of numerical models..., Applied Optics 32(36),
c     7484-7504.  This azimuthally averaged source function guarantees
c     that the Raman contribution to IRRADIANCES will be computed
c     exactly.

      a0R = fmisc(21)
!dbg      write(10,*) '### Raman exp: ',RamanEXP

c     loop over all wavelengths shorter than the current wavelength and
c     compute the Raman-scatter contribution from each
c
      do i=1,nmu
         sum1(i) = 0.0
      do iwave=1,jwave-1
c
c     Compute the second moment E2 at the current depth and incident
c     (excitation) wavelength.  
      E2zi = yinterp(i1, nspl, depth, zspl, E2spl(1,iwave))
      E2zi = exp(E2zi)
c
c     Now compute the contribution from wavelength iwave, using the 
c     Raman scattering (L&W 'absorption') coefficient (5.89) and the
c     previously discretized Raman wavelength redistribution function
      wavei = wave(iwave)
      E2term = 0.3097*0.5*E2zi*(3.0*fmu(i)*fmu(i) - 1.0)
      sum1(i) = sum1(i) + (Eozi(iwave) + E2term) *
     1    a0R*(wave0/wavei)**ramanEXP * fijraman(iwave,jwave)  
c
      end do
      end do
c
         do i=1,nmu
            Srcm(i) = Srcm(i) + sum1(i)*factor1(i)
         end do
      i = 1
      endif ! end of if(iraman.ne.0) test
c
c*****BIOLUMINESCENCE****************************************************

  999 continue
      if(ibiolum.ne.0) then
c
c     this section computes the bioluminescence source contribution
c
      nwave = imisc(7)
      if(nwave.eq.1) then
c        monochromatic run:  use exact wavelength
         wavenm = fmisc(13)
c        The generic call to the bioluminescence routine (the exact
c        routine called is specified via the GUI in Incfiles.for):
         s0 = s0bioSub(depth,wavenm)
      else
c        average the source over the current wavelength band:
         wave1 = waveb(jwave)
         wave2 = waveb(jwave+1)
         ilam1 = ifix(wave1 + 0.5)
         ilam2 = ifix(wave2 + 0.5)
         count = 0.0
         sum0 = 0.0
         do ilam=ilam1,ilam2
            count = count + 1.0
            wavenm = float(ilam)
            sum0 = sum0 + s0bioSub(depth,wavenm)
         end do
         s0 = sum0/count
      endif
c
c     add bioluminescence to the source term
         do i=1,nmu
            Srcm(i) = Srcm(i) + s0*factor1(i)
         end do

      endif  ! end of if(ibiolum.ne.0) test
c********************************************************************
c     set S+t = S-t
      do i=1,nmu
         Srcp(i) = Srcm(i)
      enddo
      i =1
      return
      end
