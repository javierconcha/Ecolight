      subroutine PAR
c
c     core routine on file PAR.f
c
c     This routine computes PAR as a function of depth (at the
c     same depths where irradiances are printed out), after
c     the entire run has been completed.  
c
c     This routine is intended for use only with multi-wavelength runs
c     that span the visible spectrum.
c
c     Change for H5.0:  PAR is defined using the scalar irradiance Eo.
c     However, it is sometimes estimated using the downwelling
c     plane irradiance Ed (e.g., this is what is done in MODIS).
c     PAR from Ed will be significantly less (typically about 3/4)
c     than PAR from Eo.  To enable this comparison, H5.0 now
c     computes both true PAR from Eo and PAR estimated from Ed.
c
      INCLUDE "DIMENS_XL.INC"
c
      Common /CPAR/fPar(0:mxz),Edpar(0:mxz),fKpar(mxz),
     1             Eoquant(0:mxz,mxwave)   !shared with EXCEL.f

      COMMON /CEospl/ nspl,zspl(mxz),Eospl(mxz,mxwave),
     1                E2spl(mxz,mxwave),Edspl(mxz,mxwave)
c     Cxcl holds arrays accumulated in routine storexcl
      COMMON /Cxcl/  nzxcl,zxcl(mxz),izxcl(mxz),
     1               axcl(mxz,mxwave,0:mxcomp),
     2               bxcl(mxz,mxwave,0:mxcomp),
     3               bbxcl(mxz,mxwave,0:mxcomp),
     4               Edxcl(0:mxz,mxwave),Euxcl(0:mxz,mxwave),
     5               Eodxcl(0:mxz,mxwave),Eouxcl(0:mxz,mxwave),
     5               Rxcl(0:mxz,mxwave),
     6               Raduxcl(0:mxz,mxwave),Radwxcl(mxwave)
c     CKxcl holds K-function arrays accumulated in routine Kfcn
      COMMON /CKxcl/ nzKxcl,zKxcl(mxz),fKdxcl(mxz,mxwave),
     1  fKuxcl(mxz,mxwave),fKoxcl(mxz,mxwave),
     2  fKnetxcl(mxz,mxwave),fKLuxcl(mxz,mxwave)
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave)
      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cirrad/ Eou(0:mxz),Eod(0:mxz),Eu(0:mxz),Ed(0:mxz), 
     1                fMUu(0:mxz),fMUd(0:mxz),fMUtot(0:mxz),R(0:mxz),
     2                E2(0:mxz)
      COMMON /CMISC/ IMISC(30),FMISC(30) 
c
c     for use in PAR:
      character Erootname*120
      common /Efilename/ Erootname
      COMMON /Ctitle/ ititle
      character ititle*120 
      real kfact, Eotemp, zk(mxz), zkmin(mxz), zkmax(mxz), Edsum
      dimension Eoband(0:mxz,mxwave), dwave(mxwave),Edbroadband(0:mxz)
c
      real wavemin,wavemax, Eo, E2
c
      logical IamEL
      external IamEL
c
      character*10 getHE5
      external getHE5
c     -------------------------------------------------------------

      nwave = imisc(7)
      ibotm = imisc(12) 
c
c     Eo at each output depth and wavelength was saved in Eospl for
c     use with inelastic scattering computations
c
c     integrate Eo over wavelength at each depth to get PAR
c
      wavemin = fmisc(26)! the min wavelength for PAR computation
      wavemax = fmisc(27) ! the max wavelength for PAR computation
      planck = 6.626e-34
      speed = 2.998e8
      fnmtom = 1.0e-9
      fmoltoumol = 1.0e6
      avagadro = 6.023e23
      factor = fnmtom*fmoltoumol/(planck*speed*avagadro)
c
c     *** Calculate in-air values
      parsum  = 0.
      do iwave=1,nwave
          dwave(iwave) = (waveb(iwave+1)- waveb(iwave))
          temp = wave(iwave) * dwave(iwave)
          Eotemp = Eouxcl(0, iwave) + Eodxcl(0, iwave)       !use Excel array since has in-air(wavel)
          Eoband(0,iwave) = Eotemp * dwave(iwave) 
          Eoquant(0,iwave) = factor * temp * Eotemp
          if(wave(iwave).ge.wavemin .and. wave(iwave).le.wavemax) then
            parsum =  parsum  + temp * Eotemp
            Edparsum = Edparsum + temp * Edxcl(0,iwave)
            Edsum = Edsum + Edxcl(0,iwave)*dwave(iwave)
          endif   
      enddo
c
      fpar(0) = factor*parsum
      Edpar(0) = factor*Edparsum
      Edbroadband(0) = Edsum

C     *** Calculate in-water values
      Do izz = 1,npirad
         parsum  = 0.0
         parsum2 = 0.0
         Edparsum = 0.0
         Edsum = 0.0
         iz = izirad(izz)
         if(izz.gt.nzKxcl) then     !protect case where EL is being run
           iz2 = 1    !dummy value to ensure we don't access invalid memory loc
         elseif(izz.eq.npirad .and. ibotm.ne.0) then    !if finite depth, use next to last depth at bottom
           iz2 = iz - 1
           zkmin(iz) = z(iz2)
           zkmax(iz) = z(iz)
         else
           iz2 = iz + 1
           zkmin(iz) = z(iz)
           zkmax(iz) = z(iz2)
         endif
         kfact = - 2.0 / ( z(iz2) - z(iz) )
         zk(iz) = 0.5*( z(iz)+z(iz2) )
         do iwave=1,nwave
          temp = wave(iwave)* dwave(iwave)
          Eo = exp(Eospl(iz,iwave))
          E2 = exp(E2spl(iz,iwave))
          Eoband(iz,iwave) = Eo * dwave(iwave)
          Eoquant(iz,iwave) = factor * temp * Eo
c         sum only over 400-700 nm to compute PAR
          if(wave(iwave).ge.wavemin .and. wave(iwave).le.wavemax) then
            parsum =  parsum  + temp * Eo
            Eo = exp(Eospl(iz2,iwave))
            parsum2 = parsum2 + temp * Eo
            Edparsum = Edparsum + temp * Edspl(iz,iwave)
            Edsum = Edsum + Edspl(iz,iwave)*dwave(iwave)
          endif
         end do
         fpar(iz) = factor*parsum
         Edpar(iz) = factor*Edparsum
         fKpar(iz) = kfact*(parsum2-parsum)/(parsum+parsum2)
         Edbroadband(iz) = Edsum
      End do
c
c     *** printout and digital file (at user requested depths)

c     open the file for scalar irradiance and PAR output
c     (this file is "Erootname.txt" in the digital directory)
!      nunitEo = 75
!      open(unit=nunitEo, file=Erootname)
!      write(nunitEo,*) trim(getHE5())," Run Title:  ", trim(ititle)
!      write(nunitEo,fmt='(2i5, 5x,"(nz, nwave)")') npirad,nwave
c
c     Band-integrated energy scalar irradiance (depth, wavelength)
!      write(nunitEo,'(/a)') "Band-integrated energy Eo [W/m^2]"
!      write(nunitEo,400) (wave(i),i=1,nwave)
      write(10,100) (wave(i),i=1,nwave)
!      write(nunitEo,501) (dwave(i),i=1,nwave)
!      write(nunitEo,502)    !depth header
c     write in-air values
!      write(nunitEo,401) (Eoband(0,i),i=1,nwave)
      write(10,402) (Eoband(0,i),i=1,nwave)
c     write in-water values
      do ip=1,npirad
         iz = izirad(ip)
         write(10,202) z(iz),(Eoband(iz,i),i=1,nwave)
!         write(nunitEo,202) z(iz),(Eoband(iz,i),i=1,nwave)
      enddo

c     Band-integrated quantum scalar irradiance (depth, wavelength)
!      write(nunitEo,'(/a)') 
!     1    "Band-integrated quantum Eo [(micromol photons)/(m^2 s)]"
!      write(nunitEo,400) (wave(i),i=1,nwave)
      write(10,200) (wave(i),i=1,nwave)
!      write(nunitEo,501) (dwave(i),i=1,nwave)
!      write(nunitEo,502)    !depth header
c     write in-air values
!      write(nunitEo,401) (Eoquant(0,i),i=1,nwave)
      write(10,402) (Eoquant(0,i),i=1,nwave)
c     write in-water values
      do ip=1,npirad
         iz = izirad(ip)
         write(10,202) z(iz),(Eoquant(iz,i),i=1,nwave)
!         write(nunitEo,202) z(iz),(Eoquant(iz,i),i=1,nwave)
      enddo
c
c
!     only calculate PAR if we have sufficient wavelengths
      IF((waveb(1).gt.wavemin).or.(waveb(nwave+1).lt.wavemax))then
          ! we do not have suff wavelengths; zero array and return
          fPAR(0)=-999.
          EdPAR(0)=-999.
          do ip=1,nzxcl
             iz = izxcl(ip)
             fPAR(iz)=-999.
             EdPAR(iz)=-999.
             fKPAR(iz)=-999.  
          end do
          write(10,110) wavemin,wavemax
          return    
      ENDIF

c     PAR(depth)
      write(10,300) wavemin,wavemax
!      write(nunitEo,299) wavemin,wavemax
c
!     PRINTOUT PAR
c     write in-air values
      write(10,301) fpar(0),Edpar(0),Edbroadband(0),
     1 fpar(0)/Edbroadband(0)
!      write(nunitEo,301) fpar(0),Edpar(0)
c     write in-water values
      do ip=1,npirad
         iz = izirad(ip)
         write(10,302) z(iz),fpar(iz),Edpar(iz),Edbroadband(iz),
     1 fpar(iz)/Edbroadband(iz)
!         write(nunitEo,302) z(iz),fpar(iz),Edpar(iz)
      enddo
!
!     PRINTOUT kPAR (in-water)
      if(IamEL()) then            !EL
        write(10,601) 
!        write(nunitEo,601) 
      else                        !HL
        write(10,600) 
!        write(nunitEo,600) 
      endif
      do ip=1,nzKxcl
         iz = izirad(ip)
         write(10,602) zkmin(iz), zkmax(iz), zk(iz), fKpar(iz)
!         write(nunitEo,303) zk(iz), fKpar(iz)
      enddo

!     close(nunitEo)
c
  299 format(/2x,'PAR as a function of depth (mumol phot/m^2/s) '/
     2 5x,'The wavelengths used for this PAR calculation',
     3' were ',f6.1,' to',f6.1,' nm'/
     4'      depth(m)     PAR(from Eo)      PAR(from Ed)')
  300 format(//2x,'PAR and broadband Ed as a function of depth'/
     2 5x,'The wavelengths used for these PAR and Ed calculations',
     3' were ',f6.1,' to',f6.1,' nm'//
     4'         depth   PAR(from Eo)      PAR(from Ed)       broadband E
     4d       PAR(Eo)/Ed'/,
     5'          (m) (mumol phot/m^2/s) (mumol phot/m^2/s)     (W/m^2)
     5     (mumol phot/s/W)'/)
  301 format(8x,"in air",1p,e15.4,3(3x,e15.4))
  302 format(f14.2,1p,e15.4,3(3x,e15.4))
  303 format(f11.3,1p,e15.4)
  100 format(//'  Energy Eo as a function of depth and wavelength'/
     1'     Units are W/m^2  [values are integrated over the wavelength
     1bands]'//' depth',(100f11.1))
  110 format(//2x,'PAR was not calculated for this run'/
     15x,'(PAR is meaningful only for runs covering the near-UV ',
     2'and visible spectrum)'/
     32x,'PAR range currently set to',f6.1,' to', f6.1,' nm'/)
  200 format(//'  Quantum Eo as a function of depth and wavelength'/
     1'     Units are (micromol photons)/(m^2 s)  [values are '
     2'integrated over the wavelength bands]'//' depth',(1000f11.1))
  202 format(f6.1,1p,(1000e11.3))
  400 format(5x,(1000f11.1))
  401 format("in air",1p,(1000e11.3))
  402 format(/"in air",1p,(1000e11.3))
  501 format(" dwave(nm)",1000(f6.1,5x))
  502 format(" depth(m)")
c
  600 FORMAT(//2x,'K-function for PAR (units of 1/meter)',
     1' (valid only when',
     1' zupper and zlower are closely spaced)'//
     2'    zupper    zlower',7X,'z',5X,
     3'KPAR(from Eo)'/)
  601 FORMAT(//2x,'K-function for PAR (units of 1/meter)'/,
     1' (valid only when',
     1' zupper and zlower are closely spaced)'//
     2'    zupper    zlower',7X,'z',5X,
     3'KPAR(from Eo, layer-averaged)'/)
  602 FORMAT(3F10.3,F10.5) 
c
      return
      end


