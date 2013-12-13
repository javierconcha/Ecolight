C     Last change:  LKS   7 Aug 2008    6:52 pm
      Subroutine inishrad(nuphas)
c 
c     core routine on file inishrad.f
c 
c     called by RADIANCE [MAIN->RADIANCE->INISHRAD]
c
c     This routine initializes subroutine RADIANCE
c
c     calls zetatoz, ztozeta, qasky, and the "ab" routine
c 
      INCLUDE "DIMENS_XL.INC"
c
      CHARACTER pftitl*120
      Character txtstring*200, txtunits*200, txticomp*2
c
      common /Cbbopt/ ibbopt(mxcomp), bbtilde(mxcomp), 
     1                BfrefPL(mxcomp), Bf0PL(mxcomp), BfmPL(mxcomp)
      integer ibbopt
      real BfrefPL, Bf0PL, BfmPL
c       flag for calculating bb shared with IOP routine (calc bb iff flag=1)
      Common /CbbCalc/ ibbCalc
c
      COMMON /Cmisc/ imisc(30),fmisc(30) 
      COMMON /Crhotau/ rho(mxmu,mxmu),tau(mxmu,mxmu),
     1                 betatP(mxmu,mxmu,mxcomp),betatM(mxmu,mxmu,mxcomp)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Ciop/ acoef(mxz,0:mxcomp),bcoef(mxz,0:mxcomp),
     1		      atten(mxz),albedo(mxz), bbcoef(mxz,0:mxcomp)
      COMMON /Cpirrad/ npirad,izirad(mxz)

      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
c
c     storage of DPF file quad info for comparison to Surfwind file
      real Udum
      integer mdum(mxcomp), iqpartdum(mxcomp), nraydum
      integer iprob
c     storage of SURFWIND values to compare with pf used in selpfbb
      COMMON /Cbbquadchk/ m,iqpart     !loadsurface, inishrad & selpfbb
c
      dimension nuphas(mxcomp)

c     temporaray local storage: 
      dimension acomp(mxcomp),bcomp(mxcomp)
      dimension astore(mxz,mxcomp),bstore(mxz,mxcomp)
      dimension bbstore(mxz,mxcomp)
      character*17 betaPMname
c
      data kall/0/
c 
      save
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     load surface file for this wavelen
      wavelen = fmisc(13)
      windspd = fmisc(15)
      call loadsurface(windspd, wavelen)
c
c.....add dynamic depth calc depths IFF iDynZ is selected
c     (this improves calc accuracy when inelastic sources and inf bottom 
c      are selected)  i(27)=iDynZ; i(12)=ibotm; i(8) = isource
      If(imisc(27).eq.1 .and. imisc(12).eq.0 .and. imisc(8).gt.0) then
        iwave = imisc(11)
        call getDynZ(iwave)    ! this sets imisc(4)
      Endif
c
c     Compute the (input) quad-averaged radiances for the sky 
c     at this wavelength
      call qasky
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nmu = imisc(1)
      nz = imisc(4)
      iop = imisc(5)
      ncomp = imisc(6)
      nwave = imisc(7) 
      wavelen = fmisc(13)
      windspd = fmisc(15)
c******************************************************
c
c     Read previously discretized quad-averaged phase functions
c     from files nuphas1, nuphas2, ... (these are the *.dpf files
c     created by the routines in the ..\discpf directory)
c     The files were opened in routine radiance
c
c     If ibbopt > 0 , call selpfbb to select phase function with given bb/b
c
c     Don't load any file if ibbopt<0 (signals CDOM or other nonscattering components)
c
      do k=1,ncomp

        IF(ibbopt(k) .eq. 1) then
c         use constant bbfract to select the dpf if ibbopt=1
          bbfractm = bbtilde(k)
          call selpfbb(bbfractm,k)
          bbtilde(k) = bbfractm
        ELSEIF(ibbopt(k) .eq. 2) then
c         use bbfract(z=0) to select the initial dpf if ibbopt=2
          bbfractm = bbtilde(k)
        ELSEIF(ibbopt(k) .eq. 3) then
c         If bb/b set by a power law, calc and load new DPF
          call Bfcalc(k, wavelen, bbtilde(k))
        ELSEIF(ibbopt(k) .eq. 0) then
c         load the data from the appropriate phase function file
          rewind nuphas(k)
c         read header records from phase function file
          read(nuphas(k),'(a)') pftitl
c         the following line is new in H4.1
          read(nuphas(k),402) Udum,mdum(k),iqpartdum(k),nraydum
          read(nuphas(k),'(13x,f10.5)') bbtilde(k)
c         read the arrays of discretized phase functions
          read(nuphas(k),fmt='(a)') betaPMname
          do i=1,nmu
             read(nuphas(k),404) (betatP(i,j,k),j=1,nmu)
          end do
          read(nuphas(k),fmt='(a)') betaPMname
          do i=1,nmu
             read(nuphas(k),404) (betatM(i,j,k),j=1,nmu)
          end do
        ELSE
c         zero out array if component is nonscattering
          do j=1,nmu
          do i = 1,nmu
            betatP(i,j,k) = 0.0
            betatM(i,j,k) = 0.0
          end do
          end do

        ENDIF     ! end ibbopt(k) control

 100  end do    ! end k loop
c
c     Check to see that the discretized phase functions have
c     the same quads as the surfwind file
      iprob = 0
      Do k = 1, ncomp
        if(ibbopt(k).ne.0) goto 101
        if(mdum(k).ne.m) then
          iprob = 1
          write(10,*) 'mdum(k),m = ',mdum(k),m
          write(10,600) pfname(k)
        endif
        if(iqpartdum(k).ne.iqpart) then
          if(iprob.eq.0) iprob = -1
          write(10,602) pfname(k)
        endif
 101  Enddo
c     Kill run if quads do not match
      If(iprob.gt.0) then
        call HERR("inishrad","Quads of DPF files do not match surfwind")  !stop run
      Endif
c 
c     set up correspondence between optical and geometrical depth
c
c     if iop.eq.1, it is assumed that the ab model takes OPTICAL depth as input
c     if iop.eq.0, it is assumed that the ab model takes GEOMETRIC depth as input
c
      if(iop.eq.1) then
c
c     Model in being run in "optical-depth mode," which requires that the
c     run be for only one wavelength
      if(nwave.ne.1) then
         write(10,1030)
        call HERR("inishrad","depth type is not compatible")  !stop run
      endif
c
c     zeta AS READ contains OPTICAL depths.
c     Compute the GEOMETRIC depths corresponding (at this wavelength) to the
c     optical depths where output is requested. 
c 
c     Since we do not need to calc bbfrac for building up z to zeta table, 
c     set flag to instruct ab routine to NOT recalc b/bb for this part
      ibbCalc = 0

      call zetatoz 

c     Reset ab routine flag to calculate bb/b for depths as needed, but no the dpf
      ibbCalc = 1
c
c
C     This strange bit of code allows the ab routine to "overwrite" the number of
C     columns of "component" printout without requiring reading extraneous input
C     from Iroot.txt.  e.g., abCase1 has 2 components (water + chl&cdom).
C     To allow the output to segregate CHL and CDOM absorption, nconc (which is
C     given in Iroot.txt as 2 since only two concentration sets need be read) is
C     reset to 3 on first call to abcase1.       (3/6/00 by lks)
      If(kall.eq.0) then
         depth = 0.0
c        insert the call to the desired "ab" routine:
         call abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
         nconc = imisc(23)
c
c       Generate the header line for component printout IFF that printout is included
C       and set up to print out NCONC columns of data
         If(imisc(9).ge.0) then
           npr = min(10,nconc)
           txtstring='  iz   Opt Depth  Geo Depth'
           txtunits ='                     (m)      (1/m)'	!protect trailing blanks
           Do iipr=1,npr
             write(txticomp, '(i2)') iipr
             ilen = lenstr(txtstring)
             ilen2 = lenstr(txtunits)
             txtstring = txtstring(1:ilen) // '   Comp' // txticomp
             if (iipr.ne.npr)
     1         txtunits  = txtunits(1:ilen2) // '    (1/m)'
           Enddo
           ilen = lenstr(txtstring)
           txtstring = txtstring(1:ilen) // '    Total'
           ilen2 = lenstr(txtunits)
           txtunits  = txtunits(1:ilen2) // '     (1/m)'
         Endif
      Endif
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do iz=1,nz
c
c        call the ab model with optical depth
         depth = zeta(iz)
c        insert the call to the desired "ab" routine:
         call abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
c
         acoef(iz,0) = atotal
         bcoef(iz,0) = btotal
         Do jj = 1,nconc
         	acoef(iz,jj) = acomp(jj)
         	bcoef(iz,jj) = bcomp(jj)
         Enddo
         atten(iz) = atotal + btotal
         albedo(iz) = btotal/(atotal + btotal)
c
c        store a, b and bb at each depth for printout below
c        bb is computed from the scattering coefs and the bbtilde
c        values from the phase function files
         bbcoef(iz,0) = 0.
         DO jj=1,nconc
            astore(iz,jj)  = acomp(jj)
            bstore(iz,jj)  = bcomp(jj)
            bbstore(iz,jj) = bcomp(jj)*bbtilde(jj)
            bbcoef(iz,jj)  = bbstore(iz,jj)
            bbcoef(iz,0)   = bbcoef(iz,0) + bbstore(iz,jj)
         end do
      end do
c     
c     Reset ab routine flag to calculate bb/b for depths as needed AND the dpf
      ibbCalc = 2

      else
C
c     Model is being run in "geometric-depth mode," which can be for
c     one or more wavelengths. 
C     On the first call, zeta AS READ contains GEOMETRIC depths in meters.
c     Move these values to z, and then compute the corresponding
c     OPTICAL depths.  Finally, define a lookup table "zetaz"
c     for converting zeta to z for use in the Riccati integrations
c 
C     This strange bit of code allows the ab routine to "overwrite" the number of
C     columns of "component" printout without requiring reading extraneous input
C     from Iroot.txt.  e.g., abCase1 has 2 components (water + chl&cdom).
C     To allow the output to segregate CHL and CDOM absorption, nconc (which is
C     given in Iroot.txt as 2 since only two concentration sets need be read) is
C     reset to 3 on first call to abcase1.       (3/6/00 by lks)
      If(kall.eq.0) then
         depth = 0.0
c        insert the call to the desired "ab" routine:
         call abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
         nconc = imisc(23)
c
c       Generate the header line for component printout IFF that printout is included
C       and set up to print out NCONC columns of data
         If(imisc(9).ge.0) then
           npr = min(10,nconc)
           txtstring='  iz   Opt Depth  Geo Depth'
           txtunits ='                     (m)      (1/m)'	!protect trailing blanks
           Do iipr=1,npr
             write(txticomp, '(i2)') iipr
             ilen = lenstr(txtstring)
             ilen2 = lenstr(txtunits)
             txtstring = txtstring(1:ilen) // '   Comp' // txticomp
             if (iipr.ne.npr)
     1          txtunits  = txtunits(1:ilen2) // '    (1/m)'
           Enddo
           ilen = lenstr(txtstring)
           txtstring = txtstring(1:ilen) // '    Total'
           ilen2 = lenstr(txtunits)
           txtunits  = txtunits(1:ilen2) // '     (1/m)'
         Endif
      Endif
c  
c     Since we do not need to calc bbfrac for building up z to zeta table, 
c     set flag to instruct ab routine to NOT recalc b/bb fore this part
      ibbCalc = 0

      call ztozeta

c     Reset ab routine flag to calculate bb/b for depths as needed, but no the dpf
	ibbCalc = 1
c
      do iz=1,nz
c
c        call the ab routine with geometric depth
         depth = z(iz)
c
c        insert the call to the desired "ab" routine:
	   call abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
c
         acoef(iz,0) = atotal
         bcoef(iz,0) = btotal
         Do jj = 1,nconc
         	acoef(iz,jj) = acomp(jj)
         	bcoef(iz,jj) = bcomp(jj)
         Enddo
         atten(iz) = atotal + btotal
         albedo(iz) = btotal/(atotal + btotal)
c
c        store a,b, and bb at each depth for printout below
c        bb is computed from the scattering coefs and the bbtilde
c        values from the phase function files
         bbcoef(iz,0) = 0.
         DO jj=1,nconc
            astore(iz,jj)  = acomp(jj)
            bstore(iz,jj)  = bcomp(jj)
            bbstore(iz,jj) = bcomp(jj)*bbtilde(jj)
            bbcoef(iz,jj)  = bbstore(iz,jj)
            bbcoef(iz,0)   = bbcoef(iz,0) + bbstore(iz,jj)
         end do
      end do
c
c     Reset ab routine flag to calculate bb/b for depths as needed AND the dpf
	ibbCalc = 2

      endif 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Print out the absorption coefficients for the first
c     10 components and the total absorption in one block,
c     do the same for the scattering coef, the backscattering
c     coef, and then print a summary block of IOP data.  
c     Printout is at the user-requested depths.
c
c
C     This strange bit of code allows the ab routine to "overwrite" the number of
C     columns of "component" printout without requiring reading extraneous input
C     from Iroot.txt.  e.g., abCase1 has 2 components (water + chl&cdom).
C     To allow the output to segregate CHL and CDOM absorption, nconc (which is
C     given in Iroot.txt as 2 since only two concentration sets need be read) is
C     reset to 3 on first call to abcase1.       (3/6/00 by lks)
      If(kall.eq.0) then
         depth = 0.0
c        insert the call to the desired "ab" routine:
	  call abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
      	  nconc = imisc(23)
c
c       Generate the header line for component printout IFF that printout is included
C       and set up to print out NCONC columns of data
	 If(imisc(9).ge.0) then
	   npr = min(10,nconc)
	   txtstring='  iz   Opt Depth  Geo Depth'
	   txtunits ='                     (m)      (1/m)'	!protect trailing blanks
	   Do iipr=1,npr
	 	write(txticomp, '(i2)') iipr
	 	ilen = lenstr(txtstring)
	 	ilen2 = lenstr(txtunits)
	 	txtstring = txtstring(1:ilen) // '   Comp' // txticomp
	 	if (iipr.ne.npr)
     1            	  txtunits  = txtunits(1:ilen2) // '    (1/m)'
	   Enddo
	   ilen = lenstr(txtstring)
	   txtstring = txtstring(1:ilen) // '    Total'
	   ilen2 = lenstr(txtunits)
          txtunits  = txtunits(1:ilen2) // '     (1/m)'
	 Endif
      Endif
c	Don't print out component a's and b's if "minimal" printout is selected
	if(imisc(9).ge.0) then
		write(10,1124) wavelen,txtstring, txtunits
		do iii=1,npirad
		   izz = izirad(iii)
		   write(10,1125) izz,zeta(izz),z(izz),
     1                 (astore(izz,i),i=1,npr),acoef(izz,0)
		end do
c
c     loop over depth again to print out scattering coefficients.
		write(10,1126) wavelen,txtstring, txtunits
		do iii=1,npirad
		   izz = izirad(iii)
		   write(10,1125) izz,zeta(izz),z(izz),
     1                 (bstore(izz,i),i=1,npr),bcoef(izz,0)
		end do
c
c     loop over depth again to print out backscattering coefficients,
		write(10,1128) wavelen,txtstring, txtunits
		do iii=1,npirad
		   izz = izirad(iii)
		   write(10,1129) izz,zeta(izz),z(izz),
     1                 (bbstore(izz,i),i=1,npr),bbcoef(izz,0)
		end do
	endif
c
c     loop over depth again to print out remaining quantities.
      write(10,1024) wavelen
      do iii=1,npirad
         izz = izirad(iii)
         write(10,1025) izz,zeta(izz),z(izz),acoef(izz,0),bcoef(izz,0),
     1                 atten(izz),albedo(izz),bbcoef(izz,0),
     2                 bbcoef(izz,0)/bcoef(izz,0)
      end do
c
c     save acoef and atten at bottom of water column (for use in
c     routine infbotm, if ibotm = 1)
      fmisc(11) = acoef(nz,0)
      fmisc(12) = atten(nz)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      kall = 1
       
      return
c 
c     Format statements 
C
  402 format (f6.2,2i3,i12 )
  403 format (16i5,i12,3i5,10i2 / 10(e12.6,1x))
  404 format(10(e12.6,1x))

 600  format(/'Error in INISHRAD: surfwind and DPF files must have the'
     1'same spatial quads!',/'The MU quads read from phase function ',
     25x,a, /10x,' do not match the surfwind file')
 601  format(/'Error in INISHRAD: surfwind and DPF files must have the'
     1'same spatial quads!',/'The PHI quads read from phase function ',
     25x,a, /10x,' do not match the surfwind file')
 602  format(/'Error in INISHRAD: surfwind and DPF files must have the'
     1'same spatial quads!',/'The quad TYPE read from phase function ',
     25x,a, /10x,' do not match the surfwind file')

 1024 format(/2x,'Summary of Inherent Optical Properties',
     3 ' at ',f6.1,' nm',
     1//'  iz   Opt Depth  Geo Depth',3X,'total a',3X,
     2'total b',4X,'total c    albedo   total bb  total bb/b'/
     3     21x,'(m)',8x,'(1/m)',5X,'(1/m)',6X,'(1/m)',15x,'(1/m)'/)
 1025 FORMAT(i4,f10.3,f11.3,2x,f10.4,1X,f9.4,3X,f8.4,f10.4,2f11.5)
 1124 format(//2x,'Absorption Coefficients of Individual Components',
     3 ' at ',f6.1,' nm',
     1 //
     $ a, /, a, /) 
 1125 format(i4,f10.3,f11.3,2x,11f9.4)
 1126 format(//2x,'Scattering Coefficients of Individual Components',
     3 ' at ',f6.1,' nm',
     1 //
     $ a, /, a, /) 
 1128 format(//2x,'Backscattering Coefficients of Individual Components',
     3 ' at ',f6.1,' nm',
     1/' (computed from the scattering coefficients and the phase ',
     2'functions)'//
     $ a, /, a, /) 
 1129 format(i4,f10.3,f11.3,2x,11f9.5)
 1030 format(//'STOP:  Attempt to run model in multiwavelength mode',
     1' with OPTICAL depth output.  MUST use GEOMETRIC depths.')
c
      END
