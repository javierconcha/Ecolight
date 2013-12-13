C     Last change:  LKS  25 May 2008    7:32 pm
      subroutine loadsurface(windspd, wavelen)
c
      real windspd, wavelen
      real refr
c
      INCLUDE "DIMENS_XL.INC"
c
C     Common blocks used or defined here
      COMMON /Cmisc/ imisc(30),fmisc(30)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      Character*120 datadir, digitdir, spreadir,
     1              phasedir, surfdir,bottdir, Pdir
      COMMON /Cdirnames/ datadir,digitdir,spreadir,
     1                   phasedir,surfdir,bottdir, Pdir
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      COMMON /Crtsurf/ raw(mxmu,mxmu),taw(mxmu,mxmu),rwa(mxmu,mxmu),
     1                 twa(mxmu,mxmu)
c     storage of SURFWIND values to compare with pf used in selpfbb
      COMMON /Cbbquadchk/ m,iqpart     !loadsurface, inishrad & selpfbb

c     temporary local storage:
      character*120 surftitl, rtname
      character char2*2, char4*4
      integer imisc1(30), nmu, nunit
      integer lensfdir, mdum,iqdum,nraydum
      real rt2(mxmu,mxmu)
      real fmisc1(30), Udum
c
      data nusrt1/15/, nusrt2/16/, nunit/10/
c
c***********************************************************************
c***********************************************************************
c
c     Call common routine that returns index of refractrion at this wavelength
      call getrefr(wavelen, windspd, refr, ireturn)
      if(ireturn.ne.0) return    !flag that no new surface file is needed
c
c     Check that this refr is valid
      call checkrefr(refr, ir1, ir2)
c
c     GET text for index of refraction
      write(char4,fmt='(i4)') ir1
c
c     CALL common routine to find bracketing Windspeed files
      lensfdir = lenstr(surfdir) 
      surfname = surfdir(1:lensfdir) // 'windlist.txt'
      call slctsurf(i1, i2, iexact, xinterp, windspd, surfname)
c
c     read the two bracketting windspeed files
c
c      ***First file
      if(i1.lt.10) then
         write(char2,fmt='(i1,1x)') i1
      else
         write(char2,fmt='(i2)') i1
      end if
      lenchar2 = lenstr(char2)
c
c     concatenate the wind-speed character string with the data 
c     directory name to create the name of the file containing the
c     surface data for the given wind speed
      lensfdir = lenstr(surfdir) 
      surfname = surfdir(1:lensfdir) // 
     1           'surfwind_' // char4 //
     2           '.' // char2(1:lenchar2)
c
c
      open(nusrt1,file=surfname, form='formatted',status='old', err=999)

      read(nusrt1,fmt='(a)') surftitl
      read(nusrt1,402) Udum,mdum,iqdum,nraydum
      read(nusrt1,403) imisc1,fmisc1
      nmu  = imisc1(1)
c
c     check for array dimensions being larger than allowed
      if(nmu.gt.mxmu) then
         write(6,fmt='(" nmu =",i3," gt mxmu =",i3)') nmu,mxmu
         write(nunit,fmt='(" nmu =",i3," gt mxmu =",i3)') nmu,mxmu
         call HERR("LOADSURFACE","increase mu LIMIT in UI")  !stop run
      endif
c
      read(nusrt1,404) (fmu(i),i=1,nmu)
      read(nusrt1,404) (bndmu(i),i=1,nmu)
      read(nusrt1,404) (omega(i),i=1,nmu)
      read(nusrt1,404) (deltmu(i),i=1,nmu)
c
      call readhat(nusrt1, nmu, raw)
      call readhat(nusrt1, nmu, taw)
      call readhat(nusrt1, nmu, rwa)
      call readhat(nusrt1, nmu, twa)
      close(nusrt1)
c
      if(iexact.eq.1) goto 250
c
c     if the windspeed does not exactly matches that of the first
c     file (i.e., if iexact = 0), then read the second bracketing
c     file and interpolate with the first
      if(i2.lt.10) then
         write(char2,fmt='(i1,1x)') i2
      else
         write(char2,fmt='(i2)') i2
      end if
      lenchar2 = lenstr(char2)
c
c     concatenate the wind-speed character string with the data 
c     directory name to create the name of the file containing the
c     surface data for the given wind speed
      lensfdir = lenstr(surfdir) 
      surfname = surfdir(1:lensfdir) // 
     1           'surfwind_' // char4 //
     2           '.' // char2(1:lenchar2)
!      write(10,*) 'opening #2: ',trim(surfname)
      open(nusrt2,file=surfname, form='formatted',status='old', err=999)

      read(nusrt2,fmt='(a)') surftitl
      read(nusrt2,402) Udum,mdum,iqdum,nraydum
      read(nusrt2,403) imisc1,fmisc1
      nmu2 = imisc1(1)
c
c     check to see if nmu's are the same (same sized arrays)
      if(nmu.ne.nmu2) then
        write(nunit,*) 'ERR: Trying to iterate between two files ',
     1              'that have different dimensions!'
         call HERR("LOADSURFACE","Surface files have different nmu")  !stop run
      endif
c
         read(nusrt2,404) (fmu(i),i=1,nmu)
         read(nusrt2,404) (bndmu(i),i=1,nmu)
         read(nusrt2,404) (omega(i),i=1,nmu)
         read(nusrt2,404) (deltmu(i),i=1,nmu)
c
c     Load and get averages
      call readhat(nusrt2, nmu, rt2)
      call gethat(nmu, xinterp, raw, rt2)

      call readhat(nusrt2, nmu, rt2)
      call gethat(nmu, xinterp, taw, rt2)

      call readhat(nusrt2, nmu, rt2)
      call gethat(nmu, xinterp, rwa, rt2)

      call readhat(nusrt2, nmu, rt2)
      call gethat(nmu, xinterp, twa, rt2)
c
c     close all files
      close (nusrt2)
c
c*********************
 250  imisc(1) = imisc1(1)
c
      fmisc(1) = fmisc1(1)
      fmisc(2) = fmisc1(2)
      fmisc(3) = fmisc1(3)
c
c     store values to compare to DPF in inishamp
      m = mdum
      iqpart = iqdum
c*********************
c
      return
  999 call nofile(nunit, 'LOADSURFACE', surfname)   !err opening file
c 
  402 format (f6.2,3i3,i12)
  403 format (16i5,i12,3i5,10i2 / 10(e12.6,1x))
  404 format (10(e12.6,1x))

      end subroutine
c
!*******************************************************
      subroutine gethat(nhat, xinterp, hat1, hat2)
      INCLUDE "DIMENS_XL.INC"
      PARAMETER (mxhat=mxmu) 
c
      integer nhat
      real xinterp
      real hat1(mxhat,mxhat), hat2(mxhat,mxhat)
c
c	**** average values
      Do j=1,nhat
      Do i=1,nhat
        hat1(i,j)= (1.0 - xinterp)*hat1(i,j) +
     1             xinterp*hat2(i,j)
      Enddo
      Enddo
      end subroutine

!*******************************************************
      subroutine readhat(nusrt, nhat, hat)
      INCLUDE "DIMENS_XL.INC"
      PARAMETER (mxhat=mxmu) 
c
      integer nusrt, nhat
      real hat(mxhat,mxhat) 

c	**** read in pairs of arrays
         read(nusrt,fmt='(a)') rtname     !unique to EL
	DO I=1,nhat
	   READ(nusrt,404) (hat(I,J),J=1,nhat)
	end do

  404 format (10(e12.6,1x))
      end subroutine
