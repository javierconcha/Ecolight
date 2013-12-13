C     Last change:  LKS  25 Oct 2013    0:44 am
!     CONTAINS routines called by loadsurface_XL that are common
!     to both EL and HL versions

      subroutine slctsurf(i1, i2, iexact, xinterp, windspd, filen)
c
c     input:
c     windspd:  desired windspeed
c     filen:    name of file containing list of avail windspeeds
c
c     ouput:
c     i1, i2:  bracketing windspeed files
c     iexact:  one file matches exactly iff iexact=1
c     xinterp: interp fraction; if=0, i1 matches windspd
c
      parameter (mxws=50)
      integer i1, i2, iexact, iws(mxws) 
      real xinterp, windspd  
      Character*120 filen
c     temporary variables
      integer nfiles
c 
      data nusrt/7/
c********************************************************************
c     open list file to read in existing files
      open(nusrt, file=filen,status='old', err=999)

      Do ii=1,mxws
        read(nusrt,end=555,fmt='(i3)') iws(ii)
        if(iws(ii).eq.-1) goto 555
      Enddo
 555  close(nusrt)
      nfiles = ii-1
c
c     locate the available windspeeds that bracket the requested windspeed
c     (use the nearest value if windspeed is less than or greater than
c     the available values)
      iexact = 0
cc
      if(windspd .lt. float(iws(1))) then
c        error: U < 0; print message and use U = 0
         write(10,787) windspd,iws(1)
  787 format(' WARNING:  The requested wind speed of U =',f8.3,
     1' is less than the minimum available speed.  U =',i3,
     2' will be used.')
         i = 1
         iexact = 1
     
      elseif(windspd .gt. float(iws(nfiles))) then
c        error: U > 15; print message and use U = 15
         write(10,788) windspd,iws(nfiles)
  788 format(' WARNING:  The requested wind speed of U =',f8.3,
     1' is greater than the maximum available speed.  U =',i3,
     2' will be used.')
         i = nfiles + 1
         iexact = 1

      else
         do i = 1,nfiles
            if(windspd .le. float(iws(i))) goto 200
         enddo
         i = nfiles + 1
      endif
cc
  200 continue
      if(i.eq.1) then
         i1 = 1
         i2 = 1
         xinterp = 1.0
         iexact = 1
      elseif(i.eq.nfiles + 1) then
         i1 = nfiles
         i2 = nfiles
         xinterp = 0.0
         iexact = 1
      else
         i1 = i - 1
         i2 = i
         xinterp = (windspd - float(iws(i1)))/float(iws(i2) - iws(i1))
c        check to see if at one of the bracketing U's
         if(xinterp .lt. 1.0e-3) then     !just use first file
           iexact = 1
         elseif(xinterp .gt.0.999) then   !just use 2nd file
           i1 = i2
           iexact = 1
         endif
      endif
c
      i1 = iws(i1)
      i2 = iws(i2)
      return

  999 call nofile(10, 'SLCTSURF', filen)
      end subroutine

!-----------------------------------------------------------------------
                 
      subroutine getrefr(wavelen, windspd, refr, ireturn)
c
c     INPUT:
      real wavelen, windspd
c     OUPUT:
      integer ireturn
      real refr
c
c     TEMPORARY VARIABLES:
      real S, T, refrOld
c
      INCLUDE "DIMENS_XL.INC"
C     Common blocks used or defined here
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
      data nunit/10/, kall/0/
      save kall,refrOld
c********************************************************************
c
      ireturn = 0   !will be set to 1 if no surfwind file needs to be loaded

c     CALL common routine to get index of refraction at this wavelength
      S = fmisc(28) 
      T = fmisc(29)
      If(fmisc(30).gt.0) then      !constant refr
        refr = fmisc(30)
      Else                        !calc refr
        refr = Surfn(wavelen, S, T)
      Endif
c
      if(kall.gt.0) then   
c       don't print on first call by initial 
        write(nunit,401) refr
c       only load surface once if const refr or same file
        If(fmisc(30).gt.0) ireturn=1  
        if(nint(1000*refr).eq.nint(1000*refrOld)) ireturn=1 
      endif
      kall = 1
      refrOld = refr
c
  401 format (//2x,'Index of refraction at this wavelength is ',f5.3)
      end subroutine
!-----------------------------------------------------------------------
      subroutine checkrefr(refr, ir1, ir2)
!     This routine checks to see if refr is contained in set
      real refr
      integer kall, ir1, ir2
c
      INCLUDE "DIMENS_XL.INC"
C     Common blocks used or defined here
      COMMON /Cmisc/ imisc(30),fmisc(30)
      Character*120 datadir, digitdir, spreadir,
     1              phasedir, surfdir,bottdir, Pdir
      COMMON /Cdirnames/ datadir,digitdir,spreadir,
     1                   phasedir,surfdir,bottdir, Pdir
c
      parameter (mxn=100)
      integer irefr(mxn), nfiles, nrefr 
      Character*120 filen
c
      logical IamEL
      external IamEL

      data kall/0/
      data nusrt/7/  
      save 
c********************************************************************
c 
      If(kall.eq.0) then
c      open list file to read in existing files
       filen = surfdir(1:lenstr(surfdir) ) // 'nlist.txt'
       open(nusrt, file=filen,status='old', err=999)
c
       Do ii=1,mxn
        read(nusrt,end=555, fmt=*) irefr(ii)
        if(irefr(ii).eq.-1) goto 555
       Enddo
 555   close(nusrt)
       nfiles = ii-1
       kall = 1
      Endif
c
c     ...subsequent calls start here
c     Demo version will only have one entry
      if(nfiles.eq.1) then
        refr = 0.001*irefr(1)
        fmisc(30) = refr
      endif
c
c     Find bracketting refrs
      nrefr =  nint(1000*refr)    !value to match
      if(nrefr.lt.irefr(1)) then
        refr = 0.001*irefr(1)    
        return      
      elseif(nrefr.gt.irefr(nfiles)) then
        refr = 0.001*irefr(nfiles)        
        return      
      else
         do i = 1,nfiles
            if(nrefr .le. irefr(i)) goto 200
         enddo
      endif
c
  200 continue
c
      !save bounds
      ir1 = irefr(i)
      ir2 = irefr(min(i+1,nfiles))
c
      if(IamEL()) then            !EL, refr must be exact
        if(nrefr-irefr(i).ne.0) then
          if(abs(nrefr-irefr(i)).ge.
     1       abs(nrefr-irefr(min(i+1,nfiles))))then
             refr = 0.001*irefr(i) 
          else
             refr = 0.001*irefr(min(i+1,nfiles)) 
             ir1 = irefr(min(i+1,nfiles))
          endif
        endif
        ir2 = ir1   !no interp in EL
      endif

      return
  999 call nofile(10, 'CHECKREFR', filen)
      end subroutine

