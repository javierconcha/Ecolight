C     Last change:  LKS  15 Aug 2013   11:50 am
      real function yinterp(i1, nvals, x, xvals, yvals)
c     This function returns the y value corresponding to the given x
c     value.  Initial guess at lower bound is i1
c
      real x, xvals(nvals), yvals(nvals), f
      integer i1, i2, i, j, nvals
!*********************************************************************
!     input QC of initial index guess
      If(nvals.eq.1) then   !only one data pt
        yinterp = yvals(1)
      ElseIf(i1.lt.1) then
        i1=1
      Elseif(i1.ge. nvals) then
        i1 = nvals - 1
      Endif

!     start searching relative to last index
      If(x.gt. xvals(i1)) then
        Do j=i1, nvals
            if(xvals(j).ge.x) goto 30
        Enddo
        j = nvals
      Else
        Do j=i1, 1, -1
            if(xvals(j).lt.x) then
              goto 20
            endif
        Enddo
        j = 0
      Endif
!
 20   j = j + 1   !exit point for decreasing search
c     find factor for interpolation
 30   i = j       !exit point for increasing search
      if(i.eq.1) then
        i1=1
        i2=2
        f=0.0
      elseif(i.gt.nvals) then
        i1=nvals-1
        i2=nvals
        f=1.0
      else
        i1 = i-1
        i2 = i
        f = (x-xvals(i1))/(xvals(i2)-xvals(i1))
      endif
!
      yinterp = (1.-f) * yvals(i1) + f * yvals(i2)
!      write(80,'(2f8.4)') x, yinterp
      end function

!_______________________________________________________________________
      Subroutine getonGrid1(itype,nheadr,nh2p, ngrid,filen, 
     1                      xref,yref,xgrid,ygrid, zgrid)
!
!     routine reads the data from file filen, then puts the y-data on
!     the xgrid via linear interpolation.
!     Returns ygrid (ngrid), and yref(xref)
!     May be called to read and fit ONE data vector (itype=0, ygrid), or
!     two data vectors with the same abscissas (itype=2, ygrid and zgrid)

!     IO vars
      integer nheadr,ngrid, nh2p    !nh2p is #header lines to print
      real  xref,yref,xgrid(ngrid),ygrid(ngrid), zgrid(*)
      CHARACTER (LEN = * ) filen

!     local variables
      real, allocatable,dimension(:) :: xdat,ydat, zdat
      integer iclose
      integer getndat
      external getndat
!*********************************************************************
!     load the data
      nudata = 40
      iclose = -1
      nDat = getndat(iclose, nudata,nheadr, nh2p, filen)  
      allocate(xDat(nDat), yDat(nDat), zDat(ndat))  !dim to #data
      call loadDat(nudata,nheadr,nDat,itype,filen, xDat,yDat,zDat)

!     put the data on our grid
      call GridDat(itype,ngrid,ndat,xref,yref,xgrid,ygrid,zgrid,
     1                   xDat, yDat, zDat)

!     clean up and exit
      deallocate(xDat, yDat, zDat)
      return

  200 call nofile(10,"GetOnGrid1",filen)
      end subroutine
!_______________________________________________________________________
      Subroutine GridDat(itype,ngrid,ndat,xref,yref,xgrid,ygrid,zgrid,
     1                   xDat, yDat, zDat)
!
!     routine puts the y-data on the xgrid via linear interpolation.
!     Returns ygrid (ngrid), and yref(xref)
!     May be called to read and fit ONE data vector (itype=0, ygrid), or
!     two data vectors with the same abscissas (itype=2, ygrid and zgrid)

!     IO vars
      integer itype, ngrid, ndat 
      real xref,yref
      real,dimension(ngrid) :: xgrid, ygrid, zgrid
      real,dimension(ndat) :: xdat,ydat, zdat

!     local variables
      integer i, i1
!*********************************************************************

!     load the y grid array
      i1 = 1
      Do i=1,ngrid
        ygrid(i) = yinterp(i1, nDat, xgrid(i), xDat, yDat)
      Enddo
      if(xRef.gt.0) then
        i1 = 1
        yref = yinterp(i1, nDat, xref, xDat, yDat)
      endif
!     load the z grid IFF itype = 2 (NOTE: currently z-var does not include a zREF)
      If(itype.eq.2) then
        i1 = 1
        Do i=1,ngrid
          zgrid(i) = yinterp(i1, nDat, xgrid(i), xDat, zDat)
        Enddo
      Endif

      return
      end subroutine

!_______________________________________________________________________
      integer function getndat(iclose,nudata,nheadr, nh2p, filen)
!
!     This function opens and counts lines of data
!     if iclose<=0, file left open; if >0 file is closed

      integer nudata, nheadr, nh2p
      CHARACTER (LEN = * ) filen
      character(len = 120) header
      integer, intent(in) :: iclose
      real x
      integer i, ndat
!*********************************************************************
      call IOshorten(nudata, filen)
      if(nh2p.gt.0) write(10,102) trim(filen)
      Do i=1,nheadr
        read(nudata,'(a)') header
        if (i.le.nh2p) write(10,103) trim(header)
      Enddo
      nDat = 0
      Do 
        read(nudata,*, end=10) x
        if (x.lt.0) goto 10
        nDat = nDat+1
      Enddo
   10 getndat = ndat
      write(10,104) ndat
      if(iclose.gt.0) close(nudata)
      return
  102 format(/5x,'The data file header records from ',a,' are:')
  103 format(8x,a)
  104 format(/,5x,i6,' records read')
      end function

!_______________________________________________________________________
      integer function getnwavac(nudata,nheadr)
!
!     This function reads nwave of open ac-9 type data file
      integer, intent(in) :: nudata, nheadr
      character(len = 2) header
      integer i
!*********************************************************************
      rewind(nudata)
      Do i=1,nheadr
        read(nudata,'(a)') header
      Enddo
      read(nudata,*) getnwavac      !nwave in an ac-9 type file
      end function

!!_______________________________________________________________________
      subroutine LoadDat(nudata,nheadr,nDat,itype,filen,xDat,yDat,zDat)
!     loads opened file

!     itype = 0:  one col of data; wavel is indep variable x
!           = 1:  one col of data; depth is indep variable x
!           = 2:  two col of data; wavel is x, depth is y, dep var is z

      integer nudata, nheadr, nDat, itype
      real xDat(nDat), yDat(nDat), zDat(nDat)
      CHARACTER (LEN = * ) filen
      character(LEN=2) header
      character(Len=80) msg
!*********************************************************************
      rewind(nudata)
      Do i=1,nheadr
        read(nudata,'(a)') header
      Enddo
      Do i=1,nDat    !load all data
        If(itype.eq.2) then
          read(nudata,fmt=*, end=20) xDat(i), yDat(i), zDat(i)
c       *Check goodness of data (must be =>0)
          call checkDat2x(yDat(i), xDat(i), zDat(i))  !z, wavel, data
        Else
          read(nudata,fmt=*, end=20) xDat(i), yDat(i)
c       *Check goodness of data (must be =>0)
        call checkDat(xDat(i), yDat(i))
        Endif
c       Data must be monotonically increasing
        If(i.gt.1) then
          If(xDat(i).le.xDat(i-1)) then
!            msg="Data abscissas not monotonically increasing." //
!     2          "  Clean up your data and rerun."
            write(10,fmt='(/10x,a,/15x,f8.2,a,f8.2,a,i8)')
     1          "Data abscissas not monotonically increasing: ",
     2          xDat(i), " < ", xDat(i-1),' for datum #:',i
            msg = "  Clean up your data and rerun."
            call HERR('LoadDat', msg)
          Endif
        Endif
      Enddo          
  20  close(nudata)

!     save info to report data Quality
      If(itype.eq.1) then         ! x var is depth
        zmin = xDat(1)
        zmax = xDat(ndat)
        wmin=0.0
        wmax = 0.0
      Else                        !x var is wavelen
        zmin=0.0
        zmax = 0.0
        wmin = xDat(1)
        wmax = xDat(ndat)
      Endif
      CALL RangeCheck(1, zmin,zmax,wmin,wmax,filen)

      RETURN
      END SUBROUTINE

!_______________________________________________________________________
      subroutine checkDat(x, y)
      real x, y
c     *Check goodness of data (must be =>0)
      IF(y.lt.0.) THEN
        write(10,*) 'CheckDat: Error in data read from file'
        write(10,'(5x,a,f6.2)') 
     1        'negative, non-physical data at abscissa value: ', x
             write(10,*) 'Datum set to zero.'
             y=0.0
      ENDIF
      RETURN
      END SUBROUTINE

!_______________________________________________________________________
      subroutine checkDat2x(x, y, z)
      real x, y, z
c     *Check goodness of data (must be =>0)
      IF(z.lt.0.) THEN
        write(10,*) 'CheckDat: Error in data read from file'
        write(10,'(5x,a,f6.2)') 
     1        'negative, non-physical data at depth: ', x,
     2        '                      and wavelength: ', y
             write(10,*) 'Datum set to zero.'
             z=0.0
      ENDIF
      RETURN
      END SUBROUTINE

!_______________________________________________________________________
      subroutine checkiWav(jwave, wavel)
      real wavel, dwav, dwav1
      integer jwave, i
      INCLUDE "DIMENS_XL.INC"
      COMMON /Cmisc/ imisc(30),fmisc(30)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
!*********************************************************************
!     Check that wavel is the same wavelen as iwave
      IF(abs(wavel-wave(jwave)).gt.1e-6) THEN
        jwave = 1
        dwav = 10000
        Do i=1,imisc(7)        !find nearest iwave   
          dwav1 = abs(wavel-wave(i))
          if(dwav1.lt.dwav) then
            dwav = dwav1
            jwave = i
          endif
        Enddo      
      ENDIF
      RETURN
      END SUBROUTINE
