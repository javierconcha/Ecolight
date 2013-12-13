C     Last change:  LKS  24 Feb 2009    4:27 pm
c
      Subroutine cpright(iounit)
!      program cpright

      integer iounit 
      character*5 version

      character*10 getHE5
      external getHE5

      include "version.inc"
!!!
!      iounit = 10
!      open(iounit, file='output.txt')
!!!
      Write(iounit,9998) getHE5(), version
      call license(iounit)
      Write(iounit,9999)

 9998 Format(/
     *'     ========================================================'/,
     *'     |                                                      |'/,
     *'     |        ',a,' Version ',a5,' (MAIN CODE)          |'/,
     *'     |       Copyright (c) 2013 by Curtis D. Mobley         |'/,
     *'     |                                                      |')

 9999 format(
     *'     |                                                      |'/,
     *'     | THIS PROGRAM IS EXPERIMENTAL AND IS LICENSED "AS IS" |'/,
     *'     | WITHOUT REPRESENTATION OF WARRANTY OF ANY KIND,      |'/,
     *'     | EITHER EXPRESS OR IMPLIED.  THE ENTIRE RISK AS TO    |'/,
     *'     | THE QUALITY AND PERFORMANCE OF HYDROLIGHT AND        |'/,
     *'     | ECOLIGHT IS WITH THE USER.  HYDROLIGHT AND ECOLIGHT  |'/,
     *'     | ARE NOT FAULT TOLERANT.                              |'/,
     *'     |                                                      |'/,
     *'     ========================================================'/)


      return
!      close(iounit)
      end
c
C--------------------------------------------------------------------------
      Subroutine date_stamp(iounit, iflag, jwave, nwave)
c
c     *****Variables for Check date and time *****
      character*10 tdate, zone,time, txttime, txtdate
      integer iv(8), iounit, iflag
      integer iv0(8), ijul, ijul0, isec, kall
      integer iv2(8), ijul2, jwave, nwave
      real nettime
      character*10 txtflag(0:1)
c
      data txtflag/'started','completed'/
      data kall/0/
      save
c
c     Calculate total runtime using julian days
c     NOTE:  will only handle crossing New Years ONCE
c
c     this routine times the wavelength and total run times
c         iflag = 0:  start of run
c         iflag = 1:  completion of run
c         iflag = 2:  initialization timing
c         iflag = 3:  incremental waveband timing (want time since last called)
c
c
c     save previous run if "incremental" time desired
      if(kall.gt.0 .and. (iflag.eq.2 .or. iflag.eq.3)) then
          Do j=1,8
              iv2(j) = iv(j)
          Enddo
          ijul2 = julianday(iv2)  !calculate the julian day of the year
      endif

c     ***Call internal system routine to get date and time
c     Optional time stamp routine for Lahey Fortran only
      call DATE_AND_TIME(tdate,time,zone,iv)
c     date: ccyymmdd  time: hhmmss.ss
c     iv:  1-ccyy, 2-mm, 3-dd, 4-UTC, 5-hh, 6-min, 7-sec, 8-msec
      txtdate = tdate(5:6)//char(47)//tdate(7:8)//char(47)//tdate(1:4)
      txttime = time(1:2)// ':' // time(3:4) // ':' //time(5:9)
c
c
c     ***Print text if start or completion of run
      if( iflag.eq.0 .or. iflag.eq.1) then
       Write(6,600) 'Run ',txtflag(iflag)(1:lenstr(txtflag(iflag))),
     1             ' on ',txtdate,' at ',txttime
       Write(iounit,600) 'Run ',txtflag(iflag)
     1    (1:lenstr(txtflag(iflag))),' on ',txtdate,' at ', txttime
      endif
c

c     ***Store start time data if call is "start time"
c     if kall=0, call to "completed time" occurred before call for start time
      if(iflag.eq.0 .or. kall.eq.0) then
        Do j=1,8
          iv0(j) = iv(j)
        Enddo
        ijul0 = julianday(iv0)    !calculate the julian day of the year
c
c     ***Calculate total runtime if call is for "completed time"
      Else

        ijul = julianday(iv)      !calculate the julian day of the year
c
        if((iv(1)-iv0(1)).ne.0) then  !we have crossed over into a new year
          idayyr = 365
          if(abs(modulo(1.0*abs(iv(1)-2000),4.0)).lt.1.e-4) idayyr = 366
          ijul = ijul + idayyr
        endif
c
c       ** Print out incremental time
        If(iflag.eq.2) then
c         Calculate total number of elapsed seconds during initialization
          isec = (((ijul-ijul2)*24+(iv(5)-iv2(5)))*60+
     1               (iv(6)-iv2(6)))*60 + (iv(7)-iv2(7))
            nettime = 0.001*(iv(8)-iv2(8)) + 1.0*isec
c
            write(6,400) nettime
            write(iounit,400) nettime
        ElseIf(iflag.eq.3) then
c         Calculate total number of elapsed seconds in waveband
          isec = (((ijul-ijul2)*24+(iv(5)-iv2(5)))*60+
     1               (iv(6)-iv2(6)))*60 + (iv(7)-iv2(7))
            nettime = 0.001*(iv(8)-iv2(8)) + 1.0*isec
c
            write(6,301) jwave, nwave, nettime
            write(iounit,300) jwave, nwave, nettime
        Else
c       ** Print out completion time
c         Calculate total number of elapsed seconds in run
          isec = (((ijul-ijul0)*24+(iv(5)-iv0(5)))*60 +
     1                (iv(6)-iv0(6)))*60 +(iv(7)-iv0(7))
            nettime = 0.001*(iv(8)-iv0(8)) + 1.0*isec
c
            write(6,200)  nettime
            write(iounit,200)  nettime
        Endif
      Endif
c  
      kall = 1
      return


 200  format(/5x,'Total (wall clock) run time = ',f12.1,' sec')
 300  format(//5x,'Waveband ',i3,' of ',i3,' completed in ',
     1       f6.1,' sec.'/)
 301  format(5x,'Waveband ',i3,' of ',i3,' completed in ',
     1       f6.1,' sec.'/)
 400  format(//5x,'Initialization of Hydrolight Completed in ',
     1       f4.1,' sec.')
 600  format(/5x, 6a,2x/)
      end
c
C--------------------------------------------------------------------------
      Integer function julianday(iv)
c     Calculate Julian Day of the year
      integer iv(8), mday(12), iday
c     iv:  1-ccyy, 2-mm, 3-dd, 4-UTC, 5-hh, 6-min, 7-sec, 8-msec
      data mday/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
c
c     February has 29 days if it is a leap year
      leap_day = 0 	! most years are not leap years
      If(modulo(1.0*iv(1),400.0).lt.1.e-4) Then
          leap_day = 1  ! years divisible by 400 are leap years
      ElseIf(modulo(1.0*iv(1),100.0).lt.1.e-4)Then
          leap_day = 0 ! other centuries are not leap years
      ElseIf(modulo(1.0*iv(1),4.0).lt.1.e-4)Then
          leap_day = 1 ! otherwise every 4th year is a leap year
      End If
      IF(leap_day .eq. 1) mday(2)=29
c
      iday = iv(3)            !iday = day of current month
      If( iv(2) .gt.1) then   !if it is not January...
        Do j=1, iv(2)-1       !...sum days over months that have passed
          iday = iday + mday(j)
        Enddo
      Endif
      julianday = iday
      return
      end

C--------------------------------------------------------------------------
c
      INTEGER FUNCTION LENSTR ( String )

C     This function returns the number of characters in a character string.
C     Counting begins at the end of the string and decrements so that
C     imbedded blanks are allowed.

      Implicit None

      Integer I
      Character String * (*)

      Do I = Len(String),1,-1
         IF ( String(I:I) .ne. ' ' ) THEN
            LENSTR = I
            Return
         END IF
      end do
c
c     The string is all blanks:
      LENSTR = 0

      Return
      End
C--------------------------------------------------------------------------

      character*10 function getHE5()
c     This routine returns either "hydrolight" or "ecolight", 
c     which is currently being run 
c
      logical IamEL
      external IamEL              !in routine IamEL_XL.f

      if(IamEL()) then            !EL
        getHE5 = "  ECOLIGHT"
      else                        !HL
        getHE5 = "HYDROLIGHT"
      endif
 
      return
      end function
C--------------------------------------------------------------------------

      subroutine license(iunit)
!    This subroutine gets and prints the license info for this copy of HE5
! Dummy arguments
      INTEGER :: Iunit
      INTENT (IN) Iunit
! Local variables
      CHARACTER(48) :: sn , user

c      DLL_IMPORT :: info

c      call info(user, sn) !call to DLL

c      WRITE (Iunit,100) user, sn
 
c100   FORMAT (
c     &        "     |     This copy of HE5.2 is ",                        &
c     &        "licensed solely to           |",/"     |   ",a48,        &
c     &        "   |",/                                                  &
c     &   '     |                                                      |'&
c     &      /"     |     S/N ",a44," |"   )
      END SUBROUTINE LICENSE
c     **************************************************

      Real Function yBounds(y, ymin, ymax)
c     This function returns the value y, unless y is outside the bounds
c     (ymin, ymax) in which case it returns the nearest bound.
      Real y, ymin, ymax

        If(y.lt.ymin) then
           yBounds = ymin
        Elseif(y.gt. ymax) then
           yBounds = ymax
        Else
           yBounds = y
        Endif

        Return
        End
