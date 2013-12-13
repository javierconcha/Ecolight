C     Last change:  LKS  21 Jan 2008    8:36 am
!-----------------------------------------------------------------
      real function get_o3_climat(jday, rlat, rlon)
!     This subroutine reads in the HL ozone climatology and returns o3

!     input variables:
      integer jday
      real rlat, rlon

!    temp variables
      integer iunit, im 
      integer cdiy(12)
      character*120 header, filenm
      real o3row(12)
      real rlatb(2), rlonb(2)

      data iunit/34/
      data cdiy/31,59,90,120,151,181,212,243,273,304,334,366/
!     ******************************************************
      filenm = '..\data\TOMS_O3.txt'
      open(unit=iunit, file=filenm, status='old',err=40)
c
c    read and print header lines
      write(10,*)
      Do i=1,3
        read(iunit, '(a120)', end=40) header
        write(10,'(5x,a)') header
      Enddo
      read(iunit, '(a30)', end=40) header
c
c     search thru table of Cummulative Days in Year for month containing jday
      imonth = jday/31 + 1  !initial low guess
      Do i=imonth,12
        if(jday.lt.cdiy(i)) goto 5
      Enddo
   5  im = i
!      print *, im, 'jday: ',jday, cdiy(i)

c     Read in data
  10  read(iunit,*) rlatb(1),rlatb(2),rlonb(1),rlonb(2), 
     1              (o3row(i),i=1,12)

      if(rlat.ge.rlatb(1).and.(rlat.lt.rlatb(2).or.rlatb(2).ge.90)) then
        do j=1, 36         !find our long bin
          if(rlon.ge.rlonb(1).and.(rlon.le.rlonb(2).or.rlonb(2).ge.180))
     1    then
            get_o3_climat = o3row(im)
!            write(10,*) 'Ozone from file: ',get_o3_climat
            close(iunit) 
            return     !we're done!
          else
            goto 10   !read next line; we're close!
          endif  
        enddo
      else
        do i=1,35
          read(iunit,*) rlatb(1)   !skip to next lat band
        enddo
        goto 10       !read next line of data
      endif

  40  write(10,*) 'ERROR in get_o3_climat:  ',
     1            trim(filenm),' not found; use 300'
      get_o3_climat = 300.0
      end
