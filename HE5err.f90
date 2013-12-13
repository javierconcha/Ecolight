!     Last change:  LKS   4 Apr 2008    5:32 pm
      subroutine nofile(nuinit, callname, filen)
!    This routine is called IFF a required data file cannot be opened

      implicit none
      integer nuinit
      character (LEN = * ) callname
      CHARACTER (LEN = * ) filen
      CHARACTER (LEN = 120 ) mssg

      mssg = 'required data file ' // TRIM(filen) // ' not found'
!
      call HERR(TRIM(callname),mssg)
      end subroutine

!_______________________________________________________________________
      subroutine HERR(subr, mssg)
      character (LEN = * ) subr
      CHARACTER (LEN = * ) mssg

      !print mssg
      WRITE(10,100) 'ERROR: ',TRIM(mssg),'STOP in routine ',TRIM(subr)
      WRITE(6,100) 'ERROR: ',TRIM(mssg),'STOP in routine ',TRIM(subr)
 100  FORMAT(/2x,2a,/2x,2a)

      !close all file connections
      call CLOSEALL

      !end run
      STOP
      end subroutine

!_______________________________________________________________________
      subroutine CLOSEALL
      INTEGER i
      !close all files
      DO i=1,100
        CLOSE(i)
      END DO

      END subroutine
