!     Last change:  LKS  20 Feb 2008    9:24 am
!   include "..\windows.f90"  -- this is just a reminder that windows.f90 must
!                             be compiled BEFORE this file since it is USEd below

    module windows_module
!     Note:  this module must occur BEFORE any subroutine that uses it.
!          Similarly, windows_h (windows.f90) must be compiled before this module.
      use WIN32MOD
      implicit none

!      integer :: GetShortPathNameA
      DLL_IMPORT GetShortPathNameA
    end module windows_module


   subroutine IOshorten(iounit, longfn)
!     this routine calls an API function to find the short (8.3) filename SHORTFN 
!     correspponding to the user-supplied 120-character LONGFN.

!     Note that LONGFN can be a full pathname, a relative pathname, or just a
!     filename (if the file is in the same directory).

!     If the file exists, lreturn will equal the length of shortfn and 
!     this routine will open the file with the unit number IOUNIT.  
!     If the file does NOT exist (lreturn = 0), the routine will report the 
!     error and stop

      use windows_module   
!     module windows_module must appear in this file before IOshorten which USEs it
!      implicit none
      character (len=120) :: longfn, longfn0      !long  filename
      character (len=120) :: shortfn     !short filename

      integer :: lenfn, lreturn, iounit

      lenfn   = 80
      lreturn = 0
      shortfn = char(0)
      longfn0  = trim(longfn(1:119)) // char(0)

      lreturn = GetShortPathNameA(carg(offset(longfn0)),carg(offset(shortfn)),carg(lenfn))

      IF(lreturn.gt.0) then
!        make sure file is opened as OLD (status) and READONLY (action) to avoid
!        any conflict opening the file if file is marked 'readonly' by Windows
        open(iounit, file=shortfn, status='old', ACTION='read')
      Else
        call nofile(10, 'IOSHORTEN', longfn)   !err opening file

      Endif

    end subroutine IOshorten


    subroutine stoprun90
!      this routine replaces the traditional f77 stop command with the <snazier>
!      f90 version (for the Windows users).  The f90 stop brings up a message
!      box.
       implicit none

       WRITE(6,*) 'Run completed successfully! Hit return to exit.'
       call CLOSEALL

   end subroutine stoprun90
