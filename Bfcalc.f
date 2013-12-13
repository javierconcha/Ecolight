C     Last change:  LKS  21 Apr 2008   11:06 am
      subroutine Bfcalc(icomp, wavenm, Bf)
!     This subroutine is called IFF the power law option is selected
!     to specify bb/b
      integer icomp
      real, intent(in) :: wavenm
      real, intent(out) :: Bf
c
      INCLUDE "DIMENS_XL.INC"
      common /Cbbopt/ ibbopt(mxcomp), bbfrac(mxcomp), 
     1                BfrefPL(mxcomp), Bf0PL(mxcomp), BfmPL(mxcomp)
      integer ibbopt
      real bbfrac, BfrefPL, Bf0PL, BfmPL
!---------------------------------------------------------------------

c-------- EVALUATE Power Law model OPTION ----------
      Bf = Bf0PL(icomp)*(BfRefPL(icomp)/wavenm)**BfmPL(icomp)
      call selpfbb(Bf,icomp)

!        write(10,'(a,f6.3)') '  Bf for this wavelength: ',Bf

      return
      end subroutine
