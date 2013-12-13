C     Last change:  LKS  25 Oct 2013    3:51 pm
!  This file contains a listing of all routines common to
!  the Hydrolight and Ecolight code
!
!  Lydia Sundman, July 2007
!
!     Common incfiles routines
!     **note:  when STND version built, incfiles_default is copied
!              to incfiles_user before compiling by makestndexe.bat
      INCLUDE "..\common\incfiles_stnd.for"
!
!     Common IOP models
      INCLUDE "..\common\abacbb.f"
      INCLUDE "..\common\abcase1.f"
!      INCLUDE "..\common\abcase1H.f"
      INCLUDE "..\common\abcase2.f"
      INCLUDE "..\common\abconst.f"
!      INCLUDE "..\common\ablayer.f"
      INCLUDE "..\common\abnewcase1.f"
      INCLUDE "..\common\abpure.f"
!
!     Common component routines
      INCLUDE "..\common\astar.f"
      INCLUDE "..\common\bstar.f"
      INCLUDE "..\common\Bfcalc.f"
      INCLUDE "..\common\chlzdata.f"
      INCLUDE "..\common\cdomdata.f"
      INCLUDE "..\common\conczdata.f"   !new in HE5.2
      INCLUDE "..\common\minzdata.f"
      INCLUDE "..\common\pureh2o.f"
!
!     Common inelastic routines
      INCLUDE "..\common\s0bdata.f"
      INCLUDE "..\common\wrfuncs.f"
!
!     Common atmospheric models
      INCLUDE "..\common\cosirrad.f"
      INCLUDE "..\common\cosrad.f"
      INCLUDE "..\common\get_03_climat.f"
      INCLUDE "..\common\hcnrad.f"
      INCLUDE "..\common\irradat.f"
      INCLUDE "..\common\sunang.f"
!
!     Common BC models
      INCLUDE "..\common\getDynZ.f"
      INCLUDE "..\common\rbottom.f"
!
!     Common processing routines
      INCLUDE "..\common\radanal.f"

!     Common IO and admin routines
      INCLUDE "..\common\admin.f"
!      INCLUDE "..\common\CIExyY.f"
      INCLUDE "..\common\Color.f"
      INCLUDE "..\common\dataQC.f"
      INCLUDE "..\common\excel.f"
      INCLUDE "..\common\initial.f"
      INCLUDE "..\common\secchi.f"
      INCLUDE "..\common\setdflts.f"
      INCLUDE "..\common\slctsurf.f"
      INCLUDE "..\common\Surfn.f"
      INCLUDE "..\common\par.f"
      INCLUDE "..\common\pntinit.f"
      INCLUDE "..\common\zetatoz.f"
      INCLUDE "..\common\ztozeta.f"
!
!     Common mathematical routines
      INCLUDE "..\common\BLAS.f"
      INCLUDE "..\common\eigenvv.f"
      INCLUDE "..\common\interp.f"     !in house linear interpolation routines
      INCLUDE "..\common\matxmat.f"
      INCLUDE "..\common\matinv.f"
      INCLUDE "..\common\ode.f"
      INCLUDE "..\common\p2aray.f"
      INCLUDE "..\common\p3aray.f"
      INCLUDE "..\common\SLAcom.f"
      INCLUDE "..\common\sort.f"
      INCLUDE "..\common\VECxMAT.f"

