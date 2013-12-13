C     Last change:  LKS   8 Dec 2011    3:17 pm
!     This is the dummy version used for compiling the STANDARD version of HL
!------------------------------------------

      subroutine user_abmodel(depth, dcompconc)
c
      INCLUDE "DIMENS_XL.INC"
c
      real depth, dcompconc(mxcomp)

      call HERR("user_abmodel","no concentrations specified")  !stop run
      end
c     -----------------------------------------

      subroutine user_abscat(depth,wavelen,ncomp,acomp,bcomp,
     1                       atotal,btotal)
c
c     This file selects the 'ab' routine to be used in the run,
c     called by routines inishamp, rhotau, zetatoz, and ztozeta.
c
      INCLUDE "DIMENS_XL.INC"
c
      integer ncomp
      real depth,wavelen,atotal,btotal
      real acomp(mxcomp),bcomp(mxcomp)

      call HERR("user_abscat","no IOP model specified")  !stop run
      end
c     -----------------------------------------

      subroutine user_acdomsub(depth, wavelen, abscdom)
c
c     This file specifies the 'aCDOM' to be used in the run,
c     by calls from the 'ab' routines, and the CDOM fluorescence routine 'shatcdom'
      real depth,wavelen, abscdom

c     CDOM is not included in this run
      abscdom = 0
      call HERR("user_acdomsub","subroutine not specified")  !stop run
      return
      end
c     -----------------------------------------

      subroutine user_achlz(depth, wavelen, achl)
c
c     This file selects how the chlorophyll profile is specified for the run,
c     for use by the 'ab' routines and the chlorophyll fluorescence
c     routine 'shatchl'
      real depth,wavelen, abscdom

c     CDOM is not included in this run
      achl = 0
      call HERR("user_achlz","subroutine not specified")  !stop run
c
      return
      end
c     -----------------------------------------

      subroutine user_skyirrad(suntheta,sunphi,eddif,eddir)
c
c     This file selects the sky irradiance model to be used in the run,
c     and is called by routine qasky.
c
      real suntheta,sunphi,eddif,eddir
      call HERR("user_skyirrad","subroutine not specified")  !stop run
c
      return
      end
c     -----------------------------------------

      subroutine user_skyrad(suntheta,sunphi,skytheta,skyphi,skyrad0)
c
c     This file selects the sky radiance model to be used in the run,
c     and is called by routine qasky.
c
      real suntheta,sunphi,skytheta,skyphi,skyrad0
c
      call HERR("user_skyrad","subroutine not specified")  !stop run
      return
      end
c     -----------------------------------------

      Real Function s0bioSub(depth,wavelen)
c
c     This file calls the bioluminescent function to be used in the run,
c     (if needed at all) by shat.f
c
      Common /Csource0/ ibiolum,ichlfl,icdomfl,iraman, ramanEXP
      real depth,wavelen, ramanEXP
      INTEGER ibiolum,ichlfl,icdomfl,iraman
c
      IF(ibiolum.eq.0) then
c      Bioluminescence is not included in this run
       s0bioSub = 0.0
      ELSEIF(ibiolum.eq.1) then
c      read a data file for bioluminescence
       s0bioSub = s0bdata(depth,wavelen)
      ELSE
c      not supported without recompiling a USER version
      call HERR("s0bioSub","Bioluminescence FUNC option not supported")  !stop run
      ENDIF
c
      return
      end
c     -----------------------------------------

