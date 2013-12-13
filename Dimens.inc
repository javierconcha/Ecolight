c     This file gives the parameter values used to specify array
c     dimensions in Hydrolight.  All other array dimensions found
c     in the code can be computed from these.  These parameter
c     statements are inserted into the code via the INCLUDE DIMENS.INC
c     statements in the various routines.  In order to change array
c     dimensions, change the appropriate values here and remake the
c     executable (i.e., recompile and relink).
c
c         mxmu is the maximum allowed value of nmu = M/2 (= m)
c            (e.g., nmu = 10 for 10 theta bands in each hemisphere)
c         mxphi is the maximum allowed value of nphi = N (= 2n)
c            (e.g., nphi = 24 for 15-degree delta phi bands)
c         mxz is the maximum number of depths at which output can be 
c            requested (e.g. 100, for 50 pairs of closely spaced depths)
c         mxwave is the maximum number of wavelength bands
c            (e.g. 90 for 5 nm bands from 350 to 800 nm)
c         mxnsky is the maximum number of parameters used by the 
c            skyspecs to store input for the sky models
c         mxcomp is the maximum number of components used to build up
c            water IOP's in the  models
c         mxnzvals is the maximum number of zeta-to-z pairs used to
c            convert between optical and geometrical depths;
c            default is mxnzvals = 7 * mxz
c
      integer mxmu,mxphi,mxz,mxwave,mxnsky,mxcomp,mxnzvals

      PARAMETER (mxz=200, mxwave=150, mxmu=10, mxphi=24, mxcomp=10)
      PARAMETER (mxnsky=10, mxnzvals=1400)
