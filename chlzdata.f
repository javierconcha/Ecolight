C     Last change:  LKS   7 Aug 2008    7:06 pm
      subroutine chlzdata(icomp, depth, chlconc)
c
c     user routine on file chlzdata.f
c
c     called by "ab" routine
c
c     This routine is called when the user has a data file of measured
c     (or simulated) Chlorophyll values at discrete depths.  [Routine
c     chlzfunc is called when the user has a continuous function that
c     defines Chl(z).]  
c
c     This routine is used when determining absorption and scattering
c     coefficients via bio-optical models, and when computing
c     chlorophyll fluorescence.
c
c     This routine is designed for reading in Chlorophyll profile data,
c     which are assumed to be on the Hydrolight "Chlorophyll-data 
c     standard format." 
c
c     When this routine is called for the first time, the data file is
c     read.
c
c     On subsequent calls, linear interpolation is used to define Chl
c     at any depth that lies between the first and last depths in the data
c     files.  Depths above (below) the first (last) data depth are given
c     the Chl value for the first (last) depth measured. 
c
c     INPUT:
c        depth: the depth in meters where the Chl concentration is requested
c
c     OUTPUT:
c        chlconc:  the chl concentration in mg chl / m^3 at the input depth
c
      INCLUDE "DIMENS_XL.INC"

      common /Cmisc/ imisc(30),fmisc(30)
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
c
c     local variables
      integer iop, nudata
c     dummy variables (IO)
      integer icomp
      real depth, chlconc
c
!     local variables
      integer kall, itype, nheadr, nh2p, i1, iclose
      real, allocatable,dimension(:) :: zdat, chldat
      real dum
      integer getndat
      external getndat
      data kall/0/, itype/1/,nheadr/10/,nh2p/2/
c     Note on dynamic storage:  MUST save values between subroutine calls:
      save kall, i1, nDat, zDat, chlDat
c
c-----  Begin Initialization on first call  ---------------------------
      if(kall.eq.0) then
c     initialization on first call
         write(10,100) datafiles(4)
c
         iop = imisc(5)
         if(iop.ne.0) then
c           stop the run if it is using optical depth; chlzdata expects
c           geometric depth as input
            write(10,104) iop
            call HERR("chlzdata","depth type is not compatible")  !stop run
         endif
c
c       open and read the CDOM data file
        nudata = 40
        iclose = -1
        nDat = getndat(iclose, nudata,nheadr, nh2p, datafiles(4)) 
        allocate(zDat(nDat),chlDat(ndat))  !dim to #data
        call loadDat(nudata,nheadr,nDat,itype, datafiles(4),
     1               zDat,chlDat,dum)
        kall = 1
        i1 = 1
      endif
c-----  End of initialization on first call----------------------------
c
c-----  Subsequent calls start here:
c     chlconc from interp routine:
      chlconc =  yinterp(i1, nDat, depth, zDat, chlDat)
      return
c
c     Formats:
  100 format(/5x,
     1'The Chlorophyll concentration is obtained from measured',
     2' values read from the file:'/
     38x,a,/8x,
     4'(Linear interpolation is used to define Chl(z) ',
     5'at all depths)')
  104 format(//'Error in sub chlzdata:  called with optical depth: iop
     1= ',i2) 
      end
