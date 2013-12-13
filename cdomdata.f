C     Last change:  LKS   7 Aug 2008    7:05 pm
      subroutine cdomdata(icomp, depth, cdomconc)
c
c     user routine on file cdomdata.f		
c
c     called by "acdomsub" routine in incfiles.for
c
c     This routine is called when the user has a data file of measured
c     (or simulated) CDOM absorption values at discrete depths and at a
c     known reference wavelength.  
c
c     This routine is used when determining absorption via bio-optical 
c	models, and when computing cdom fluorescence.
c
c     This routine is designed for reading in CDOM profile data,
c     which are assumed to be on the Hydrolight "CDOM-data 
c     standard format."
c
c     Linear Interpolation is used to define CDOM at any
c     depth that lies between the first and last depths in the data
c     file.  Depths above (below) the first (last) data depth are given
c     the CDOM value for the first (last) depth measured. 
c
c     INPUT:
c        depth: the depth in meters where the CDOM absorption is requested
c
c     OUTPUT:
c        cdomconc:  the CDOM absorption in 1/m at the input depth and ref waven
c
      INCLUDE "DIMENS_XL.INC"

      common /Cmisc/ imisc(30),fmisc(30)
c
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
c
c     local variables
      integer iop, nudata
c
c     dummy variables (IO)
      integer icomp
      real depth, cdomconc
c
!     local variables
      integer kall, itype, nheadr, nh2p, i1, iclose
      real, allocatable,dimension(:) :: zdat, cdomdat
      real dum
      integer getndat
      external getndat
      data kall/0/, itype/1/,nheadr/10/,nh2p/2/
c     Note on dynamic storage:  MUST save values between subroutine calls:
      save kall, i1, nDat, zDat, cdomDat
c
c-----  Begin Initialization on first call  ---------------------------
      if(kall.eq.0) then
c     initialization on first call
        write(10,100) datafiles(5)	!CDOM is third concentration in model
c
        iop = imisc(5)
        if(iop.ne.0) then
c          stop the run if it is using optical depth; cdomdata expects
c          geometric depth as input
           write(10,104) iop
           call HERR("cdomdata","depth type is not compatible")  !stop run
        endif
c
c       open and read the CDOM data file
        nudata = 40
        iclose = -1
        nDat = getndat(iclose, nudata,nheadr, nh2p, datafiles(5)) 
        allocate(zDat(nDat),cdomDat(ndat))  !dim to #data
        call loadDat(nudata,nheadr,nDat,itype,datafiles(5), 
     1               zDat,cdomDat,dum)
        kall = 1
        i1 = 1
      endif
c-----  End of initialization on first call----------------------------
c
c-----  Subsequent calls start here:
c     cdomconc from interp routine:
      cdomconc = yinterp(i1, nDat, depth, zDat, cdomDat)
      return
c
c     Formats:
  100 format(//5x,
     1'The CDOM absorption at a reference wavelength is obtained ',
     2' from measured '/5x,
     2'values read from the file:'/
     38x,a,/)
  104 format(//'Error in sub cdomdata:  called with optical depth: ',
     1'iop = ',i2) 
      end
