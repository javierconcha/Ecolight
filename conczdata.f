C     Last change:  LKS  25 Oct 2013    4:06 pm
      subroutine conczdata(icomp, depth, conc)
c
c     user routine on file conczdata.f
c
c     called by "ab" routine
c
c     This routine is called when the user has a data file of measured
c     (or simulated) concentration values at discrete depths.   
c
c     This routine is designed for reading in concentration profile data,
c     which are assumed to be on the Hydrolight "Concentration data 
c     standard format" such as for chlorophyll profiles. 
c
c     When this routine is called for the first time, the data file is
c     read.
c
c     On subsequent calls, linear interpolation is used to define the conc
c     at any depth that lies between the first and last depths in the data
c     files.  Depths above (below) the first (last) data depth are given
c     the value for the first (last) depth measured. 
c
c     INPUT:
c        depth: the depth in meters where the concentration is requested
c
c     OUTPUT:
c        conc:  the concentration in mg / m^3 at the input depth
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
!     User adjustable max # depths in concentration data file
      parameter (MaxZdat = 500)

      real,dimension(MaxZdat,mxcomp) :: zdatS, concdatS
      integer,dimension(mxcomp) :: nDat

c     local variables  
      integer iop, nudata
      
c     dummy variables (IO)
      integer icomp
      real depth, conc
c
!     local variables
      integer kall, itype, nheadr, nh2p, iclose, jcomp
      integer kallS(mxcomp), i1(mxcomp)
      real, allocatable,dimension(:) :: zdat, concdat
      real dum
      integer getndat
      external getndat
      data kall/0/, itype/1/,nheadr/10/,nh2p/2/

c     Note on dynamic storage:  MUST save values between subroutine calls:
      save kall, kallS, i1, nDat, zDatS, chlDatS
c
c-----  Begin Initialization on first call  ---------------------------
      if(kall.eq.0) then
        ncomp = imisc(6)
        do i=1,ncomp
          kallS(i)=0
        enddo
        kall = 1
      endif

      if(kallS(icomp).eq.0) then
c        initialization on first call
         jcomp = icomp + 6  !index into datafile since component files start at 7
         write(10,100) icomp, datafiles(jcomp)
c
         iop = imisc(5)
         if(iop.ne.0) then
c           stop the run if it is using optical depth; expects
c           geometric depth as input
            write(10,104) iop
            call HERR("compzdata","depth type is not compatible")  !stop run
         endif
c
c       open and read the CDOM data file
        nudata = 40
        iclose = -1
        nDat(icomp)=getndat(iclose,nudata,nheadr,nh2p,datafiles(jcomp)) 
        allocate(zDat(nDat(icomp)),concDat(ndat(icomp)))  !dim to #data
        call loadDat(nudata,nheadr,nDat(icomp),itype, datafiles(jcomp),
     1               zDat,concDat,dum)
        If(nDat(icomp).gt.maxZDat) then
          write(10,105) nDat(icomp), maxZDat
          call HERR("compzdata",
     1         "Cleanup CONC data or increase MaxZDat")  !stop run
        Endif
c
c       store data into permanent memory and deallocate
        zDatS(1:nDat(icomp),icomp) = zDat
        concDatS(1:nDat(icomp),icomp) = concDat

        kallS(icomp) = 1
        i1(icomp) = 1
        conc  = yinterp(i1(icomp), nDat, depth, zDat, concDat)
        deallocate(zDat,concDat)
        return
      endif
c-----  End of initialization on first call----------------------------
c
c-----  Subsequent calls start here:
      allocate(zDat(nDat(icomp)),concDat(nDat(icomp))) 
      zDat = zDatS(1:nDat(icomp),icomp) 
      concDat = concDatS(1:nDat(icomp),icomp)

c     conc from interp routine:
      conc  = yinterp(i1(icomp), nDat(icomp), depth, zDat, concDat)
      deallocate(zDat,concDat)
      return
c
c     Formats:
  100 format(/5x,
     1'The concentration for component',i3,/5x,
     1'is obtained from measured values read from the file:'/
     38x,a,/8x,
     4'(Linear interpolation is used to define conc(z) ',
     5'at all depths)')
  104 format(//'Error in sub conczdata:  called with optical depth: iop
     1= ',i2) 
  105 format(//'Error in sub conczdata:  'i6,' data read but only ', 
     1i6,'allocated in MaxZDat') 
      end
