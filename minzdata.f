C     Last change:  LKS   7 Aug 2008    7:04 pm
      subroutine minzdata(imincomp, depth, fminconc)
c
c     user routine on file minzdata.f
c
c     called by "ab" routines
c
c     This routine is called when the user has a data file of measured
c     (or simulated) mineral particle concentrations at discrete depths.
c     [Routine minzfunc is called when the user has a continuous function that
c     defines conc(z).]  
c
c     This routine is used when determining absorption and scattering
c     coefficients via geo-optical models
c
c     This routine is designed for reading in mineral particle profile data,
c     which are assumed to be on the Hydrolight "concentration-data 
c     standard format."  (Users can redefine that format by changing
c     the statments in this routine that read the data file.)
c
c     This routine opens and reads the file "minzdatafile", where 
c     "minzdatafile" is the path and filename as specified when running 
c     the front end program.
c
c     When this routine is called for the first time, the data file is
c     read.
c
c     On subsequent calls, linear interpolation is used to define conc 
c     at any depthbetween the first and last depths in the data
c     files.  Depths above (below) the first (last) data depth are given
c     the value for the first (last) depth measured. 
c
c     INPUT:
c        depth: the depth in meters where the concentration is requested
c
c     OUTPUT:
c        fminconc:  the mineral concentration in gm / m^3 at the input depth
c
      INCLUDE "DIMENS_XL.INC"
c
      integer imincomp
      real depth, fminconc
!
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      common /Cmisc/ imisc(30),fmisc(30)
c
!     local variables
      integer kall, itype, nheadr, nh2p, i1, iclose, iop
      real, allocatable,dimension(:) :: zdat, minDat
      real dum
      integer getndat
      external getndat
      data kall/0/, itype/1/,nheadr/10/,nh2p/2/
c     Note on dynamic storage:  MUST save values between subroutine calls:
      save kall, i1, nDat, zDat, minDat
c
c-----  Begin Initialization on first call  ---------------------------
      if(kall.eq.0) then
c     initialization on first call
         write(10,100) datafiles(imincomp+6)
c
         iop = imisc(5)
         if(iop.ne.0) then
c           stop the run if it is using optical depth; minzdata expects
c           geometric depth as input
            write(10,104) iop
            write(6,104) iop
            call HERR("Minzdata","cannot use optical depth")  !stop run
         endif
c
c       open and read the CDOM data file
        nudata = 40
        iclose = -1
        nDat = getndat(iclose,nudata,nheadr,nh2p,datafiles(imincomp+6)) 
        allocate(zDat(nDat),minDat(ndat))  !dim to #data
        call loadDat(nudata,nheadr,nDat,itype,datafiles(imincomp+6),
     1               zDat, minDat, dum)
        kall = 1
        i1 = 1
      endif
c-----  End of initialization on first call----------------------------

c-----  Subsequent calls start here:
c     fminconc from interp routine:
      fminconc =  yinterp(i1, nDat, depth, zDat, minDat)
      return
c
c     Formats:
  100 format(//5x,
     1'The mineral particle concentration is obtained from a file of',
     2' measured values read from a file named:'/
     38x,a120,/5x,
     4'(Linear interpolation is used to define conc(z) at all',
     5' depths)')
  104 format(//'Error in sub minzdata:  called with optical depth: iop
     1= ',i2) 
      end
