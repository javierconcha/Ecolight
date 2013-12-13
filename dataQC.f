C     Last change:  LKS   7 Aug 2008    7:26 pm

      subroutine RangeCheck(iopt, zmin, zmax, wmin, wmax, datafile)
c     Checks range of data; stores needed info to print 
!     NOTE:  a z or w value =0 signifies z or w not abscissa for data

!     iopt = 0 means data outside bounds will be set = 0 
!            1 means data outside bounds will be set = value at bounds
c
      INCLUDE "DIMENS_XL.INC"
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      COMMON /Cmisc/ imisc(30),fmisc(30)

c    INPUT:
      real zmin,zmax,wmin,wmax   !range limits of data file
      real zRmin,zRmax,wRmin,wRmax   !range limits of HE RUN
      character*120 datafile
C     OUPUT variables shared with RangeWarn via COMMON block
      COMMON /CRANGEW/ nkall, ierr, zdat, wdat, dfwarn, ibc
      integer nkall, ierr(2*mxcomp+8, 4), ibc(2*mxcomp+8)
      real zdat(2*mxcomp+8, 2), wdat(2*mxcomp+8, 2)
      character*120 dfwarn(2*mxcomp+8)

c      temp variables:
      integer i, kall, ierrT, nz, nwave
c
      data kall/0/
      save kall 
c     *************************************************************
c
      nkall = kall
      kall = kall + 1
c
c      get RUN range limits
      nz = imisc(4)
      nwave = imisc(7)
      zRmin = 0.0
      zRmax = zeta(nz)
      wRmin = wave(1)
      wRmax = wave(nwave)
c     re-init err code flags
      Do i=1,4
        ierr(kall,i) = 0
      Enddo
c
c     Check data limits
      If(zmin.gt.zRmin) ierr(kall,1) = 1
      If(zmax.gt.0 .and. zmax.lt.zRmax) ierr(kall,2) = 1
      If(wmin.gt.wRmin) ierr(kall,3) = 1
      If(wmax.gt.0 .and. wmax.lt.wRmax) ierr(kall,4) = 1
c     ierr(0) flags whether there are any limit errors
      ierrT = ierr(kall,1) + ierr(kall,2) + 
     1        ierr(kall,3) + ierr(kall,4)
c
c     Do we need ANY warnings for this file?
      If(ierrT.eq.0) then
        kall = kall - 1     !on next call, overwrite vars
        return
      Else
        dfwarn(kall) = trim(datafile)
        ibc(kall) = iopt
        zdat(kall,1) = zmin        
        zdat(kall,2) = zmax        
        wdat(kall,1) = wmin        
        wdat(kall,2) = wmax        
        nkall = kall
      Endif

      end subroutine

!_______________________________________________________________________

      subroutine RangeWarn(nuout)
c     Checks range of data; prints warnings if data does not cover full
c     range for run
!     NOTE:  a z or w value =0 signifies z or w not abscissa for data
!     iopt = 0 means data outside bounds will be set = 0 
!            1 means data outside bounds will be set = value at bounds
c
      INCLUDE "DIMENS_XL.INC"

c    INPUT:
      integer nuout              
c
C     INPUT variables shared with RangeWarn via COMMON block
      COMMON /CRANGEW/ nkall, ierr, zdat, wdat, dfwarn, ibc
      integer nkall, ierr(2*mxcomp+8, 4), ibc(2*mxcomp+8)
      real zdat(2*mxcomp+8, 2), wdat(2*mxcomp+8, 2)
      character*120 dfwarn(2*mxcomp+8)
c     *************************************************************
c
      if(nkall.eq.0) return   !no warnings needed

c     print header
      write(nuout, 200) nkall
c
      Do kall = 1, nkall
        write(10,100) trim(dfwarn(kall))
        If(ibc(kall).eq.1) then
         If(ierr(kall,1).gt.0) write(nuout,101)zdat(kall,1),zdat(kall,1)
         If(ierr(kall,2).gt.0) write(nuout,102)zdat(kall,2),zdat(kall,2)
         If(ierr(kall,3).gt.0) write(nuout,103)wdat(kall,1),wdat(kall,1)
         If(ierr(kall,4).gt.0) write(nuout,104)wdat(kall,2),wdat(kall,2)
        Else
         If(ierr(kall,1).gt.0) write(nuout,301)zdat(kall,1)
         If(ierr(kall,2).gt.0) write(nuout,302)zdat(kall,2)
         If(ierr(kall,3).gt.0) write(nuout,303)wdat(kall,1)
         If(ierr(kall,4).gt.0) write(nuout,304)wdat(kall,2)
        Endif
      Enddo
c
c     print summary
      write(nuout,201)  
      return
      
  200 format(//2x,"***** WARNING MESSAGES FOR INPUT DATA FILES *****",
     1       /5x,"NUMBER OF FILES WITH WARNINGS:  ",i3)
  100 format(/5x, "Data file: ",a,
     1       /8x," does not extend over the full ",
     1       "range of depths and wavelengths selected for this run."/)
c
  101 format(8x,"For depths less than ",f5.1," m the value at ",
     1       f5.1, " m will be used")
  102 format(8x,"For depths greater than ",f5.1," m the value at ",
     1       f5.1, " m will be used")
  103 format(8x,"For wavelengths less than ",f5.1," nm the value at ",
     1       f5.1," nm will be used")
  104 format(8x,"For wavelengths greater than ",f5.1," nm the value at ",
     1       f5.1," nm will be used")
c
  301 format(8x,"For depths less than ",f5.1," value set to zero")
  302 format(8x,"For depths greater than ",f5.1," value set to zero")
  303 format(8x,"For wavelengths less than ",f5.1," value set to zero")
  304 format(8x,"For wavelengths greater than ",f5.1,
     1          " value set to zero")
      !finish up this warning section with
  201 format(/5x,"HYDROLIGHT/ECOLIGHT OUTPUTS MAY BE UNPHYSICAL OR ",
     1      "INACCURATE"/5x," FOR DEPTHS AND WAVELENGTHS OUTSIDE ",
     2      "THE RANGE OF THE INPUT DATA."//2x,
     3      "***** END WARNING MESSAGES FOR INPUT DATA FILES *****",
     4      //)

      end subroutine

