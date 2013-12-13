C     Last change:  LKS   7 Aug 2008    6:55 pm
      subroutine PNTBC
c     Prints out bottom BC info to Proot file at beginning of run
c     Called by initial
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
      real depthbot
c********************************************************************
c
      iop = imisc(5)
      ibotm = imisc(12)
c
c     zeta still contains the depths as read in, whether optical or 
c     geometrical
      depthbot = zeta(izirad(npirad))
c
      if(ibotm.eq.0 .and. iop.eq.1) write(10,1020) depthbot
      if(ibotm.eq.0 .and. iop.eq.0) write(10,1021) depthbot
      if(ibotm.ge.1) then
        write(10,1031) depthbot
c       initialize the BRRF routine: 
        call pntBRDFbotm
        R=rbottom(ibotm,0.0)   !init to generate msgs      
      endif 
c
      return
c
 1020 format(/5x,'The bottom boundary is an infinitely deep,',
     1' homogeneous water body below optical depth',f7.2)
 1021 format(/5x,'The bottom boundary is an infinitely deep,',
     1' homogeneous water body below depth',f7.2,' m')
 1031 format(/5x,'The bottom boundary is an opaque reflecting',
     1' surface at depth',f7.2,' m')
      end subroutine
!-----------------------------------------------------------------------

      subroutine PNTDPF(k)
c     prints out the DPF info on initialization of IOP model
c
c     INPUT:    
      integer k  != component index\

c     COMMON BLOCKS:    
      INCLUDE "DIMENS_XL.INC"
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      common /Cbbopt/ ibbopt(mxcomp), bbfrac(mxcomp), 
     1                BfrefPL(mxcomp), Bf0PL(mxcomp), BfmPL(mxcomp)
      integer ibbopt
      real bbfrac, BfrefPL, Bf0PL, BfmPL
c       flag for calculating bb shared with IOP routine (calc bb iff flag=1)
      Common /CbbCalc/ ibbCalc

C     TEMPORARY LOCAL VARIABLES:
      integer nuphas
      character*120 pftitl
c********************************************************************

      IF(ibbopt(k) .eq. 1) then
c       use constant bbfract to select the dpf if ibbopt=1
        write(10,551) k,bbfrac(k)
      ELSEIF(abs(ibbopt(k)) .eq. 2) then
        write(10,553) k,trim(datafiles(1)),trim(datafiles(3))
      ELSEIF(ibbopt(k) .eq. 3) then
        write(10,554) k,Bf0PL(k),BfrefPL(k),BfmPL(k)

      ELSEIF(ibbopt(k) .eq. 0) then
c       read header records from phase function file
        nuphas = 49+k
        OPEN(nuphas,file=pfname(k),status='old', err=999)
        read(nuphas,'(a)') pftitl
        close(nuphas)
        write(10,552) k,trim(pftitl)
      ENDIF     ! end ibbopt(k) control
c
      return

  999 call nofile(nuphas, 'PNTDPF', pfname(k))

  551 format(/5x,'The phase function for component',i3,
     1' will be chosen to have a bb/b ratio of ',f8.4)
  552 format(/5x,'The phase function for component',i3,
     1' comes from a file titled'/5x,a)
  553 format(/5x,'The phase function for component',i3,
     1' will be chosen to have a bb/b ratio ',
     2/5x,'as specified by ac and bb data files: ',
     3/8x,a,'  and',/8x,a)
  554 FORMAT(/5x,'The phase function for component',
     1 i2, ' will be chosen to have a',
     2 /8x,'bb/b ratio calculated by a Power Law model:',
     2 /10x, 'bb/b = ',f6.4,' * (',f4.0,'/wavenm)^(',f6.4,')')
      end subroutine
