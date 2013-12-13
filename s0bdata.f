C     Last change:  LKS   7 Aug 2008    5:33 pm
      function S0bdata(depth,wavenm)
c
!     wavenm not used in this call anymore !! values read from fmisc() array

c     called by S0bioSub <- from Shatbiol
c
c     calls routines to read, sort, and average data found in file abacbb.f
c
c     This routine is called when the user has a data file of measured
c     (or simulated)bioluminescence values at discrete depths and wavel.  
c
c     This bioluminescence data is assumed to be on the Hydrolight "ac-9
c     standard format." 
c
c     When this routine is called for the first time, the data file is
c     read and processed as follows:
c     1) the depth and corresponding bioluminescence value are read for
c	   each depth where bioluminescence was measured
c     2) the discrete-depth bioluminescence values are linearly 
c        interpolated onto a 1nm resolution grid on which sources are evaluated.
c
c     The Bioluminescent source is set to ZERO outside the provided data. 
c
c     INPUT:
c        depth
c        wavelength
c     OUTPUT:
c        S0bdata  the bioluminescence in W/(m^3 nm) 
c
      INCLUDE "DIMENS_XL.INC"
c
      Common /CSourceBdata/  S0datafile    !read in initial
      character*120 S0datafile
      real S0bdata
c
      real, allocatable,dimension(:) :: ydat,zdat,xgrid,ygrid,zgrid
      real, allocatable,dimension(:,:) :: S0grid
      Integer nzdata, nwavdat
      real wmin, wmax
      real, allocatable,dimension(:) :: zdata, wavdat
      real, allocatable,dimension(:,:) :: S0data
c
      Common /CS0xcl/ S0pnt(mxwave,mxz)
      real S0pnt
c
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      common /cmisc/ imisc(30),fmisc(30)
!
      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
!
      integer getndat, getnwavac
      external getndat, getnwavac
c
      data i1/1/, kall/0/, iclose/0/, nheadr/10/,nh2p/2/
c     Note on dynamic storage:  MUST save values between subroutine calls:
      save  kall, S0grid, wmin, wmax, zdata
c
c-----  Begin Initialization on first call  ---------------------------
c
c     On the first call to S0bdata, read in the bioluminescence data on its standard
c     format and fit a 1-D spline to the discrete data.  This spline
c     will be used in the subsequent calls to get bioluminescence at any depth.
      if(kall.eq.0) then
c
c     initialization on first call
         write(10,100) S0datafile
c
         iop = imisc(5)
         if(iop.ne.0) then
c           stop the run if it is using optical depth; S0bdata expects
c           geometric depth as input
            write(10,104) iop
            call HERR("S0bdata","depth type is not compatible")  !stop run
         endif
c
c     open and read the bioluminescence data file, using Filtered ac9 routine
      nudata = 42
      nzdata = getndat(iclose,nudata,nheadr,nh2p, S0datafile)  - 1 
      nwavdat = getnwavac(nudata,nheadr)
      allocate(zdata(nzdata), wavdat(nwavdat))
      allocate(S0data(nzdata, nwavdat))

      call readac9F(nudata,nzdata,nwavdat,
     1              zdata,wavdat,S0data)
c     generate warnings if data range doesn't cover run range
      CALL RangeCheck(0, zdata(1),zdata(nzdata), 
     1                wavdat(1), wavdat(nwavdat), S0datafile )

c     Linear interpolate data onto run wavelength grid!
      itype = 1  !(only fit one column of data; depth is indep var)
      xRef = 0.0
      nwave = imisc(7)
      allocate(yDat(nwavdat),zDat(nwavdat))    !input
      allocate(xgrid(nwave),ygrid(nwave),zgrid(nwave))  !intermediates
      allocate(S0grid(nzdata, nwave))                   !output

      xgrid(1:nwave) = wave(1:nwave)
      Do i=1,nzdata
        yDat = S0data(i,:)
        Call GridDat(itype,nwave,nwavdat,xref,yref,
     1               xgrid, ygrid,zgrid,
     1               wavdat,yDat,zDat)
        S0grid(i,:) = ygrid(:)
      Enddo
      wmin = wavdat(1)
      wmax = wavdat(nwavdat)
      deallocate(xgrid,ygrid,zgrid)  !intermediates
      deallocate(yDat,zDat,wavdat,S0data)    !input
!
!     Store and print out the S0 grid to Proot
      write(10,200) (z(izirad(iiz)),iiz=1,npirad)
      Do i=1,nwave
      do iiz=1,npirad
          depth = z(izirad(iiz))
          S0pnt(i,iiz)=yinterp(i1, nzdata, depth, zdata, S0grid(1,i))
        Enddo
        waven = wave(i)
        write(10,201) waven, (S0pnt(i,j),j=1,npirad)
      Enddo
      i1 = 1   !reset
      kall = 1
      endif
c
c-----  End of initialization on first call----------------------------
c
c-----  Subsequent calls start here:
c
      jwave = imisc(11)
      wavenm = fmisc(13)

c     If the requested depth is outside the bounds of the original 
c     data, bioluminescence set to zero.
      if((depth.lt.zdata(1)).or.(depth.gt.zdata(nzdata)) ) then
          S0bdata = 0.0
          return
      endif
c
c     If the requested wavelen is outside the bounds of the original 
c     data, bioluminescence set to zero.
      if((wavenm.lt.wmin).or.(wavenm.gt.wmax) ) then
          S0bdata = 0.0
          return
      endif
c
      S0bdata = yinterp(i1, nzdata, depth, zdata, S0grid(1,jwave) )
      return
c     Formats:
  100 format(/5x,
     1'The bioluminescence as a function of depth and wavelength ',
     2'is obtained from the file:'/
     38x,a)
  104 format(//'Error in sub S0bdata:  called with optical depth: ',
     1'iop = ',i2) 
 200  format(//8x,'Bioluminescent source So tabulated on depth and ',
     1'wavelength grid [W/(m3 nm)]: ',
     2 //8x,'wavel/depth',/16x,1000f12.1)
 201  format(8x,f8.1,1p1000E12.4)
      end
