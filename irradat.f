C     Last change:  LKS  15 Dec 2010    5:07 pm
      subroutine irradat(itype, Ed1, Ed2, filen)
c
c	called by qasky -> skyirrad (in incfiles.for)
c
c     This routine opens the data file specified by the user in the frontend
c     and reads in all the incidient irradiance from the file with the format:
c     10 header records; 
c     wavelength   Edtotal          !if itype = 0, or
c     wavelength   Eddif    Eddir   !if itype = 2
c
c     If itype=0, RADTRANX is called to partition the total irradiance
c     into direct and diffuse components.
c
      INCLUDE "DIMENS_XL.INC"

!     Declare argument list 
      integer itype
      real Ed1, Ed2
      CHARACTER (LEN = 120 ) filen

!     COMMON blocks
      common /cmisc/ imisc(30),fmisc(30)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 

!     Declare local arrays
      real, allocatable ::  wavdat(:),Edat(:)       !data as read in
      real, allocatable ::  zdat(:)         !optional 2nd array of data iff itype=2
      real Edgrid(mxwave)                  !gridded Edat
      real Ed2grid(mxwave)                 !opt gridded zdat

!     local variables needed for calls to data readers
      integer iclose, nheadr, nh2p, navg, i1, wmin, wmax, it
      real edd, edz
      integer getndat
      external getndat
c
      data kall/0/, iclose/0/, nheadr/10/,nh2p/3/, i1/1/

      save kall, Edgrid, Ed2grid

c---------------------------------------
      nwave = imisc(7)
c
      If(kall.eq.0) then
c        write(10,100) trim(filen)
        nudata = 45
        ndata = getndat(iclose,nudata,nheadr,nh2p,trim(filen))
        allocate(wavdat(ndata),Edat(ndata),zdat(ndata))
        it = abs(itype)  
        call LoadDat(nudata,nheadr,ndata,it,filen,wavDat,Edat,zDat)
        if(itype.eq.0) zDat(:)=0.0    !fill array if not used   

!       If only one data point read, assume it is a 1nm bandwidth lidar source
        If(ndata.eq.1) then
          del = 1.0e-6
c         copy data pt to pt #3
          wavdat(3) = wavdat(1)
          edd = Edat(1)
          edz = zdat(1)
c         generate extra data points at +/-0.5nm to make our tophat
          wavdat(1) =  wavdat(3) - 0.5 - del
          wavdat(2) =  wavdat(3) - 0.5 
          wavdat(4) =  wavdat(3) + 0.5 - del
          wavdat(5) =  wavdat(3) + 0.5 
          Edat(1) = 0.0
          zdat(1) = 0.0
          Edat(5) = 0.0
          zdat(5) = 0.0
          Do i=2,4
            Edat(i) = edd
            zdat(i) = edz
          Enddo
          ndata = 5
        Endif

!       get Ed on the run grid
        if(itype.eq.0) then
          write(10,110) trim(filen)
        elseif(itype.gt.0) then
          write(10,111) trim(filen)
        else
          write(10,112) trim(filen)
        endif
        Do jwave=1,nwave
c         If the requested wavelength is outside the bounds of the original
c         wavelength data, then the value at the nearest wavelength is used.
          if(waveb(jwave+1).le.wavdat(1)) then
            Edgrid(jwave)  = Edat(1)
            Ed2grid(jwave) = zdat(1)
          elseif(waveb(jwave).ge.wavdat(ndata)) then
            Edgrid(jwave)  = Edat(ndata)
            Ed2grid(jwave) = zdat(ndata)
          else
            waven = wave(jwave)
c           Find the upper and lower bounds of datapoints CONTAINED within 
c           the waveband of interest, interp on 1nm res and average over bounds
            edd = 0.0
            edz = 0.0
            del = waveb(jwave+1)-waveb(jwave)   !our range
            navg=max(10,int(del) ) !at least 10 steps or 1nm
            del = del/float(navg)
            Do ii = 1,navg   !avg over 1nm interp or less
              waven = waveb(jwave)+ (float(ii)-0.5)*del     !all pts within wbounds
              edd = edd + yinterp(i1, ndata, waven, wavDat, Edat)
              edz = edz + yinterp(i1, ndata, waven, wavDat, zdat)
            Enddo
            Edgrid(jwave) = edd/float(navg)
            Ed2grid(jwave) = edz/float(navg)
          endif
          if(itype.eq.0) then
            write(10,120) waveb(jwave),waveb(jwave+1),wave(jwave),
     1                Edgrid(jwave)
          else
            write(10,120) waveb(jwave),waveb(jwave+1),wave(jwave),
     1                Edgrid(jwave),Ed2grid(jwave)
          endif
        ENDDO
        deallocate(wavdat, Edat, zdat)
        kall = 1
      Endif
c-----  End of initialization on first call----------------------------

c-----  Subsequent calls start here:
c     Get current wavelength
      jwave = imisc(11)
c     get Ed from gridded values
      Ed1 = Edgrid(jwave)
      Ed2 = Ed2grid(jwave)

      return

c     Formats:
c  100 format(5x,
c     1'The incident total (direct + diffuse) sky irradiance Ed is '
c     2'obtained from user data file: ',a,/8x)
  110 format(//5x,
     1'The incident total (direct + diffuse) sky irradiance is '
     2'obtained from user data file: ',//8x,a,//5x,
     3'RADTRANX will be called to partition the total into direct'
     4' and diffuse parts.'//
     5/8x,'Input values interpolated and averaged over each ',
     6'run waveband:',
     7//8x,' wave_min    wave_max    wave_center     Ed_total',
     8 /8x,'    (nm)        (nm)         (nm)       (W/m^2/nm)',/)
  111 format(//5x,
     1'The incident direct and diffuse sky irradiances are '
     2'obtained from user data file: ',//8x,a,//5x,
     3'RADTRANX will NOT be called.'//
     4/8x,'Input values interpolated and averaged over each ',
     5'run waveband:',
     6//8x,' wave_min    wave_max    wave_center     Ed_direct  ',
     7'   Ed_diffuse',
     8 /8x,'    (nm)        (nm)         (nm)       (W/m^2/nm)  ',
     9'   (W/m^2/nm)',/)
  112 format(//5x,
     1'The incident total (direct + diffuse) sky irradiance and '
     2'direct fraction (Ed_direct/Ed_total) are obtained from '
     3'user data file: ',//8x,a,//5x,'RADTRANX will NOT be called.'
     4/8x,'Input values interpolated and averaged over each ',
     5'run waveband:',
     6//8x,' wave_min    wave_max    wave_center     Ed_total',
     7'    Ed_dir/Ed_tot',
     8 /8x,'    (nm)        (nm)         (nm)       (W/m^2/nm)',
     9'    (nondimen)',/)
  120 format(8x,f8.1, 4x,f8.1, 4x,f10.2, 4x, 1pE12.5, 3x,E12.5)
      end