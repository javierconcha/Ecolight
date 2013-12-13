C     Last change:  LKS   8 Dec 2011    3:17 pm
c     File excel.f contains the routines that write selected 
c     Hydrolight output to files for postprocessing in Excel 
c     spreadsheets.
c
c     Users can change these routines to alter the data that are 
c     saved in a Hydrolight run for Excel postprocessing.
c
c     This file contains routines
c
c        WRTXCELS (called by RADANAL; creates the Sroot.txt files)
c        WRTXCELM (called by MAIN; creates the Mroot.txt files)
c        STORXCEL (called by RADANAL)
c
c----------------------------------------------------------------------
c
      SUBROUTINE WRTXCELS
c
c     (WRTXCELS = WRite To eXCEL for Single-wavelength format)
c
c     This subroutine writes selected data to an ascii file for later 
c     analysis in a spreadsheet.  The format corresponds to what is
c     expected by the Hydrolight Excel macro singlewl.xcl.  This format 
c     is intended for use with single-wavelength runs.  However, an 
c     option (variable iwrtss1, set in routine SETDFLTS) allows for 
c     the same information to be written for each wavelength in turn.
c
c     All character strings are written enclosed in double quotes "..." 
c     for recognition by Excel
c
c     In-water output is at the user-requested output depths
C
      INCLUDE "DIMENS_XL.INC"
c
c     Cxcl holds arrays accumulated in routine storexcl
      COMMON /Cxcl/  nzxcl,zxcl(mxz),izxcl(mxz),
     1               axcl(mxz,mxwave,0:mxcomp),
     2               bxcl(mxz,mxwave,0:mxcomp),
     3               bbxcl(mxz,mxwave,0:mxcomp),
     4               Edxcl(0:mxz,mxwave),Euxcl(0:mxz,mxwave),
     5               Eodxcl(0:mxz,mxwave),Eouxcl(0:mxz,mxwave),
     5               Rxcl(0:mxz,mxwave),
     6               Raduxcl(0:mxz,mxwave),Radwxcl(mxwave)
c     CKxcl holds arrays accumulated in routine Kfcn
      COMMON /CKxcl/ nzKxcl,zKxcl(mxz),fKdxcl(mxz,mxwave),
     1  fKuxcl(mxz,mxwave),fKoxcl(mxz,mxwave),
     2  fKnetxcl(mxz,mxwave),fKLuxcl(mxz,mxwave)
c
      COMMON /Cirrad/ Eou(0:mxz),Eod(0:mxz),Eu(0:mxz),Ed(0:mxz), 
     1                fMUu(0:mxz),fMUd(0:mxz),fMUtot(0:mxz),R(0:mxz),
     2                E2(0:mxz)
      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cpkfcn/ npkfcn,izkfcn(mxz) 
      COMMON /Cfrstcls/ iabscat,iqasky,iradamps,iradian, iradxcl
      COMMON /Cmisc/ imisc(30),fmisc(30)
      COMMON /Ctitle/ ititle
      Character*120 ititle
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      integer ioD,lenD
c
      character*10 getHE5
      external getHE5
c     --------------------------------------------------------------------
c
c     iwrtss1 = 0 does not write this file
c             = 1 writes this file only for the first wavelength
c             = 2 writes data for all wavelengths; all information
c                 is repeated for wavelengths after the first
c
      nwave = imisc(7)
      iwvl = imisc(11)
      iwrtss1 = imisc(20)
      wavenm = fmisc(13)
c
      nuxcel1 = 13
c
      if(iwrtss1.eq.0 .or. (iwrtss1.eq.1 .and. iwvl.gt.1)) return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Open the file and write the header records
c
c     open file to append data after 1st call (at first wavelength):
c
      if (iradxcl.eq.1) then
         iradxcl = 0
         ioD = 0
         open(nuxcel1,file=Srootname)
      else
c       CM: changed access='append' to position='append' for Lahey:
        open(nuxcel1,file=Srootname,position='append',iostat=ioD,err=66)
      endif
c
c     If cannot append file, open a new file after alerting user.
c     The cause of this rare error is not enough memory allocated to the
c     stack.  The fix for this it to increase the allocation by editing
c     the automake.fig file and adding/icreasing the -stack 500000 flag on
c     the linker line.
 66   IF(ioD.gt.0) then
        lenD = lenstr(Srootname)
        write(10,300) nuxcel1,Srootname(1:lenD)
        write(10,302)
        Srootname = Srootname(1:lenD-4) // 'X.txt'
        lenD = lenstr(Srootname)
        write(10,301) Srootname(1:lenD)
        close(nuxcel1)
        open(nuxcel1,file=Srootname)
      Endif
c
  300 format(/2x,'There has been an unexpected error while attempting',
     1' to open the appending file: unit =',i5,/5x,a)
  301 format(/2x,'A new file will be opened for writing out the',
     1' remaining data, named: ',/5x,a)
  302 format(/2x,'Hydrolight was unable to append this output file.',/
     12x,'The cause of this rare error is not enough memory allocated',
     2' to the stack.',/2x,'The fix for this it to increase the',
     3' allocation by editting the automake.fig file ',/2x,
     4'and adding/icreasing the -stack 500000 flag on the linker line.',
     5//2x,'To salvage this run, a new file will be opened for output.')
c
c     write file header records
c
      if(iwrtss1.eq.1) then
         numwave = 1
      else
         numwave = nwave
      endif

      write(nuxcel1,100) trim(getHE5()), ititle(1:lenstr(ititle))
  100 format('"',a,' Run Title:  ',a,'"')
      write(nuxcel1,101) wavenm
  101 format('"','Wavelength =',f6.1,' nm"')
      nsheets = 4
      write(nuxcel1,102) numwave,nsheets
  102 format(i5,5x,i5,5x,'"(num of wavelengths, num of sheets)"') 
c     blank lines for possible later use
      write(nuxcel1,'(a)') '"line 4: reserved for future use"'
      write(nuxcel1,'(a)') '"line 5: reserved for future use"'
c
c     wavelength block header
    
      write(nuxcel1,106) wavenm,iwvl
  106 format(f6.1,i5,5x,'"(wavelength and wavelength band)"')
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Write the blocks of data that will become worksheets in the
c     Excel workbook.  The first block of data written here will
c     become the _rightmost_ worksheet (on the worksheet tabs at the
c     bottom of the workbook), and the last block of data will
c     become the _leftmost_ worksheet.
c
c     Worksheet 4:  K functions
c
      nrows = nzKxcl
      ncols = 6
      write(nuxcel1,410) nrows,ncols,
     1  '"(sheet title, num rows, num cols)"'
  410 format('"K funcs"',5x,i3,5x,i3,5x,a)
      write(nuxcel1,412)
  412 format(' "depth"',1x,'"Kd"',1x,'"Ku"',1x,'"Ko"',1x,'"Knet"',
     1   1x,'"KLu"')
      write(nuxcel1,414)
  414 format(' "(m)"',1x,'"(1/m)"',1x,'"(1/m)"',1x,'"(1/m)"',1x,
     1'"(1/m)"',1x,'"(1/m)"')

      do i=1,nzKxcl
      write(nuxcel1,416) zKxcl(i),fKdxcl(i,iwvl),fKuxcl(i,iwvl),
     1  fKoxcl(i,iwvl),fKnetxcl(i,iwvl),fKLuxcl(i,iwvl)
  416 format(f8.3,1p,5e13.4)
      end do
c
c     Worksheet 3:  AOPs in air and in water
c
      nrows = nzxcl + 1
      ncols = 6
      write(nuxcel1,310) nrows,ncols,
     1  '"(sheet title, num rows, num cols)"'
  310 format('"AOPs"',5x,i3,5x,i3,5x,a)
      write(nuxcel1,312)
  312 format(' "depth"',1x,'"Lu/Ed"',1x,'"R=Eu/Ed"',1x,'"Mud"',
     1   1x,'"Muu"',1x,'"Mu"')
      write(nuxcel1,314)
  314 format(' "(m)"',1x,'"(1/sr)"',1x,'" "',1x,'" "',1x,'" "',1x,'" "')
c
      write(nuxcel1,315) Raduxcl(0,iwvl)/Edxcl(0,iwvl),
     1   Rxcl(0,iwvl),fMUd(0),fMUu(0),fMUtot(0)
  315 format(' "in air"',1p,e13.4,0p,4f9.4)
      do i=1,nzxcl
c     must use explicit indexing for mu's, since they were not stored before
      iz = izxcl(i)
      write(nuxcel1,316) zxcl(i),Raduxcl(i,iwvl)/Edxcl(i,iwvl),
     1Rxcl(i,iwvl),fMUd(iz),fMUu(iz),fMUtot(iz)
  316 format(f8.3,1p,e13.4,0p,4f9.4)
      end do
c
c     Worksheet 2:  irradiances and Lu in air and in water
c
      nrows = nzxcl + 1
      ncols = 6
      write(nuxcel1,210) nrows,ncols,
     1  '"(sheet title, num rows, num cols)"'
  210 format('"IRRAD"',5x,i3,5x,i3,5x,a)
      write(nuxcel1,212)
  212 format(' "depth"',1x,'"Ed"',1x,'"Eu"',1x,'"Eo"',1x,'"Lu"',
     1       1x,'"Lu/Ed"')
      write(nuxcel1,214)
  214 format(' "(m)"',1x,'"(W/m^2 nm)"',1x,'"(W/m^2 nm)"',1x,
     1   '"(W/m^2 nm)"',1x,'"(W/m^2 sr nm)"',1x,'"(1/sr)"')
c
      write(nuxcel1,215) Edxcl(0,iwvl),Euxcl(0,iwvl),
     1  Eouxcl(0,iwvl)+Eodxcl(0,iwvl),Raduxcl(0,iwvl),
     2  Raduxcl(0,iwvl)/Edxcl(0,iwvl)
  215 format(' "in air"',1p,5e13.4)
      do i=1,nzxcl
      write(nuxcel1,216) zxcl(i),Edxcl(i,iwvl),Euxcl(i,iwvl),
     1  Eouxcl(i,iwvl)+Eodxcl(i,iwvl),Raduxcl(i,iwvl),
     2  Raduxcl(i,iwvl)/Edxcl(i,iwvl)
  216 format(f8.3,1p,5e13.4)
      end do
c
c     Worksheet 1:  IOP values
c
      nrows = nzxcl
      ncols = 7			!modified 2/21/00
      write(nuxcel1,110) nrows,ncols,
     1  '"(sheet title, num rows, num cols)"'
  110 format('"IOPs"',5x,i3,5x,i3,5x,a)
      write(nuxcel1,112)
  112 format(' "depth"',1x,'"a"',1x,'"b"',1x,'"c"',1x,'"omega0"',1x,
     1'"bb"', 1x, '"bb/b"')
      write(nuxcel1,114)
  114 format(' "(m)"',1x,'"(1/m)"',1x,'"(1/m)"',1x,'"(1/m)"',1x,'" "',
     11x,'"(1/m)"')
c
c     IOP's at output depths
c
      do i=1,nzxcl
      acoef = axcl(i,iwvl,0)
      bcoef = bxcl(i,iwvl,0)
      ccoef = acoef + bcoef
      albss = bcoef/ccoef
      bbcoef = bbxcl(i,iwvl,0)
      bbob = bbcoef / bcoef		!added 2/21/00
      write(nuxcel1,116) zxcl(i),acoef,bcoef,ccoef,albss,bbcoef,bbob	!modified 2/21/00
  116 format(f8.3,4f12.4,2f12.5)
      enddo
      close(nuxcel1)
c
      return
      end
c
c----------------------------------------------------------------------
c
      SUBROUTINE WRTXCELM
c
c     WRTXCELM = WRite To eXCEL for Multiple wavelengths
c
c     This subroutine is called once at the end of a multi-wavelength run
c     to write selected data to an ascii file for later 
c     analysis in a spreadsheet.  The format corresponds to what is
c     expected by the Hydrolight Excel macro multiwl.xcl.  This format
c     is intended for use with multi-wavelength runs, where 2-D plots
c     are to be made as a function of wavelength for various depths, or
c     where 3-D plots are to be made as a function of depth and wavelength.   
c
c     All character strings are written enclosed in double quotes "..." 
c     for recognition by Excel
c
c     In-water output is at the user-requested depths
c 
c     The arrays written out are accumulated in routines "storexcl" and
c     "Kfcn" as Hydrolight cycles over the wavelengths of a multi-
c     wavelength run
c
      INCLUDE "DIMENS_XL.INC"
c
C     The following common block array was written for Worksheet 5,
C	to store and printout the above-water direct and diffuse Ed's
C     at each wavelength (block shared with qasky.f only)
      common /CErik/ EdifOut(mxwave), EdirOut(mxwave)
c
      Common /CPAR/fPar(0:mxz),Edpar(0:mxz),fKpar(mxz),
     1             Eoquant(0:mxz,mxwave)   !shared with PAR.f
c
      Common /CS0xcl/ S0pnt(mxwave,mxz)    !shared with S0bdata
      Common /Csource0/ ibiolum,ichlfl,icdomfl,iraman, ramanEXP
      real S0pnt, ramanEXP
      integer ibiolum,ichlfl,icdomfl,iraman
c
c     Cxcl holds arrays accumulated in routine storexcl
      COMMON /Cxcl/  nzxcl,zxcl(mxz),izxcl(mxz),
     1               axcl(mxz,mxwave,0:mxcomp),
     2               bxcl(mxz,mxwave,0:mxcomp),
     3               bbxcl(mxz,mxwave,0:mxcomp),
     4               Edxcl(0:mxz,mxwave),Euxcl(0:mxz,mxwave),
     5               Eodxcl(0:mxz,mxwave),Eouxcl(0:mxz,mxwave),
     5               Rxcl(0:mxz,mxwave),
     6               Raduxcl(0:mxz,mxwave),Radwxcl(mxwave)
c     CKxcl holds K-function arrays accumulated in routine Kfcn
      COMMON /CKxcl/ nzKxcl,zKxcl(mxz),fKdxcl(mxz,mxwave),
     1  fKuxcl(mxz,mxwave),fKoxcl(mxz,mxwave),
     2  fKnetxcl(mxz,mxwave),fKLuxcl(mxz,mxwave)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      COMMON /Cmisc/ imisc(30),fmisc(30)
      COMMON /Ctitle/ ititle
      Character*120 ititle
c
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
c
      logical IamEL
      external IamEL
c
      character*10 getHE5
      external getHE5
c     --------------------------------------------------------------------
c
      nwave = imisc(7)
      nconc = imisc(23)
c      wavenm = fmisc(13)
c
      nuxcelM = 12
c
c     open the file
c
      open(nuxcelM,file=Mrootname)
c
c     write file header records
c
      write(nuxcelM,100) trim(getHE5()), ititle(1:lenstr(ititle))
  100 format('"',a,' Run Title:  ',a,'"')

c	set the number of sheets to be loaded into the spreadsheet
      nsheets = 20
      if(ibiolum.gt.0) nsheets = nsheets + 1

      write(nuxcelM,102) nwave,nsheets
  102 format(i5,5x,i5,5x,'"(num of wavelengths, num of sheets)"') 
      write(nuxcelM,103)
  103 format('"wavelength (nm)"',1x,'"depth (m)"')
      write(nuxcelM,'(a)') '"line 5: reserved for future use"'
c
c     Write the blocks of data that will become worksheets in the
c     Excel workbook.  The first block of data written here will
c     become the _rightmost_ worksheet (on the worksheet tabs at the
c     bottom of the workbook), and the last block of data will
c     become the _leftmost_ worksheet.
c     --------------------------------------------------
c
c     [optional] Worksheet 21:  S0data
      if(ibiolum.gt.0) then
        ncols = nzxcl+1
        write(nuxcelM,2110) nwave,ncols
 2110   format('"S0"',1x,'"S0 W/(m^3 nm)"',5x,i3,5x,i3,5x,a)
        write(nuxcelM,*)'" "'
        write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
        do j=1,nwave
          write(nuxcelM,114) wave(j),(S0pnt(j,i),i=1,nzxcl)
        end do
      endif

c     Worksheet 20:  KPAR     !note ipirad depths and xcl depths are the same
      ncols = 2
      write(nuxcelM,1810) nzKxcl,ncols
      if(waveb(1).le.fmisc(26) .and. waveb(nwave+1).ge.fmisc(27)) then
        write(nuxcelM,*)'" "'
      else
        write(nuxcelM,*)'"KPAR NOT CALCULATED FOR THIS RUN"'
      endif
 1810 format('"KPAR"',1x,'"KPAR (mumol phot/m^2 s)"',5x,i3,5x,i3,5x,a)
      if(IamEL()) then            !EL
        write(nuxcelM,1812)
      else                        !HL
        write(nuxcelM,1813)
      endif
 1812 format('" " "depth" "K_PAR (from Eo, layer-averaged)"')
 1813 format('" " "depth" "K_PAR (from Eo)"')
      do i=1,nzKxcl
        iz = izxcl(i)
        write(nuxcelM,1714) zKxcl(i), fKpar(iz) 
      enddo
c
c     Worksheet 19:  PAR     !note ipirad depths and xcl depths are the same
      ncols = 3
      write(nuxcelM,1710) nzxcl+1,ncols
      if(waveb(1).le.fmisc(26) .and. waveb(nwave+1).ge.fmisc(27)) then
        write(nuxcelM,*)'" "'
      else
        write(nuxcelM,*)'"PAR NOT CALCULATED FOR THIS RUN"'
      endif
 1710 format('"PAR"',1x,'"PAR (mumol phot/m^2 s)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,1712)
 1712 format('" " "depth" "PAR_Eo"  "PAR_Ed" ')
      write(nuxcelM,1713) fPar(0),Edpar(0)
 1713 format('" " "in air"', 1p2E15.4) 
      do i=1,nzxcl
         iz = izxcl(i)
        write(nuxcelM,1714) zxcl(i),fPar(iz),Edpar(iz) 
      enddo
 1714 format(f9.3,1p,2E15.4)
c

c     Worksheet 18:  diffuse attenuation KLu
c
      ncols = nzKxcl + 1
      write(nuxcelM,1610) nwave,ncols
      write(nuxcelM,*)'" "'
 1610 format('"KLu"',1x,'"KLu (1/m)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,112) (zKxcl(i),i=1,nzKxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(fKLuxcl(i,j),i=1,nzKxcl)
      end do
c

c     Worksheet 17:  diffuse attenuation Ku
c
      ncols = nzKxcl + 1
      write(nuxcelM,1510) nwave,ncols
      write(nuxcelM,*)'" "'
 1510 format('"Ku"',1x,'"Ku (1/m)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,112) (zKxcl(i),i=1,nzKxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(fKuxcl(i,j),i=1,nzKxcl)
      end do
c

c     Worksheet 16:  diffuse attenuation Kd
c
      ncols = nzKxcl + 1
      write(nuxcelM,1410) nwave,ncols
      write(nuxcelM,*)'" "'
 1410 format('"Kd"',1x,'"Kd (1/m)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,112) (zKxcl(i),i=1,nzKxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(fKdxcl(i,j),i=1,nzKxcl)
      end do
c
     
c     Worksheet 15:  Scalar irradiance Eo
c
      ncols = nzxcl + 2
      write(nuxcelM,1310) nwave,ncols
      write(nuxcelM,*)'" "'
 1310 format('"Eo_quantum"',1x,'"Eo [microEinst/(m^2 s)]"',5x,
     1       i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(Eoquant(izxcl(i),j),i=0,nzxcl)
      end do

c     Worksheet 14:  Remote-sensing reflectance Rrs and related quantities
c                   in the air
c
      ncols = 5
      write(nuxcelM,1210) nwave,ncols
      write(nuxcelM,*)'" "'
 1210 format('"Rrs"',1x,'"Rrs (1/sr)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,1212)
 1212 format('" " "in air"',1x,'"Rrs"',1x,'"Ed"',1x,'"Lw"',1x,'"Lu"')
      do j=1,nwave
      Rrs = Radwxcl(j)/Edxcl(0,j)
      write(nuxcelM,114) wave(j),Rrs,Edxcl(0,j),Radwxcl(j),
     1  Raduxcl(0,j)
      end do
c
c     Worksheet 13:  Irradiance reflectance R
c
      ncols = nzxcl + 2
      write(nuxcelM,1110) nwave,ncols
      write(nuxcelM,*)'" "'
 1110 format('"R"',1x,'"R = Eu/Ed"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(Rxcl(i,j),i=0,nzxcl)
      end do
c
c     Worksheet 12:  Lu/Ed
c
      ncols = nzxcl + 2
      write(nuxcelM,1010) nwave,ncols
      write(nuxcelM,*)'" "'
 1010 format('"Lu over Ed"',1x,'"Lu/Ed (1/sr)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
        write(nuxcelM,114) wave(j),(Raduxcl(i,j)/Edxcl(i,j),i=0,nzxcl)
      enddo
c
c
c     Worksheet 11: Upwelling plane irradiance Eu
c
      ncols = nzxcl + 2
      write(nuxcelM,910) nwave,ncols
      write(nuxcelM,*)'" "'
  910 format('"Eu"',1x,'"Eu (W/m^2 nm)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(Euxcl(i,j),i=0,nzxcl)
      end do
c
c
c     Worksheet 10:  Downwelling plane irradiance Ed
c
      ncols = nzxcl + 2
      write(nuxcelM,810) nwave,ncols
      write(nuxcelM,*)'" "'
  810 format('"Ed"',1x,'"Ed (W/m^2 nm)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(Edxcl(i,j),i=0,nzxcl)
      end do
c
c
c     Worksheet 9:  Scalar irradiance Eo
c
      ncols = nzxcl + 2
      write(nuxcelM,710) nwave,ncols
      write(nuxcelM,*)'" "'
  710 format('"Eo"',1x,'"Eo (W/m^2 nm)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(Eouxcl(i,j)+Eodxcl(i,j),i=0,nzxcl)
      end do
c
c
c
c     Worksheet 8:  Scalar irradiance Eou
c
      ncols = nzxcl + 2
      write(nuxcelM,711) nwave,ncols
      write(nuxcelM,*)'" "'
  711 format('"Eou"',1x,'"Eou (W/m^2 nm)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(Eouxcl(i,j),i=0,nzxcl)
      end do
c
c
c
c     Worksheet 7:  Scalar irradiance Eod
c
      ncols = nzxcl + 2
      write(nuxcelM,712) nwave,ncols
      write(nuxcelM,*)'" "'
  712 format('"Eod"',1x,'"Eod (W/m^2 nm)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(Eodxcl(i,j),i=0,nzxcl)
      end do
c
c
c     Worksheet 6:  Lu
c
      ncols = nzxcl + 2
      write(nuxcelM,620) nwave,ncols
      write(nuxcelM,*)'" "'
  620 format('"Lu"',1x,'"Lu (W/m^2 nm sr)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,612) (zxcl(i),i=1,nzxcl)
  612 format('" " "wavel" "in air"',1000f9.3)
      do j=1,nwave
        write(nuxcelM,114) wave(j),(Raduxcl(i,j),i=0,nzxcl)
      enddo
c
c
c     Workheet 5:  Direct and Diffuse Ed just above the sea surface
c
      ncols = 4
      write(nuxcelM,510) nwave,ncols
  510 format('"Ed_in_air"',1x,'"Ed in air (W/m^2 nm)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,*)'" "'
      write(nuxcelM,512)
  512 format('" " "wavel" "Ed_diffuse" "Ed_direct" "Ed_total"')
      do j=1,nwave
  514   format(f6.1,1p,3e12.4)
        write(nuxcelM,514) wave(j), EdifOut(j), EdirOut(j),
     1                     EdifOut(j)+EdirOut(j)
      end do
c
c     Worksheet 4:  backscattering ratio bb/b
c
      ncols = nzxcl + 1
      write(nuxcelM,410) nwave + nconc*(nwave+3),ncols
      write(nuxcelM,*)'" "'
  410 format('"bb fraction"',1x,'"backscat ratio bb/b"',5x,i3,5x,i3,
     1       5x,a)
      write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
       write(nuxcelM,114) wave(j),(bbxcl(i,j,0)/bxcl(i,j,0),i=1,nzxcl)
      end do
c
      Do ii=1,nconc
        write(nuxcelM,*)'" "'
        write(nuxcelM,411) ii
  411  format('" " "bb/b ratio for component ', i2,'"')
       write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
       do j=1,nwave
         write(nuxcelM,114) wave(j),(bbxcl(i,j,ii)/bxcl(i,j,ii),
     1                               i=1,nzxcl)
       end do
      Enddo

c
c     Worksheet 3:  backscattering coef bb
c
      ncols = nzxcl + 1
      write(nuxcelM,310) nwave + nconc*(nwave+3),ncols
      write(nuxcelM,*)'" "'
  310 format('"bb"',1x,'"backscat coef b (1/m)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(bbxcl(i,j,0),i=1,nzxcl)
      end do
c
      Do ii=1,nconc
       write(nuxcelM,*)'" "'
       write(nuxcelM,311) ii
  311  format('" " "bb (1/m) for component ', i2,'"')
       write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
       do j=1,nwave
         write(nuxcelM,114) wave(j),(bbxcl(i,j,ii),i=1,nzxcl)
       end do
      Enddo

c     Worksheet 2:  scattering coef b
c
      ncols = nzxcl + 1
      write(nuxcelM,210) nwave + nconc*(nwave+3),ncols
      write(nuxcelM,*)'" "'
  210 format('"b"',1x,'"scat coef b (1/m)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(bxcl(i,j,0),i=1,nzxcl)
      end do
c
      Do ii=1,nconc
       write(nuxcelM,*)'" "'
       write(nuxcelM,211) ii
  211  format('" " "b (1/m) for component ', i2,'"')
       write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
       do j=1,nwave
         write(nuxcelM,114) wave(j),(bxcl(i,j,ii),i=1,nzxcl)
       end do
      Enddo
c
c     Workheet 1:  absorption coef a
c
      ncols = nzxcl + 1
      write(nuxcelM,110) nwave + nconc*(nwave+3),ncols
      write(nuxcelM,*)'" "'
  110 format('"a"',1x,'"abs coef a (1/m)"',5x,i3,5x,i3,5x,a)
      write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
  112 format('" " "wavelen"',1000f9.3)
      do j=1,nwave
      write(nuxcelM,114) wave(j),(axcl(i,j,0),i=1,nzxcl)
  114 format(f6.1,1p,1000(1x,e12.4))
      end do
c
      Do ii=1,nconc
       write(nuxcelM,*)'" "'
       write(nuxcelM,111) ii
  111  format('" " "a (1/m) for component ', i2,'"')
       write(nuxcelM,112) (zxcl(i),i=1,nzxcl)
       do j=1,nwave
         write(nuxcelM,114) wave(j),(axcl(i,j,ii),i=1,nzxcl)
       end do
      Enddo
c
c
      close(nuxcelM)
c
      return
      end

