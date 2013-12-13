      SUBROUTINE SAVEDATA
c
c     core routine on file savedata.f
c
c     called by RADANAL [MAIN->RADANAL->SAVEDATA]
c
c     This routine writes the EcoLight Droot.txt files.

c     This subroutine writes all pertinent data to a file for later
c     graphical or other analysis.  For convenience, a complete set
c     of data is written for each wavelength, even though the grid
c     data are independent of the wavelength.
c
c     Array elements are written in the order that is most efficient
c     for reading in by IDL graphics routines (see routine READALL.pro
c     in the IDL directory).
c
c     data format changed form 10e13.5 to 10es15.5e3 Aug 2013 for H v5.2
c     to allow for 3 numbers in the exponent, e.g., 1.23456E-100
c
      INCLUDE "DIMENS_XL.INC"
c
      integer ioD, lenD
c
      COMMON /Ciop/ acoef(mxz,0:mxcomp),bcoef(mxz,0:mxcomp),
     1		      atten(mxz),albedo(mxz), bbcoef(mxz,0:mxcomp)
c     radiance arrays are stored as follows:
c        RADsky(mu) = (old RAD0Pa) = L+(z=a;mu) = total incident sky radiance
c        RADwla(mu) = Lw(mu,z=a) = water-leaving radiance
c        RADrsa(mu) = (old RAD0Ma) = Lreflsky(mu,z=a) = surface-reflected sky rad
c        RADupa(mu) = (old RADMa) = Lu(mu,z=a) = total upwelling rad in air
c           note: RADupa = RADwla + RADrsa, i.e., Lu = Lw + Lreflsky in air
c        RADupz(mu,z) =(old RADMz) = L-(mu,z) = total upwelling rad in water = Lu(mu,z)
c        RADdnz(mu,z) = (old RADPz) = L+(mu,z) = total downwelling rad in water = Ld(mu,z)
      COMMON /Crad/ RADsky(mxmu),RADwla(mxmu),RADrsa(mxmu),RADupa(mxmu),
     1              RADupz(mxmu,mxz),RADdnz(mxmu,mxz)
      COMMON /Cirrad/ Eou(0:mxz),Eod(0:mxz),Eu(0:mxz),Ed(0:mxz), 
     1                fMUu(0:mxz),fMUd(0:mxz),fMUtot(0:mxz),R(0:mxz),
     2                E2(0:mxz)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      COMMON /Csky/ iskyflag,skyspecs(mxnsky)
      COMMON /Cfrstcls/ iabscat,iqasky,iradamps,iradian, iradxcl
      COMMON /Cmisc/ imisc(30),fmisc(30)
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
      COMMON /Ctitle/ ititle
      Character*120 ititle
      COMMON /Cpirrad/ npirad,izirad(mxz)
c
      character*10 getHE5
      external getHE5
c     ----------------------------------------------------------------
c
      nmu = imisc(1)
!      nz = imisc(4)        !this is the last computed depth (dynZ changes)
      nz = izirad(npirad)   !this is the last requested output depth computed
      nwave = imisc(7)
      iwvl = imisc(11)
      nurad = imisc(16)
      wavenm = fmisc(13)
      nconc = imisc(23)
c
c     open file to append data after 1st call:
c
      if (iradian.eq.1) then
c
         iradian = 0
	   ioD = 0
         open(nurad,file=Drootname, status = 'unknown')
      else
         open(nurad,file=Drootname,position='append', iostat=ioD,err=66)
      endif
c
c     If cannot append file, open a new file after alerting user.
c     The cause of this rare error is not enough memory allocated to the
c     stack.  The fix for this it to increase the allocation by editting
c     the automake.fig file and adding/icreasing the -stack 500000 flag on
c     the linker line.
 66   IF(ioD.gt.0) then
		lenD = lenstr(Drootname)
		write(10,300) nurad,Drootname(1:lenD)
          write(10,302)
		Drootname = Drootname(1:lenD-4) // 'X.txt'
		lenD = lenstr(Drootname)
		write(10,301) Drootname(1:lenD)
		close(nurad)
		open(nurad,file=Drootname)
	Endif
c
c     grid information
c
 67   write(nurad,fmt='(a)') trim(getHE5()) // " Run Title: " 
     1                       // ititle(1:lenstr(ititle))
      write(nurad,fmt='(a)') "NOTE: This file is formatted for reading b
     1y IDL routine read_E_Dfile.pro, v5.2"

      write(nurad,fmt='(" nmu, nz, nwave, nconc =",4i5)')
     1 nmu,nz,nwave, nconc
      write(nurad,fmt='(" wavelength band",i3,
     1 "   nominal wavelength =",f8.2)') iwvl, wavenm
      write(nurad,fmt='(a)') 'imisc'
      write(nurad,102) (imisc(i),i=1,30)
      write(nurad,fmt='(a)') 'fmisc'
      write(nurad,104) (fmisc(i),i=1,30)
      write(nurad,fmt='(a)') 'fmu (quad-center mu values)'
      write(nurad,104) (fmu(i),i=1,nmu)
      write(nurad,fmt='(a)') 'zeta (optical depths of output)'
      write(nurad,104) (zeta(i),i=1,nz)
      write(nurad,fmt='(a)') 'z (geometric depths of output)'
      write(nurad,104) (z(i),i=1,nz)
      write(nurad,fmt='(a)') 'bndmu (quad boundary mu values)'
      write(nurad,104) (bndmu(i),i=1,nmu)
      write(nurad,fmt='(a)') 'omega (quad solid angles)'
      write(nurad,104) (omega(i),i=1,nmu)
      write(nurad,fmt='(a)') 'wave (wavelength band-center values)'
      write(nurad,104) (wave(i),i=1,nwave)
      write(nurad,fmt='(a)') 'waveb (wavelength band-boundary values)'
      write(nurad,104) (waveb(i),i=1,nwave+1)
c
c     IOP's at output depths
c
      write(nurad,fmt='(a)') 'acoef (Total absorption coefficient)'
      write(nurad,104) (acoef(i,0),i=1,nz)
      Do icomp = 1,nconc
      	  write(nurad,fmt='(a,i2,a)') 'acoef (for component ',icomp,')'
      	  write(nurad,104) (acoef(i,icomp),i=1,nz)
      Enddo
      write(nurad,fmt='(a)') 'bcoef (Total scattering coefficient)'
      write(nurad,104) (bcoef(i,0),i=1,nz)
      Do icomp = 1,nconc
      	  write(nurad,fmt='(a,i2,a)') 'bcoef (for component ',icomp,')'
      	  write(nurad,104) (bcoef(i,icomp),i=1,nz)
      Enddo
      write(nurad,fmt='(a)') 'bbcoef (backscattering coefficient)'
      write(nurad,104) (bbcoef(i,0),i=1,nz)
      Do icomp = 1,nconc
      	  write(nurad,fmt='(a,i2,a)') 'bbcoef (for component ',icomp,')'
      	  write(nurad,104) (bbcoef(i,icomp),i=1,nz)
      Enddo
      write(nurad,fmt='(a)') 'atten (beam attenuation coefficient)'
      write(nurad,104) (atten(i),i=1,nz)
      write(nurad,fmt='(a)') 'albedo (albedo of single scattering)'
      write(nurad,104) (albedo(i),i=1,nz)
c
c     AOP's in air and at output depths
c
      write(nurad,fmt='(a)') 'Eou (upwelling scalar irradiance)'
      write(nurad,104) (Eou(i),i=0,nz)
      write(nurad,fmt='(a)') 'Eod (downwelling scalar irradiance)'
      write(nurad,104) (Eod(i),i=0,nz)
!      write(nurad,fmt='(a)') 'Eo (scalar irradiance)'
!      write(nurad,104) (Eod(i)+Eou(i),i=0,nz)
      write(nurad,fmt='(a)') 'Eu (upwelling plane irradiance)'
      write(nurad,104) (Eu(i),i=0,nz)
      write(nurad,fmt='(a)') 'Ed (downwelling plane irradiance)'
      write(nurad,104) (Ed(i),i=0,nz)
      write(nurad,fmt='(a)') 'fMUu (upwelling average cosine)'
      write(nurad,104) (fmuu(i),i=0,nz)
      write(nurad,fmt='(a)') 'fMUd (downwelling average cosine)'
      write(nurad,104) (fmud(i),i=0,nz)
      write(nurad,fmt='(a)') 'fMUtot (total average cosine)'
      write(nurad,104) (fMUtot(i),i=0,nz)
      write(nurad,fmt='(a)') 'R (irradiance reflectance)'
      write(nurad,104) (R(i),i=0,nz)
c
c     radiances:
c
      write(nurad,fmt='(a)') 'RADsky (Lsky = incident sky radiance))'
      write(nurad,104) (RADsky(i),i=1,nmu)
      write(nurad,fmt='(a)') 'RADwla (Lw = water-leaving radiances))'
      write(nurad,104) (RADwla(i),i=1,nmu)
      write(nurad,fmt='(a)') 'RADrsa (Lrefl = surface-reflected sky radi
     1ance))'
      write(nurad,104) (RADrsa(i),i=1,nmu)
      write(nurad,fmt='(a)')'RADupz (Lu(theta,z) = upward radiances in w
     1ater)'
      write(nurad,104) ((RADupz(i,k),i=1,nmu),k=1,nz)
      write(nurad,fmt='(a)')'RADdnz (Ld(theta,z) = downward radiances in
     1 water)'
      write(nurad,104) ((RADdnz(i,k),i=1,nmu),k=1,nz)

      close(nurad)
      return
c
  102 format(15i5)
  104 format(10es15.5e3)
  300 format(/2x,'There has been an unexpected error while attempting',
     1' to open the appending file: unit =',i5,/5x,a)
  301 format(/2x,'A new file will be opened for writing out the ',
     1'remaining data, named: ',/5x,a)
  302 format(/2x,'Hydrolight was unable to append this output file.',
     1/2x,'The cause of this rare error is not enough memory allocated',
     2' to the stack.',/2x,'The fix for this it to increase the',
     3' allocation by editting the automake.fig file ',/2x,
     4'and adding/icreasing the -stack 500000 flag on the linker line.',
     5//2x'To salvage this run, a new file will be opened for output.')

      end
