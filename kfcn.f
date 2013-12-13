C     Last change:  CDM 05 Dec 2008
      subroutine kfcn
c 
c     core routine on file kfcn.f   
c
c     called by RADANAL [MAIN->RADANAL->KFCN]
c 
c     This routine computes the K-functions associated with the scalar
c     and plane irradiances.  The functions are computed as rates 
c     of change with respect to geometric depth. 
c 
c+++++WARNING:  Any pair of depths z(iz) and z(iz+1) can be used to 
c               compute layer-averaged K's.  These layer-averaged K's will 
c               be good estimates of the local K's at the layer
c               midpoint only if the depths are very closely spaced. 

c     Changed for HE5: use finite differences of ln(E), rather than E,
c                      as shown in the HE5 UG section 5.7
c 
      INCLUDE "DIMENS_XL.INC"
c 
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cpkfcn/ ipkfcn,izkfcn(mxz) 
      COMMON /Cirrad/ Eou(0:mxz),Eod(0:mxz),Eu(0:mxz),Ed(0:mxz), 
     1                fMUu(0:mxz),fMUd(0:mxz),fMUtot(0:mxz),R(0:mxz),
     2                E2(0:mxz)
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
      COMMON /Cmisc/ imisc(30),fmisc(30)
      COMMON /CKxcl/ nzKxcl,zKxcl(mxz),fKdxcl(mxz,mxwave),
     1  fKuxcl(mxz,mxwave),fKoxcl(mxz,mxwave),
     2  fKnetxcl(mxz,mxwave),fKLuxcl(mxz,mxwave)
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)
c
c     declare temp vars
      integer nmu, iwvl, jwave, izc
      real wavelen
c
      nmu = imisc(1)
      iwvl = imisc(11)
      wavelen = fmisc(13)

c     printout down to the the max depth where quantities were
c     computed by solving the RTE.  This depth was determined in varzmax.
      jwave = imisc(11)
      izc = min(indexz(jwave),ipkfcn+1)
c 
c     geometric-depth K functions
c
      if(imisc(9).ge.0) then
		write(10,400) wavelen
	else
		write(10,500) wavelen
	endif
c
      do iiz=1,izc-1
        iz = izkfcn(iiz)
        c = -1.0/(z(iz+1) - z(iz)) 
        zmid = 0.5*(z(iz+1) + z(iz))
      
        fKou = c*alog(Eou(iz+1)/Eou(iz)) 
        fKod = c*alog(Eod(iz+1)/Eod(iz)) 
        fKo = c*alog((Eou(iz+1)+Eod(iz+1))/(Eou(iz) + Eod(iz)))
        fKu = c*alog(Eu(iz+1)/Eu(iz)) 
        fKd = c*alog(Ed(iz+1)/Ed(iz))
        fKnet = c*alog((Ed(iz+1) - Eu(iz+1))/(Ed(iz) - Eu(iz)))
        fKLu = c*alog(radupz(nmu,iz+1) / radupz(nmu,iz))     
c
      if(imisc(9).ge.0) then
	  write(10,302) z(iz),z(iz+1),zmid,fKou,fKod,fKo,fKu,fKd,
     1                   fKnet,fKLu
	else
	  write(10,302) z(iz),z(iz+1),zmid,fKo,fKu,fKd,fKnet,fKLu
	endif
c
c     save selected K-functions for the Excel spreadsheet routines
c
      zKxcl(iiz) = zmid
      fKdxcl(iiz,iwvl) = fKd
      fKuxcl(iiz,iwvl) = fKu
      fKoxcl(iiz,iwvl) = fKo
      fKnetxcl(iiz,iwvl) = fKnet
      fKLuxcl(iiz,iwvl) = fKLu
c
      end do
c
      return
c 
  400 FORMAT(///2x,'LAYER-AVERAGE K-functions (units of 1/meter)',
     1' at ',f6.1,' nm'//
     2'    zupper    zlower',6X,'zmid',4X,
     3'Kou(z)    Kod(z)    Ko(z)     Ku(z)     Kd(z)    Knet(z)',
     4'    KLu(z)'/)
  500 FORMAT(//2x,'LAYER-AVERAGE K-functions (units of 1/meter)',
     1' at ',f6.1,' nm'//
     2'    zupper    zlower',6X,'zmid',4X,
     3'Ko(z)     Ku(z)     Kd(z)    Knet(z)    KLu(z)'/)
  302 FORMAT(3F10.3,7F10.5) 
      END 
