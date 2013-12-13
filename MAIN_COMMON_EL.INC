c     ALL COMMON BLOCKS MUST BE DEFINED HERE !!!!
c
c     Establish all common blocks used throughout the main code:
c
      COMMON /CRTS/ Rzw(mxmu,mxmu,mxz),Twz(mxmu,mxmu,mxz),
     1              Rzb(mxmu,mxmu,mxz),Sptwz(mxmu,mxz),Smtbz(mxmu,mxz)
      COMMON /CsourceE/ SrcM(mxmu),SrcP(mxmu)
      Common /Csource0/ ibiolum,ichlfl,icdomfl,iraman, ramanEXP
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      COMMON /Crhotau/ rho(mxmu,mxmu),tau(mxmu,mxmu),
     1                 betatP(mxmu,mxmu,mxcomp),betatM(mxmu,mxmu,mxcomp)
      COMMON /CBOTBC/ Rmb(mxmu,mxmu)
      COMMON /Crtsurf/ raw(mxmu,mxmu),taw(mxmu,mxmu),rwa(mxmu,mxmu),
     1                 twa(mxmu,mxmu)
      COMMON /Ciop/ acoef(mxz,0:mxcomp),bcoef(mxz,0:mxcomp),
     1		      atten(mxz),albedo(mxz), bbcoef(mxz,0:mxcomp)
c
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

c       depth index 0 is used to store values in air (z=a)
c       depth 1 is z=0 in water (z=w)  ; E2 is used in Raman calculations
      COMMON /Cirrad/ Eou(0:mxz),Eod(0:mxz),Eu(0:mxz),Ed(0:mxz), 
     1                fMUu(0:mxz),fMUd(0:mxz),fMUtot(0:mxz),R(0:mxz),
     2                E2(0:mxz)

      COMMON /CMISC/  IMISC(30),FMISC(30) 
c    
      COMMON /Cfrstcls/ iabscat,iqasky,iradian,iradanal,iradxcl
      COMMON /Csky/ iskyflag,skyspecs(mxnsky)

      COMMON /Ctitle/ ititle
      Character ititle*120

      Character*120 datadir, digitdir, spreadir,
     1              phasedir, surfdir,bottdir, Pdir
      COMMON /Cdirnames/ datadir,digitdir,spreadir,
     1                   phasedir,surfdir,bottdir, Pdir

      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
C
C******NOTE:  About the array datafiles in common block Cfilenames
C               filenames are stored with the following indices:
C        (0): abPureWater  (2 options; set in setdefaults)
C        (1): ac9 data file
C        (2): Hydroscat data file
C        (3): chl(z) data file
C        (4): CDOM data file
C        (5): bottom reflectance data file
C        (6 to ncomp+5):  [optional] data files for each component 
C
C     The filenames stored in slots (1 thru ncomp+5) are set in initial
C     The filename for the pure water model (0) is set in setdefaults
c
c     for use in PAR:
      character Erootname*120
      common /Efilename/ Erootname
c
      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cpkfcn/ npkfcn,izkfcn(mxz)
c
c     common blocks for communication with the Gregg & Carder sky
c     irradiance model (routine GCEd, RADTRANX):
      parameter (nlt=701)
      common /cgcirr1/ iein,iblw,jday,rlon,rlat,the,GMThr,pres,am,rh,
     1                 wv,wsm,ws,vi,ro3
      common /cgcirr2/ Eddirgc(nlt),Eddifgc(nlt),Edtotgc(nlt)
          character*32 gcifile
      COMMON/gcifname/gcifile
c
c     Arrays for 1-D spline fitting and interpolation routines
      parameter(mxkx=4)  ! kx = 4 for cubic polynomials
      Integer kx
      COMMON /CEospl/ nspl,zspl(mxz),Eospl(mxz,mxwave),
     1                E2spl(mxz,mxwave),Edspl(mxz,mxwave)
      COMMON /CEospl2/kx, inbv(mxwave),  Tx(mxz+mxkx, mxwave),
     1			 Bsplin(mxz, mxwave),
     2                inbv2(mxwave), Tx2(mxz+mxkx, mxwave),
     3			 Bsplin2(mxz, mxwave)
c
c     Cxcl holds arrays accumulated in routine storexcl
      COMMON /Cxcl/  nzxcl,zxcl(mxz),izxcl(mxz),
     1               axcl(mxz,mxwave,0:mxcomp),
     2			bxcl(mxz,mxwave,0:mxcomp),
     3         	bbxcl(mxz,mxwave,0:mxcomp),
     4               Edxcl(0:mxz,mxwave),Euxcl(0:mxz,mxwave),
     5               Eodxcl(0:mxz,mxwave),Eouxcl(0:mxz,mxwave),
     5               Rxcl(0:mxz,mxwave),
     6               Raduxcl(0:mxz,mxwave),Radwxcl(mxwave)
c     CKxcl holds arrays accumulated in routine Kfcn
      COMMON /CKxcl/ nzKxcl,zKxcl(mxz),fKdxcl(mxz,mxwave),
     1  fKuxcl(mxz,mxwave),fKoxcl(mxz,mxwave),
     2  fKnetxcl(mxz,mxwave),fKLuxcl(mxz,mxwave)
c
c     Common blocks that are new in H4.1
      integer indexchl
      common /CindexChl/ indexchl

c     Common blocks new to Ecolight 3.0
      COMMON /Cbbquadchk/ m,iqpart     !shared w/ inishrad & selpfbb
c
c     Common blocks associated with depth optimization
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)
!!!
      character*120 Lrootname, S0datafile
      common /Cradfile/ Lrootname
      Common /CSourceBdata/  S0datafile    !read in initial, used in s0bdata
c
      integer ibbopt, itype
      real bbfrac, BfrefPL, Bf0PL, BfmPL
      Integer iastropt(mxcomp+2), ibstropt(mxcomp)
      Real astarRef(mxcomp+2), astar0(mxcomp+2)
      Real asgamma(mxcomp+2)
      Real bstarRef(mxcomp), bstar0(mxcomp)
      Real CompN(mxcomp), PLm(mxcomp)
      Real GAMm(mxcomp),GAMi(mxcomp)
      Character*120 astarfile(mxcomp+2), bstarfile(mxcomp+2)
      common /Cconc/  itype(mxcomp)
      common /Cbbopt/ ibbopt(mxcomp), bbfrac(mxcomp),
     1                BfrefPL(mxcomp), Bf0PL(mxcomp), BfmPL(mxcomp)
      Common /Castar/ iastropt, astarRef, astar0, asgamma
      Common /Castar2/ astarfile
      Common /Cbstar/ ibstropt, bstarRef, bstar0, CompN, Plm,
     1                GAMm, GAMi
      Common /Cbstar2/ bstarfile
      Common /Cabac9n/ nac9files
      INTEGER iparam
      COMMON /CmodelOpts/ iparam(0:5)
      REAL compconc
      COMMON /CConstConc/ compconc(mxcomp)
c     flag for calculating bb shared with ab routine (calc bb iff flag=1)
      Common /CbbCalc/ ibbCalc
c
C     OUPUT variables shared with RangeWarn via COMMON block
      integer nkall, ierr(2*mxcomp+8, 4), ibc(2*mxcomp+8)
      real zdat(2*mxcomp+8, 2), wdat(2*mxcomp+8, 2)
      character*120 dfwarn(2*mxcomp+8)
      COMMON /CRANGEW/ nkall, ierr, zdat, wdat, dfwarn, ibc
c
C     The following common block shared between excelM and qasky.f
      common /CErik/ EdifOut(mxwave), EdirOut(mxwave)
      Common /CPAR/fPar(0:mxz),Edpar(0:mxz),fKpar(mxz),
     1             Eoquant(0:mxz,mxwave)   !shared with PAR.f
c
c     Arrays for 1-D spline fitting and interpolation ZTOZETA
      Integer kxz, inbvz
      Real Txz,Bsplinz
      common /Czetazspl/ nzvalsz, kxz, inbvz,
     1                   Txz(mxnzvals+mxkx),Bsplinz(mxnzvals)
c
c     Bioluminescent Source printout array
      real S0pnt(mxwave,mxz)
      Common /CS0xcl/ S0pnt


