C     Last change:  LKS   1 Nov 2007    1:01 pm
        subroutine FillWave
c       This subroutine fills in the scalar irradiances for printout and for use
c       in calculating PAR when only alternating wavebands were solved
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave)
      COMMON /CEospl/ nspl,zspl(mxz),Eospl(mxz,mxwave),
     1                E2spl(mxz,mxwave),Edspl(mxz,mxwave)
      COMMON /CMISC/ IMISC(30),FMISC(30)

      nwave = imisc(7)
      nwskip = imisc(26)

      do iwave=2, nwave-1, nwskip
        do izz=1, npirad
          iz = izirad(izz)
          Eospl(iz,iwave) = 0.5*(Eospl(iz,iwave-1) + Eospl(iz,iwave+1))
        enddo
      enddo
c
      end
