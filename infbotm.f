C     Last change:  LKS   7 Aug 2008    6:26 pm
      SUBROUTINE infbotm
C 
C     core routine on infbotm.f
C
C     called by BOTMBC  [MAIN->RADIANCE->BOTMBC->INFBOTM]
C
C     calls routines RHOTAU and MATxMAT, and canned routines:
C     sgeev, spsort, and sgefs in files eigenvv.f, sort.f, and matinv.f
C 
c     This routine sets up the bottom boundary condition for the case of
c     infinitely deep water.  Asymptotic values of various quantities are
c     also found.
c
C     The routine first sets up and solves the eigenmatrix problem 
c     K*E = E*kappa (eq. 9.50) as described in Light & Water Section 9.3
C
C     The submatrices EP = E(+) and EM = E(-) are extracted, and
C     R(infinity,L) = -E(-) * E(+)inverse (eq. 9.76) is computed
c     for use as the bottom boundary reflectance for infinitely deep 
c     water.
C 
C     The asymptotic radiance distribution is computed using eq (9.82), 
c     and associated quantities are also found using eq (9.86)
c
      INCLUDE "DIMENS_XL.INC"
      PARAMETER(mxmu2=2*mxmu)
C 
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Crhotau/ rho(mxmu,mxmu),tau(mxmu,mxmu),
     1                 betatP(mxmu,mxmu,mxcomp),betatM(mxmu,mxmu,mxcomp)
      COMMON/CBOTBC/ rmb(mxmu,mxmu)
      COMMON/CMISC/ IMISC(30),fmisc(30) 
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)

c     temporary local storage:
      DIMENSION ip(mxmu2),evals(mxmu2),eigv(mxmu),esort(mxmu2) 
      dimension ak(mxmu2,mxmu2),emk(mxmu2,mxmu2),fullf(mxmu2,mxmu2)
      dimension em(mxmu,mxmu),ep(mxmu,mxmu),epnv(mxmu,mxmu)
      dimension iperm(mxmu2),ip2(mxmu2),etemp(mxmu2)
      dimension thet(mxmu),rminf(mxmu),rpinf(mxmu)
c
c     arrays for use with LAPACK routines sgeev and sgefs:
      parameter(mxwork=1000)
      dimension evimag(mxmu2),evecr(mxmu2,mxmu2),work(mxwork)
      dimension iwork(mxmu2)
c
      nmu = IMISC(1)
c      nz = IMISC(4)
      jwave = imisc(11)
      radeg = fmisc(3)
c
c     when depth optimization is used, apply the bottom boundary
c     at the wavelength-dependent max computation depth determined in varzmax

      nz = indexz(jwave)
c
      call rhotau(zeta(nz),z(nz))
c 
      m = nmu
      m2 = 2*m
c 
      acoef = fmisc(11)
      atten = fmisc(12)
c
C     Initialize the K matrix, using (9.33)
C
      DO I=1,M2
        DO J=1,M2
          IF(I.LE.M) THEN 
            IF(J.LE.M) ak(I,J) = -tau(I,J)
            IF(M.LT.J) ak(I,J) = rho(I,J-M)
          ELSE
            IF(J.LE.M) ak(I,J) = -rho(I-M,J)
            IF(M.LT.J) ak(I,J) = tau(I-M,J-M)
          ENDIF 
          temp = temp +abs(ak(I,J))
        end do
      end do
C 
C     Find the eigenvalues and eigenvectors of K
C
cccccccccccc  Canned routine calls  ccccccccccccccccccccccccccc
c
c     Use public-code routine SGEEV from LAPACK (Linear Algebra PACKage)
c     to find the eigenvalues and eigenvectors of a real, general matrix.
c     SGEEV and related routines are contained on files EIGENVV.f,
c     BLAS.f, and SLACOM.f
c
      call sgeev('N','V',m2,ak,mxmu2,evals,evimag,vleft,1,evecr,
     1  mxmu2,work,mxwork,info)
c 
      if(work(1).gt.mxwork) print*,' sub infbotm:  increase mxwork to ',
     2work(1),' for maximum efficiency in sgeev'     
c
c**********************************************************************
C
C     Sort the eigenvalues in evals (smallest to largest) and generate
c     a corresponding permutation array, iperm, for use in
c     sorting the eigenvectors
c
      kflag = 1
      call spsort(evals,m2,iperm,kflag,ier)
c
c     Sort evals according the the index array iperm
      Do ii=1,m2
         etemp(ii) = evals(iperm(ii))
      end do
      do ii=1,m2
         esort(ii) = etemp(ii)
      end do
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     re-order: positive (smallest to largest), then negative (smallest
c     magnitude to largest magnitude); same reorder for permutations
c
      do i=1,m
         ip2(i) = iperm(i+m)
         ip2(i+m) = iperm(m+1-i)
         evals(i) = esort(i+m)
         eigv(i) = evals(i)
         evals(i+m) = esort(m+1-i)
      end do
c
c     ip2 shows how to reorder eigenvalues to get the original order.
c     We need to reorder the original eigenvectors to match the rearranged
c     eigenvalues, so exchange permutation index and permutation value:
c
      do i=1,m2
         ip(ip2(i)) = i
      end do
C 
C     Define real, ordered eigenvector matrix EMK 
C 
      DO J=1,M2
        JJ = IP(J)
          DO I=1,M2
c           for use with LAPACK:
            emk(i,jj) = evecr(i,j)
          end do
      end do
c
C     Extract the submatrices EP = E(+) and EM = -E(-) via (9.65)
C 
      DO I=1,M
        DO J=1,M
         em(i,j) = -emk(i+m,j)
         ep(i,j) = emk(i,j)
c         eptemp(i,j) = ep(i,j)
        end do
      end do
C 
C     Invert E(+) and compute R(infinity), using (9.76)
C
cccccccccc  Canned routine calls for matrix inversion  ccccccccccccc
c
      DO ii = 1,m
         DO jj = 1,m
            epnv(ii,jj) = 0.
         end do
         epnv(ii,ii) = 1.
      end do
c
      itask = 1
      call sgefs(ep,mxmu,m,epnv(1,1),itask,ind,work,iwork)
      if(ind.lt.0) 
     1   write(10,*)'Error on first call to SGEFS from infbotm',IND
      itask = 2
      do j = 2,m
        call sgefs(ep,mxmu,m,epnv(1,j),itask,ind,work,iwork)
        if(ind.lt.0) 
     1   write(10,*)'Error on second call to SGEFS from infbotm',IND
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call matxmat(em,epnv,m,m,m,mxmu,mxmu,rmb,mxmu)
c
c     invert the full eigenvector matrix E = emk to get F = fullf
c
c     correct possible sign errors in paired eigenvectors
      do jj=m+1,m2
        do ii=1,m
          emk(ii,jj) = emk(ii+m,jj-m)
          emk(ii+m,jj) = emk(ii,jj-m)
        end do
      end do
c
cccccccccc  Canned routine calls  ccccccccccccccccccccccccccccccc
c
      DO ii = 1,m2
         DO jj = 1,m2
            fullf(ii,jj) = 0.
         end do
         fullf(ii,ii) = 1.
      end do
c
      itask = 1
      call sgefs(emk,mxmu2,m2,fullf(1,1),itask,ind,work,iwork)
      if(ind.lt.0) 
     1   write(10,*)'Error on third call to SGEFS from infbotm',IND
c
      itask = 2
      do j = 2,m2
        call sgefs(emk,mxmu2,m2,fullf(1,j),itask,ind,work,iwork)
        if(ind.lt.0) 
     1   write(10,*)'Error on fourth call to SGEFS from infbotm',IND
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Get the asymptotic radiance distribution from the first row of F.
c     Note that that Rad(+) [the downward irradiances] are obtained 
c     from F(+) (= cols 1 to nmu of F) and that Rad(-) [the upward
c     radiances] are obtained from F(-) (= cols nmu+1 to 2*nmu of F).
c 
c     normalize to the nadir radiance
c
      anorm = 1.0/fullf(1,nmu)
      if(imisc(9).ge.0) write(10,655)
c
c           get exact theta values at cell centers for printout
      do i=1,nmu
         if(i.eq.1) then
            thet(i) = 0.5*(90.0 + radeg*acos(bndmu(1)))
         elseif(i.eq.nmu) then
            thet(i) = 0.
         else
            thet(i) = 0.5*radeg*(acos(bndmu(i-1)) + acos(bndmu(i)))
         endif
      end do
c
      do i=nmu,1,-1
         rminf(i) = anorm*fullf(1,i+nmu)
         rpinf(i) = anorm*fullf(1,i)
         if(imisc(9).ge.0) write(10,657) i,thet(i),rpinf(i)
      end do
      do i=1,nmu
         if(imisc(9).ge.0) write(10,657) i,180.-thet(i),rminf(i)
      end do
c
C     Use the asymptotic radiance distribution to get other asymptotic
C     quantities 
C 
C     accumulate irradiance sums (Ed(inf), etc.)
      Eou = 0.
      Eod = 0.
      Eu = 0.
      Ed = 0.
      DO I=1,nmu
         dmu = deltmu(I) 
         Eou = Eou + rminf(I)*dmu
         Eod = Eod + rpinf(I)*dmu
         Eu = Eu + rminf(I)*fmu(I)*dmu 
         Ed = Ed + rpinf(I)*fmu(I)*dmu
      end do
c
c     mean cosines
      fmuuinf = Eu/Eou 
      fmudinf = Ed/Eod 
      fmuinf = (Ed - Eu)/(Eod + Eou)
c
C     kinf by (9.83)
      fkinf = atten*eigv(1) 
c
C     Rinf by the asymptotic limit of (5.67)
      rinf = (fkinf - acoef/fmudinf)/(fkinf + acoef/fmuuinf)
      rinf = Eu/Ed
c
c     absorption check by asymptotic limit of (5.69)
      absinf = fkinf*fmuinf
c
      if(imisc(9).ge.0) write(10,672) fmuuinf,fmudinf,fmuinf,rinf,
     1					fkinf,absinf
c
      RETURN
C 
C     FORMATS 
C 
  655 format(//2x,'The shape of the asymptotic radiance distribution',
     1' for source-free water is given by'//
     2'    u    theta     Linf(theta)')
  657 FORMAT(I5,f9.2,1P,E15.4)
  672 FORMAT(//2x,'Asymptotic source-free water AOP values are'/
     1'        MUu(inf) =',F7.4/
     2'        MUd(inf) =',F7.4/
     3'         MU(inf) =',F7.4/
     4'        R(m,inf) =',f7.4/
     5'          k(inf) =',f7.4,' (1/m)'/
     6'   absorp(check) =',f7.4,' = k(inf)*MU(inf)')  
c 
      END 

