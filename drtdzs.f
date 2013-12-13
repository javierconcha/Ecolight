C     Last change:  LKS   9 Dec 2011   10:44 am
      SUBROUTINE drtdzs(zetanow,RT,deriv)
c
c     core routine on file drtdzs.f 
c
c     called by ODE [MAIN->RADIANCE->RICCATI->ODE->DRTDZS]
c
c     calls rhotau, sumsrct (in sources.f), and eval1D (in dataintp.f)
c
c     This subroutine evaluates deriv = d(RT)/d(zeta) at zeta = zetanow
c     (the right hand side of 8.74-8.85) for use by the canned routine
c     (ode) that solves the Riccati equation system.  Internal source
c     terms are included if isource .ne. 0, and omitted if isource = 0.
c
c     For the downward integration sweep, arrays Rzw, etc are stored
c     in the linear array RT as follows (for a given zeta value): 
C 
c     Rzw(I,J)  IS RT(I + (J-1)*nmu) 
c     Twz(I,J)  IS RT(I + (J-1)*nmu + nmu*nmu)
c     Sptwz(i) is RT(i + 2*nmu*nmu)
C
      INCLUDE "DIMENS_XL.INC"
      parameter(mxeqn=2*mxmu*mxmu + mxmu)
c
c     lookup table for z to zeta grid
      Common /Cztozeta/ nzvals, zetavals(mxnzvals),zvals(mxnzvals)
      integer nzvals
c
      dimension RT(mxeqn),deriv(mxeqn)

      COMMON /Crhotau/ rho(mxmu,mxmu),tau(mxmu,mxmu),
     1                 betatP(mxmu,mxmu,mxcomp),betatM(mxmu,mxmu,mxcomp)
      COMMON /CsourceE/ SrcM(mxmu),SrcP(mxmu)
      Common /Csource0/ ibiolum,ichlfl,icdomfl,iraman, RamanEXP
      COMMON/CMISC/ imisc(30),fmisc(30)
c 
c     temporary local storage:
      dimension work(mxmu,mxmu)
c
      External yinterp 
c
!     iznow is the last known index into the z to zeta array
      integer iznow
      data iznow/1/
      save iznow
!*********************************************************************
c
      nmu = imisc(1)
      iop = imisc(5)
      isource = imisc(8)
      isweep = imisc(13) 
      nsq = nmu*nmu
c
      if(iop.ne.1) then
c        If the ab routine is to be called with geometric depth, use the
c        splines determined in ztozeta to compute the geometric depth znow
c        corresponding to the current optical depth zetanow.  (znow is
c        not needed if the ab routine is to be called with optical depth.)
c        (depth calls are made in routines rhotau and sumsrc
c
c       call to linear interp routine that returns 
c       index and linear interp coef for znow
        znow = yinterp(iznow, nzvals, zetanow, zetavals, zvals)
        fmisc(18) = znow
      end if
c
c----------------------------------------------------------------
c 
c     Determine rho and tau at the current zeta value
c     NOTE:  rhotau calls abscat; the total a and b values returned
c     from that call are then used in routine sumsrc.
C
      call rhotau(zetanow,znow)
c
c     Determine the source terms at the current zeta value
      if(isource.ne.0) call sumsrc(zetanow,znow)
c
c---------------------------------------------------------------
c
      if(isweep.eq.1) then
c     Compute dRT/dzeta as defined for the downward integration sweep
C 
c     Compute work = tau + rho*Rzw (used in 8.74, 8.75 and 8.78)
      DO I=1,nmu
      DO J=1,nmu
      work(I,J) = tau(I,J)
         DO K=1,nmu
         work(I,J) = work(I,J) + rho(I,K)*RT(K + (J-1)*nmu)
         end do
      end do
      end do
c
c     Compute d(Rzw)/dzeta BY EQ. 8.74 
C 
      DO I=1,nmu
      DO J=1,nmu
      temp1 = 0.
      temp2 = 0.
         DO K=1,nmu
         temp1 = temp1 + RT(I + (K-1)*nmu)*work(K,J)
         temp2 = temp2 + tau(I,K)*RT(K + (J-1)*nmu)
         end do
      deriv(I + (J-1)*nmu) = rho(I,J) + temp1 + temp2
      end do
      end do
C 
c     Compute d(Twz)/dzeta BY EQ. 8.75 
C 
      DO I=1,nmu
      DO J=1,nmu
      temp1 = 0.
         DO K=1,nmu
         temp1 = temp1 + RT(I + (K-1)*nmu + nsq)*work(K,J)
         end do
      deriv(I + (J-1)*nmu + nsq) = temp1
      end do
      end do
C
      if(isource.ne.0) then
c     include source terms:
c     Compute d(Sptwz)/dzeta by EQ. 8.78
C
      DO J=1,nmu
         temp1 = 0.
         temp2 = 0.
            DO K=1,nmu
               temp1 = temp1 + RT(k + 2*nsq)*work(K,J)
               temp2 = temp2 + SrcM(K)*RT(K + (J-1)*nmu)
            end do
         deriv(j + 2*nsq) = SrcP(J) + temp1 + temp2
      end do

      endif  ! end isource test
C
      endif  ! end isweep = 1
C
c--------------------------------------------------------------------
c
      if(isweep.eq.2) then
c
c     upward integration sweep to solve Eqs. 8.80 and 8.84
c
c     For the upward sweeps, arrays Rzb and Smtbz are stored in
c     the linear array RT as follows (for a given zeta value):
c
c     Rzb(I,J) is RT(I + (J-1)*nmu)
c     Smtbz(i) is RT(i + nmu*nmu)
C
c     Compute dRT/dzeta as defined for the upward integration sweep
C 
c     Compute work = tau + rho*Rzb (used in 8.80 and 8.84)
      DO I=1,nmu
      DO J=1,nmu
         work(I,J) = tau(I,J)
         DO K=1,nmu
            work(I,J) = work(I,J) + rho(I,K)*RT(K + (J-1)*nmu)
         end do
      end do
      end do
c
c     Compute d(Rzb)/dzeta BY EQ. 8.80 
C 
      DO I=1,nmu
      DO J=1,nmu
      temp1 = 0.
      temp2 = 0.
         DO K=1,nmu
         temp1 = temp1 + RT(I + (K-1)*nmu)*work(K,J)
         temp2 = temp2 + tau(I,K)*RT(K + (J-1)*nmu)
         end do
      deriv(I + (J-1)*nmu) = -rho(I,J) - temp1 - temp2
      end do
      end do
c
      if(isource.ne.0) then
c     Compute d(Smtbz)/dzeta by 8.84
c
      do j=1,nmu
      temp1 = 0.
      temp2 = 0.
         do k=1,nmu
         temp1 = temp1 + RT(k+nsq)*work(k,j)
         temp2 = temp2 + SrcP(k)*RT(k + (j-1)*nmu)
         end do
      deriv(j + nsq) = -SrcM(j) - temp1 - temp2
      end do
      endif  ! end isource test
c
      endif  ! end isweep = 2
c
      RETURN
      END 

