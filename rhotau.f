C     Last change:  LKS   1 Nov 2007    1:05 pm
      subroutine rhotau(zetanow,znow)
c 
c     core routine on file rhotau.f
c 
C     called by INFBOTM  [MAIN->RADIANCE->BOTMBC->INFBOTM->RHOTAU]
C     called by DRTDZS   [MAIN->RADIANCE->RICCATI->ODE->DRTDZS->RHOTAU]
c
c     This routine calls the "ab" routine
c
c     This routine computes the local reflectances and
c     transmittances rho and tau (for the current depth zetanow)
c     using Eq. 8.42 with the quad (band) averaged betatilde+/-
c     = betatP/M, which were computed in discpf
c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /Crhotau/ rho(mxmu,mxmu),tau(mxmu,mxmu),
     1                 betatP(mxmu,mxmu,mxcomp),betatM(mxmu,mxmu,mxcomp)
      COMMON /Cmisc/ imisc(30),fmisc(30) 
c
c     temporary local storage:
      dimension bcomp(mxcomp), acomp(mxcomp)
      dimension totlpp(mxmu,mxmu),totlpm(mxmu,mxmu)
C
      nmu = imisc(1)
      iop = imisc(5)
      ncomp = imisc(6)
      wavelen = fmisc(13)
c
c-------------------------------------------------------------------
c
c     Get the total phase functions betatilde+ = totlpp and
c     betatilde- = totlpm, and the albedo, at optical depth zetanow
c
c     Get the component scattering coefficients
c
c     Call the abscat routine with either geometric (the usual case)
c     or optical depth, as is appropriate for the run.
c     NOTE:  the total a and b values computed in this abscat call
c     will be used in the Source routines (called by drtdzs)
c
      depth = znow
      if(iop.eq.1) depth = zetanow
c
c     insert the call to the desired "ab" routine:
      call abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
c
      fmisc(19) = atotal
      fmisc(20) = btotal
c
      albedo = btotal/(atotal + btotal)
c
c     Sum the component phase functions, weighted by the scattering
c     coefficients as in Eq. (3.13), to get the total quad-averaged
c     phase function at depth zetanow for use in Eq. (8.34)
c
      do j=1,nmu
      do i=1,nmu
      temppp = 0.
      temppm = 0.
        if(btotal.ge.1.0e-6) then
          do k=1,ncomp
            temppp = temppp + (bcomp(k)/btotal)*betatP(i,j,k)
            temppm = temppm + (bcomp(k)/btotal)*betatM(i,j,k)
          end do
        end if
      totlpp(i,j) = temppp
      totlpm(i,j) = temppm
      end do
      end do
c
c-----------------------------------------------------------------
c
c     Compute rho and tau at this depth
c
      do ir=1,nmu
       do iu = 1,nmu
         rho(ir,iu) = albedo*totlpm(ir,iu)/fmu(iu)
         if(ir.eq.iu) then
            delt = 1.
         else
            delt = 0.
         endif
         tau(ir,iu) = (albedo*totlpp(ir,iu) - delt)/fmu(iu)
       end do   ! iu loop
      end do   ! ir loop
c
      RETURN
      END 

