C     Last change:  LKS   1 Nov 2007    8:18 pm
      SUBROUTINE riccati
c 
c     core routine on file riccati.f
c
c     called by RADIANCE [MAIN->RADIANCE->RICCATI]
c
c     calls ODE (and passes DRTDZS to ODE)
c 
c     This version uses the SODE package to solve the ODE system.
c     Modified 2/26/03 by lks to use extra integration intervals to 
c     assure good solutions of the ODE system.
c
c
c     This routine solves for the standard operators Rzw = R(zeta,w),
c     Twz = T(w,zeta), Sptwz = S+t(w,zeta), etc by integrating
c     Eqs (8.74), (8.75), and (8.78) in a downward sweep with
c     initial values of R(w,w) = 0, T(w,w) = 0, S+t(w,w) = 0, etc
c     as given by Eqs. (8.72).
c
c     Values of Rzb = R(zeta,b), Smtbz = S-t(b,zeta), etc.
c     are obtained by integrating Eqs. (8.80) and (8.84) in an upward
c     sweep from zeta = m to zeta = w, with values for the lower
c     boundary S[m,b] built into the initial conditions as shown in Eq.
c     (8.94).
c
c     If isource.ne.0, the Riccati integration includes
c     internal source terms (the sum of true sources and effective 
c     sources from inelastic scatter).  The source terms and the 
c     associated equations are omitted if isource = 0, in order to 
c     decrease the run time when no sources are present.
c
c     The lower boundary is assumed to be a reflecting surface at
c     depth zetz = m, with no radiance coming up through the lower
c     boundary.  This surface can described either an opaque bottom
c     or an infinitely deep layer of water below depth zeta = m.
c     The lower boundary S[m,b] is also assumed to be source free.
c     These lower boundary assumptions simplify the integration
c     scheme outlined in Fig 8.2.
c
c     For the downward integration sweep, arrays Rzw, etc are stored
c     in the linear array RT as follows (for a given zeta value): 
c 
c     Rzw(I,J)  IS RT(I + (J-1)*nmu) 
c     Twz(I,J)  IS RT(I + (J-1)*nmu + nmu*nmu)
c     Sptwz(i) is RT(i + 2*nmu*nmu)
c
      INCLUDE "DIMENS_XL.INC"
      PARAMETER(mxeqn = 2*mxmu*mxmu + mxmu)
c
      COMMON /CRTS/ Rzw(mxmu,mxmu,mxz),Twz(mxmu,mxmu,mxz),
     1              Rzb(mxmu,mxmu,mxz),Sptwz(mxmu,mxz),Smtbz(mxmu,mxz)
      COMMON /CBOTBC/ Rmb(mxmu,mxmu)
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)
      COMMON /CMISC/ imisc(30),fmisc(30)
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)
c
      DIMENSION RT(mxeqn)
c     local storage for SODE:
      DIMENSION work(100+21*mxeqn),iwork(5)
C 
C     Subroutine drtdzs evaluates the rhs of Eqs. (8.74), (8.75), etc. 
c     for use by the ODE solver.
c 
      EXTERNAL drtdzs
C
      nmu = imisc(1)
      isource = imisc(8)
      jwave = imisc(11)
      nz = indexz(jwave)
      relerrs = fmisc(6)
      abserrs = fmisc(7) 
      nmu2 = nmu*nmu 
C 
c     ------------------------------------------------------------
c
c     Begin integration of (8.74), (8.75) and (8.78) in a downward
c     sweep from zeta = w to zeta = m.
c
C     Initialize the arrays at zeta = w using (8.72) 
C 
      if(isource.eq.0) then
         neqns = 2*nmu2
      else
         neqns = 2*nmu2 + 2*nmu
      endif
c
      DO j=1,nmu
         Sptwz(j,1) = 0.       ! will not be used if isource = 0
         RT(j + 2*nmu2) = 0.
            DO i=1,nmu
               Rzw(I,J,1) = 0.
               RT(I+(J-1)*nmu) = 0.
               delt = 0.
               IF(I.EQ.J) delt = 1.
               Twz(I,J,1) = delt
               RT(I+(J-1)*nmu+nmu2) = delt
           end do
      end do
C 
      zetastrt = zeta(1) 
c     set flag to "downward sweep" (= 1) for drtdzs:
      imisc(13) = 1
c
C     Integrate (8.74), etc to find R(zeta,w), etc at each zeta level 
C 
      DO iz=2,nz
c
ccccccccccc  Call canned routine to solve the ODE system  ccccc
c
      relerr = relerrs
      abserr = abserrs

clks  add extra evaluations iff zetas are too widely spaced
c
      dzeta = zeta(iz) - zeta(iz - 1) 
      if( dzeta  .lt. fmisc(26) ) then 
         nsteps = 1
      else
         nsteps = int(dzeta/fmisc(26)) +  1
         dzeta = dzeta / nsteps
      endif
      
      Do istep = 1, nsteps
c       special case if last step; make sure zetaend is our output depth
        If(istep.eq.nsteps) then
          zetaend = zeta(iz) 
        Else
          zetaend = zetastrt + dzeta
        Endif
       
c     iflag is set to -1 to prevent ode (a predictor-corrector
c     algorithm) from extrapolating past the integration interval 
c     (e.g. to below the bottom on the downward sweep or to above 
c     the surface on the upward sweeps), which might cause problems 
c     in the calls to the abscat routine (even though the over-
c     extrapolation would be corrected to end up at the proper 
c     point).
c
        iflag = -1

        call ode(drtdzs,neqns,RT,zetastrt,zetaend,relerr,abserr,
     1           iflag,work,iwork)

      Enddo		!isteps
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C 
C     Save the solution at zeta = zetaend 
C 
         DO J=1,nmu
         Sptwz(j,iz) = RT(j + 2*nmu2)
            DO I=1,nmu
               Rzw(I,J,iz) = RT(I+(J-1)*nmu)
               Twz(I,J,iz) = RT(I + (J-1)*nmu + nmu2)
            end do   ! i loop
         end do      ! j loop

      end do         ! iz loop

c
c     ============= End of downward integration sweep ==================
c
c     Begin integration of (8.80) and (8.84) in an upward integration
c     sweeps from zeta = m to zeta = w.
c     The bottom boundary conditions are incorporated into the
c     initial values at zeta = m. 
c     Subroutine BOTMBC already has computed the needed value of
c     Rmb = R(m,b) for the current bottom type.
c
c     For the upward sweep, arrays Rzb and Smtbz are stored in
c     the linear array RT as follows (for a given zeta value):
c
C     Rzb(I,J) IS RT(I + (J-1)*nmu)
c     Smtbz(i) is RT(i + nmu*nmu)
c
C     Integrate 8.80 from m to w to find R(zeta,b) and integrate
c     8.84 to get S-t(b,zeta), at each zeta level
C 
C     Initialize at zeta = m with Rmb = R(m,b), using 8.94
C 
      DO J=1,nmu
      Smtbz(j,nz) = 0.
      RT(j + nmu2) = 0.
         DO I=1,nmu
            Rzb(I,J,nz) = Rmb(I,J)
            RT(I+(J-1)*nmu) = Rmb(I,J)
         end do
      end do
C 
      zetastrt = zeta(nz)
c
c     set flag to "upward sweep" for drtdzs:
      imisc(13) = 2
c
C     Integrate 
c
      DO iz=1,nz-1

      izrev = nz-iz
c
ccccccccccc  Call canned routine to solve the ODE system  ccccc
c
      relerr = relerrs
      abserr = abserrs

clks  add extra evaluations iff zetas are too widely spaced
c
      dzeta = zeta(izrev+1) - zeta(izrev)
      if( dzeta  .lt. fmisc(26) ) then 
         nsteps = 1
      else
         nsteps = int(dzeta/fmisc(26)) +  1
         dzeta = dzeta / nsteps
      endif
      
      Do istep = 1, nsteps
c       special case if last step; make sure zetaend is our output depth
        If(istep.eq.nsteps) then
          zetaend = zeta(izrev) 
        Else
          zetaend = zetastrt - dzeta
        Endif
        
c     iflag is set to -1 to prevent ode (a predictor-corrector
c     algorithm) from extrapolating past the integration interval 
c     (e.g. to below the bottom on the downward sweep or to above 
c     the surface on the upward sweeps), which might cause problems 
c     in the calls to the abscat routine (even though the over-
c     extrapolation would be corrected to end up at the proper 
c     point).
c
      iflag = -1

        call ode(drtdzs,neqns,RT,zetastrt,zetaend,relerr,abserr,
     1           iflag,work,iwork)

      Enddo		!isteps

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     Save the solution at zetaend 
         DO J=1,nmu
         Smtbz(j,izrev) = RT(j + nmu2)
            DO I=1,nmu
               Rzb(I,J,izrev) = RT(I+(J-1)*nmu)
            end do
         end do
      end do   ! iz loop
c
c     End of upward sweep -----------------------
c
 999  RETURN
      END
