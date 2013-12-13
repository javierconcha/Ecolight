C     Last change:  LKS   9 Dec 2011   10:45 am
      subroutine RADw
c 
c     core routine on file RADw.f
c
c     called by RADIANCE [MAIN->RADIANCE->RADW]
c
c     calls canned routines MATxMAT and sgefs (in matinv.f)
c 
c     This routine computes the radiances L+(w) = RADdnz(z=0) and
c     and L-(w) = RADupz(z=0) using Eqs (8.98) and (8.102), respectively.
c     It is via these equations that the air-water surface boundary
c     effects are incorporated into the final solution.
c     Source terms are omitted if isource = 0.
C 
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /CRTS/ Rzw(mxmu,mxmu,mxz),Twz(mxmu,mxmu,mxz),
     1              Rzb(mxmu,mxmu,mxz),Sptwz(mxmu,mxz),Smtbz(mxmu,mxz)
      COMMON /Crtsurf/ raw(mxmu,mxmu),taw(mxmu,mxmu),rwa(mxmu,mxmu),
     1                 twa(mxmu,mxmu)
C
c     radiance arrays are stored as follows:
c        RADsky(mu) = total incident sky radiance
c        RADwla(mu) = Lw(mu,z=a) = water-leaving radiance
c        RADrsa(mu) = Lreflsky(mu,z=a) = surface-reflected sky rad
c        RADupa(mu) = Lu(mu,z=a) = total upwelling rad in air
c           note: RADupa = RADwla + RADrsa, i.e., Lu = Lw + Lreflsky in air
c        RADupz(mu,z) = Lup(mu,z) = total upwelling rad in water
c        RADdnz(mu,z) = Ldn(mu,z) = total downwelling rad in water
      COMMON /Crad/ RADsky(mxmu),RADwla(mxmu),RADrsa(mxmu),RADupa(mxmu),
     1              RADupz(mxmu,mxz),RADdnz(mxmu,mxz)

      COMMON /CMISC/ imisc(30),fmisc(30)
c
c     temporary local storage:
      Dimension temp1(mxmu,mxmu),temp2(mxmu,mxmu),temp3(mxmu),
     1 temp4(mxmu),temp5(mxmu)

c     work arrays for matrix inversion routine sgefs
      dimension work(mxmu),iwork(mxmu)
c
      nmu = imisc(1)
      isource = imisc(8) 
C
c     Begin evaluation of Eq. (8.98) -----------------------------
c
C     compute temp1 = R(w,b) * r(w,a)
c        Note: Rwb(mu,mu) = Rzb(mu,mu,z=0)
C
      call MATxMAT(Rzb(1,1,1),rwa,nmu,nmu,nmu,mxmu,mxmu,temp1,mxmu)
C
C     Compute the inverse of temp1 = I - R(w,b)*r(w,a)
C
      do j=1,nmu
         do i=1,nmu
             if(i.eq.j) then
               delt = 1.0
             else
               delt = 0.0
             endif
             temp1(i,j) = delt - temp1(i,j)
         end do
      end do
c
cccccccccc  Canned routine calls for matrix inversion  ccccccccccccccc
c
c     The matrix temp1 is inverted by first obtaining an LU
c     decomposition, and the backsubstituting column by column
c     on the identity matrix
c 
      DO i = 1,nmu
         DO j = 1,nmu
            temp2(i,j) = 0.0
         end do
         temp2(i,i) = 1.0
      end do
c     temp 2 now contains the identity matrix
c
c    obtain the LU decomposition of temp1 (itask = 1) and then 
c    backsubstitute for the first column vector
c
      itask = 1
      call sgefs(temp1,mxmu,nmu,temp2(1,1),itask,ind,work,iwork)
	if(ind.lt.0)
     1   write(10,*)'Error on first call to SGEFS from RADW',IND
c
c     temp1 now contains something related to the LU decomposition.
c     use the existing LU decomposition (itask = 2) and backsubstitute 
c     for the remaining column vectors to obtain the full inverse
c
      itask = 2
      do j = 2,nmu
         call sgefs(temp1,mxmu,nmu,temp2(1,j),itask,ind,work,iwork)
	if(ind.lt.0)
     1   write(10,*)'Error on second call to SGEFS from RADW',IND
      end do
c
c     temp2 now contains the inverse of I - R(w,b)*r(w,a)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     compute temp3 = RADsky*t(a,w)
c
      call VECxMAT(RADsky,taw,nmu,nmu,mxmu,temp3)
c
      if(isource.ne.0) then
c
c        include internal sources:  compute temp4 = S-t(b,w) * r(w,a)
c        Note: S-t(b,w) = S-t(b,z=0) = Smtbz(mu,1)
c
         call MATxMAT(Smtbz(1,1),rwa,1,nmu,nmu,1,mxmu,temp4,1)
         factor = 1.0
c
      else
c        no internal sources:
         factor = 0.0
      endif
c
      do i=1,nmu
         temp4(i) = temp3(i) + factor*temp4(i)
      end do
c
c     compute RADdnz(w) by Eq. (8.98)
c
      call VECxMAT(temp4,temp2,nmu,nmu,mxmu,temp5)
c     copy the row vector into the depth = 0 column of RADdnz
      do i=1,nmu
         RADdnz(i,1) = temp5(i)
      enddo
c
c     Begin evaluation of Eq. (8.102) -------------------------------
c
c     first term:
c
c     compute temp4 = RADsky*t(a,w)*(inverse)
c     temp2 contains the inverse of I - R(w,b)*r(w,a) from above
c     temp3 contains RADsky*t(a,w) from above
c
      call VECxMAT(temp3,temp2,nmu,nmu,mxmu,temp4)
c
      call VECxMAT(temp4,Rzb(1,1,1),nmu,nmu,mxmu,temp5)
c
      if(isource.ne.0) then
c        include internal sources:
c
c        compute r(w,a) * R(w,b)
c
      call MATxMAT(rwa,Rzb(1,1,1),nmu,nmu,nmu,mxmu,mxmu,temp1,mxmu)
c
C     Compute the inverse of temp1 = I - r(w,a) * R(w,b)
C
      do j=1,nmu
         do i=1,nmu
            if(i.eq.j) then
               delt = 1.0
            else
               delt = 0.0
            endif
         temp1(i,j) = delt - temp1(i,j)
         end do
      end do
c
cccccccccc  Canned routine calls for matrix inversion  ccccccccccccccc
c
      DO i = 1,nmu
         DO j = 1,nmu
            temp2(i,j) = 0.0
         end do
      temp2(i,i) = 1.0
      end do
c
      itask = 1
      call sgefs(temp1,mxmu,nmu,temp2(1,1),itask,ind,work,iwork)
	if(ind.lt.0)
     1   write(10,*)'Error on third call to SGEFS from RADW',IND
c
      itask = 2
      do j = 2,nmu
         call sgefs(temp1,mxmu,nmu,temp2(1,j),itask,ind,work,iwork)
	if(ind.lt.0)
     1   write(10,*)'Error on fourth call to SGEFS from RADW',IND
      end do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        second term:  S-t(b,w) * (inverse)
c
         call MATxMAT(Smtbz(1,1),temp2,1,nmu,nmu,1,mxmu,temp4,1)
         factor = 1.0
c
      else
c        no internal sources:
         factor = 0.0
      endif   ! end isource test
c
c     Compute RADupz(w) by Eq. (8.102)
c
      do i=1,nmu
         RADupz(i,1) = factor*temp4(i) + temp5(i)
      end do
c
      return
      END 

