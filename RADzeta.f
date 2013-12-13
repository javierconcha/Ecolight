C     Last change:  LKS  17 Apr 2008    7:27 pm
      SUBROUTINE RADzeta
C 
C     core routine on file RADzeta.f
C
c     called by RADIANCES [MAIN->RADIANCE->RADZETA]
C
C     calls canned routines P2ARAY, SGEFS (in matinv.f), and MATxMAT
C 
c     This routine computes the radiances L+(zeta) = RADdnz
c     and L-(zeta) = RADupz at all INTERIOR depths w .lt. zeta .le. m,
c     using Eqs. (8.105) and (8.106).  Source terms are omitted if
c     isource = 0.
C 
      INCLUDE "DIMENS_XL.INC"
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
c
      COMMON /CRTS/ Rzw(mxmu,mxmu,mxz),Twz(mxmu,mxmu,mxz),
     1              Rzb(mxmu,mxmu,mxz),Sptwz(mxmu,mxz),Smtbz(mxmu,mxz)
      COMMON /CMISC/ imisc(30),fmisc(30)
      Common /Cvarz/ indexz(0:mxwave),zopt(mxwave),zFPAR(mxwave)
c
c     temporary local storage:
      DIMENSION Twzb(mxmu,mxmu),temp1(mxmu,mxmu),temp2(mxmu,mxmu),
     1  temp3(mxmu), temp4(mxmu),temp5(mxmu)
      DIMENSION Raddnw(mxmu)
c     work arrays for matrix inversion routine sgefs
      dimension work(mxmu),iwork(mxmu)
C 
      nmu = imisc(1)
      isource = imisc(8)
c     the bottom depth was determined in varzmax
      jwave = imisc(11)
      nz = indexz(jwave)
C
C     Compute the radiances at each interior zeta level
      DO iz=2,nz
c
c     Begin evaluation of Eq. (8.105) ---------------------------------
c  
C     Compute Twzb = T(w,zeta,b) of Eq. (8.103)
C 
C     Compute temp1 = I - R(z,b) * R(z,w)
      DO i=1,nmu
         DO j=1,nmu
           sum = 0.
             DO k=1,nmu
               sum = sum + Rzb(i,k,iz)*Rzw(k,j,iz)
             end do
           delt = 0.
           IF(I.EQ.J) delt = 1.
           temp1(i,j) = delt - sum
         end do
      end do
c
c     Invert I - R(z,b) * R(z,w)
c
cccccccccc  Canned routine calls  ccccccccccccccccccccccc
c
      DO i = 1,nmu
        DO j = 1,nmu
          temp2(i,j) = 0.
        end do
        temp2(i,i) = 1.0
      end do
c
      itask = 1
      call sgefs(temp1,mxmu,nmu,temp2(1,1),itask,ind,work,iwork)
	if(ind.lt.0) 
     1   write(10,'(/2x,a,i4)')
     2   'Error on first call to SGEFS from RADzeta',IND
c
      if(ind.eq.-10) then
        call HERR("RADzeta","SGEFS returned ind = -10")  !stop run
	endif
c
      itask = 2
      do j = 2,nmu
         call sgefs(temp1,mxmu,nmu,temp2(1,j),itask,ind,work,iwork)
	if(ind.lt.0) 
     1   write(10,*)'Error on second call to SGEFS from RADzeta',IND
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     temp2 now contains the inverse of I - R(z,b)*R(z,w)
c
c     Compute T(w,zeta,b) = T(w,zeta) * (Inverse)
      DO i=1,nmu
         DO j=1,nmu
         sum = 0.0
            DO k=1,nmu
            sum = sum + Twz(i,k,iz)*temp2(k,j)
            end do
         Twzb(i,j) = sum
         end do
      end do
c
c     first term in (8.105):
      do i=1,nmu
         RADdnw(i) = RADdnz(i,1)
      enddo
c
      call VECxMAT(RADdnw,Twzb,nmu,nmu,mxmu,temp3)
c
      if(isource.ne.0) then
c
c     include the source terms:  evaluate the second term in (8.105):
c
      do j=1,nmu
        sum = SPtwz(j,iz)
          do k=1,nmu
            sum = sum + SMtbz(k,iz)*Rzw(k,j,iz)
          end do
        temp4(j) = sum
      end do
c
      call VECxMAT(temp4,temp2,nmu,nmu,mxmu,temp5)
      factor = 1.0
c
      else
c     no source terms:
      factor = 0.0
      endif
c
      do i=1,nmu
         RADdnz(i,iz) = temp3(i) + factor*temp5(i)
      end do
c
c     Begin evaluation of Eq. (8.106) ---------------------------------
c
c     first term (using previously computed Twzb)
c
      do i=1,nmu
        do j=1,nmu
          sum = 0.0
             do k=1,nmu
               sum = sum + Twzb(i,k)*Rzb(k,j,iz)
             end do
          temp2(i,j) = sum
        end do
      end do
c
      call VECxMAT(RADdnw,temp2,nmu,nmu,mxmu,temp3)
c
      if(isource.ne.0) then
c
c     include internal sources:  evaluate the second term
c
c     compute temp1 = I - R(z,w) * R(z,b)
c
      do i=1,nmu
        do j=1,nmu
          sum = 0.0
            do k=1,nmu
              sum = sum + Rzw(i,k,iz)*Rzb(k,j,iz)
            end do
          delt = 0.0
          if(i.eq.j) delt = 1.0
          temp1(i,j) = delt - sum
        end do
      end do
c
c     invert I - R(z,w) * R(z,b)
c
cccccccccc  Canned routine calls  ccccccccccccccccccccccc
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
c
	if(ind.lt.0) 
     1   write(10,*)'Error on third call to SGEFS from RADzeta',IND

      itask = 2
      do j = 2,nmu
         call sgefs(temp1,mxmu,nmu,temp2(1,j),itask,ind,work,iwork)
	if(ind.lt.0) 
     1   write(10,*)'Error on fourth call to SGEFS from RADzeta',IND
      end do
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     temp2 now contains the inverse if I - R(z,w)*R(z,b)
c
      do j=1,nmu
        sum = SMtbz(j,iz)
        do k=1,nmu
          sum = sum + SPtwz(k,iz)*Rzb(k,j,iz)
        end do
        temp4(j) = sum
      end do
c
      call VECxMAT(temp4,temp2,nmu,nmu,mxmu,temp5)
      factor = 1.0
c
      else
c       no source terms:
        factor = 0.0
      endif
c
      do i=1,nmu
         RADupz(i,iz) = temp3(i) + factor*temp5(i)
      end do
C  
      END DO  ! end iz loop
C 
      RETURN
      END
