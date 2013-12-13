C     Last change:  LKS  25 May 2008    7:58 am
      subroutine getDynZ(iwave)
!     This function returns the max depth that should be made for
!     wavelen index iwave
!     Used when inelastic scatt and inf bottom selected

      INCLUDE "DIMENS_XL.INC"
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               zgeo(mxz),zeta(mxz)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 
      COMMON /Cmisc/ imisc(30),fmisc(30)
c     temporary local storage:
      real acomp(mxcomp),bcomp(mxcomp)
!
      real,allocatable,dimension(:) :: ztemp, zmax
      integer,allocatable,dimension(:) :: izmax  
      integer nizmax(mxwave)
      real zb, wavel, atotal,btotal
      integer i, i0, i1, nwave, ncomp    
      save nizmax
!-----------------------------------------------------------------------
      If(iwave.eq.0 ) then  !initialize
        nwave = imisc(7)
        ncomp = imisc(6)
        nz = imisc(4)
        allocate(zmax(nwave),ztemp(nwave), izmax(nwave))
        zb = zeta(imisc(4))   !requested final depth, at init still in zeta array 
        !   calc the max depth for each wavelength
        Do i = 1, nwave  !nwave
          wavel = wave(i)
          imisc(11) = i
!         insert the call to the desired "ab" routine:
          call abscat(zb,wavel,ncomp,acomp,bcomp,atotal,btotal)
          ztemp(i) = zb + 20./(atotal + btotal)           !add 20 dtau
        Enddo
        !  set max needed for each iwave equal to the biggest zmax at higher wavel
        i1 = 1
        niz = 0
 10     z1 = ztemp(i1)
        i0 = i1 + 1
        Do i=i0, nwave   !find largest zb between wavel i0 and nwave
          if(ztemp(i).ge.z1) then
            z1 = ztemp(i)
            i1 = i
          elseif(ztemp(i)-z1 .ge. -1.0e-1) then  !force depths to be at least 0.1m apart
            i1 = i
          endif
        Enddo
        niz = niz + 1
        zmax(niz) = z1
        Do i=i0-1,i1
          izmax(i) = niz
        Enddo
        If(i1.eq.nwave-1) then
          niz = niz + 1 
          izmax(nwave) = niz
          zmax(niz) = ztemp(nwave)
        ElseIF(i1.ne.nwave) then
          i1 = i1 + 1
          goto 10
        Endif
        deallocate(ztemp)

!       check to see if we have enough z memory to store ALL wavelen depths
        nzAvail = mxz - nz
        If(nzAvail.ge. niz+4) then  !we have plenty!
!        put a few depths between zb and first zcalc_max
         Do i=1,4
           zgeo(nz+i) = zgeo(nz) + i*0.2*(zmax(niz)-zgeo(nz)) 
         Enddo
         nz = nz + 4

        ElseIf(nzAvail.le.0) then    !cannot use DynZ
         write(10,120) mxz, nz+niz+4
         Do i=1, nwave
           nizmax(i) = nz
         Enddo
         deallocate(zmax, izmax)
         imisc(27) = 0 !turn OFF DynZ
         return   !cannot use DynZ
        
        Else    !need to budget our dynamic depths; space evenly
         write(10,110) nzAvail, niz+4, nz+niz+4
         dz = (zmax(1) - zb) / float(nzAvail)   !spacing size
         Do j = 2, nzAvail
           zlim = zmax(1) - dz * (j-1)
           Do i = 1, niz
             if(zmax(i).le.zlim) then
               izmax(i) = j
             endif
           Enddo
           zmax(j) = zlim
         Enddo
         niz = nzAvail
        Endif

!        store each needed wavelen calc depth in ascending order
        Do i=niz, 1, -1
          zgeo(nz + i) = zmax(niz - i + 1)
        Enddo
        Do i=1, nwave
          nizmax(i) = nz + niz - izmax(i) + 1
        Enddo
        imisc(4) = nz + niz     !update nz    
        deallocate(zmax, izmax)

!     subsequent calls start here....
      ELSE
        imisc(4) = nizmax(iwave)    !override nz to solve
        write(10,100) zgeo(imisc(4)), wave(iwave)
      ENDIF
      return
 100  format(/2x,'The infinitely deep bottom boundary condition is ',
     1 'applied at ',f5.1,' m at ',f6.1,' nm'/)
 110  format(/5x,' NOTICE:  consider increasing',
     1   ' the output depth LIMIT in the UI',/15x,
     2   'to make the dynamic depth option run more efficiently ',
     3   '(faster)',/15x,
     3   'nz Available: ',i4,5x,'nz Requested: ',i4,5x,
     4   'Optimum mxz for this run: ',i4)
 120  format(/5x,' NOTICE:  Dynamic Depth cannot be used because all ',
     1   'output depths are in use.',/15x,'Consider increasing',
     2   ' the output depth LIMIT in the UI',/15x,
     3   'to allow the dynamic depth option to be used',/15x,
     4   'Current mxz: ',i4,5x
     5   'Optimum mxz for this run: ',i4)
      end subroutine 
