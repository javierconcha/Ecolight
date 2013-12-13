C     Last change:  LKS  25 May 2008    7:27 pm
      subroutine initprintout(idbug)

c
      INCLUDE "DIMENS_XL.INC"
c
      COMMON /Cpirrad/ npirad,izirad(mxz)
      COMMON /Cpkfcn/ npkfcn,izkfcn(mxz) 
      COMMON /Cmisc/ imisc(30),fmisc(30)

      nz = imisc(4)

cccccc  SET PRINTOUT PARAMETERS  cccccccccccccccccccccccccccccccccccc
c     check for array dimensions being larger than allowed
      if(nz.gt.mxz) then
         write(6,fmt='(" nz =",i3," gt mxz =",i3)') nz,mxz/2
         write(10,fmt='(" nz =",i3," gt mxz =",i3)') nz,mxz/2
         call HERR("initprintout",
     1        "increase LIMIT of number of ouput depths in UI")  !stop run
      endif
c
c     Set parameters for how much printout is to be given on the
c     standard printout file.
c
c     The defaults as set below give output that is typically of 
c     interest to oceanographers.  These values can be changed to 
c     get more or less printout.
c
c.....IRRADIANCE and IOP printout.  DEBUG some more !!!!
c
c      use ipirad = 0 for NO printout of the irradiances or IOP values
c                 = 1  for printout only of the values IN AIR (at depth z = a)
c                      (useful for remote-sensing studies)
c                 = 2 for printout at the user-selected depths (the odd
c                     depths)
c                 = 3 for printout at all depths
c
c     DEFAULT:  print irradiances and IOPs at the user-requested depths:     
	If(idbug .eq. -1) then  
		ipirad = 0		!don't printout irrad if selected "minimal" output
	Else
		ipirad = 2		!give "standard" output
	Endif
      ipirad = 2
c
cc      if(ipirad.eq.0) then
c        no printout of irradiances or IOPs
cc         npirad = 0			<-- this will crash pntgrid
cc      else...
	if(ipirad.eq.1) then
c        printout only of the air values (irradiances only)
         npirad = 1
      elseif(ipirad.eq.2 .or. ipirad.eq.0) then
c        printout at all user-requested depths
      do iz=1,nz
         izirad(iz) = iz
      end do
      npirad = nz
      else
         print*,' invalid selection for irradiance printout: ipirad = ',
     1            ipirad
      endif 
c
c.....IRRADIANCE K-FUNCTION printout.
c
c      use ipkfcn = 0 for NO printout of the irradiance K-functions
c                 = 1 for printout of layer averaged K functions, computed
c                     from each pair of the user-selected depths

c     DEFAULT:  print IRRADIANCE K-functions at the user-selected depths
c
      ipkfcn = 1
c
      if(ipkfcn.eq.0) then
c        no printout
         npkfcn = 0
      elseif(ipkfcn.eq.1) then
c        printout at user-selected depths
         npkfcn = 0
         do iz=1,nz-1
            npkfcn = npkfcn + 1
            izkfcn(npkfcn) = iz
         end do
      else
         print*,' invalid selection for K-function printout: ipkfcn = ',
     1            ipkfcn
      endif
c
      return
      end subroutine

!-----------------------------------------------------------------------
      subroutine pntBRDFbotm
!     Prints the appropriate BRDF msg on initialization
!     in Ecolight, bottom is always assumed lambertian 
c
      write(10,100)
      
  100 format(/5x,'The bottom bidirectional reflectance distribution',
     1' function (BRDF) is Lambertian:'/8x,
     2"BRDF(mu',phi',mu,phi) = R/pi")
      end     
