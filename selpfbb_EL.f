C     Last change:  LKS   8 Nov 2008    9:37 am
       subroutine selpfbb(bbfract, komp)
c
c      Ecolight version, 10/06 lks
c
! 2006 Feb 08 Marcos Montes (MJM) [Comments modified 2006 June 27.]
!
!            This is a "drop-in" replacement for the distributed selpfbb.
!            The old version needs to be saved in a separate folder. After that,
!            put this version into the maincode folder of the HydroLight
!            distribution.
!            
!            This version has been tested with several versions of Hydrolight 4.2.
!            This has not (yet) been tested with Hydrolight 4.1.
!
!            NOTE: The default bb/b interval could easily be dropped 
!            (and should be based on my testing) to <=0.0001 after placing
!            this file in HydroLight's maincode directory.
!  
!            Notes are included to comment/uncomment certain code sections
!            if using Fortran-77 (f77) instead of a Fortran-90 or better
!            (f90+) compliant compiler. Comments can either be a "c" or a "!" 
!            in the first column of the line. Follow the directions in the code
!            as to what needs to be commented and what needs to be uncommented. 
!
!            Commented betaM1, betaM2 declarations, not needed.
!            Added calculations of betaP1, betaP2, but only in debug section.
!            Saved i1, i2: check to see if bbfract is in that interval first.
!            Rearranged search for correct bracketing bbvals(). More elegant now.
!            bbvals doesn't need to be a common. Just save it. 
!            More consistent error checking on files. Now stops only if 
!            all files fail test. Failing test includes not opening or 
!            being wrong angular grid. Files that don't open or have wrong grid 
!            are not included in list to be sorted, and their bbvals are not 
!            recorded. This yields stronger code. 
!
!            Opens are now only performed on the the very first call. 
!            All the Fournier-Forand phase functions are read in at once into 
!            an array that is allocated on the fly.
!            It is sorted into the correct order. Now, aside from the first call, 
!            all that is done is an interpolation.
!            We're at least 10, to more than 100 times faster, depends on a lot 
!            of input parameters, so hard to quantify.
!            NOTE: If you don't want to use f90's allocate, then just
!            pre-allocate the array as usual. Need to have enough locations
!            (outer dimension) for the number of Fournier-Fourand files. 
!            That could be carried in dimens.inc. However, using the LF95 compiler, 
!            this works as is and there is no need to change the code. 
!            Added variable mjm_init (trivial); initialized nfiles 
!            in declaration which makes it "saved" under F90 and later.
!            Also, moved the check for consistent quads to the same location,
!            so the check is only performed once. 
!
!
!      Original version by C. Mobley & L. Sundman, their code date was 2001 Sep. 6.  
c      core routine on file selpfbb.f 
c
c      called by INISHAMP
c
c      calls: LENSTR
c 
c      Given the backscatter fraction bbfract (0.0 to 1.0), this routine
c      (1) reads the list of phase function files that are available for
c      use in fitting backscatter fractions (the list is on file
c      ..\data\phasefun\pfbblist.txt)
c      (2) locates the two *.dpf files that have phase functions with
c      backscatter fractions bracketing the given bbfract
c      (3) linearly interpolates to generate a phase function that
c      has the desired backscatter fraction; this phase function is 
c      stored in the betatP and betatM arrays for component number komp
c
      INCLUDE "DIMENS_XL.INC"
      PARAMETER (mxL=mxphi/2, mxbeta=mxmu*(mxL+1))
      
!      parameter (mxbb=50) ! uncomment for f77, comment for f90
      integer,save::mxbb ! comment for f77, uncomment for f90      
!     NOTE: mxbb controls the maximum number of phase functions.
!     For f77, it is fixed at whatever it is defined above. For 
!     f90+, it is dynamically determined, and limited by the maximum
!     default integer, usually (2**15)-1 .     
c
      integer komp, nuphas
c
      COMMON /Crhotau/ rho(mxmu,mxmu),tau(mxmu,mxmu),
     1                 betatP(mxmu,mxmu,mxcomp),betatM(mxmu,mxmu,mxcomp)
      COMMON /Cmisc/ imisc(30),fmisc(30) 
      COMMON /Cgrid/ fmu(mxmu),bndmu(mxmu),omega(mxmu),deltmu(mxmu),
     1               z(mxz),zeta(mxz)

      Character*120 datadir, digitdir, spreadir,
     1              phasedir, surfdir,bottdir, Pdir
      COMMON /Cdirnames/ datadir,digitdir,spreadir,
     1                   phasedir,surfdir,bottdir, Pdir
c
c     storage of SURFWIND values to compare with pf (from inishamp)
      COMMON /Cbbquadchk/ m,iqpart     !loadsurface, inishrad & selpfbb
c 
! Fortran-90+ version: uncomment next three lines for f90+
      character (len=120),dimension(:),allocatable::pfbbname
      real,dimension(:),allocatable::bbvals
      integer,dimension(:),allocatable::iperm
!     For f77 only, comment above three lines and uncomment next two lines. 
!      dimension pfbbname(mxbb),bbvals(mxbb),iperm(mxbb)
!      character pfbbname*120
       character pftitl*120,filename*120
                 
c
      dimension betaP1(mxmu,mxbeta),!MJM betaM1(mxmu,mxbeta),
     1          betaP2(mxmu,mxbeta) !MJM ,betam2(mxmu,mxbeta)
      
      ! MJM - initialize so we don't do this a lot of times
      ! Not F77 compliant!
      integer:: mjm_init=0, nfiles=-1,i_mjm,good_nfiles=-1,iqqq
      ! F77 compliant (but why use F77?) - comment one line above
      ! and uncomment three lines below
      ! integer mjm_init, nfiles, i_mjm_good, good_nfiles, iqqq
      ! data mjm_init/0/, nfiles/-1/, good_nfiles/-1/ 
      ! save mjm_init, nfiles,good_nfiles
     
      integer::i1=1, i2=2 ! uncomment for f90+, comment for f77
      !integer i1, i2 ! uncomment for f77, comment for f90+
      !data i1/1/, i2/2/ ! uncomment for f77, comment for f90+
      !save i1, i2 ! uncomment for f77, comment for f90
      save bbvals
      !
      ! Not F77 compliant!
      real, allocatable, save,dimension(:,:,:)::betaMall,betaPall
      ! F77 version - comment one line above, uncomment 2 lines below
      ! real betamall(mxmu,mxbeta,mxpfbbf),betaPall(mxmu,mxbeta,mxbb)
      ! save betaMall, betaPall
      
      ! Next line is for F90+; for f77 only, comment next line....
      character (len=120),dimension(:),allocatable::good_pfbbname 
      !F77 version, comment line above and uncomment next two lines
      !character*120 good_pfbbname !uncomment for F77
      !dimension good_pfbbname(mxbb) ! uncomment for F77
      ! MJM  - so pfbbname() will have all the possible file names,
      ! including ones that don't open, or have wrong angular grid. 
      ! nfiles is its length. good_pfbbname is the subset that both opens
      ! and has the correct angular grid. good_nfiles is its length.
      real bbfuse
      character (len=120)::junk ! Comment for f77, uncomment for F90+
      integer::ios ! Comment for f77, uncomment for F90+

      nmu = imisc(1)
      kcol = nmu   !EL version only
      wavelen = fmisc(13)
      depth   = fmisc(18)
c
c     set unit number for the phase function for this component
c     Note:  both files will be opened sequentially w/ this unit number
      nuphas = 49 + komp
c
      idbug = imisc(9)
	If(idbug.gt.0) then  
c		Inform the user of the phase function selected, iff added printout requested
		write(10,551) komp,bbfract, depth, wavelen
	Endif
      
      ! Only do the listing if first entry, i.e., mjm_init=0
      if (mjm_init.eq.0) then
c
c     read the list of *.dpf files to be used in matching bbfract
         lenpd = lenstr(phasedir)
         filename = phasedir(1:lenpd) // 'pfbblist.txt'
         len = lenstr(filename)
         open(nuphas,file=filename(1:len), err=666)
! MJM First, count the number of entries there are. This requires using
! a very very large upper limit. In practice there are probably less than 
! 100 files. There are really elegant methods of coding this in strictly
! f90+. For the f90+ versions, this determines mxbb, and lets one have large
! lists of files in the pfbblist.txt file. 
! MJM: IF f77 comment next 19 lines
         ios=0
         i=0
         do
            read(nuphas,'(a)',iostat=ios)junk
            if (ios.ne.0) exit ! exits this loop
            i=i+1
         enddo
         rewind(nuphas) ! go back to start of file
         mxbb=i
         if (mxbb.le.0) then
           write(10,*)'selpfbb: no valid filenames in ',filename(1:len)
           call HERR("SELPFBB","no valid filenames in pfbblist.txt")  !stop run
         endif
         if (allocated(pfbbname)) deallocate(pfbbname)
         if (allocated(bbvals)) deallocate(bbvals)
         if (allocated(iperm)) deallocate(iperm)
         allocate(pfbbname(mxbb),bbvals(mxbb),iperm(mxbb))
 ! MJM: If f77, comment 19 lines immediately above
 
 ! This section is used by both f77 and f90+ versions
 ! For f77, mxbb is a parameter. For f90, is just a variable.
         do i=1,mxbb
            read(nuphas,'(a)',end=100) pfbbname(i)
            lenpfbb = lenstr(pfbbname(i))
            if(lenpfbb.eq.0) go to 100
         end do
         i=mxbb + 1
         
  100    nfiles = i - 1
         close(nuphas)

         ! If f77, comment next two lines
         if (allocated(good_pfbbname)) deallocate(good_pfbbname)
         allocate(good_pfbbname(nfiles))
         ! If f77, comment above two lines
c
c     open the *.dpf files one at a time, and read their bbfract values
c     from the second record
c
 
         iqqq=0 ! counter for "good" files that open
         do i=1,nfiles
!!            filename=''
            filename = phasedir(1:lenpd) // pfbbname(i)
            len = lenstr(filename)
            open(nuphas,file=filename(1:len), err=105)
            read(nuphas,'(a)') pftitl
c        the following line is new in H4.1
         read(nuphas,402) Udum,mdum,iqpartdum,nraydum
! File name check moved from end of file to here by MJM
c          Check to see that the discretized phase functions have
c           the same quads as the surfwind file
c
!           Ideally, phase functions with different discretizations should be kept in 
!           different directories. - MJM
            iprob = 0
	      if(mdum.ne.m) then
	         iprob = 1
		       write(10,600) filename
	      endif
!!	      if(ndum.ne.n) then
!!	         iprob = 1
!!		       write(10,601) filename
!!	      endif
	      if(iqpartdum.ne.iqpart) then
	         if(iprob.eq.0) iprob = -1
		       write(10,602) filename
	      endif
c
            If(iprob.gt.0) then ! bad file, will not use it.
                write(10,*)'Quads of DPF file do not match surfwind.'
                write(10,*)'Ignoring ',filename
                write(10,*)'Will attempt to read and use other FF files'
	      else ! it is good, so read bbvals, save name
                iqqq=iqqq+1
                good_pfbbname(iqqq)=pfbbname(i)
                read(nuphas,'(13x,f10.5)') bbvals(iqqq)
            endif
            close(nuphas)
  105    end do
         
         good_nfiles=iqqq
         if (good_nfiles.eq.0) then 
           call HERR("SELPFBB","no valid FF DPF files")  !stop run
         elseif (good_nfiles.eq.1) then
           call HERR("SELPFBB","only 1 valid DPF file (cannot interp)")  !stop run
         endif
c
c     sort the bb values into ascending order
c
         kflag = 2
         call spsort(bbvals,good_nfiles,iperm,kflag,ier)

         ! If using a strictly F77 compiler, comment out next 4 lines
         if (allocated(betaMall)) then ! Paranoia, really
            deallocate(betaMall)
         endif
         allocate(betaMall(nmu,kcol,good_nfiles))
         
         ! If using a strictly f77 compiler, comment out next 4 lines
         if (allocated(betaPall)) then ! Paranoia again
            deallocate(betaPall)
         endif
         allocate(betaPall(nmu,kcol,good_nfiles))


! MJM  Now we read the files in the correct order, starting at the lowest, 
! MJM  going to the highest
!      Permutations (i.e., iperm()) are not needed after this loop.
         do i=1,good_nfiles
            filename = phasedir(1:lenpd) // good_pfbbname(iperm(i)) !sorted now
            len = lenstr(filename)
            open(nuphas,file=filename(1:len))
            read(nuphas,'(a)') pftitl
            read(nuphas,402) Udum,mdum,iqpartdum,nraydum
            read(nuphas,'(13x,f10.5)') bbvals(i) ! automatically sorted now
            !print*,filename(1:len),bbvals(i)
            read(nuphas,'(a)') pftitl            !read header, only in EL DPF files
 !     print *, i, bbvals(i), trim(pftitl)
            do j=1,kcol
               read(nuphas,404) (betaPall(j, i_mjm,i),i_mjm=1,nmu)
            end do
            read(nuphas,'(a)') pftitl            !read header, only in EL DPF files
!      print *, i, bbvals(i), trim(pftitl)
            do j=1,kcol
               read(nuphas,404) (betaMall(j, i_mjm,i),i_mjm=1,nmu)
            end do
            close(nuphas)
       
         enddo
!         print*,'Done reading the files'
!         print*,'good_nfiles=',good_nfiles
         mjm_init=1
      endif ! initialization; mjm_init test above 
      
c     print a warning if the requested bb/b value is outside the range
c     of available data
      if(bbfract.lt.bbvals(1)) then
           write(10,110) bbfract, komp, depth, bbvals(1)
      elseif(bbfract.gt.bbvals(nfiles)) then
           write(10,110) bbfract, komp, depth, bbvals(nfiles)
      endif
  110 Format(/2x,'*** warning ***  Requested bb/b value of ',f5.3,
     1       ' for component # ',i4, ' at depth ', f6.1, ' m ',
     1       /5x' is outside of the range of available values.',
     2       /5x,'The phase function with a bb/b value of ',f5.3,
     3       ' will be used instead.')
c
c     locate the available bb values that bracket the requested bbfract
c     (use the nearest bb value if bbfract is less than or greater than
c     the available values)
c
      ! First, see if bbfract is in the current interval (pretty likely)
      ! This saves having to locate it again.
      if ((bbfract.ge.bbvals(i1)).and.(bbfract.le.bbvals(i2))) then
         bbfuse=bbfract
      else
!MJM  Since the range is restricted, no special cases are needed!
!MJM This restricts the range, makes loops below simpler
!MJM no matter what, range is restricted to be within the bbvals array range
         bbfuse=min(max(bbfract,bbvals(1)),bbvals(good_nfiles))
         loop: do i=2,good_nfiles 
            if(bbfuse.le.bbvals(i)) exit loop
         end do loop
         i1 = i - 1
         i2 = i
!MJM Range is restricted. So, lowest value is bbvals(1). When bbuse==bbvals(1),
!MJM i1=1, i2=1, xinterp=0, correct value obtained below.
!MJM Largest allowed value is bbvals(good_nfiles). When bbfuse==bbvals(good_nfiles),
!MJM then i=good_nfiles, i1=good_nfiles-1, i2=good_nfiles, xinterp=1.0,
!MJM which yields correct expresion below. 
      endif
      xinterp = (bbfuse - bbvals(i1))/(bbvals(i2) - bbvals(i1))
!      print*,'i1, i2, bbfract, bbfuse, xinterp=',i1, i2, bbfract, 
!     1  bbfuse, xinterp
      
c
c     interpolate to define the final phase function 
c
c      call p3aray(betatP,nmu,nmu,komp,mxmu,mxbeta,2,'betatP at start')

!MJM betaP1, betap2, betam1, betam2 are not assigned any more. Use the 
!MJM allocated array directly. i1 and i2 refer to the array not the permuted 
!MJM array. Permutations are not needed anymore since we read the array
!MJM in the correct order.
!MJM NOT F77 compliant, but noone should be using F77!!!! 
      betatP(:nmu,:kcol,komp)= (1.0 - xinterp)*betaPall(:,:,i1) +
     1                        xinterp*betaPall(:,:,i2)
      betatM(:nmu,:kcol,komp) = (1.0 - xinterp)*betaMall(:,:,i1)+
     1                        xinterp*betaMall(:,:,i2)
!MJM If you really need to use F77, then comment the above  4 lines, and 
!MJM uncomment the next 8 lines lines below
!      do iqi=1,kcol
!         do jqj=1,nmu
!            betatP(jqj,iqi,komp)= (1.0 - xinterp)*betaPall(jqj,iqi,i1) +
!     1                        xinterp*betaPall(jqj,iqi,i2)
!            betatM(jqj,iqi,komp)= (1.0 - xinterp)*betaMall(jqj,iqi,i1)+
!     1                        xinterp*betaMall(jqj,iqi,i2)
!         enddo
!      enddo

      idbug = 0  !lks
      if(idbug.ne.0) then
        betaP1(:nmu,:kcol)=betaPall(:,:,i1) !MJM
        call p2aray(betaP1,nmu,nmu,mxmu,2,'betaP1')
        betaP2(:nmu,:kcol)=betaPall(:,:,i2) !MJM
        call p2aray(betaP2,nmu,nmu,mxmu,2,'betaP2')
        call p3aray(betatP,nmu,nmu,komp,mxmu,mxbeta,2,'betatP interp')
c     Compute the check on the total scattering by the last equation
c     on page 386 (the quad-averaged form of Eq. 3.8).  Compute the
c     backscatter fraction from the betatM array and the polar cap
c     incident direction. 
c 
        write(6,208)
  208   FORMAT(//2x,'Checksums on interpolated',
     1' band-averaged phase function:'
     1//' ',2x,'r    Sum (=1)'/)

        do ir=1,nmu
          sump = 0.0
          bbtilde = 0.0
          do iu=1,nmu
            factr = omega(iu)/omega(ir)
            sump = sump + (betatP(ir,iu,komp)+betatM(ir,iu,komp))*factr 
            if(ir.eq.nmu) bbtilde = bbtilde+betatM(ir,iu,komp)*factr 
          end do
          write(6,210) ir,sump
  210     FORMAT(1x,I3,F11.5)
        end do


        write(10,212) bbtilde
  212   format(/2x,'The backscatter fraction is bbtilde =',f8.5)
        bbfract = bbtilde
        call HERR("SELPFBB","Normal stop in DEBUG mode")  !stop run
      endif
c
      return
  666 call nofile(10,"SELPFBB", filename )   !err opening file
c
 402  format (f6.2,2i3,i12 )
 404  format(10(e12.6,1x))
 551  format(/2x,'The phase function for component',i3,
     1' will be chosen to have a bb/b ratio of ',f8.4, /5x,
     2' at a depth of ', f8.3,' and wavelength ', f8.1, ' nm')
 600	format(/'Error in SELPFBB: surfwind and DPF files must have the'
     1'same spatial quads!',/'The MU quads read from phase function ',
     25x,a, /10x,' do not match the surfwind file')
 601	format(/'Error in SELPFBB: surfwind and DPF files must have the'
     1'same spatial quads!',/'The PHI quads read from phase function ',
     25x,a, /10x,' do not match the surfwind file')
 602	format(/'Error in SELPFBB: surfwind and DPF files must have the'
     1'same spatial quads!',/'The quad TYPE read from phase function ',
     25x,a, /10x,' do not match the surfwind file')

      end
