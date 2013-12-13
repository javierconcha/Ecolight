C     Last change:  LKS   7 Aug 2008    7:08 pm
      function astar(icomp, wavenm)
c
c     user routine on file astar.f   
c
c     called by SHAT (in shat.f), ABCASE1, and ABCASE2
c
c     calls LENSTR, and interpolation routines (in dataintrp.f).
c
c     This function defines the specific absorption coefficient,
c     a*(wavelength) in mg/m^3, for use by routine shat and the
c     "ab" routines.
c
c     The a* specifications are:
c         iastropt = -1:  Subroutine called for aCDOM (a* routine not called)
c                    -2:  aCDOM solved as a percentage of aw+achl in "ab" routine (a* not called)
c		    0:  user-supplied data file read to get a* values
c                     1:  Pope & Fry absorption model used (pure water only)
c                     2:  Smith and Baker absorption model used (pure water only)
c                     3:  Prieur-Sathyenranath-Morel model used (Chlorophyll only)
c                     4:  Exponential model used (CDOM only)
c                     5:  Prieur-Sathyendranath-Morel model used (CDOM only)
c                     6:  Constant: independant of wavelength
c
c         The other a* params are only used when iastropt=4 (exp model)
c         astarRef: reference wavelength for exp model
c         astar0:   astar at the reference wavelength
c         asgamma:  exponential decay constant
c
c
c     For "datafile" or PSM options, on the first call, this routine
c     reads in the data file which contains the a* data, and 
c     puts the data on the run wavel grid.  Subsequent calls return
c     astar at the given wavelength.
c
c     For the PSM option, the values returned are the ac values used in the
c     Prieur-Sathy-Morel model, rather than the true a* (which for their model
c     depends on the Chlorophyll concentration)
c
c     For the "exponential" model option (CDOM), the routine returns a*
c
c     This routine should NOT be called for the pure water options (1 and 2)
c     since a and b for water are handled by pureh2o.f (and there isn't an a* for water)
c
      INCLUDE "DIMENS_XL.INC"
c
c     Common blocks containing a* specifications for all components
      Common /Castar/ iastropt, astarRef, astar0, asgamma
      Common /Castar2/ astarfile
!
      common/cmisc/ imisc(30),fmisc(30)
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave) 

      Integer iastropt(mxcomp+2)
      Real astarRef(mxcomp+2), astar0(mxcomp+2)
      Real asgamma(mxcomp+2)
      Character*120 astarfile(mxcomp+2)
c
      Character*120 datadir, digitdir, spreadir,
     1              phasedir, surfdir,bottdir, Pdir
      COMMON /Cdirnames/ datadir,digitdir,spreadir,
     1                   phasedir,surfdir,bottdir, Pdir

      real wavenm

!     local
      integer kall(mxcomp), itype, nheadr, nh2p
      real astargrid(mxwave,mxcomp+2), a0(mxcomp), zgrid(1)
      real, allocatable,dimension(:) :: agrid
!
      data kall/mxcomp*0/,waveold/0/
      data itype/0/, nheadr/10/, nh2p/0/
!     itype=0 means abscissa is wavelen; nheadr is # header lines in file 
!     and nh2p is the number of header lines to be printed out
c
c     Note on dynamic storage:  MUST save data values between subroutine 
c     calls:
      save
c
      if(kall(icomp).eq.0) then
c     initialization on first call
c
        If((iastropt(icomp).eq.0).or.(iastropt(icomp).eq.3)) then       !data file will be read
c-------- INITIALIZE DATA FILE OPTION ----------
          If(iastropt(icomp).eq.3) then   !Prieur-Sathyenranath-Morel model OPTION
            write(10,103) icomp
          Else
            write(10,100) icomp, astarfile(icomp)
          Endif
          ngrid = imisc(7)
          allocate(agrid(ngrid))
          call getonGrid1(itype, nheadr,nh2p, ngrid,astarfile(icomp), 
     1                    astarRef(icomp),a0(icomp),wave,agrid, zgrid)
          astargrid(:ngrid, icomp) = agrid(:ngrid)
          deallocate(agrid)
        Elseif(iastropt(icomp).eq.4) then
c-------- INITIALIZE Exponential model OPTION ----------
          write(10,104) icomp,astar0(icomp),asgamma(icomp),
     1                  astarRef(icomp)
c
        Elseif(iastropt(icomp).eq.6) then
c-------- INITIALIZE Constant OPTION ----------
          write(10,106) icomp,astar0(icomp)
c
        Elseif(iastropt(icomp).lt.-2 .or. iastropt(icomp).gt.6) then
          write(10,*) 'ASTAR: component #',icomp
          write(10,*) '       Invalid option', iastropt(icomp)
          astar = 1.0
          return
        Endif
        kall(icomp) = 1
      endif
c-----End of initialization on first call.
c
c     Subsequent calls start here:
      IF (waveold.ne.wavenm) THEN
        waveold = wavenm
      END IF
c
      If((iastropt(icomp) .eq. 0).or.(iastropt(icomp) .eq. 3)) then       !data file will be read
c-------- EVALUATE DATA FILE OPTION ----------
        If(abs(wavenm-astarRef(icomp)).le.1e-4) then
          astar = a0(icomp)
        Else
          jwave = imisc(11)     !index to current wavelen
          call checkiWav(jwave, wavenm)
          astar = astargrid(jwave,icomp)
        Endif
      Elseif(iastropt(icomp).eq.4) then
c-------- EVALUATE Exponential model OPTION ----------
          astar = astar0(icomp)*exp(-asgamma(icomp)*
     1            (wavenm-astarRef(icomp)))
c
      Elseif(iastropt(icomp).eq.6) then
c-------- EVALUATE Constant OPTION ----------
          astar = astar0(icomp)
      Else
          astar= 1.0
      ENDIF

      return
 100  FORMAT(/5x,'Specific Absorption for component',
     1 i2, ' read in from file: ',a)
 103  FORMAT(/5x,'Absorption for component',
     1 i2, ' calculated from the Prieur-Sathyenranath-Morel model')
 104  FORMAT(/5x,'Specific Absorption for component',
     1 i2, ' calculated from an Exponential model:',/10x,
     2 'astar(waven) = ',f6.4,'* exp[ -',f6.4,'*(waven - ',f4.0,')]')
 106  FORMAT(/5x,'Specific Absorption for component',
     1 i2, 'assumed constant:',/10x, 'astar = ',f6.4)
      end
