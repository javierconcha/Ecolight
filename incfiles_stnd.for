C     Last change:  LKS  15 Dec 2010    3:21 pm
!------------------------------------------
c     Set run params which control how incfiles_std handles run
c      iparam(1) = IOP model
c      iparam(2) = sky model
c      iparam(3) = irradat (iff=1, get incident irrad from datafile)
c      iparam(4) = chl component index
c      iparam(5) = CDOM component index
!------------------------------------------
      subroutine abmodel(depth, dcompconc)
c
c     This file selects the concentration routine to be used in the run,
c     in the provided as well as user-developed 'ab' 
c     routines
c
c     The following models are currently supported:
c     abConst
c     abCase1
c     abCase2
c     abacbb
c     abnewcase1
c
      INCLUDE "DIMENS_XL.INC"
c
      real depth, conc, dcompconc(mxcomp)
      INTEGER icomp
c
      common /Cconc/  itype(mxcomp)
      COMMON /CConstConc/ compconc(mxcomp)
      INTEGER iparam
      COMMON /CmodelOpts/ iparam(0:5)
      Common /Castar/ iastropt, astarRef, astar0, asgamma
      COMMON /Cmisc/ imisc(30),fmisc(30)
c
c     Specific ab variables
      Integer iastropt(mxcomp+2), itype
      Real astarRef(mxcomp+2), astar0(mxcomp+2)
      Real asgamma(mxcomp+2)
c
      wavenm = fmisc(13)  !some CDOM models need wavelength
      ncomp = imisc(6)
      icompchl = iparam(4)
      icompcdom = iparam(5)
c
      IF (iparam(0).eq.1) THEN
!       CALL USER model instead if HL was recompiled for this run
!       (this occurs if using user-def IOP routine, biolum incl, or
!       concentrations set in user-defined subroutines
        call user_abmodel(depth, dcompconc)
      ELSE
c     Component 1 is always pure water in the STANDARD IOP models.
c     No calc needed.
c
c     Define other components as needed
       DO icomp=2, ncomp
        If (icomp.eq.icompcdom .And. iastropt(icomp) .eq. -2) Then
c         absorption due to CDOM will be calculated within routine acdomsub using
c         acdom = 0.2*(achl(440))*exp(-.014(wave-440))
c         [return a value of one here which will not be used]
          call acdomSub(depth, wavenm, conc)
          dcompconc(icomp) = conc
        ElseIf (itype(icomp) .eq. 0) Then
c         absorption for this component is constant with depth
          dcompconc(icomp) =  compconc(icomp)
        ELSEIF (itype(icomp) .EQ. 2) then  !read datafile
c         this block only used if case1, NEWcase1 or case2 IOP model used
          IF(icomp.eq.icompchl) then
            call chlzdata(icomp, depth, dcompconc(icomp))
          elseIF(icomp.eq.icompcdom) then
            call cdomdata(icomp, depth, dcompconc(icomp))
          elseIF(icomp.eq. 4) then
            call minzdata(icomp, depth, dcompconc(icomp))
          else
            WRITE(10,*)'stop in abmodel for comp#: ', icomp
            call HERR("abmodel","concentration routine not specified")  !stop run
          endif
        ELSE
          WRITE(10,*)'stop in abmodel for comp#: ', icomp
          call HERR("abmodel","OPTION not supported")  !stop run
        END IF
       END DO
      END IF
      return
      end
c     -----------------------------------------

      subroutine abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
c
c     This file selects the 'ab' routine to be used in the run,
c     called by routines inishamp, rhotau, zetatoz, and ztozeta.
c
      INCLUDE "DIMENS_XL.INC"
c
      integer ncomp
      real depth,wavelen,atotal,btotal
      real acomp(mxcomp),bcomp(mxcomp)
c
      INTEGER iparam
      COMMON /CmodelOpts/ iparam(0:5)
c
c     iparam(1) indicates which IOP model has been specified
c
      IF (iparam(1).eq.0) THEN
        call abconst(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
      ELSEIF (iparam(1).eq.1) THEN
        call abcase1(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
      ELSEIF (iparam(1).eq.2) THEN
        call abcase2(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
      ELSEIF (iparam(1).eq.3) THEN
        call abacbb(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
      ELSEIF (iparam(1).eq.4) THEN
        call abnewcase1(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
      ELSE
        call user_abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
      END IF
c
      return
      end
c     -----------------------------------------

      subroutine acdomsub(depth, wavelen, abscdom)
c
c     This file specifies the 'aCDOM' to be used in the run,
c     by calls from the 'ab' routines, and the CDOM fluorescence routine 'shatcdom'
c
      INCLUDE "DIMENS_XL.INC"
c
      INTEGER iparam
      COMMON /CmodelOpts/ iparam(0:5)
c
c     iparam(1) indicates which IOP model has been specified
c
      integer ncomp, icomp
      real depth,wavelen, abscdom, cdomconc, atotal,btotal
      real acomp(mxcomp),bcomp(mxcomp), achlRef, waveRef
c
c     flag for calculating bb shared with inishamp and abac9 routines
c     (calc bb iff flag=1)"
      Common /CbbCalc/ ibbCalc
      Common /Castar/ iastropt, astarRef, astar0, asgamma
      Common /Cabac9n/ nac9files
c
      common /Cconc/  itype(mxcomp)
      COMMON /CConstConc/ compconc(mxcomp)
c
c     Specific ab variables
      Integer iastropt(mxcomp+2), ibbCalc, itype,  nac9files, kall
      Real astarRef(mxcomp+2), astar0(mxcomp+2)
      Real asgamma(mxcomp+2)
c
      data kall/0/
      save kall

      icompcdom = iparam(5)

      IF(icompcdom.eq.0 .and. iparam(1).ne.1) Then
c       CDOM is not included in this run
        abscdom = 0
      Else
        IF (kall.eq.0) THEN  !printout messages
          IF(itype(icompcdom).eq.0) then		!constant with depth
            write(10,301) compconc(icompcdom)
          Elseif(iastropt(icompcdom).eq.-2) then
            write(10,303)
          Elseif(itype(icompcdom).eq.1) then
            write(10,304)
          Endif
          kall = 1
          If (iparam(1).eq.4 .or. (iparam(1).eq.3
     1        .And. nAC9files.eq.2)) return  !abnewcase1 and abacbb call once for msg
        END IF
        If ( iparam(1).eq.4 .or.                           !abnewcase1
     1      (iparam(1).eq.3 .And. nAC9files.eq.2) ) Then     !we have aCDOM via filtered ac9 data
c         aCDOM will be taken from the filtered ac-9 datafile
c         or from the new CASE 1 model, as appropriate
          ncomp = 3
          ibbCalc = 0
          call abscat(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
          abscdom = acomp(icompcdom)
        ElseIf (iparam(1).eq.1) Then    !abcase1 (classic)
c         aCDOM is given in terms of achl at 440nm
          waveRef = 440.0
          call  achlz(depth, waveRef, achlREF)
          abscdom = 0.2*achlREF*exp(-0.014*(wavelen - 440.))
        ElseIf (iastropt(icompcdom) .eq. -2) Then
c         aCDOM is given in terms of achl at a Reference wavelength"
          waveRef =  astarRef(icompcdom)
          call  achlz(depth, waveRef, achlREF)
          abscdom = astar0(icompcdom)
     1            * achlRef * exp( - asgamma(icompcdom)
     2            * (wavelen - waveRef ) )
        ELSEIF(itype(icompcdom) .eq. 1) then
c         Call user defined subroutine to give aCDOM as a function
c          of depth and wavelength
          call user_acdomsub(depth, wavelen, abscdom)
        Else
c         Get concentration of CDOM
          If (itype(icompcdom) .eq. 0) Then
c           CDOM concentration is constant with depth
            cdomconc = compconc(icompcdom)
          ELSEIF(itype(icompcdom) .eq. 2) then
c           Call subroutine to get aCDOM(z,waven) from datafile
            call cdomdata(icompcdom, depth, cdomconc)
          End If
c         Use function ASTAR to get wavelength dependence of aCDOM
          icomp = icompcdom
          abscdom = cdomconc * astar(icomp, wavelen)
        End If
      End If
c
      return
 301  format(/5x,'The CDOM absorption at each wavelength is ',
     1'constant with depth with a value of ',F8.2, ' (1/m)' )
 303  format(/5x,'CDOM absorption will covary with chlorophyll ',
     1       'according to a Case 1 model')
 304  format(/5x,'A subroutine will be called to provide aCDOM as',
     1       ' a function of depth and wavelength')
      end
c     -----------------------------------------

      subroutine achlz(depth, wavelen, achl)
c
c     This file selects how the chlorophyll profile is specified for the run,
c     for use by the 'ab' routines and the chlorophyll fluorescence
c     routine 'shatchl'
c
      INCLUDE "DIMENS_XL.INC"
c
      INTEGER iparam
      COMMON /CmodelOpts/ iparam(0:5)
      real acomp(mxcomp),bcomp(mxcomp), atotal,btotal
      real depth, wavelen, achl, chlconc, ac440, waveRef, ac
c
c     flag for calculating bb shared with inishamp and abac9 routines (calc bb iff flag=1)"
      Common /CbbCalc/ ibbCalc
      Common /Castar/ iastropt, astarRef, astar0, asgamma
      Common /Cabac9n/ nac9files
      common /Cconc/  itype(mxcomp)
      COMMON /CConstConc/ compconc(mxcomp)
c     Specific ab variables
      Integer iastropt(mxcomp+2), ibbCalc, itype,  nac9files, kall
      Real astarRef(mxcomp+2), astar0(mxcomp+2)
      Real asgamma(mxcomp+2)
c
      data kall/0/
      save kall

      icompchl = iparam(4)

      If (icompchl .eq. 0) Then
c       Chlorophyll is not included in this run
        achl = 0
      ElseIf (itype(icompchl) .eq. 1) Then
c       a subroutine will be called to give chl concentration
        call user_achlz(depth, wavelen, achl)
      ElseIF (iparam(1).eq.4 .and. kall.ne.0) THEN
        call abnewcase1(depth,wavelen,ncomp,acomp,bcomp,atotal,btotal)
        achl = acomp(2)
      Else
!       calc astar
        If (iastropt(icompchl) .eq. 3) Then !Prieur-Satheydranath achl
c         Specify component index for chorophyll for astar routines
          If(kall.eq.0) then
            waveRef = 440.0
            ac440 = astar(icompchl, waveRef)
          Endif
        End If
!
!       Set concentration information
        If (itype(icompchl) .eq. 0) Then
c         Chlorophyll is constant with depth
          chlconc = compconc(icompchl)
          IF(kall.eq.0) then
            WRITE(10,201) chlconc
          Endif
        Else
c         a subroutine will be called to read chl concentration
c         from the data file " & datafile(icompchl)
          call chlzdata(icompchl, depth, chlconc)
        End If
c       Set a* call and calculate a_chl
        If (iastropt(icompchl) .eq. 3) Then !Prieur-Satheydranath achl
          ac = astar(icompchl, wavelen)/ac440
          achl = 0.06*ac*chlconc**(0.65)
        Else
          achl = chlconc*astar(icompchl, wavelen)
        End If
      End If
c
      kall = 1
      return
 201  format(/5x,'The chlorophyll concentration is ',
     1'constant with depth with a value of ',F8.2, ' (mg/m^3)')
      end
c     -----------------------------------------

      subroutine skyirrad(suntheta,sunphi,eddif,eddir)
c
c     This file selects the sky irradiance model to be used in the run,
c     and is called by routine qasky.
c
      real suntheta,sunphi,eddif,eddir
c
      INCLUDE "DIMENS_XL.INC"
c
      INTEGER iparam
      COMMON /CmodelOpts/ iparam(0:5)
c
      COMMON /Cmisc/ imisc(30),fmisc(30)
      Character surfname*120,pfname*120,
     1           Drootname*120,Srootname*120,Mrootname*120,
     2           datafiles*120
      COMMON /Cfilenames/ surfname,pfname(mxcomp),
     1                    Drootname,Srootname,Mrootname,
     2                    datafiles(0:7+mxcomp)
c
      data kall/0/
      save kall
c
      ncomp = imisc(6)  !used iff irrad from datafile
c
      IF (iparam(2) .eq. 0) THEN
        call cosirrad(suntheta,sunphi,eddif,eddir)
      ELSE
        IF(iparam(3).EQ.-1) THEN
c         call irradat to read the  datafile for {Eddir, Eddif}
          itype = 2
          call irradat(itype, eddir, eddif, datafiles(ncomp + 7) )
        ELSEIF(iparam(3).EQ.-2) THEN
c         call irradat to read the  datafile for {Edtot, fracdir}
          itype = -2
          call irradat(itype, Edtot, eddir, datafiles(ncomp + 7) )
          eddif = MAX(Edtot * (1.-eddir), 0.0)
          eddir = MAX(Edtot * eddir,0.0)
        ELSEIF(iparam(3).eq.1) THEN
c         call irradat to read the  datafile for {Edtot}
          itype = 0
          call irradat(itype, Edtot, eddir, datafiles(ncomp + 7) )
c         Routine RADTRANX provides the RADTRAN direct:diffuse ratio
          write(10,531)
          call RADTRANX(suntheta,sunphi,eddif,eddir)
          IF(kall.eq.0) then !return after printing msgs
            kall = 1
            return
          Endif
          dirfrac = eddir/(eddir+eddif)
c         Use dirfrac from RADTRAN and the total Ed from the data file
c         to calculate the direct and diffuse parts
          eddir = dirfrac * Edtot
          eddif = Edtot   - eddir
          write(10,532) dirfrac
          write(10,533) Eddif,Eddir,Edtot
        ELSE
!         call RADTRAN semi-empirical model to get eddir and eddif
          call RADTRANX(suntheta,sunphi,eddif,eddir)
        ENDIF
      END IF

      return

  531 format(//5x,'Gregg and Carder values used to get the ratio of ',
     1'diffuse to direct irradiance'/)
  532 format(5x,'Ratio of direct to total Irradiance ',
     1'as calculated from the Gregg and Carder model:',f6.3)
  533 format(//5x,'Sky spectral irradiances used by Hydrolight for ',
     1'this run (from datafile and G&C ratio):'/10x,'Ed(diffuse) =',
     2 1pe11.3,5x,'Ed(direct) =',e11.3,5x,'Ed(total) =',e11.3)
      end
c     -----------------------------------------

      subroutine skyrad(suntheta,sunphi,skytheta,skyphi,skyrad0)
c
c     This file selects the sky radiance model to be used in the run,
c     and is called by routine qasky.
c
      real suntheta,sunphi,skytheta,skyphi,skyrad0
c
c
      INTEGER iparam
      COMMON /CmodelOpts/ iparam(0:5)
c
      IF (iparam(2) .eq. 0) THEN
        call cosrad(suntheta,sunphi,skytheta,skyphi,skyrad0)
      ELSEIF (iparam(2) .eq. 1) THEN
        call hcnrad(suntheta,sunphi,skytheta,skyphi,skyrad0)
      ELSE
        call user_skyrad(suntheta,sunphi,skytheta,skyphi,skyrad0)
      END IF
      return
      end
c     -----------------------------------------

