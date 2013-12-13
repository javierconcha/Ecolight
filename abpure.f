C     Last change:  LKS  31 Jan 2008    8:10 pm
      subroutine abpure(z,wavenm,ncomp, acomp,bcomp,atotal,btotal)
c
c     user routine on file abpure.f     
c
c     called by INISHAMP, BEAMc1, BEAMcz, and RHOTAU
c
c     This routine returns the absorption and scattering
c     coefficients (in 1/meter) for pure water.  Routine pureH20.f
c     is called to obtain the a and b values.
c
      dimension acomp(*),bcomp(*)
      common/cmisc/ imisc(30),fmisc(30)
c
      data kall/0/
c
      save
c
      if(kall.eq.0) then
         write(10,100)
         if(imisc(6).ne.1) then
            write(10,106) imisc(6)
            call HERR("abpure","# components is not compatible")  !stop run
            endif
        kall = 1 
      endif
c
      call pureH2O(wavenm, aw,bw,cw)
c
      acomp(1) = aw
      bcomp(1) = bw
      atotal = acomp(1)
      btotal = bcomp(1)
c
      return
c
 100  format(/5x,'The IOP routine "abpure" is being used:'//5x,
     1'Absorption and scattering coefficients for pure sea water are',
     2' used.')
 106  format(//'Error in sub abpure:  ncomp =',
     1i3,' (should = 1 for pure water)')

      end









