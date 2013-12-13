C     Last change:  LKS  31 Jan 2008    8:04 pm
      subroutine abconst(z,wavel,ncomp, acomp,bcomp,atotal,btotal)
c
c     user routine on file abconst.f     
c
c     called by INISHAMP, BEAMc1, BEAMcz, and RHOTAU
c
c     This routine returns the same absorption and scattering
c     coefficients (in 1/meter) for all depths and wavelengths.
C     Depth z can be either optical or geometric depth.
c     This routine is useful for single-wavelength runs in which
c     the a and b values are specified via "areset" and "breset".
c     The default is to return a and b for pure water at 500 nm.
c  
      dimension acomp(*),bcomp(*)
      common/cmisc/ imisc(30),fmisc(30)
c
      save kall,abs,scat
      data kall/0/
c
      if(kall.eq.0) then
         write(10,100)
c         ncomp = 1
         areset = fmisc(8)
         breset = fmisc(9)
         if(areset.lt.0.0) then
            write(10,103) areset
         call HERR("abconst","absorption negative")  !stop run
         else
            abs = areset
            write(10,104) areset
         endif
         if(breset.lt.0.0) then
            write(10,105) breset
         call HERR("abconst","scattering negative")  !stop run
         else
            scat = breset
            write(10,106) breset
         endif
         if(ncomp.ne.imisc(6)) then
            write(10,102) imisc(6)
         call HERR("abconst","# components is incompatable ")  !stop run
         endif
c--------- DPF printout
         CALL PNTDPF(1)
         kall = 1
      endif
c
      acomp(1) = abs
      bcomp(1) = scat
      atotal = acomp(1)
      btotal = bcomp(1)
c
      return
c
  100 format(/5x,'The IOP routine "abconst" is being used.'//5x,
     1'Constant absorption and scattering coefficients are used:')
  102 format(//'Error in subroutine "abconst":  ncomp =',
     1i3,' (ncomp must = 1)')
  103 format(//'Error in input:  absorption coef = ',1pe12.3)
  104 format(/'     The absorption coefficient is a =',f8.4,' 1/m')
  105 format(//'Error in input:  scattering coef = ',1pe12.3)
  106 format(/'     The scattering coefficient is b =',f8.4,' 1/m')    

      end









