C     Last change:  LKS   7 Aug 2008    7:25 pm
      Real Function Surfn(wl, S, T)
!     this function returns the index of refraction of sea water using Eq. 3 of
!     Quan and Fry, Appl Optics 34(18), 3477-3480.
!     T is in deg C; S is in PSU (ppt); wl is wavelength in nm
c
      real wl, S, T      !in
c
c     temporary parameters
      real n0, n1, n2, n3, n4, n5, n6, n7, n8, n9
c
      n0 = 1.31405
      n1 = 1.779e-4
      n2 = -1.05e-6
      n3 = 1.6e-8
      n4 = -2.02e-6
      n5 = 15.868
      n6 = 0.01155
      n7 = -0.00423
      n8 = -4382.0
      n9 = 1.1455e6

      Surfn = n0 + (n1 + n2*T + n3*T**2)*S + n4*T**2 + 
     1       (n5 + n6*S + n7*T)/wl + n8*wl**(-2) + n9*wl**(-3)

      return
      end function
