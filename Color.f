C     Last change: CM 11/18/2010 3:35:00 PM
      subroutine Color
c
c     core routine on file Color.f
c
c     Written by CDM, 11 Nov 04, originally named CIExyY.f
c     Modified 05 Oct 2010 to add Forel-Ule numbers to the output and renamed to Color.f
c     Modified 18 Nov 2010 to fix bug in Forel-Ule numbers for extreme CIE colors
c
c     This routine first computes the 1931 CIE chromaticity coordinates
c     (x,y) for various radiometric quantities after a HydroLight or EcoLight
c     run has been completed.
c     It is intended for use only with multi-wavelength runs
c     that span the visible spectrum from 400 to 700 nm (at least;
c     for best results, run from 350 to 780 nm).
c
c     CIE calculations are based on Light and Water, Section 2.5
c
c     The Forel-Ule number is computed using the CIE (x,y) values given for
c     the Forel-Ule colors in Wernand and van der Woerd, 2010. "
c     Spectral analysis of the Forel-Ule ocean color comparator scale."
c     J. Euro. Opt. Soc. 5, 10014S.
c     The closest Forel-Ule numbers of the nadir-viewing Rrs = Lw/Ed and
c     above-surface Lu are computed.  These are usually the same, but can
c     differ if Lu includes a significant amount of surface-reflected sun
c     or sky radiance (sun and sky glint).
c
      INCLUDE "DIMENS_XL.INC"
c
c     Cxcl holds arrays accumulated in routine storexcl
      COMMON /Cxcl/  nzxcl,zxcl(mxz),izxcl(mxz),
     1               axcl(mxz,mxwave,0:mxcomp),
     2               bxcl(mxz,mxwave,0:mxcomp),
     3               bbxcl(mxz,mxwave,0:mxcomp),
     4               Edxcl(0:mxz,mxwave),Euxcl(0:mxz,mxwave),
     5               Eodxcl(0:mxz,mxwave),Eouxcl(0:mxz,mxwave),
     5               Rxcl(0:mxz,mxwave),
     6               Raduxcl(0:mxz,mxwave),Radwxcl(mxwave)
c
      COMMON /Cwave/ wave(mxwave),waveb(mxwave+1),fijchl(mxwave,mxwave),
     1               fijcdom(mxwave,mxwave),fijraman(mxwave,mxwave)

      COMMON /CMISC/ IMISC(30),FMISC(30)

c     The CIE 1931 tristimulus functions
c     Ref: L&W Table 2.4 (typo: zbar(380) = 0.0065)
      integer, parameter :: nCIE=43
      real, parameter :: fKm = 683.0
      real :: CIEwave(nCIE),CIExbar(nCIE),CIEybar(nCIE),CIEzbar(nCIE)
c     CIE tristimulus functions for use with current wavelengths
      real, allocatable,dimension(:) :: xbarwave, ybarwave, zbarwave
c     CIE coordinates of the Forel-Ule numbers
      real :: FUx(23),FUy(23)
      integer :: FU_Rrs,FU_Lu,FU_Lw

      data CIEwave   /360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0,
     1  430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0, 500.0, 510.0,
     2  520.0, 530.0, 540.0, 550.0, 560.0, 570.0, 580.0, 590.0, 600.0,
     3  610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0, 680.0, 690.0,
     4  700.0, 710.0, 720.0, 730.0, 740.0, 750.0, 760.0, 770.0, 780.0/

      data CIExbar /  0.0001,0.0004,0.0014,0.0042,0.0143,0.0435,0.1344,
     1  0.2839,0.3483,0.3362,0.2908,0.1954,0.0956,0.0320,0.0049,0.0093,
     2  0.0633,0.1655,0.2904,0.4334,0.5945,0.7621,0.9163,1.0263,1.0622,
     3  1.0026,0.8544,0.6424,0.4479,0.2835,0.1649,0.0874,0.0467,0.0227,
     4  0.0114,0.0058,0.0028,0.0014,0.0007,0.0003,0.0002,0.0001,0.0000/

      data CIEybar /  0.0000,0.0000,0.0000,0.0001,0.0004,0.0012,0.0040,
     1  0.0116,0.0230,0.0380,0.0600,0.0910,0.1390,0.2080,0.3230,0.5030,
     2  0.7100,0.8620,0.9540,0.9950,0.9950,0.9520,0.8700,0.7570,0.6310,
     3  0.5030,0.3810,0.2650,0.1750,0.1070,0.0610,0.0320,0.0170,0.0082,
     4  0.0041,0.0021,0.0011,0.0005,0.0003,0.0001,0.0001,0.0000,0.0000/

      data CIEzbar /  0.0006,0.0019,0.0065,0.0201,0.0679,0.2074,0.6456,
     1  1.3856,1.7471,1.7721,1.6692,1.2876,0.8130,0.4652,0.2720,0.1582,
     2  0.0782,0.0422,0.0203,0.0087,0.0039,0.0021,0.0017,0.0011,0.0008,
     3  0.0003,0.0002,0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     4  0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000/

c     CIE x and y coordinates of the Forel-Ule numbers 1 to 21 (from the
c     spreadsheet of M. Wernand)
c     Phony FU =0 appended at CIE(0.176,0.005) to correct bug at very blue colors
c     Phony FU =22 appended at CIE(0.735,0.265) to correct bug at very red colors
c     The FU number is now 1 less that the FUx array index
c      data FUx /0.176,0.1887,0.1961,0.2130,0.2293,0.2422,0.2631,0.2903,
c     1  0.3108,0.3373,0.3631,0.3883,0.3944,0.3972,0.4040,0.4098,0.4181,
c     2  0.4275,0.4398,0.4527,0.4625,0.4730,0.735/
c      data FUy /0.005,0.1614,0.1941,0.2548,0.3010,0.3315,0.3735,0.4154,
c     1  0.4389,0.4630,0.4798,0.4905,0.4883,0.4864,0.4821,0.4778,0.4724,
c     2  0.4660,0.4581,0.4485,0.4404,0.4285,0.265/

c     Revised CIE x and y coordinates of the Forel-Ule numbers 1 to 21 (from the
c     16 May 2013 email of M. Wernand)
c     Phony FU =0 appended at CIE(0.176,0.005) to correct bug at very blue colors
c     Phony FU =22 appended at CIE(0.735,0.265) to correct bug at very red colors
c     The FU number is now 1 less that the FUx array index
      data FUx /0.176,0.191,0.199,0.210,0.227,0.246,0.266,0.291,
     1  0.315,0.337,0.363,0.386,0.402,0.416,0.431,0.446,0.461,
     2  0.475,0.489,0.503,0.516,0.528,0.735/
      data FUy /0.005,0.167,0.200,0.240,0.288,0.335,0.376,0.412,
     1  0.440,0.462,0.476,0.487,0.481,0.474,0.466,0.458,0.449,
     2  0.441,0.433,0.425,0.416,0.408,0.265/

      nwave = imisc(7)
      allocate(xbarwave(nwave), ybarwave(nwave), zbarwave(nwave))
      i1 = 1
      Do i=1,nwave
        wl = wave(i)
        if(wl.le.CIEwave(1) .or. wl.gt.CIEwave(nCIE)) then
           xbarwave(i) = 0.0
           ybarwave(i) = 0.0
           zbarwave(i) = 0.0
        else
          xbarwave(i) = yinterp(i1, nCIE, wl, CIEwave, CIExbar)
          ybarwave(i) = yinterp(i1, nCIE, wl, CIEwave, CIEybar)
          zbarwave(i) = yinterp(i1, nCIE, wl, CIEwave, CIEzbar)
       endif
      Enddo
c
c     compute X,Y,Z values for various radiometric quantities, using
c     L&W Eq. (2.10)
c
c     nadir-viewing, total upwelling radiance in air
c     nadir-viewing, water-leaving radiance
c     total Eu and Ed in air
c      
      xLusum = 0.0
      yLUsum = 0.0
      zLusum = 0.0
      xLwsum = 0.0
      yLwsum = 0.0
      zLwsum = 0.0
      xEusum = 0.0
      yEusum = 0.0
      zEusum = 0.0
      xEdsum = 0.0
      yEdsum = 0.0
      zEdsum = 0.0
      xRrssum = 0.0
      yRrssum = 0.0
      zRrssum = 0.0
        do iwave=1,nwave
        deltawave = waveb(iwave+1)-waveb(iwave)

        xLusum = xLusum + raduxcl(0,iwave)*xbarwave(iwave)*deltawave
        yLusum = yLusum + raduxcl(0,iwave)*ybarwave(iwave)*deltawave
        zLusum = zLusum + raduxcl(0,iwave)*zbarwave(iwave)*deltawave

        xLwsum = xLwsum + radwxcl(iwave)*xbarwave(iwave)*deltawave
        yLwsum = yLwsum + radwxcl(iwave)*ybarwave(iwave)*deltawave
        zLwsum = zLwsum + radwxcl(iwave)*zbarwave(iwave)*deltawave

        xEusum = xEusum + Euxcl(0,iwave)*xbarwave(iwave)*deltawave
        yEusum = yEusum + Euxcl(0,iwave)*ybarwave(iwave)*deltawave
        zEusum = zEusum + Euxcl(0,iwave)*zbarwave(iwave)*deltawave

        xEdsum = xEdsum + Edxcl(0,iwave)*xbarwave(iwave)*deltawave
        yEdsum = yEdsum + Edxcl(0,iwave)*ybarwave(iwave)*deltawave
        zEdsum = zEdsum + Edxcl(0,iwave)*zbarwave(iwave)*deltawave

        xRrssum = xRrssum +
     1 (radwxcl(iwave)/Edxcl(0,iwave))*xbarwave(iwave)*deltawave
        yRrssum = yRrssum +
     1 (radwxcl(iwave)/Edxcl(0,iwave))*ybarwave(iwave)*deltawave
        zRrssum = zRrssum +
     1 (radwxcl(iwave)/Edxcl(0,iwave))*zbarwave(iwave)*deltawave

        end do

c     compute (x,y,Y) from L&W Eq. (2.11)
      XYZLu = xLusum + yLusum + zLusum
      capYLu = fKm*yLusum
      xLu = xLusum/XYZLu
      yLu = yLusum/XYZLu
      xLuair = xLu ! save for Forel-Ule calculations
      yLuair = yLu

      XYZLw = xLwsum + yLwsum + zLwsum
      capYLw = fKm*yLwsum
      xLw = xLwsum/XYZLw
      yLw = yLwsum/XYZLw

      XYZEu = xEusum + yEusum + zEusum
      capYEu = fKm*yEusum
      xEu = xEusum/XYZEu
      yEu = yEusum/XYZEu

      XYZEd = xEdsum + yEdsum + zEdsum
      capYEd = fKm*yEdsum
      xEd = xEdsum/XYZEd
      yEd = yEdsum/XYZEd

      XYZRrs = xRrssum + yRrssum + zRrssum
      xRrs = xRrssum/XYZRrs
      yRrs = yRrssum/XYZRrs

      write(10,200)
  200 format(//
     1' CIE 1931 chromaticity coordinates (x,y) and (il)luminance Y '
     2'for selected quantities.  Output is (x,y,Y).',//,5x,
     3'Y is illuminance in lm/m^2 for Ed and Eu.  ',
     4'Y is luminance in lm/(m^2 sr) for Lu and Lw.  ',
     5'Y is not defined for Rrs',//5x,
     6'  depth',13x,'Ed',25x,'Eu',25x,'Lu',25x,'Lw',24x,'Rrs',/)

      write(10,202) xEd,yEd,capYEd,xEu,yEu,capYEu,xLU,yLu,capYLu,
     1 xLw,yLw,capYLw,xRrs,yRrs,0.0
  202 format(5x,' in air',5('  (',0pf5.3,', ',0pf5.3,','1pe10.3,')'))

c     in-water, depth-dependent quantities

      do iz=1,nzxcl
      xLusum = 0.0
      yLUsum = 0.0
      zLusum = 0.0
      xEdsum = 0.0
      yEdsum = 0.0
      zEdsum = 0.0
      xEusum = 0.0
      yEusum = 0.0
      zEusum = 0.0
        do iwave=1,nwave
        deltawave = waveb(iwave+1)-waveb(iwave)

        xLusum = xLusum + raduxcl(iz,iwave)*xbarwave(iwave)*deltawave
        yLusum = yLusum + raduxcl(iz,iwave)*ybarwave(iwave)*deltawave
        zLusum = zLusum + raduxcl(iz,iwave)*zbarwave(iwave)*deltawave

        xEusum = xEusum + Euxcl(iz,iwave)*xbarwave(iwave)*deltawave
        yEusum = yEusum + Euxcl(iz,iwave)*ybarwave(iwave)*deltawave
        zEusum = zEusum + Euxcl(iz,iwave)*zbarwave(iwave)*deltawave

        xEdsum = xEdsum + Edxcl(iz,iwave)*xbarwave(iwave)*deltawave
        yEdsum = yEdsum + Edxcl(iz,iwave)*ybarwave(iwave)*deltawave
        zEdsum = zEdsum + Edxcl(iz,iwave)*zbarwave(iwave)*deltawave

        end do

c     compute (x,y,Y) from L&W Eq. (2.11)
      XYZLu = xLusum + yLusum + zLusum
      capYLu = fKm*yLusum
      xLu = xLusum/XYZLu
      yLu = yLusum/XYZLu

      XYZEu = xEusum + yEusum + zEusum
      capYEu = fKm*yEusum
      xEu = xEusum/XYZEu
      yEu = yEusum/XYZEu

      XYZEd = xEdsum + yEdsum + zEdsum
      capYEd = fKm*yEdsum
      xEd = xEdsum/XYZEd
      yEd = yEdsum/XYZEd

      write(10,204)zxcl(iz),xEd,yEd,capYEd,xEu,yEu,capYEu,xLu,yLu,capYLu
  204 format(5x,f7.2,3('  (',0pf5.3,', ',0pf5.3,','1pe10.3,')'))

      enddo

c     Compute the Forel-Ule number for CIE (x,y) = (xRrs, yRrs)

c       Must check each pair of FU points to find the closest FU number
c       to the line from the CIE white point (1/3,1/3) through the
c       (xRrs,yRrs) point.  This is done by finding the intersection of the
c       line from (1/3,1,3) through (xRrs,yRrs) and the line connecting
c       each pair of FU points.
c       Need to compute only the x coordinate of the intersection, since
c       the FU CIE(x,y) values are monotonically increasing from FU = 1 to 21.
c       However, there can be double solutions for the line intersections
c       for some waters, so search from low to high FU numbers for CIEx values
c       less than 0.3333 (left of the white point), and search from high to low
c       for CIEx greater than 0.3333.

        if(xRrs .le. 0.3333) then
          i1 = 1
          i2 = 22
          i3 = 1
        else
          i1 = 22
          i2 = 1
          i3 = -1
        endif
        
        do i=i1,i2,i3
          FUratio = (Fuy(i+1) - FUy(i))/(FUx(i+1) - FUx(i))
          x = ( (xRrs - 0.3333)*(FUy(i) - FUx(i)*FUratio) +
     1       0.3333*(yRrs - 0.3333) -0.3333*(xRrs - 0.3333) )/
     2       (yRrs - 0.3333 - (xRrs - 0.3333)*FUratio)

          ! check to see if this x is between FU points i and i+1
          if(x.gt.FUx(i) .and. x.le.FUx(i+1)) then ! correct FU pair found
            ! check to see which point is closest
            if (x - FUx(i) .le. FUx(i+1) - x) then
              FU_Rrs = i-1 ! FU number is one less that array index
            else
              FU_Rrs = i
            endif
            exit ! solution found; done with the do loop
          endif
        end do
       if(FU_Rrs .lt. 1) FU_Rrs = 1
       if(FU_Rrs .gt. 21) FU_Rrs = 21


      write(10,206) FU_Rrs
  206 Format(//' The Forel-Ule number of the nadir-viewing Rrs is FU ='
     1 ,i3)

c     Compute the Forel-Ule number for for CIE (x,y) = (xLuair, yLuair)

        if(xLuair .le. 0.3333) then
          i1 = 1
          i2 = 22
          i3 = 1
        else
          i1 = 22
          i2 = 1
          i3 = -1
        endif

        do i=i1,i2,i3
          FUratio = (Fuy(i+1) - FUy(i))/(FUx(i+1) - FUx(i))
          x = ( (xLuair - 0.3333)*(FUy(i) - FUx(i)*FUratio) +
     1       0.3333*(yLuair - 0.3333) -0.3333*(xLuair - 0.3333) )/
     2       (yLuair - 0.3333 - (xLuair - 0.3333)*FUratio)

          ! check to see if this x is between FU points i and i+1
          if(x.gt.FUx(i) .and. x.le.FUx(i+1)) then ! correct FU pair found
            ! check to see which point is closest
            if (x - FUx(i) .le. FUx(i+1) - x) then
              FU_Lu = i-1
            else
              FU_Lu = i
            endif
            exit
          endif
        end do
       if(FU_Lu .lt. 1) FU_Lu = 1
       if(FU_Lu .gt. 21) FU_Lu = 21


      write(10,208) FU_Lu
  208 Format(/
     1' The Forel-Ule number of the nadir-viewing, in-air Lu is FU =',
     2 i3)

c     Compute the Forel-Ule number for for CIE (x,y) = (xLw, yLw)

        if(xLw .le. 0.3333) then
          i1 = 1
          i2 = 22
          i3 = 1
        else
          i1 = 22
          i2 = 1
          i3 = -1
        endif

        do i=i1,i2,i3
          FUratio = (Fuy(i+1) - FUy(i))/(FUx(i+1) - FUx(i))
          x = ( (xLw - 0.3333)*(FUy(i) - FUx(i)*FUratio) +
     1       0.3333*(yLw - 0.3333) -0.3333*(xLw - 0.3333) )/
     2       (yLw - 0.3333 - (xLw - 0.3333)*FUratio)

          ! check to see if this x is between FU points i and i+1
          if(x.gt.FUx(i) .and. x.le.FUx(i+1)) then ! correct FU pair found
            ! check to see which point is closest
            if (x - FUx(i) .le. FUx(i+1) - x) then
              FU_Lw = i-1
            else
              FU_Lw = i
            endif
            exit
          endif
        end do
       if(FU_Lw .lt. 1) FU_Lw = 1
       if(FU_Lw .gt. 21) FU_Lw = 21

      write(10,210) FU_Lw
  210 Format(/
     1' The Forel-Ule number of the nadir-viewing Lw is FU =',
     2 i3)

      return
      end


