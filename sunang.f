C     Last change:  LKS  11 Aug 2008    5:18 pm
c**********************************************************************
c
      subroutine sunang(iday,hr,rad,xlon,ylat,sunz,suna)
c
c     Called by qasky_XL
c
c     Computes sun azimuth and zenith angles for a given
c     Julian day, GMT, latitude, and longitude.  This program
c     is from the NMFS ELAS computer code.  Modified for
c     standard coordinates (W long. negative), to correct
c     for dateline problem, and to correct coefficients (taken
c     from Iqbal, 1983, An Introduction to Solar Radiation).
c     Watson Gregg, Research and Data Systems, Corp.
c     Further modified to remove unused input (CDM, March 98)
c
c     INPUT:
c     iday = day of year ('julian' day)
c     hr = GMT in hours (e.g. 21.5 for 21 hours, 30 minutes GMT)
c     rad = radian to degrees conversion factor
c     ylat = latitude of pixel
c     xlon = longitude of pixel
c     OUTPUT:
c     sunz = solar zenith angle in degrees
c     suna = solar azimuth angle in degrees
c
c     internal variables:
c     sdec = solar declination angle in degrees
c     thez = theta zero orbital position in degrees
c     tc = time correction
c     xha = solar hour angle in degrees
c
c  Compute solar declination angle
      thez = 360.0*(iday-1)/365.0
      rthez = thez/rad
      sdec = 0.396372-22.91327*cos(rthez) + 4.02543*sin(rthez)
     *   - 0.387205*cos(2.0*rthez) + 0.051967*sin(2.0*rthez)
     *   - 0.154527*cos(3.0*rthez) + 0.084798*sin(3.0*rthez)
      rsdec = sdec/rad
c
c  Time correction for solar hour angle, and solar hour angle
      tc = 0.004297 + 0.107029*cos(rthez) - 1.837877*sin(rthez)
     *   - 0.837378*cos(2.0*rthez) - 2.342824*sin(2.0*rthez)
      xha = (hr-12.0)*15.0 + xlon + tc
      if (xha .gt.  180.0)xha = xha-360.0
      if (xha .lt. -180.0)xha = xha+360.0
      rlat = ylat/rad
c      rlon = xlon/rad
      rha = xha/rad
c
c  Sun zenith
      costmp = sin(rlat)*sin(rsdec) +
     *       cos(rlat)*cos(rsdec)*cos(rha)
      eps = abs(costmp)
      if (eps .gt. 1.1)then
       write(10,*)'Error in acos argument in sun zenith'
       write(10,*) costmp
      else if (eps .gt. 1.0)then
       if (costmp .gt. 0.0)costmp=1.0
       if (costmp .lt. 0.0)costmp=-1.0
      endif
      rsunz = acos(costmp)
c  Sun azimuth
      sna=cos(rsdec)*sin(rha)
      csa=sin(rlat)*cos(rha)*cos(rsdec)-sin(rsdec)*cos(rlat)
      rsuna=(atan2(sna,csa)+4.0*atan(1.0))
c
c  Convert to degrees
      suna = rsuna * rad
      sunz = rsunz*rad
c
      return
      end

