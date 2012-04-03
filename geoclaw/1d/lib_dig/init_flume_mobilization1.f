c Dummy routine
c supply functions for initialization
c eg: z=hinit(x), z=huinit(x) etc.

c=======================================================================
      function hinit(x)
c
c     function to return the depth, given x
c     this is for the flume mobilization/failure experiments
c=======================================================================

      implicit none

      double precision x,hinit,h,deg2rad
      double precision alpha1,alpha2,alpha3,Depth
      double precision D1,D2,D3,D4,x0,x1,x2,x3,x4

      deg2rad = 3.14159d0/180.d0
      alpha1 = 30.d0*deg2rad
      alpha2 = 10.d0*deg2rad
      alpha3 = 19.d0*deg2rad
      Depth = 0.65

      D1 = Depth/tan(alpha1)
      D2 = 4.d0
      D3 = Depth/tan(alpha2)
      D4 = Depth/tan(alpha3)

      x3 = 0.d0 !gate
      x2 = x3 - D3 ! failure surface
      x1 = x2 - D2 ! back slope
      x0 = x1 - D1 ! start of debris
      x4 = x3 + D4 ! front ramp

      if (x.ge.x3.or.x.le.x0) then
         h = 0.d0
      else
         if (x.ge.x0.and.x.lt.x1) then
            h = (x-x0)*tan(alpha1)
         elseif (x.ge.x1.and.x.lt.x2) then
            h = Depth
         elseif (x.ge.x2.and.x.lt.x3) then
            h = Depth - (x-x2)*tan(alpha2)
         endif
      endif
      hinit = h

      return
      end

c=======================================================================
      function theta_init(x)
c
c     function to return the depth, given x, for the flume hopper
c     uses parameters set in digparams
c=======================================================================

      implicit none

      double precision x,hinit,h,deg2rad
      double precision alpha1,alpha2,alpha3,Depth
      double precision D1,D2,D3,D4,x0,x1,x2,x3,x4

      deg2rad = 3.14159d0/180.d0
      alpha1 = 30.d0
      alpha2 = 10.d0
      alpha3 = 19.d0
      Depth = 0.65

      D1 = Depth/tan(alpha1)
      D2 = 4.d0
      D3 = Depth/tan(alpha2)
      D4 = Depth/tan(alpha3)

      x3 = 0.d0 !gate
      x2 = x3 - D3 ! failure surface
      x1 = x2 - D2 ! back slope
      x0 = x1 - D1 ! start of debris
      x4 = x3 + D4 ! front ramp
      if (x.ge.x4.or.x.le.x2) then
         theta = 31.d0
      else
         if (x.ge.x2.and.x.lt.x3) then
            theta = 31.d0 -alpha2
         elseif (x.ge.x3.and.x.lt.x4) then
            theta = 31.d0 + alpha3
         endif
      endif
      theta_init = theta*deg2rad

      return
      end