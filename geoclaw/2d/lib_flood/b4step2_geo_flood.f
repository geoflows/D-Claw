c     ============================================
      subroutine b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,t,dt,maux,aux)
c     ============================================
c
c     # called before each call to step
c     # use to set time-dependent aux arrays or perform other tasks.
c
c     This particular routine sets negative values of q(i,j,1) to zero, 
c     as well as the corresponding q(i,j,m) for m=1,meqn.  
c     This is for problems where q(i,j,1) is a depth.
c     This should occur only because of rounding error.

c     Also calls movetopo if topography might be moving.

      implicit double precision (a-h,o-z)

      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)


c=====================Parameters===========================================
      include "geo.i"
      include "topo.i"


c     # check for NANs in solution:
      call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,1)

c     # check for h < 0 and reset to zero
c     # check for h < drytolerance and reset to zero, also adjusting bathy
c     # also set hu = hv = 0 in all these cells

      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
          if (q(i,j,1).lt.drytolerance) then
             do m=1,meqn
                q(i,j,m)=0.d0
                enddo
c             if (q(i,j,1).lt.0.d0) then
c                q(i,j,1)=0.d0
c                endif 
             endif
        enddo
      enddo


      write(26,*) 'B4STEP2: t, mdtopo: ', t,mdtopo
      if (mdtopo.gt.0) then
          call movetopo(maxmx,maxmy,mbc,mx,my,
     &         xlower,ylower,dx,dy,t,dt,maux,aux,dtopowork)
       endif

      return
      end
