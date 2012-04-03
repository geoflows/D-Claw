c     ============================================
      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
c     ============================================
c
c     # set auxiliary arrays
c
      implicit none
      double precision aux(1-mbc:maxmx+mbc, *)

      integer maxmx,mbc,mx,maux
      double precision xlower,dx

      integer m,i
      double precision xcell,theta_init,theta_init_flume,topo_init
      double precision phi_init

      include "digparamsdec.i"
      include "digparamscommon.i"

      do i=1-mbc,mx+mbc
         do m=1,maux
            aux(i,m)=0.d0
         enddo
         xcell = xlower + (i-0.5d0)*dx

         if (init_htype.eq.0) then
             aux(i,i_theta) = theta_init_flume(xcell)
             aux(i,i_topo) = 0.d0
             aux(i,i_phi) =  phi_bed
         else
             aux(i,i_theta) = theta_init(xcell)
             aux(i,i_topo) = topo_init(xcell)
             aux(i,i_phi) = phi_init(xcell,phi_bed,phi_int)
         endif
      enddo

      return
      end

c=======================================================================
      function theta_init_flume(xcell)
c
c     function to return the flume geometry, given x, for the flume hopper
c     uses parameters set in digparams
c=======================================================================
      implicit none

      double precision theta_init_flume
      double precision D2,theta,xcell

      include "digparamsdec.i"
      include "digparamscommon.i"

      D2 = flumelen + flumerad*(theta1-theta2)

      if (xcell.lt.flumelen) then
          theta=theta1
      elseif (xcell.gt.D2) then
          theta=theta2
      else
          theta=theta1-(xcell-flumelen)/flumerad
      endif
      theta_init_flume = theta
      return
      end
