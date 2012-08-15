c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set dambreak initial conditions for q.
c     # Includes dry bed if used with high topography
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)
c
       ! set user-defined inital conditions for q
       eta_left = .55d0
c       eta_left=.25d0
       eta_right= .15d0
       damloc=.3d0
c
       do 150 i=1,mx
	     xcell = xlower + (i-0.5d0)*dx
c             call mapc2p(xcell,xp)
              xp=xcell
             if (xp .lt. damloc) then
                 q(i,1) = eta_left - aux(i,1)
             else
                 q(i,1) = eta_right - aux(i,1)
             endif

             if (q(i,1).le.0.d0) then
                q(i,1)=0.d0
             endif
c
             q(i,2) = 0.d0
c 
  150        continue
c
      return
      end
