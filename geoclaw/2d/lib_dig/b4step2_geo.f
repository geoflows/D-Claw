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

      use geoclaw_module
      use topo_module
      use dtopo_module
      use digclaw_module

      implicit double precision (a-h,o-z)

      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)


c=====================Parameters===========================================
      gmod = grav

c     # check for NANs in solution:

      call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,1)

c     # check for h < 0 and reset to zero
c     # check for h < drytolerance
c     # set hu = hv = 0 in all these cells

      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
           theta = 0.d0
           if (bed_normal.eq.1) theta=aux(i,j,i_theta)
           call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),q(i,j,5)
     &                     ,u,v,sv,theta)
           q(i,j,6) = min(q(i,j,6),q(i,j,1))
           q(i,j,6) = max(q(i,j,6),0.0)
        enddo
      enddo

c      write(26,*) 'B4STEP2: t, num_dtopo: ', t,num_dtopo
      do i=1,num_dtopo
          call movetopo(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,t,dt,maux,aux,
     &      dtopowork(i0dtopo(i):i0dtopo(i)+mdtopo(i)-1),
     &      xlowdtopo(i),ylowdtopo(i),xhidtopo(i),yhidtopo(i),
     &      t0dtopo(i),tfdtopo(i),dxdtopo(i),dydtopo(i),dtdtopo(i),
     &      mxdtopo(i),mydtopo(i),mtdtopo(i),mdtopo(i),
     &      minleveldtopo(i),maxleveldtopo(i),topoaltered(i))
      enddo

c======find factor of safety ratios===================================
      call calc_taudir(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &                     q,maux,aux)


c=============mobilization =============================================
      !write(*,*) 'b4step2, t:',t,init_pmin_ratio
      if (p_initialized.eq.1) then
         return
      endif

      if (init_ptype.lt.1) then
         return
      endif

      select case (init_ptype)
         !1,2 set to failure pressure in qinit

         case(3:4)
            !raise pressure
            if (t.gt.max(init_ptf,init_ptf2)) then
               p_initialized = 1
               return
            endif

            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                  if (q(i,j,1).le.drytolerance) cycle
                  theta = 0.d0
                  p_ratioij = init_pmin_ratio
                  if (bed_normal.eq.1) then
                     theta=aux(i,j,i_theta)
                     gmod = grav*dcos(theta)
                     p_ratioij = init_pmin_ratio*
     &                 (1.0 + aux(i,j,1)/q(i,j,1)) -aux(i,j,1)/q(i,j,1)
                  endif
                  call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),
     &                           q(i,j,4),q(i,j,5),u,v,sv,theta)

                  rho = sv*rho_s + (1.0-sv)*rho_f
                  pfail = p_ratioij*rho*gmod*q(i,j,1)
                  q(i,j,5) = pfail - abs(pfail)*(init_ptf - t)/init_ptf

               enddo
            enddo

      end select


      return
      end
