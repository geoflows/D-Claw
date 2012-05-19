

c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c      # Set initial sea level flat unless mqinitfiles>0, in which case
c      # an initial perturbation of the q(i,j,iqinit) is specified and has
c      # been strored in qinitwork.


      use geoclaw_module
      use qinit_module
      use digclaw_module

      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      gmod = grav

      ! depth set for sealevel. For debris flows can be used for a lake
      ! q to zero, before reading mqinit files is, so this remains.
      do i=1-mbc,mx+mbc
         x = xlower + (i-0.5d0)*dx
         do j=1-mbc,my+mbc
             y = ylower + (j-0.5d0)*dy
             q(i,j,1) = dmax1(0.d0,sealevel-aux(i,j,1))
             do m=2,meqn
               q(i,j,m)=0.d0
             enddo
         enddo
      enddo
c
      xhigher = xlower + (mx-0.5d0)*dx
      yhigher = ylower + (my-0.5d0)*dy
c
      do mf =1,mqinitfiles

         if ((xlower.le.xhiqinit(mf).and.xhigher.ge.xlowqinit(mf)).and.
     &      (ylower.le.yhiqinit(mf).and.yhigher.ge.ylowqinit(mf))) then

            xintlow = dmax1(xlower,xlowqinit(mf))
            xinthi  = dmin1(xhigher,xhiqinit(mf))
            istart  = min(1,int(0.5 + (xintlow-xlower)/dx))
            iend    = max(mx,int(1.0 + (xinthi-xlower)/dx))

            yintlow = dmax1(ylower,ylowqinit(mf))
            yinthi  = dmin1(yhigher,yhiqinit(mf))
            jstart  = int(0.5 + (yintlow-ylower)/dy)
            jend    = max(my,int(1.0 + (yinthi-ylower)/dy))

            do i=istart,iend
               x = xlower + (i-0.5d0)*dx
               xim = x - 0.5d0*dx
               xip = x + 0.5d0*dx
               do j=jstart,jend
                  y = ylower + (j-0.5d0)*dy
                  yjm = y - 0.5d0*dy
                  yjp = y + 0.5d0*dy

                  if (xip.gt.xlowqinit(mf).and.xim.lt.xhiqinit(mf)
     &               .and.yjp.gt.ylowqinit(mf)
     &               .and.yjm.lt.yhiqinit(mf)) then

                     xipc=dmin1(xip,xhiqinit(mf))
                     ximc=dmax1(xim,xlowqinit(mf))
                     xc=0.5d0*(xipc+ximc)

                     yjpc=dmin1(yjp,yhiqinit(mf))
                     yjmc=dmax1(yjm,ylowqinit(mf))
                     yc=0.5d0*(yjmc+yjpc)

                     dq = topointegral(ximc,xc,xipc,yjmc,yc,yjpc,
     &                  xlowqinit(mf),ylowqinit(mf),
     &                  dxqinit(mf),dyqinit(mf),
     &                  mxqinit(mf),myqinit(mf),
     &                   qinitwork(i0qinit(mf):i0qinit(mf)+mqinit(mf)-1)
     &                     ,1)
                     dq=dq/((xipc-ximc)*(yjpc-yjmc)*aux(i,j,2))

                     if (iqinit(mf).le.meqn) then
                        q(i,j,iqinit(mf)) = q(i,j,iqinit(mf)) + dq
                     else
                        q(i,j,1) = dmax1(dq-aux(i,j,1),0.d0)
                     endif
c
                  endif
c
               enddo
            enddo
         endif
      enddo

      do i=1-mbc,mx+mbc
         do j=1-mbc,my+mbc
            q(i,j,4) = m0*q(i,j,1)
         enddo
      enddo

c=============== Pressure initialization for Mobilization Modeling======

      if (init_ptype.gt.0) then
         call calc_pmin(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &                     q,maux,aux)
      endif

      if (init_ptype.eq.1) then
         do i=1-mbc,mx+mbc
            do j=1-mbc,my+mbc
               if (bed_normal.eq.1) gmod = grav*dcos(aux(i,j,i_theta))
               q(i,j,5) = init_pmin_ratio*rho_f*gmod*q(i,j,1)
            enddo
         enddo
         p_initialized = 1
      elseif (init_ptype.eq.0) then
         do i=1-mbc,mx+mbc
            do j=1-mbc,my+mbc
               if (bed_normal.eq.1) gmod = grav*dcos(aux(i,j,i_theta))
               q(i,j,5) = rho_f*gmod*q(i,j,1)
            enddo
         enddo
         p_initialized = 1
      endif


      return
      end
