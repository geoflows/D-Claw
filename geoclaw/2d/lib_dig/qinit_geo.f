

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
     &     (ylower.le.yhiqinit(mf).and.yhigher.ge.ylowqinit(mf))) then

            xintlow = dmax1(xlower,xlowqinit(mf))
            xinthi  = dmin1(xhigher,xhiqinit(mf))
            istart  = max(1-mbc,int(0.5 + (xintlow-xlower)/dx)-mbc)
            iend    = min(mx+mbc,int(1.0 + (xinthi-xlower)/dx)+mbc)

            yintlow = dmax1(ylower,ylowqinit(mf))
            yinthi  = dmin1(yhigher,yhiqinit(mf))
            jstart  = max(1-mbc,int(0.5 + (yintlow-ylower)/dy)-mbc)
            jend    = min(my+mbc,int(1.0 + (yinthi-ylower)/dy)+mbc)

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

      initm = 0
      do mf =1,mqinitfiles
         if (iqinit(mf).eq.4) initm = 1
      enddo

      do i=1-mbc,mx+mbc
         do j=1-mbc,my+mbc
               if (initm.eq.0) then
                  q(i,j,4) = m0*q(i,j,1)
               else
                  q(i,j,4) = q(i,j,1)*q(i,j,4)
               endif
               if (q(i,j,1).le.drytolerance) then
                  do m = 1,meqn
                     q(i,j,m) = 0.d0
                  enddo
               endif
         enddo
      enddo

c=============== Pressure initialization for Mobilization Modeling======

      select case (init_ptype)
         case (-1)
            !p should already be 0 or set by qinit file
            p_initialized = 1
         case (0)
            !set to hydrostatic
            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                 if (bed_normal.eq.1) gmod = grav*dcos(aux(i,j,i_theta))
                 q(i,j,5) = rho_f*gmod*q(i,j,1)
               enddo
            enddo
            p_initialized = 1
         case(1)

           call calc_pmin(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                    dx,dy,q,maux,aux)

            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                  p_ratioij = init_pmin_ratio
                  if (q(i,j,1).le.drytolerance) then
                     q(i,j,5) = init_pmin_ratio*rho_f*gmod*q(i,j,1)
                     cycle
                  endif
                  if (bed_normal.eq.1) then
                     gmod = grav*dcos(aux(i,j,i_theta))
                     p_ratioij = max(0.0,init_pmin_ratio
     &                   + (init_pmin_ratio - 1.0)*aux(i,j,1)/q(i,j,1))
                  endif


                 q(i,j,5) = p_ratioij*rho_f*gmod*q(i,j,1)
               enddo
            enddo
            p_initialized = 1
         case(2)
            !p will be updated in b4step2
            call calc_pmin(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                        dx,dy,q,maux,aux)
            do i=1-mbc,mx+mbc
               do j=1-mbc,my+mbc
                 p_ratioij = init_pmin_ratio
                  if (q(i,j,1).le.drytolerance) then
                     q(i,j,5) = init_pmin_ratio*rho_f*gmod*q(i,j,1)
                     cycle
                  endif
                 if (bed_normal.eq.1) then
                     gmod = grav*dcos(aux(i,j,i_theta))
                     p_ratioij = init_pmin_ratio
     &                   + (init_pmin_ratio - 1.0)*aux(i,j,1)/q(i,j,1)
                 endif
                     pfail = p_ratioij*rho_f*gmod*q(i,j,1)
                     q(i,j,5) = pfail - abs(pfail)
               enddo
            enddo
            p_initialized = 0
      end select
      !write(*,*) 'qinit:init,dx,dy:',init_pmin_ratio,dx,dy

      do i=1,mx
         do j=1,my
            call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),
     &                       q(i,j,4),q(i,j,5),u,v,sv,aux(i,j,i_theta))
         enddo
      enddo
c===============set factor of safety====================================
      call calc_fs(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,
     &                     q,maux,aux)

      return
      end
