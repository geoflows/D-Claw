c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlow,ylow,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays
c
c     # aux(i,j,1) = Z(x,y) topography
c     #                     (negative below sea level for topoymetry)
c
c     # If icoordsys=2 then lat-lon coordinates on the sphere and
c     #    aux(i,j,2) = area ratio (capacity function -- set mcapa = 2)
c     #    aux(i,j,3) = length ratio for edge
c

      use geoclaw_module
      use topo_module
      use digclaw_module
      use auxinit_module

      implicit double precision (a-h,o-z)

      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      include "call.i"

      if (icoordsys.eq.2) then
         if (mcapa .ne. 2 .or. maux.lt.3) then
            write(6,*) 'ERROR in setaux:  for icoordsys=2'
            write(6,*) '     need mcapa = 2 and maux >= 3'
            write(6,*) '     have mcapa = ',mcapa,'  maux = ',maux
            stop
            endif
         endif

      do j=1-mbc,my+mbc
         ycell = ylow +(j-0.5d0)*dy
         yjm = ylow +(j-1.d0)*dy
         yjp = ylow + j*dy

         do i=1-mbc,mx+mbc
            xcell= xlow + (i- 0.5d0)*dx
            xim = xlow + (i - 1.d0)*dx
            xip = xlow + i*dx

            if (icoordsys.eq.2) then
c           # for lat-lon grid on sphere:
               deg2rad = pi/180.d0
               aux(i,j,2)= deg2rad*Rearth**2*
     &               (sin(yjp*deg2rad)-sin(yjm*deg2rad))/dy
               aux(i,j,3)= yjm*deg2rad
            else
               aux(i,j,2) = 1.d0
               aux(i,j,3) = 1.d0
            endif

            if (mtopofiles.gt.0) then
               topoint=0.d0
               call cellgridintegrate(topoint,xim,xcell,xip,yjm,ycell,
     &           yjp,xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo,
     &           mxtopo,mytopo,mtopo,i0topo,mtopoorder,
     &           mtopofiles,mtoposize,topowork)
               aux(i,j,1) = topoint/(dx*dy*aux(i,j,2))

            else
               aux(i,j,1) = 0.d0
c               # or set-up your own topo
            endif
c        #------- zero aux variables that will be set by files
            do mf = 1,mauxinitfiles
               aux(i,j,iauxinit(mf)) = 0.d0
            enddo

         enddo
      enddo


c     --------------integrate auxinit files if they exist---------------
      xhigher = xlow + (mx-0.5d0)*dx
      yhigher = ylow + (my-0.5d0)*dy

      do mf =1,mauxinitfiles
         if ((xlow.le.xhiauxinit(mf).and.
     &                  xhigher.ge.xlowauxinit(mf)).and.
     &                  (ylow.le.yhiauxinit(mf).and.
     &                  yhigher.ge.ylowauxinit(mf))) then

            xintlow = dmax1(xlow,xlowauxinit(mf))
            xinthi  = dmin1(xhigher,xhiauxinit(mf))
            istart  = min(1,int(0.5 + (xintlow-xlow)/dx))
            iend    = max(mx,int(1.0 + (xinthi-xlow)/dx))

            yintlow = dmax1(ylow,ylowauxinit(mf))
            yinthi  = dmin1(yhigher,yhiauxinit(mf))
            jstart  = int(0.5 + (yintlow-ylow)/dy)
            jend    = max(my,int(1.0 + (yinthi-ylow)/dy))

            do i=istart,iend
               x = xlow + (i-0.5d0)*dx
               xim = x - 0.5d0*dx
               xip = x + 0.5d0*dx
               do j=jstart,jend
                  y = ylow + (j-0.5d0)*dy
                  yjm = y - 0.5d0*dy
                  yjp = y + 0.5d0*dy

                  if (xip.gt.xlowauxinit(mf)
     &                     .and.xim.lt.xhiauxinit(mf)
     &                     .and.yjp.gt.ylowauxinit(mf)
     &                     .and.yjm.lt.yhiauxinit(mf)) then

                     xipc=dmin1(xip,xhiauxinit(mf))
                     ximc=dmax1(xim,xlowauxinit(mf))
                     xc=0.5d0*(xipc+ximc)

                     yjpc=dmin1(yjp,yhiauxinit(mf))
                     yjmc=dmax1(yjm,ylowauxinit(mf))
                     yc=0.5d0*(yjmc+yjpc)

                     daux =topointegral(ximc,xc,xipc,yjmc,yc,yjpc,
     &                        xlowauxinit(mf),ylowauxinit(mf),
     &                        dxauxinit(mf),dyauxinit(mf),
     &                        mxauxinit(mf),myauxinit(mf),
     &                        auxinitwork
     &                      (i0auxinit(mf):i0auxinit(mf)+mauxinit(mf)-1)
     &                        ,1)
                     daux=daux/((xipc-ximc)*(yjpc-yjmc)*aux(i,j,2))
                     aux(i,j,iauxinit(mf)) = aux(i,j,iauxinit(mf))+daux
               !write(*,*) 'theta,i,mx,j,my',aux(i,j,i_theta),i,mx,j,my
                  endif
               enddo
            enddo
         endif
      enddo

      do mf = 1,mauxinitfiles
         if (iauxinit(mf).eq.i_phi) return
      enddo

      do j=1-mbc,my+mbc
         do i=1-mbc,mx+mbc
            aux(i,j,i_phi) = phi_bed
         enddo
      enddo


      return

c     -----------------------------------------------------------------
c     # output aux array for debugging:
      open(23, file='fort.aux',status='unknown',form='formatted')
      write(23,*) 'Setting aux arrays'
      write(23,*) ' '
      do i=1,mx
        do j=1,my
           write(23,231) i,j,(aux(i,j,m),m=1,maux)
           enddo
        enddo
 231  format(2i4,4d15.3)
      close(23)

      return
      end
