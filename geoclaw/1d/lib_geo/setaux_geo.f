c     ============================================
      subroutine setaux(maxmx,mbc,mx,xlower,dx,maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # aux(i,1) = b in cell i, bottom topography "hillinbasin"
c
c     
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc, *)

      ! set user-defined topography as aux(i,1)
      do 50 i=1,mx 
         xcell = xlower + (i-0.5d0)*dx
         aux(i,1)=1.8d0*dexp(-1000.d0*xcell**2)
     &         + .5d0*1.0d0*dexp(-50.d0*(xcell-0.5d0)**2)
     &         + 1.8d0*dexp(-1000.d0*(xcell-1.d0)**2)


         write(17,701) aux(i,1)
 701       format(e22.12)


 50   continue
     


      return
      end
