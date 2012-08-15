      subroutine novaccum1(fp,fm,density,dtdx,meqn,meqnlim,tol)

       !novaccum prevents the appearance of negative states
       !by limiting the interface fluxes, fp,fm
       !at the right, left, of a cell
       !mass is actually the density that you will allow to be lost
       ! ie. q(i,1) or some percentage of that.

       ! meqnlim: limit components, 1 to meqnlim


      implicit none

      !passed in
      integer meqn, meqnlim
      double precision density,dtdx
      double precision fp(meqn)
      double precision fm(meqn)
      double precision tol

       !locals
      integer m
      double precision outdensity,phi

       !limit interfaces fluxes to preserve positivity. Consider qnew, not qold in limiting
       ! test if relimiting is necessary for cell i
     
      outdensity =dtdx*(dmax1(0.d0,fp(1))-dmin1(0.d0,fm(1)))

      if (outdensity.gt.0.d0) then
         ! find ratio of mass to total loss of mass. 
         ! if ratio is less than one, flux must be limited
         phi = dmin1(1.d0,density/outdensity) 
         phi = dmax1(0.d0,phi)
	else ! cell is safe
	   phi = 1.d0
      endif

      ! limit fluxes in case cell is not safe
      if (phi.lt.1.d0) then
c         if (phi.lt.tol) phi = 0.d0 
         if (fp(1).gt.0.d0) then
            do m=1,meqnlim
               fp(m)=phi*fp(m)
            enddo
         endif
                         
         if (fm(1).lt.0.d0) then
            do m=1,meqnlim
               fm(m)=phi*fm(m)
            enddo 
         endif

      endif

      return
      end subroutine 
