
c
c
c =========================================================
      subroutine src1(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
c =========================================================
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
c
c     integrate a source term for friction if ifrictiontype>0.

      common /geo/ grav,abs_tol 
      common /fluxlimiting/ relimit1,relimit2
      common /sourceterms/ frictioncoeff,ifrictiontype

      g=grav
      coeff = frictioncoeff
      tol = 1.d-30  !# to prevent divide by zero in gamma

      if (frictioncoeff.eq.0.d0 .or. ifrictiontype.eq.0) return

      if (ifrictiontype.eq.1) then ! integrate source term based on Manning formula
        do i=0,mx+1
           h=q(i,1)
           hu=q(i,2)
           if (h.lt.tol) then !set momentum to zero
                q(i,2)=0.d0
           else !integrate friciton
                gamma= dsqrt(hu**2)*(g*coeff**2)/(h**(7/3))
                dgamma=1.d0 + dt*gamma
                q(i,2)= q(i,2)/dgamma
           endif
        enddo
      endif 

      return
      end
