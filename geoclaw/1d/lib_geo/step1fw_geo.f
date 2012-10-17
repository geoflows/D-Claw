c
c
c ===================================================================
      subroutine step1(maxmx,meqn,mwaves,mbc,mx,q,aux,dx,dt,
     &              method,mthlim,cfl,f,fwave,s,amdq,apdq,dtdx)
!
!     # Take one time step, updating q.
!
!     ----------------------------------------------------------------------
!     # step1fw is a modified version of step1 to use fwave instead of wave.
!     # A modified Riemann solver rp1 must be used in conjunction with this
!     # routine, which returns fwave's instead of wave's.
!     # See http://amath.washington.edu/~claw/fwave.html
!
!     # Limiters are applied to the fwave's, and the only significant
!     # modification of this code is in the "do 110" loop, for the 
!     # second order corrections.
!     ----------------------------------------------------------------------
!
!     method(1) = 1   ==>  Godunov method
!     method(1) = 2   ==>  Slope limiter method
!     mthlim(p)  controls what limiter is used in the pth family
!
!
!     amdq, apdq, fwave, s, and f are used locally:
!
!     amdq(1-mbc:maxmx+mbc, meqn) = left-going flux-differences
!     apdq(1-mbc:maxmx+mbc, meqn) = right-going flux-differences
!        e.g. amdq(i,m) = m'th component of A^- \Delta q from i'th Riemann
!                         problem (between cells i-1 and i).
!
!     fwave(1-mbc:maxmx+mbc, meqn, mwaves) = waves from solution of
!                                           Riemann problems,
!            fwave(i,m,mw) = mth component of jump in f across
!                           wave in family mw in Riemann problem between
!                           states i-1 and i.
!
!     s(1-mbc:maxmx+mbc, mwaves) = wave speeds,
!            s(i,mw) = speed of wave in family mw in Riemann problem between
!                      states i-1 and i.
!
!     f(1-mbc:maxmx+mbc, meqn) = correction fluxes for second order method
!            f(i,m) = mth component of flux at left edge of ith cell 
!     --------------------------------------------------------------------
!     --------------------------------------------------------------------
!	  This version step1fw_geo.f  accumulates the updates
!       into interface fluxes.  It then relimits the interface fluxes so that
!       non-negativity of depth is maintained for depth averaged flow equations
!       where q(i,1) is the depth. 
!      
!       To enable this relimiting set relimit1=.true. and relimit2=.true.
!       To only relimit the second order correction fluxes, set relimit2=.true.

!       if relimit1 = relimit2 = .false. then it should behave the same as standard clawpack.
!
!       last modified 1/29/09, David George.
!     -----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension    q(1-mbc:maxmx+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, *)
      dimension    f(1-mbc:maxmx+mbc, meqn)
      dimension    s(1-mbc:maxmx+mbc, mwaves)
      dimension fwave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension amdq(1-mbc:maxmx+mbc, meqn)
      dimension apdq(1-mbc:maxmx+mbc, meqn)
      dimension dtdx(1-mbc:maxmx+mbc)

      dimension method(7),mthlim(mwaves)
      dimension fp(2), fm(2), fi(2)
      
      logical limit,relimit1,relimit2


      common /geo/ grav,abs_tol 
      common /fluxlimiting/ relimit1,relimit2
c
c     # check if any limiters are used:
      limit = .false.

      tol=abs_tol

      do 5 mw=1,mwaves
	 if (mthlim(mw) .gt. 0) limit = .true.
   5     continue
c
      mcapa = method(6)
      do 10 i=1-mbc,mx+mbc
	 if (mcapa.gt.0) then
	     if (aux(i,mcapa) .le. 0.d0) then
		write(6,*) 'Error -- capa must be positive'
		stop
		endif
             dtdx(i) = dt / (dx*aux(i,mcapa))
	    else
             dtdx(i) = dt/dx
	    endif
   10	 continue
c
c
c
c     # solve Riemann problem at each interface 
c     -----------------------------------------
c      write(*,*) 'rp1'
      call rp1(maxmx,meqn,mwaves,mbc,mx,q,q,aux,aux,fwave,s,amdq,apdq)

c     # Modify q for Godunov update:
c     # Note this may not correspond to a conservative flux-differencing
c     # for equations not in conservation form.  It is conservative if
c     # amdq + apdq = f(q(i)) - f(q(i-1)).

c==================================Geoclaw===========================================
c     # if relimitall=.true. the fluctuations, amdq and apdq are converted
c     # to true interface fluxes. 
c     # ordinarily amdq = f(q_{i-1/2}^-)-f(q(i-1)) and apdq = f(q(i))-f(q_{i-1/2}^+)
c     # relimitall makes amdq --> f(q_{i-1/2}^-) and apdq --> f(q_{i-1/2}^+)
c      write(*,*) 'relimit'
      if (relimit1) then
         do i=1,mx+1
            ! convert to interface fluxes
            ! fi is the flux in the cell center
            fi(1)=q(i,2)
            fi(2)= 0.d0
            if (q(i,1).gt.abs_tol) then 
                fi(2) = 0.5d0*grav*(q(i,1))**2 + q(i,2)**2/q(i,1)
            endif
            ! note that fp refers to the right of cell i and fm the left
            do m=1,2
               fp(m) = amdq(i+1,m) + fi(m)
               fm(m) = fi(m) - apdq(i,m)
            enddo

            ! test if relimiting is necessary for cell i
            outdensity = dtdx(i)*
     &          (dmax1(0.d0,fp(1))-dmin1(0.d0,fm(1)))
            density = q(i,1)

            if (outdensity.gt.density) then ! limit interface fluxes

               ! call subroutine novaccum1(fp,fm,density,dtdx,meqn,meqnlim,tol)
               call novaccum1(fp,fm,density,dtdx(i),meqn,meqn,tol)

               !Note: only one of apdq or amdq at an interface is limited
               !      since apdq(1)=amdq(1), the mass flux at an interface
               !      and so apdq and amdq remove mass from only one cell.
               !      so limiting will not happen twice for i and i+1.

               ! determine fluctuations from newly limited interface fluxes.
               do m=1,2
                  sumi = apdq(i,m)+amdq(i,m)
                  apdq(i,m) = fi(m)-fm(m)
                  amdq(i,m) = sumi - apdq(i,m)
                  sumip = apdq(i+1,m) + amdq(i+1,m)
                  amdq(i+1,m) = fp(m) - fi(m)
                  apdq(i+1,m) = sumip - amdq(i+1,m)
               enddo
            endif
         enddo
      endif
c============================================================================


      do 40 i=1,mx+1
         do 39 m=1,meqn
            q(i,m) = q(i,m) - dtdx(i)*apdq(i,m)
            q(i-1,m) = q(i-1,m) - dtdx(i-1)*amdq(i,m)
   39       continue
   40    continue

c     # compute maximum wave speed:
      cfl = 0.d0
      do 50 mw=1,mwaves
	 do 45 i=1,mx+1
c          # if s>0 use dtdx(i) to compute CFL,
c          # if s<0 use dtdx(i-1) to compute CFL:
	   cfl = dmax1(cfl, dtdx(i)*s(i,mw), -dtdx(i-1)*s(i,mw))
   45      continue
   50    continue
c
      if (method(2) .eq. 1) go to 900
c
c     # compute correction fluxes for second order q_{xx} terms:
c     ----------------------------------------------------------
c
      do 100 m = 1, meqn
            do 100 i = 1-mbc, mx+mbc
               f(i,m) = 0.d0
  100          continue
c
c     # apply limiter to waves:
      if (limit) then
         call limiter(maxmx,meqn,mwaves,mbc,mx,fwave,s,mthlim)
         endif

c
      do 120 i=1,mx+1
	 do 120 m=1,meqn
	    do 110 mw=1,mwaves
	       dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
	       f(i,m) = f(i,m) + 0.5d0 * dsign(1.d0,s(i,mw))
     &		   * (1.d0 - dabs(s(i,mw))*dtdxave) * fwave(i,m,mw)
c

  110          continue
  120       continue
c
c
  140 continue
c
c     # update q by differencing correction fluxes 
c     ============================================
c
c     # (Note:  Godunov update has already been performed above)


!     relimit the correction fluxes to ensure they maintain positivity of mass      
      if (relimit2) then
         do i=1,mx+1
            ! test if relimiting is necessary for cell i
            outdensity=dtdx(i)*(dmax1(0.d0,f(i+1,1))-dmin1(0.d0,f(i,1)))
            density = q(i,1)
            if (outdensity.gt.density) then ! limit interface fluxes
               ! call subroutine novaccum1(fp,fm,density,dtdx,meqn,meqnlim,tol)
               call novaccum1(f(i+1,1:meqn),f(i,1:meqn),
     &                       density,dtdx(i),meqn,meqn,tol)
            endif
         enddo
      endif

c
      do 150 m=1,meqn
	 do 150 i=1,mx
	    q(i,m) = q(i,m) - dtdx(i) * (f(i+1,m) - f(i,m))
  150       continue
c
  900 continue
      return
      end
