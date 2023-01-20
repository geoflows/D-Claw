

   !=========================================================
      subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
   !=========================================================
      use geoclaw_module
      use digclaw_module

      implicit none

      !i/o
      double precision :: q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      double precision :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
      double precision :: xlower,ylower,dx,dy,t,dt
      integer :: maxmx,maxmy,meqn,mbc,mx,my,maux

      !local
      double precision :: gmod,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi,dti
      double precision :: D,tau,sigbed,kperm,compress,pm,coeff,tol
      double precision :: m_eqn
      double precision :: zeta,p_hydro,p_litho,p_eq,gamma,dgamma
      double precision :: krate,krate_m,krate_h
      double precision :: vnorm,hvnorm,theta,dtheta,w,taucf,fsphi,hvnorm0
      double precision :: shear,sigebar,pmtanh01,rho_fp,seg
      double precision :: b_xx,b_yy,b_xy,chi,beta
      double precision :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      double precision :: vlow,m2,vreg,slopebound
      double precision :: b_eroded,b_remaining,dtcoeff
      double precision :: k1,k2,k3,k4
      integer :: i,j,ii,jj,jjend,icount
      logical :: ent

      !source fountain
      integer :: numCellsX,numCellsY,numCells
      integer :: numCellsHalfX,numCellsHalfY,numCellsHalf
      integer :: srcI,srcJ,srcJlo,srcIlo,srcIhi,srcJhi
      double precision :: src_xloc, src_yloc
      double precision :: s_xloclo,s_yloclo
      double precision :: s_xlochi,s_ylochi
      double precision :: s_angle,s_h,s_m0,s_q,s_Qp,s_qx,s_qy
      double precision :: s_slope,s_slope_x,s_slope_y,s_tend
      double precision :: s_vel,s_velx,s_vely,s_Vtot,s_xloc,s_yloc
      double precision :: srcXComp,srcYComp,x,y

      ! analytic calculation of m when overshoots meqn.
      double precision :: deltat,dt_remaining
      logical :: overshoot
      integer :: i_fast

      ! level awareness
      double precision :: dxmin, dymin
      common /comfine/ dxmin,dymin

      double precision, allocatable :: moll(:,:)

      ! check for NANs in solution:
      call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

      gmod=grav
      coeff = coeffmanning
      tol = drytolerance !# to prevent divide by zero in gamma
      !write(*,*) 'src:init,value',p_initialized,init_pmin_ratio
      if (entrainment>0) then
         ent = .true.
      else
         ent = .false.
      endif

      do i=1-mbc+1,mx+mbc-1
         do j=1-mbc+1,my+mbc-1
            theta = 0.d0
            dtheta = 0.d0
            if (bed_normal==1) then
               theta = aux(i,j,i_theta)
               gmod = grav*cos(theta)
               dtheta = -(aux(i+1,j,i_theta) - theta)/dx
            endif

            h = q(i,j,1)

            ! do not integrate if p not initialized or if h<drytol
            if (h<=drytolerance) cycle
            if (p_initialized==0) cycle

            ! loop to get full dt completed, even if estimated dt < dt.
            ! use dti = estimated_dt*0.9 which ensures not going over m_eqn

            i_fast = 1
            dt_remaining = dt

            do while (dt_remaining>0.0)

              !call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),q(i,j,5),u,v,m,theta)
              h = q(i,j,1)
              if (h<=drytolerance) cycle
              hu = q(i,j,2)
              hv = q(i,j,3)
              hm = q(i,j,4)
              p =  q(i,j,5)
              phi = aux(i,j,i_phi)
              pm = q(i,j,6)/h
              pm = max(0.0,pm)
              pm = min(1.0,pm)
              fsphi = aux(i,j,i_fsphi)

              call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
              call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm,m_eqn)

              ! handle any changes to seg, rho_fp, and pmtanh01 associted with segregation.
              if (dabs(alpha_seg-1.0)<1.d-6) then
                 seg = 0.0
                 rho_fp = rho_f
                 pmtanh01=0.0
              else
                 seg = 1.0
                 call calc_pmtanh(pm,seg,pmtanh01)
                 rho_fp = max(0.d0,(1.0-pmtanh01))*rho_f
              endif
              !pmtanh01 = seg*(0.5*(tanh(20.0*(pm-0.80))+1.0))
              !pmtanh01 = seg*(0.5*(tanh(40.0*(pm-0.90))+1.0))

!             calculate amount of time to reach m_eqn. This is equivalent to log product approach
!              dm/dt = -Dm/h
!              m(t) = m0 * exp(-Dt/h)
!              for m(deltat) = m_eqn = m0 * exp(-D*deltat/h)
!              deltat = -(h/D) * log(m_eqn/m0)

              deltat = -(h/D) * log(m_eqn/m)
              dti = min(dt_remaining, 0.9*deltat)

              if (dti<0.000001) then
                dti = min(0.0001, dt_remaining)
              endif

                ! if dti < 0 this indicates that D has switched signs before m
                ! has reached m_eqn. This is an indication that the basal pressure.
                ! has relaxed to hydrostatic faster than m has relaxed to m_eqn.
                ! At this point in time, the equilibrium values for h,hu,hv,hm,pb
                ! have reached their short timescale equilibria and the integration
                ! of src terms is complete.

              overshoot = dti<dt_remaining

              if (overshoot) then
                write(*,*) "iteration = ", i_fast
                write(*,*) "dti = ", dti
                write(*,*) "dx =", dx
                write(*,*) "m, m/meqn at t=0            = ", m, m/m_eqn
              endif



              !integrate momentum source term
              !tau = max(tau*(1.0-fsphi),0.0)

              vnorm = sqrt(u**2.0 + v**2.0)
              hvnorm = sqrt(hu**2.0 + hv**2.0)
              hvnorm0 = hvnorm

              !integrate friction
              hvnorm = dmax1(0.d0,hvnorm - dti*tau/rho)
              hvnorm = hvnorm*exp(-(1.d0-m)*2.0d0*mu*dti/(rho*h**2.0))
              !hvnorm = hvnorm*exp(-(1.d0-m)*2.0d0*0.1*dti/(rho*h**2.0))
              if (hvnorm<1.e-16) hvnorm = 0.0

              if (hvnorm>0.0.and.curvature==1) then
                 b_xx=(aux(i+1,j,1)-2.d0*aux(i,j,1)+aux(i-1,j,1))/(dx**2)
                 b_yy=(aux(i,j+1,1)-2.d0*aux(i,j,1)+aux(i,j-1,1))/(dy**2)
                 b_xy=(aux(i+1,j+1,1)-aux(i-1,j+1,1) -aux(i+1,j-1,1)+aux(i-1,j-1,1))/(4.0*dx*dy)
                 chi = (u**2*b_xx + v**2*b_yy + 2.0*u*v*b_xy)/gmod
                 chi = max(chi,-1.0)
                 taucf = chi*tau
                 hvnorm = dmax1(0.d0,hvnorm - dti*taucf/rho)
                 taucf = u**2.0*dtheta*tau/gmod
                 hvnorm = dmax1(0.d0,hvnorm - dti*taucf/rho)
              endif

              if (hvnorm0>0.0) then
                 hu = hvnorm*hu/hvnorm0
                 hv = hvnorm*hv/hvnorm0
              endif

              call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
              call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm,m_eqn)

              vnorm = sqrt(u**2.0 + v**2.0)

              ! KRB not clear this is used. I think this is now done in auxeval
              !integrate shear-induced dilatancy
              sigebar = rho*gmod*h - p + sigma_0
              shear = 2.0*vnorm/h
              krate = 1.5*shear*m*tanpsi/alpha
              sigebar = sigebar*exp(krate*dti)
              !p = rho*gmod*h + sigma_0 - sigebar
              ! end of KRB not clear this is used. I think this is now done in auxeval

              ! start of pressure.
              ! dp/dt = k1 (p-k2) + k3
              ! with k1 = krate = (- zeta * 2 * kperm) /(h * mu)
              ! k2 = p_equilibrium = h*rho_fp*gmod
              ! k3 = 3.0*vnorm*tanpsi/(h*compress)
              ! this also has an analytical solution.
              ! https://www.wolframalpha.com/input?i=dx%2Fdt+%3D+k1*%28x-k2%29+-+k3
              ! p(t) = k4 exp(k1*t) + k3/k1 + k2
              ! p(t=0) = p0 = k4 + k3/k1 + k2
              ! k4 = p0 - (k3/k1) - k2
              ! if k3 = 0
              ! k4 = p0 - k2
              ! p(t) = k4 * exp(k1*t) + k2
              ! which translates to the standard eqn.
              ! p = p_eq + (p-p_eq)*exp(krate*dti)

              k3 = 0.0
              if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
                 k3 = 3.0*vnorm*tanpsi/(h*compress)
              endif

              !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
              !call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm,m_eqn)

              !second integrate pressure relaxation
              !if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
              !   zeta = 3.d0/(compress*h*2.0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
              !else
              !   zeta = (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
              !endif
              zeta = ((m*(sigbed +  sigma_0))/alpha)*3.d0/(h*2.0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
              krate = -zeta*2.0*kperm/(h*max(mu,1.d-16))
              p_hydro = h*rho_fp*gmod
              p_litho = (rho_s*m + (1.d0-m)*rho_fp)*gmod*h

              !if (abs(compress*krate)>0.0) then
              !   p_eq = p_hydro + 3.0*vnorm*tanpsi/(compress*h*krate)
              !else
              !   p_eq = p_hydro
              !endif
              !if (abs(pm-.5)>.49) then
              !pmtanh01 = 0.5*(tanh(20.0*(pm-0.80))+1.0)
              p_eq = p_hydro !*(1.0-pmtanh01)

              k1 = krate
              k2  = p_eq
              k4 = p - (k3/k1) - k2

              !p_eq = max(p_eq,0.0)
              !p_eq = min(p_eq,p_litho)

              p = (k4 * exp(k1*dti)) + (k3/k1) + k2

              ! should aux get updated again here? TODO?
              ! originally it was, but now that there is a fast timescale loop?
              !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
              !call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm,m_eqn)

!              if (overshoot) then
!                write(*,*) ""
!                write(*,*) "adjusted dti for meqn"
!                write(*,*)  "i_fast", i_fast
!                write(*,*) "dt_remaining", dt_remaining
!                write(*,*) "deltat", deltat
!                write(*,*)  "D=", D
!                write(*,*) "t0: h, m, hm: ", h, m, h*m
!                write(*,*) "t0: meqn, m/meqn: ", m_eqn, m/m_eqn
!              endif

              krate = D*(rho-rho_fp)/rho
              hu = hu*exp(dti*krate/h)
              hv = hv*exp(dti*krate/h)
              hm = hm*exp(-dti*D*rho_fp/(h*rho))
              h = h + krate*dti

!              if (overshoot) then
!                write(*,*) "t1: h, m, hm: ", h, hm/h, hm
!                write(*,*) "t1: meqn, m/m_eqn: ", m_eqn, hm/h/m_eqn
!              endif

              ! call admissibleq and auxeval after updating all parts of src
              call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
              call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm,m_eqn)

!              if (overshoot) then
!                write(*,*) "Updated D: ", D
!                write(*,*) "Updated m_eqn ", m_eqn
!                write(*,*) ""
!              endif

              ! decrease dt_remaining
              dt_remaining = dt_remaining - dti
              i_fast = i_fast + 1

            enddo ! end fast timescale integration.

100         continue

            ! set dti to dt for remainder of src2 routine.
            dti = dt
            !======================mass entrainment===========================

            vnorm = sqrt(u**2.0 + v**2.0)
            vlow = 0.1d0

            if (ent.and.vnorm.gt.vlow.and.(aux(i,j,i_theta)>0.0)) then
               b_x = (aux(i+1,j,1)+q(i+1,j,7)-aux(i-1,j,1)-q(i-1,j,7))/(2.d0*dx)
               b_y = (aux(i,j+1,1)+q(i,j+1,7)-aux(i,j-1,1)-q(i,j-1,7))/(2.d0*dy)
               dbdv = (u*b_x+v*b_y)/vnorm
               slopebound = 1.e10
               b_eroded = q(i,j,7)
               if (dbdv<slopebound.and.b_eroded<aux(i,j,i_theta)) then
                  b_remaining = aux(i,j,i_theta)-b_eroded
                  m2 = 0.6d0
                  rho2 = m2*2700.d0 + (1.d0-m2)*1000.d0
                  beta2 = 0.66d0
                  t1bot = beta2*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d-2))
                  !write(*,*) '------------'
                  !write(*,*) 'vu',t1bot
                  beta = 1.0-m!tanh(10.d0*m) !tan(1.5*p/(rho*gmod*h))/14.0
                  gamma= rho*beta2*(vnorm**2)*(beta*gmod*coeff**2)/(tanh(h+1.d-2)**(1.0/3.0))
                  !write(*,*) 'gamma', gamma
                  t1bot = t1bot + gamma
                  t1bot = t1bot + tau!+p*tan(phi)
                  !write(*,*) 'tau',tau
                  t2top = min(t1bot,(1.d0-beta*entrainment_rate)*(tau))
                  !write(*,*) 't2top',t2top
                  prat = p/(rho*h)
                  !dh = dti*(t1bot-t2top)/(beta2*tanh(vnorm+1.d-2)*rho2)
                  vreg = ((vnorm-vlow)**2/((vnorm-vlow)**2+1.d0))
                  dtcoeff = entrainment_rate*dti*vreg/(beta2*(vnorm+vlow)*rho2)
                  !dh = dtcoeff*t1bot/(1.d0 + dtcoeff*tan(phi))
                  dh = dtcoeff*(t1bot-t2top)
                  dh = entrainment_rate*dti*(t1bot-t2top)/(rho2*beta2*vnorm)
                  !write(*,*) 'dh',dh
                  !write(*,*) 'dh/dti', dh/dti
                  dh = min(dh,b_remaining)
                  h = h + dh
                  hm = hm + dh*m2
                  q(i,j,7) = q(i,j,7) + dh

                  if (hm/h.gt.0.8.or.hm/h.lt.0.3) then
                    write(*,*) "src2_entr call admissibleq(h,hu,hv,hm,p,u,v,m,theta) m=", hm/h
                  endif

                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm,m_eqn)
                  p = prat*rho*h
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
               endif
            endif
            !===================================================================

            q(i,j,1) = h
            q(i,j,2) = hu
            q(i,j,3) = hv
            q(i,j,4) = hm
            q(i,j,5) = p
            q(i,j,6) = pm*h

         enddo
      enddo

      !mollification
      if (.false.) then !mollification ?
      if (.not.allocated(moll)) then
         allocate(moll(mx,my))
      endif
      do i=1,mx
         do j=1,my
         if (q(i,j,1)<=drytolerance) cycle
            p = 0.0
            icount = 0
            do ii=-1,1
               do jj=-1,1
                  !cell weights
                  if ((abs(ii)+abs(jj))==0) then
                     w = 0.5
                  elseif((abs(ii)+abs(jj))==1) then
                     w = 0.3/4.0
                  else
                     w = 0.2/4.0
                  endif
                  if (q(i+ii,j+jj,1)>drytolerance) then
                     p = p + w*q(i+ii,j+jj,5)/q(i+ii,j+jj,1)
                  else
                     p = p + w*q(i,j,5)/q(i,j,1)
                     icount = icount + 1
                  endif
               enddo
            enddo
            moll(i,j) = p !/icount
         enddo
      enddo
      do i=1,mx
         do j=1,my
            if (q(i,j,1)<=drytolerance) cycle
            q(i,j,5) = moll(i,j)*q(i,j,1)
         enddo
      enddo
      if (allocated(moll)) then
         deallocate(moll)
      endif
      endif

      ! Manning friction------------------------------------------------
      if (ifriction==0) return
      if (coeffmanning>0.d0.and.frictiondepth>0.d0) then
         do i=1,mx
            do j=1,my

               if (bed_normal==1) gmod = grav*cos(aux(i,j,i_theta))
               h=q(i,j,1)
               if (h<=frictiondepth) then
                  !# apply friction source term only in shallower water
                  hu=q(i,j,2)
                  hv=q(i,j,3)
                  hm = q(i,j,4)
                  p =  q(i,j,5)
                  phi = aux(i,j,i_phi)
                  theta = aux(i,j,i_theta)

                  if (hm/h.gt.0.8.or.hm/h.lt.0.3) then
                    write(*,*) "src2_man call admissibleq(h,hu,hv,hm,p,u,v,m,theta) m=", hm/h
                  endif

                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  if (h<tol) cycle
                  pm = q(i,j,6)/h
                  pm = max(0.0,pm)
                  pm = min(1.0,pm)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm,m_eqn)

                  if (h.lt.tol) then
                     q(i,j,1)=0.d0
                     q(i,j,2)=0.d0
                     q(i,j,3)=0.d0
                  else
                     beta = 1.0-m !tan(1.5*p/(rho*gmod*h))/14.0
                     gamma= beta*dsqrt(hu**2 + hv**2)*(gmod*coeff**2)/(h**(7.d0/3.d0))
                     dgamma=1.d0 + dt*gamma
                     q(i,j,2)= q(i,j,2)/dgamma
                     q(i,j,3)= q(i,j,3)/dgamma
                  endif
               endif
            enddo
         enddo
      endif
      ! ----------------------------------------------------------------


      !------------------------------
      ! In domain source fountain
      ! RPJ: 3/7/22
      ! KRB: 11/30/2022
      ! Specifically for post-fire debris flow hazard assessment
      ! assumes triangle hydrograph as several locations
      ! Peak volume and concentration given
      !------------------------------

      if (src_fountain_active .eqv. .TRUE.) then
         if (src_ftn_num .eq. 0) then
            src_fountain_active = .FALSE.
         endif
         if (t .gt. src_ftn_end_time) then
            src_fountain_active = .FALSE.
         endif

         do ii = 1,src_ftn_num
            ! only add material at finest level
            if (dx .gt. dxmin) then
               cycle
            endif

            !write(*,*) "fountain valid"

            ! calculate discharge to determine if
            ! end time has occured.
            s_Vtot = src_ftn_vtot(ii)
            s_Qp = 0.1 * (s_Vtot ** 0.833)
            if (s_Qp .lt. 0.0001) then
               cycle
            endif
            s_tend = 2.0 * s_Vtot / s_Qp

            if (t .gt. s_tend) then
               cycle
            endif

            ! get number of cells for a box of dimensions [numCellsX, numCellsY]
            ! centered at the location indicated by
            ! [src_ftn_xloc[ii], src_ftn_yloc[ii]]
            numCellsX = max(1, nint(src_ftn_length/dx))
            numCellsY = max(1, nint(src_ftn_length/dy))
            numCellsHalfX = int(numCellsX/2) ! num half cells as int
            numCellsHalfY = int(numCellsY/2) ! num half cells as int

            ! get source location extent
            src_xloc = src_ftn_xloc(ii)
            src_yloc = src_ftn_yloc(ii)
!            s_xloclo = src_ftn_xloc(ii) - numCellsHalfX * dx
!            s_yloclo = src_ftn_yloc(ii) - numCellsHalfY * dy

            ! get source center index location I/J
            srcI = int((src_xloc-xlower)/dx)
            srcJ = int((src_yloc-ylower)/dy)

            !write(*,*) xlower, xlower+mx*dx, ylower,ylower+my*dy

            !write(*,*) srcI,srcJ

            ! determine if source center is on this grid
            if(srcI.lt.1 .or. srcJ.lt.1 .or. srcI.gt.mx .or. srcJ.gt.my) then
               cycle
            endif

            ! if source is on this grid, add material.

            ! get index location of edges, snapped to grid extent.
            srcIlo = max(1, srcI-numCellsHalfX)
            srcJlo = max(1, srcJ-numCellsHalfY)
            srcIhi = min(srcIlo+numCellsX-1, mx)
            srcJhi = min(srcJlo+numCellsY-1, my)

            ! calculate total number of input cells.
            numCells = (srcIhi - srcIlo + 1) * (srcJhi - srcJlo + 1)

            !write(*,*) srcIlo,srcIhi,srcJlo,srcJhi,numCells,numCellsX,numCellsY

            ! calculate local slope to use for velocity.
            s_slope_x = (aux(srcI+1,srcJ,1)-aux(srcI-1,srcJ,1))/(2.d0*dx)
            s_slope_y = (aux(srcI,srcJ+1,1)-aux(srcI,srcJ-1,1))/(2.d0*dy)
            !s_slope = sqrt((s_slope_x * srcXComp)**2 + (s_slope_y * srcYComp)**2)
            s_slope = sqrt((s_slope_x)**2 + (s_slope_y)**2)
            if(s_slope .lt. 0.0001) then
               ! if no slope use a value for the component calculations (which will still be 0 in each direction)
               s_slope = 0.1 ! about 6 degrees.
            endif
            ! get x and y components of slope.

            srcXComp = s_slope_x/s_slope
            srcYComp = s_slope_y/s_slope

            ! Calculate discharge
            ! triangle function with peak centered at 0.5 * s_end
            s_q = 2.0*s_Qp*(t/s_tend-0.5)*sign(1.0d0,(0.5-t/s_tend))+s_Qp

            ! if current discharge is very small, don't add.
            if (s_q .lt. 1e-4) then
               cycle
            endif

            ! calculate x and y components of velocity.
            s_vel = abs(2.1 * s_Qp**0.33 * s_slope**0.33) !Rickenmann Eq 21
            s_velx = s_vel * srcXComp
            s_vely = s_vel * srcYComp

            ! get m0
            s_m0 = src_ftn_m0(ii)


            !write(*,*) "adding fountain: ", srcIlo,srcIhi,srcJlo,srcJhi,numCells
            ! add volume equally across all grid cells
            do i=srcIlo,srcIhi
              do j=srcJlo,srcJhi
                 q(i,j,1) = q(i,j,1) + s_q*dt/(numCells*dx*dy)
                 q(i,j,2) = q(i,j,2) + s_q*dt/(numCells*dx*dy)*s_velx
                 q(i,j,3) = q(i,j,3) + s_q*dt/(numCells*dx*dy)*s_vely
                 q(i,j,4) = q(i,j,4) + s_q*dt/(numCells*dx*dy)*s_m0


                 h = q(i,j,1)
                 hm = q(i,j,4)
                 if (hm/h.gt.0.8.or.hm/h.lt.0.3) then
                   write(*,*) "src2_3 call admissibleq(h,hu,hv,hm,p,u,v,m,theta) m=", hm/h
                 endif

                 call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),q(i,j,5),u,v,m,theta)
                 call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm,m_eqn)

                 !  Calc RHO and set q5 to hydrostatic?
                 !p_hydro = h*rho_fp*gmod
                 !q(i,j,5) = p_hydro

                 !call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),q(i,j,5),u,v,m,theta)

              end do
            end do

         enddo ! source number

      endif ! source fountain active

      return
      end

! the following routines used to implement logprod were implemented after
! the cython routines in scipy.special.

! Both D-Claw and Scipi use BDS-3 license.

!=========================================================
  subroutine logprod(z, tol, w_z)
!=========================================================

  implicit none

  double precision :: z,w_z,tol
  double precision :: w,ew,wewz,wn

  integer ::  i

  ! initial guess
  if (abs(z + exp(-1.0)).lt.0.3) then
      call lambertw_branchpt(z, w)
  elseif (z.lt.1.5) then
      ! Empirically determined decision boundary where the Pade
      ! approximation is more accurate.
      call lambertw_pade0(z, w)
  else
      call lambertw_asy(z, 0, w)
  endif

  do i=0,100
      ew = exp(-w)
      wewz = w - z*ew
      wn = w - wewz/(w + 1 - (w + 2)*wewz/(2*w + 2))
      if (abs(wn - w).lt.(tol*abs(wn))) then
          w_z = wn
          return
      else
          w = wn
      endif
  enddo
  w_z = w
  return
  end

  !=========================================================
    subroutine lambertw_branchpt(z, w_z)
  !=========================================================

    real, dimension(3) :: coeffs
    double precision :: z, w_z
    double precision :: p,e

    e = 2.71828

    coeffs(1) = -1.0/3.0
    coeffs(2) = 1.0
    coeffs(3) = -1.0

    p = sqrt(2.0*(e * z + 1.0))
    call cevalpoly(coeffs, 2, p, w_z)

    return
    end

    !=========================================================
      subroutine lambertw_pade0(z, w_z)
    !=========================================================

      real, dimension(3) :: num
      real, dimension(3) :: denom

      double precision :: z,w_z
      double precision :: numerator, denominator

      num(1) = 12.85106382978723404255
      num(2) = 12.34042553191489361902
      num(3) = 1.0

      denom(1) = 32.53191489361702127660
      denom(2) = 14.34042553191489361702
      denom(3) = 1.0

      call cevalpoly(num, 2, z, numerator )
      call cevalpoly(denom, 2, z, denominator)
      w_z = z*numerator/denominator
      return
      end

      !=========================================================
        subroutine lambertw_asy(z, k, w_z)
      !=========================================================
      double precision :: z,w_z, q
      integer :: k
      double precision :: pi

      pi=4.d0*atan(1.d0)

      w = log(z) + 2.0*pi*real(k) !* 1j
      w_z = w - log(w)
      return
      end

      !=========================================================
        subroutine cevalpoly(coeffs, degree, z, out)
      !=========================================================

      real, dimension(3) :: coeffs ! all useage has size 3.

      double precision :: z,out
      integer :: degree
      double precision :: r,s,a,b,tmp

      a = coeffs(1)
      b = coeffs(2)
      r = 2*z
      s = z*z

      do j=2,degree
          tmp = b
          call fma(-s, a, coeffs(j), b)
          call fma(r,a,tmp,a)
      enddo

      out = z*a+b
      return
      end

      !=========================================================
        subroutine fma(x,y,z,out)
      !=========================================================

      out = (x*y) + z
      return
      end