

   !=========================================================
      subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
   !=========================================================
      use geoclaw_module
      use digclaw_module

      implicit none


      !i/o
      real(kind=8) :: q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      real(kind=8) :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
      real(kind=8) :: xlower,ylower,dx,dy,t,dt
      integer :: maxmx,maxmy,meqn,mbc,mx,my,maux

      !local
      real(kind=8) :: gacc,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi,dti,gz,gx,dtk
      real(kind=8) :: D,tau,sigbed,kperm,compress,pm,dtremaining
      real(kind=8) :: zeta,gamma,dgamma,alphainv
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,w,hvnorm0
      real(kind=8) :: shear,sigebar,pmtanh01,rho_fp,seg
      real(kind=8) :: b_xx,b_yy,b_xy,chi,beta
      real(kind=8) :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      real(kind=8) :: vlow,m2,vreg,slopebound
      real(kind=8) :: b_eroded,b_remaining,dtcoeff
      real(kind=8) :: rhoh,sig_0,sig_eff,m_eq,m_0,m_eq_0,sig_eff_0
      integer :: i,j,ii,jj,jjend,itercount,itercountmax,quad0,quad1
      logical :: ent,debug

      !source fountain
      integer :: numCellsX,numCellsY,numCells
      integer :: numCellsHalfX,numCellsHalfY,numCellsHalf
      integer :: srcI,srcJ,srcJlo,srcIlo,srcIhi,srcJhi
      real(kind=8) :: src_xloc, src_yloc
      real(kind=8) :: s_xloclo,s_yloclo
      real(kind=8) :: s_xlochi,s_ylochi
      real(kind=8) :: s_m0,s_q,s_Qp
      real(kind=8) :: s_tend,s_Vtot,s_xloc,s_yloc
      real(kind=8) :: srcXComp,srcYComp,x,y

      ! level awareness
      real(kind=8) :: dxmin, dymin
      common /comfine/ dxmin,dymin

      debug=.true.

      ! check for NANs in solution:
      call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)
      
      gz = grav  !needed later for bed-normal direction gravity
      gx = 0.d0
      theta=0.d0 

      !do i=1-mbc+1,mx+mbc-1
         !do j=1-mbc+1,my+mbc-1

      do i=1,mx
         do j=1,my

            h = q(i,j,1)
            hu = q(i,j,2)
            hv = q(i,j,3)
            hm = q(i,j,4)
            p =  q(i,j,5)
            phi = aux(i,j,i_phi)
            pm = q(i,j,6)/h
            pm = max(0.0d0,pm)
            pm = min(1.0d0,pm)

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            if (h<=drytolerance) cycle

            !Manning friction
            if ((ifriction/=0).and.(coeffmanning>0.d0)) then
               if (h<=frictiondepth) then
                  beta = 1.0d0-m
                  gamma= beta*sqrt(hu**2 + hv**2)*(gz*coeffmanning**2)/(h**(7.d0/3.d0))
                  dgamma=1.d0 + dt*gamma
                  hu= hu/dgamma
                  hv= hv/dgamma
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
               endif
            endif

            !modified gravity: bed-normal weight and acceleration
            if (bed_normal==1) then
               theta = aux(i,j,6)
               gz = grav*cos(theta)
               gx = grav*sin(theta)
            endif
            if (curvature==1) then
               b_xx=(aux(i+1,j,1)-2.d0*aux(i,j,1)+aux(i-1,j,1))/(dx**2)
               b_yy=(aux(i,j+1,1)-2.d0*aux(i,j,1)+aux(i,j-1,1))/(dy**2)
               b_xy=(aux(i+1,j+1,1)-aux(i-1,j+1,1) -aux(i+1,j-1,1)+aux(i-1,j-1,1))/(4.0*dx*dy)
               dtheta = -(aux(i+1,j,6) - theta)/dx
               gacc = max(u**2*b_xx + v**2*b_yy + 2.0*u*v*b_xy + u**2*dtheta,0.d0)!max:currently only consider enhancement not reduction of gz (ie. basin not a hump)
               gz = gz + gacc
            endif

            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
            rhoh = h*rho !this is invariant in src and always >0 below
            
            !integrate momentum source term
            hvnorm0 = sqrt(hu**2 + hv**2)
            vnorm = hvnorm0/h

            if (hvnorm0>0.d0) then
               !integrate dynamic friction !DIG: TO DO - move dynamic friction to Riemann solver
               vnorm = max(0.d0,vnorm - dt*tau/rhoh) !exact solution for Coulomb friction
               vnorm = vnorm*exp(-(1.d0-m)*2.0d0*mu*dt/(h*rhoh)) !exact solution (prior to h change) for effective viscous friction
               ! velocity determined, calculate directions etc. from vnorm
               hvnorm = h*vnorm
               hu = hvnorm*hu/hvnorm0 + gx*h*dt !gx=0 unless bed-normal !DIG: last term should ultimately be in Riemann solver
               hv = hvnorm*hv/hvnorm0
               u = hu/h
               v = hv/h
               vnorm = sqrt(u**2 + v**2)
               ! velocity now constant for remainder of src2. hu,hv adjusted due to change in h
            endif

            if (p_initialized==0) cycle !# dlg: I can't remember why we don't just return in this case but don't want to change now

            !call admissibleq(h,hu,hv,hm,p,u,v,m,theta) ! dlg: don't see why this is needed again
            ! u,v change need new aux vals
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            !integrate m and p, keep rhoh ------------------------------------------------------------
           
            !explicit integration
            dtremaining = dt
            itercountmax=100
            itercount=0
            !if (debug.and.m>0.99) then
            !      write(*,*) 'WARNING SRC2: h,hm,m', h,hm,m
            !endif

            do while (dtremaining>1.d-99)
               call mp_update_FE_4quad(dtremaining,h,u,v,m,p,rhoh,gz,dtk)
               dtremaining = dtremaining-dtk
               itercount = itercount + 1
               if (dtk==0.d0.and.(.not.debug)) exit
               if (h<drytolerance) exit
               if (itercount>=itercountmax) then
                  exit
               endif
               if (debug.and.(itercount>50)) then
                  write(*,*) 'WARNING SRC2: update iteration'
                  write(*,*) 'itercount,dt,dtremaining:',itercount,dt,dtremaining
               endif
            enddo
            
            hu = h*u
            hv = h*v
            hm = h*m
            call qfix(h,hu,hv,hm,p,u,v,m,rho,gz)
            if (h<=drytolerance) then
               cycle
            endif
            !write(*,*) 'itercount', itercount,dt,dtremaining
            !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
            !--------------------------------------------------------------------------------------------
            
            if (dabs(alpha_seg-1.0)<1.d-6) then
               seg = 0.0d0
               rho_fp = rho_f
               pmtanh01=0.0d0
            else
               seg = 1.0d0
               call calc_pmtanh(pm,seg,pmtanh01)
               rho_fp = max(0.d0,(1.0d0-pmtanh01))*rho_f
            endif

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            !======================mass entrainment===========================

            if (entrainment>0) then
               ent = .true.
            else
               ent = .false.
            endif

            vnorm = sqrt(u**2 + v**2)
            vlow = 0.1d0

            if (ent.and.vnorm.gt.vlow.and.(aux(i,j,i_theta)>0.d0)) then
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
                  t1bot = beta2*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d0-2.d0))
                  !write(*,*) '------------'
                  !write(*,*) 'vu',t1bot
                  beta = 1.0d0-m!tanh(10.d0*m) !tan(1.5d0*p/(rho*gmod*h))/14.0d0
                  !! DIG not used. segfault.gamma= rho*beta2*(vnorm**2)*(beta*gmod*coeff**2)/(tanh(h+1.d0-2.d0)**(1.0d0/3.0d0))
                  !write(*,*) 'gamma', gamma
                  !! DIG not used. segfault.t1bot = t1bot + gamma
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

                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
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

      ! ----------------------------------------------------------------


      !------------------------------
      ! In domain source fountain
      ! RPJ: 3/7/22
      ! KRB: 11/30/2022
      ! Specifically for post-fire debris flow hazard assessment
      ! assumes triangle hydrograph as several locations
      ! Peak volume and concentration given
      !------------------------------
      if (.false.) then
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

            ! Calculate discharge
            ! triangle function with peak centered at 0.5 * s_end
            s_q = 2.0*s_Qp*(t/s_tend-0.5)*sign(1.0d0,(0.5-t/s_tend))+s_Qp


            ! if current discharge is very small, don't add.
            if (s_q .lt. 1e-4) then
               cycle
            endif

            ! get m0
            s_m0 = src_ftn_m0(ii)

            !write(*,*) "adding fountain: ", srcIlo,srcIhi,srcJlo,srcJhi,numCells
            ! add volume equally across all grid cells
            do i=srcIlo,srcIhi
              do j=srcJlo,srcJhi
                 q(i,j,1) = q(i,j,1) + s_q*dt/(numCells*dx*dy)
                 q(i,j,4) = q(i,j,4) + s_q*dt/(numCells*dx*dy)*s_m0

              end do
            end do

            h = q(i,j,1)
            hm = q(i,j,4)

            call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),q(i,j,5),u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
         enddo ! source number

      endif ! source fountain active
      endif
      return
      end subroutine src2

      subroutine mp_update_FE_4quad(dt,h,u,v,m,p,rhoh,gz,dtk)
      !====================================================================
      ! subroutine mp_update_FE_4quad: integrate dp_exc/dt,dm/dt by a hybrid 
      ! explicit integration that depends on the initial quadrant of phase space.
      ! Basic idea which conforms to phase space vector field: for each quadrant 
      ! of phase space (divided by p=p_eq(m_c) and m=m_c) there is only one open 
      ! interior (physically admissible) boundary that can be crossed. These boundaries 
      ! are such that a physically admissible solution must proceed clockwise (m horizontal 
      ! axis, p vertical) if it changes quadrants. 1 (UL)-> 2(UR)-> 3 (LR)-> 4(LL). 
      ! An open boundary in 1 quadrant is a closed boundary in then next, ie. m=m_c 
      ! belongs to quads 2 & 4, p=p_eq belongs to 1 & 3. Explicit solution to 
      ! non-homogeneous exponential ode d/dt(m,p_exc) is used with dtk<=dt  such that only 
      ! the open/admissible boundary can be crossed in a given substep. For a given 
      ! substep dt, either (a) sol remains interior to quadrant (dtk=dt), (b) sol. 
      ! reaches physically inadmissible boundary dtk<=dt (c) solution reaches open 
      ! quadrant boundary (dtk<=dt) that belongs to next quadrant for next substep.
      ! Quadrants 1 & 3 contain the p and m nullclines. Substeps are taken based on 
      ! the initial position wrt these nullclines, because solution is most prone 
      ! to oscillations if a substep overshoots a nullcline significantly. Solution 
      ! converges to p-nullcline envelope when it should, except when above and below
      ! the p-nullcline on left and right of m_eq respectively,the most problematic regions.
      ! m-nullcline is simply p=p_eq(m). p-nullcline is complicated but = m-nullcline if vnorm = 0. 
      ! Otherwise rises toward (above?) lithostatic pressure in left-half plane and 0 (negative?) 
      ! in right-half plane as v increases and kperm decreases.
      !====================================================================

         use digclaw_module, only: rho_f,rho_s,sigma_0,mu,alpha,setvars,qfix,qfix_cmass,phi_bed,m_crit,delta
         use geoclaw_module, only: grav,drytolerance
   
         implicit none
   
         !i/o
         real(kind=8), intent(inout) :: h,m,p
         real(kind=8), intent(in)  :: u,v,rhoh,dt
         real(kind=8), intent(in)  :: gz
         real(kind=8), intent(out) :: dtk

         !local
         real(kind=8) :: h0,p0,m_0,p_eq0,p_exc0,sig_eff,sig_0,vnorm,m_eq
         real(kind=8) :: kappa,S,rho,rho0,tanpsi,D,tau,sigbed,kperm,phi
         real(kind=8) :: km,kp,alphainv,c_d,p_exc,dtr,dtm,dtp,dts,p_excm,p_exc_ave
         real(kind=8) :: km0,kp0,alphainv0,c_d0
         real(kind=8) :: hu,hv,hm
         real(kind=8) :: rho_c,h_c,p_eq_c,m_c,rho_c1,h_c1,p_eq_c1,m_c1,m_lower,m_upper
         real(kind=8) :: m_c0,p_eq_c0,sig_c,Nd,Nn,normc,shear,convtol,m_eq1,sig_eff1
         integer :: debugloop,iter,itermax,quad0,quad1
         logical :: outquad,debug
   
         debugloop = 0
         debug = .true.
         outquad = .true.
         phi = phi_bed
         vnorm = sqrt(u**2 + v**2)
   
         !explicit integration (hybrid FE and explicit exponential solution)---------------------------------------------------
         ! q1 = q0 + dtk*f(q0)
         call setvars(h,u,v,m,p,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)
         h0 = h
         m_0 = m
         rho0 = m_0*(rho_s-rho_f)+rho_f
         p0 = p
         p_eq0 = rho_f*gz*h0
         p_exc0 = p0 - p_eq0
         
         dtk = 0.d0
         dtr = dt
         
         quad0=0
         quad1=0
         ! if at critical point dq/dt = 0
         if ((abs(p_exc0/gz*h0)<1.d-12).and.(m==0.d0.or.abs(m-m_eq)<1.d-12)) then
            dtk = dtr
            return
         endif

         if (m==0.d0) then
            p_exc = p_exc0*exp(-kp0*dtr)
            call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
            dtk = dtr
            return
         endif

         !determine coefficients for update
         ! dm/dt = (2*k*rho^2/mu(rhoh)^2)*p_exc*m
         ! dp_eq/dt = -plamda*p_eq + c_d see George & Iverson 2014
         km0 = ((2.d0*kperm*rho0**2)/(mu*rhoh**2))
         c_d0 = -3.d0*vnorm*(alphainv*rho0/(rhoh))*tanpsi
         kp0 = (kperm/(h0*mu))*(3.d0*alphainv*rho0/rhoh - 1.5d0*rho_f*gz*h0*(rho0-rho_f)/rhoh) !(kp>0)
         if (debug.and.kp0<0.d0) then
            write(*,*) '------------SRC WARNING: kp0<0 ---------->>>>>>>>'
            write(*,*) 'kp0<0:', kp0
            write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
         endif
         kp0 = max(kp0,1.d0) !shouldn't happen but prevent small rounding error
         if (c_d0==0.d0) then
            p_exc = p_exc0*exp(-kp0*dtr)
            m = m_0*exp(km0*p_exc*dtr)
            call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
            dtk = dtr
            return
         endif

         !find critical point (mc,peq_c) which defines quadrant origin
         !requires fixed-point iteration, should converge quickly
         m_c = m_crit
         rho_c = m_c*(rho_s-rho_f) + rho_f
         h_c = rhoh/rho_c
         p_eq_c = rho_f*gz*h_c
         shear = 2.d0*vnorm/h0
         if (shear>0.d0) then
            convtol = 1.d-16
            itermax = 100
            do iter=1,itermax
               m_c0 = m_c
               p_eq_c0 = p_eq_c
               rho_c = m_c*(rho_s-rho_f) + rho_f
               h_c = rhoh/rho_c
               p_eq_c = rho_f*gz*h_c
               sig_c = max(rhoh*gz-p_eq_c,0.d0)
               Nd = rho_s*(shear*delta)**2 + sig_c
               Nn = mu*shear
               m_c = m_crit*(sqrt(Nd)/(sqrt(Nd)+sqrt(Nn)))
               normc = (1.d2*(m_c - m_c0))**2 + (1.d-3*(p_eq_c-p_eq_c0))**2
               if (normc<convtol) then
                  exit
               endif
               if (iter>10.and.debug) then
                  write(*,*) '--------------------------------------------------'
                  write(*,*) 'SRC2 WARNING: fixed point iterations:', iter 
                  write(*,*) 'norm for convergence:', normc
               endif
            enddo
            if ((iter>=itermax).and.(normc>1.d0)) then
                  m_c = m_crit
                  rho_c = m_c*(rho_s-rho_f) + rho_f
                  h_c = rhoh/rho_c
                  p_eq_c = rho_f*gz*h_c
            endif
            
         endif

         !determine quadrant of initial solution in state space

         if ((m<m_c).and.(p>=p_eq_c)) then
            quad0=1
         elseif ((m>=m_c).and.(p>p_eq_c)) then
            quad0=2
         elseif  ((m>m_c).and.(p<=p_eq_c)) then
            quad0=3
         elseif ((m<=m_c).and.(p<p_eq_c)) then
            quad0=4
         endif
         
         select case (quad0)
         
         case(1) !UL quadrant, variable material and p_exc
            debugloop = 1000
            !statically loose
            !integration should only cross right boundary (m=m_c)
            !p can increase or decrease if below/above p-nullcline (which is left of m=m_eq, above p=p_eq)
            !m can increase or decrease if above/below m-nullcline (p=p_eq)
            if (p_exc0<-1.d-3*rhoh*gz) then !p increasing, m decreasing, below m-nullcline
               debugloop=debugloop + 100
               if (c_d0>0.d0) then !don't exceed p_exc = 0 until next substep
                  dtp = -(1.d0/kp0)*log(-c_d0/(kp0*p_exc0-c_d0))
               else
                  dtp = dtr
               endif
               dts = min(dtr,dtp)
               !integrate dp_exc/dt = -kp0 p_exc + c_d with coefficients at t=0.
               dts = max(dts,0.d0) 
               p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
               !m = m_0*exp(dts*km0*0.5d0*(p_exc0+p_exc))    
               m = m_0*exp(dts*km0*p_exc)
            elseif ((-kp0*p_exc0+c_d0)>0.d0) then !p_exc is increasing/below nullcline
               debugloop=debugloop + 200
               !bound time step 
               ! to max allowed (lithostatic,rhoh*g) if nullcline exceeds it
               if ((c_d/kp)>(rhoh*gz-p_eq_c)) then
                  debugloop=debugloop + 10
                  dtp = min(dtr,-(1.d0/kp0)*log((-kp0*rhoh*gz+kp0*p_eq_c+c_d)/(-kp0*p_exc0+c_d0)))
               !bound time step so not to exceed nullcline 
               else
                  debugloop=debugloop + 20
                  dtp = dtr
               endif
               p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
               p_exc_ave = 0.5d0*(p_exc+p_exc0)
               if (p_exc_ave<=0.d0) then !should only occur if p_exc0<0 and dt small
                  debugloop=debugloop + 1
                  dtm = dtr
                  m = m_0*exp(km*0.5d0*(p_exc+p_exc0)*dts)
               else
                  debugloop=debugloop + 2
                  m_upper = m_eq - max(0.d0,(1.d0/3.d0)*kp0*h0*(p_exc_ave)/(vnorm*alphainv))
                  dtm = min(dtp,(m_eq-m_0)/(km0*p_exc_ave*m_0))
                  !dts = min(dts,(m_eq-m_0)/(km0*p_exc*m_0)) !use new/higher p_exc
                  !m = m_0 + km0*m_0*0.5d0*(p_exc+p_exc0)*dts
               endif
               m = m_0 + km0*m_0*p_exc_ave*dtm
               dts = dtm
            else !-kp0*p_exc0+c_d0)<0
               debugloop=debugloop + 300
               !p_exc is decreasing/above nullcline or static on nullcline
               ! m strictly increasing
               if (c_d0>0.d0) then !p_exc remains in quad1 for any dt 
                  debugloop=debugloop + 10
                  ! only m can take solution beyond nullcline
                  dtp = dtr
                  !p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
                  !dts = min(dtp,(m_c-m_0)/(km0*0.5d0*(p_exc0+p_exc)*m_0))
                  !m = m_0 + km0*m_0*0.5d0*(p_exc0+p_exc)*dts
                  dtm = (m_c-m_0)/(km0*p_exc0*m_0)
                  dts = min(dtm,dtp)
                  dts = max(dts,0.d0)
                  m = m_0 + km0*m_0*p_exc0*dts
                  p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
               else !p_exc is to the right of m_eq
                   debugloop=debugloop + 20
                  ! true solution should remain right of m_eq
                  ! FE gives lowest slope of dp/dm
                  dts = min(dtr,(m_c-m_0)/(km0*p_exc0*m_0))
                  dtp = min(dtr,-(1.d0/kp0)*log((c_d0)/(-kp0*p_exc0+c_d0)))
                  !if (dtp<dts) write(*,*) 'ERROR: QUAD 1: WRONG DIRECTION ACROSS m=m_eq'
                  dts = min(dts,dtp)
                  dts = max(dts,0.d0)
                  m = m_0 + km0*m_0*p_exc0*dts
                  p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
               endif
            endif
            call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
            dtk = dts
            dtr = dtr-dts

         case(2) !UR quadrant, dense or equilibrium material, p>p_eq_c>p_eq
            !intertially contracting
            !integration can only cross lower boundary (p=p_eq_c>p_eq)
            !m>m_c>m_eq
            !p is strictly decreasing, m strictly increasing
            debugloop = 2000
            dtp = min(dtr,-(1.d0/kp0)*log((c_d0)/(-kp0*p_exc0+c_d0))) 
            dtm = min(dtp,(1.d0-m_0)/(km0*p_exc0*m_0))
            p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
            !m = m_0 + km0*m_0*0.5d0*(p_exc0+p_exc)*min(dts,dtm)
            m = m_0 + km0*m_0*p_exc*dtm
            dts = max(dtp,dtm)
            call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
            dtk = dts
            dtr = dtr-dts

         case(3) !LR quadrant, dilative material, p<=p_eq
            !integration can only cross left boundary (m=m_c)
            debugloop = 3000
            if (p_exc0>1.d-3*rhoh*gz) then !p decreasing, m increasing
               debugloop=debugloop + 100
               !don't descend below p_exc = 0 until next substep
               !if dm/dt hits m=1.d0, follow boundary
               dtp = min(dtr,-(1.d0/kp0)*log(c_d0/(-kp0*p_exc0+c_d0)))
               dtm = min(dtp,(1.d0-m_0)/(km0*p_exc0*m_0))
               !integrate dp_exc/dt = -kp0 p_exc + c_d with coefficients at t=0.
               p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
               m = m_0 + km0*m_0*0.5d0*(p_exc0+p_exc)*dtm
               dts = dtp
            elseif ((-kp0*p_exc0+c_d0)<0.d0) then !p_exc is decreasing above/right of nullcline
               debugloop=debugloop + 200
               if (c_d0/kp0<-p_eq0) then !assymptotic limit of p<0
                  debugloop=debugloop + 10
                  dtp = min(dtr,-(1.d0/kp0)*log((c_d0+kp0*p_eq0)/(-kp0*p_exc0+c_d0)))
               else
                  debugloop=debugloop + 20
                  dtp = dtr
               endif 
               !p_exc decrease bounded by p=0 or nullcline if dtp>=dtr
               p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
               p_exc_ave = 0.5d0*(p_exc+p_exc0)
               if (p_exc_ave>0.d0) then 
                  debugloop=debugloop + 1
                  !p_exc still above m nullcline (p=p_eq) even though below p_eq_c, m increasing
                  !happens only if still near p_exc=0.
                  dtm = min(dtr,(1.d0-m_0)/(km0*p_exc_ave*m_0))
               else
                  debugloop=debugloop + 2  
                  m_lower = m_eq + max(0.d0,(1.d0/3.d0)*kp0*h0*(-p_exc_ave)/(vnorm*alphainv))
                  dtm = min(dtr,(m_lower-m_0)/(km0*p_exc_ave*m_0))
               endif
               !dts = min(dts,dtr)
               !m = m_0 + km0*m_0*0.5d0*(p_exc0+p_exc)*dts
               m = m_0 + km0*m_0*p_exc_ave*dtm
               dts = dtm
            else !p_exc is increasing below/left of nullcline
               debugloop=debugloop + 300
               if (c_d0<0.d0) then !between nullcline and m_eq
                  debugloop=debugloop + 10
                  !no limit needed for pressure timestep (p_exc0 < p_exc-->c_d/kp <0 ), limit based on m>m_c
                  dtp = dtr
                  dtm = (m_c-m_0)/(km0*p_exc0*m_0)
                  dts = min(dtp,dtm)
                  dts = max(dts,0.d0)
                  p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
                  m = m_0 + km0*m_0*p_exc0*dts
               else !left of m_eq, p_exc0 < p_exc-->c_d/kp > 0, bound by 0 with dtp
                  debugloop=debugloop + 20
                  dtm = (m_c-m_0)/(km0*p_exc0*m_0)
                  dtp = -(1.d0/kp0)*log(-c_d0/(kp0*p_exc0-c_d0))
                  dts = min(dtp,dtm)
                  dts = min(dtr,dts)
                  dts = max(dts,0.d0)
                  !integrate dp_exc/dt = -kp0 p_exc + c_d with coefficients at t=0.
                  p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
                  !m = m_0 + km0*m_0*0.5d0*(p_exc0+p_exc)*dts
                  m = m_0 + km0*m_0*p_exc0*dts
               endif
            endif
            call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
            dtk = dts
            dtr = dtr-dts
            
         case(4) !LL quadrant
            ! m is strictly decreasing, p strictly increasing
            debugloop = 4000
            !don't exceed p_exc = 0 until next substep
            dtp = min(dtr,-(1.d0/kp0)*log(-c_d0/(kp0*p_exc0-c_d0)))
            dtm = dtp
            !integrate dp_exc/dt = -kp0 p_exc + c_d with coefficients at t=0.
            p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
            m = m_0*exp(dtm*km0*0.5d0*(p_exc0+p_exc))
            !m = m_0*exp(dts*km0*p_exc)
            dts = dtp
            call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
            dtk = dts
            dtr = dtr-dts
         end select

         if (debug) then
            !call setvars(h,u,v,m,p,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)
            if ((m<m_c).and.(p>=p_eq_c)) then
               quad1=1
            elseif ((m>=m_c).and.(p>p_eq_c)) then
               quad1=2
            elseif  ((m>m_c).and.(p<=p_eq_c)) then
               quad1=3
            elseif ((m<=m_c).and.(p<p_eq_c)) then
               quad1=4
            endif
         endif

         if (debug) then
            call setvars(h,u,v,m,p,gz,rho,kperm,alphainv,sig_0,sig_eff1,m_eq1,tanpsi,tau)
            
            if (dtk/dtr<1.d-6.and.quad0==quad1) then
               write(*,*) '---------------SRC WARNING: SMALL TIMESTEP---------->>>>>>'
               write(*,*) 'dtk,dtr:', dtk,dtr
               write(*,*) 'dtm,dtp:', dtm,dtp
               write(*,*) 'debugloop:', debugloop
               write(*,*) ' quad0,quad1', quad0,quad1
               write(*,*) 'm_0, m_c:', m_0, m_c
               write(*,*) 'm_0 - m_c:', m_0 - m_c
               write(*,*) 'm_0 - m_eq0:', m_0 - m_eq
               write(*,*) 'meq0 - m_c:', m_eq - m_c
               write(*,*) 'm1, m_eq1', m, m_eq1
               write(*,*) 'm1 - m_c:', m - m_c
               write(*,*) 'm1 - m_eq1:', m - m_eq1
               write(*,*) 'meq1 - m_c:', m_eq1 - m_c
               write(*,*) 'sig_eff0,sig_eff1:', sig_eff,sig_eff1
               write(*,*) 'p0,p1:', p0,p
               write(*,*) 'p_exc0,p_exc:', p_exc0,p_exc
               write(*,*) 'h0,h_c,diff', h0,h_c,h0-h_c
               write(*,*)  'p,p_eq0,diff', p0,p_eq0,p-p_eq0
               write(*,*)  'p,p_eq_c,diff', p0,p_eq_c,p-p_eq_c
               write(*,*)  'p_eq0,p_eq_c,diff', p_eq0,p_eq_c,p_eq0-p_eq_c
               !stop
               write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
               !stop
            endif

            if (quad1-quad0==-1) then
               write(*,*) '------------SRC WARNING: COUNTERCLOCKWISE ---------->>>>>>>>'
               write(*,*) ' quad0,quad1',quad0,quad1
               write(*,*) 'debugloop:',debugloop
               write(*,*) 'm_0,m,m_eq,m_c',m_0,m,m_eq,m_c
               write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
               !stop
            endif

            if (dtm<0.d0.or.dtp<0.d0) then
               write(*,*) '------------SRC WARNING: negative dt ---------->>>>>>>>'
               write(*,*) 'dtk,dtr:', dtk,dtr
               write(*,*) 'dtm,dtp:', dtm,dtp
               write(*,*) ' quad0,quad1',quad0,quad1
               write(*,*) 'debugloop:',debugloop
               write(*,*) 'm_0,m,m_eq,m_c',m_0,m,m_eq,m_c
               write(*,*) 'p_exc0,-1.d-3*rhoh*gz', p_exc0,-1.d-3*rhoh*gz
               write(*,*) 'kp0,c_d0,-kp0*p_exc0+c_d0',kp0,c_d0,-kp0*p_exc0+c_d0
               write(*,*) 'h0,h_c,diff', h0,h_c,h0-h_c
               write(*,*)  'p,p_eq0,diff', p0,p_eq0,p-p_eq0
               write(*,*)  'p,p_eq_c,diff', p0,p_eq_c,p-p_eq_c
               write(*,*)  'p_eq0,p_eq_c,diff', p_eq0,p_eq_c,p_eq0-p_eq_c
               write(*,*)  'if small rounding error use max(dt,0)'
               write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
               stop
            endif
         endif
         
         return
         end subroutine mp_update_FE_4quad

   !====================================================================
   ! subroutine mp_update_FEexp: integrate dp/dt,dm/dt by a hybrid 
   ! FE and exponential integration
   !====================================================================

      subroutine mp_update_FEexp(dt,h,u,v,m,p,rhoh,gz,dtk)

      use digclaw_module, only: rho_f,rho_s,sigma_0,mu,alpha,setvars,qfix,qfix_cmass
      use geoclaw_module, only: grav,drytolerance

      implicit none

      !i/o
      real(kind=8), intent(inout) :: h,m,p
      real(kind=8), intent(in)  :: u,v,rhoh,dt
      real(kind=8), intent(in)  :: gz
      real(kind=8), intent(out) :: dtk

      !local
      real(kind=8) :: h0,p0,m_0,p_eq0,p_exc0,sig_eff,sig_0,vnorm,m_eq
      real(kind=8) :: kappa,S,rho,rho0,tanpsi,D,tau,sigbed,kperm,phi
      real(kind=8) :: mkrate,plambda,alphainv,c_dil,p_exc,dtm,dtp


      vnorm = sqrt(u**2 + v**2)

      !explicit integration (hybrid FE and explicit exponential solution)---------------------------------------------------
      ! q1 = q0 + dt*f(q0)
      call setvars(h,u,v,m,p,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)
      h0 = h
      m_0 = m
      rho0 = m_0*(rho_s-rho_f)+rho_f
      p0 = p
      p_eq0 = rho_f*gz*h0
      p_exc0 = p0 - p_eq0

      !determine coefficients for update
      ! dm/dt = (2*k*rho^2/mu(rhoh)^2)*p_exc*m
      ! dp_eq/dt = -plamda*p_eq + c_dil see George & Iverson 2014
      mkrate = ((2.d0*kperm*rho0**2)/(mu*rhoh**2))*p_exc0
      c_dil = -3.d0*vnorm*(alphainv*rho0/(rhoh))*tanpsi
      plambda = (2.d0*kperm/(h0*mu))*(((6.d0*alphainv*rho0)/(4.d0*rhoh)) &
               - ((3.d0*rho_f*gz*h0*(rho-rho_f))/(4.d0*rhoh))) !should be always >= 0.d0, but fix if small rounding error below
      
      !calculate stable dt and update 
      dtk = dt
      dtm = dt
      dtp = dt
      if (mkrate*m_0>0.d0) then
         dtm = min(dt,max(1.d0-m_0,0.d0)/(mkrate*m_0))
      endif
      if (c_dil>0.d0) then
         dtp = min(dt,max(rhoh*gz - p0,0.d0)/c_dil)
      elseif (c_dil<0.d0) then
         dtp = min(dt,max(p0,0.d0)/abs(c_dil))
      endif

      !if (dtm==0.d0) then !m_0 = 1.d0 and mkrate >0, adjust pressure if possible
      !   if (dtp>0.d0) then !integrate pressure for dtp=dtk
      !      dtk = dtp
      !      p_exc = p_exc0 + dtp*c_dil
      !      p_exc = p_exc*exp(-max(plambda,0.d0)*dtp)
      !   else !pressure is zero or lithostatic. Only physical response is relaxation
      !        !(flow wants to dilate and p=0 or contract and p=rhogh)
      !      dtk = dt
      !      p_exc = p_exc*exp(-max(plambda,0.d0)*dtk)
      !   endif
      !elseif (dtp==0.d0) then
      !   dtk = dtm
      !   if (mkrate<=0.d0) then !integrate exponential
      !      m = m_0*exp(mkrate*dtk)
      !   else 
      !      m = min(m_0 + dtk*mkrate*m_0,1.d0)
      !   endif
      !   p_exc = p_exc*exp(-max(plambda,0.d0)*dtk)
      !else !dtm,dtp>0
      !   dtk = min(dtm,dtp)
      !   if (mkrate<=0.d0) then !integrate exponential
      !      m = m_0*exp(mkrate*dtk)
      !   else 
      !      m = min(m_0 + dtk*mkrate*m_0,1.d0)
      !   endif
      !   p_exc = p_exc0 + dtk*c_dil
      !   p_exc = p_exc*exp(-max(plambda,0.d0)*dtk)
      !endif

      dtk = min(dtm,dtp)
      if (mkrate<=0.d0) then !integrate exponential
            m = m_0*exp(mkrate*dtk)
      else 
            m = min(m_0 + dtk*mkrate*m_0,1.d0)
      endif
      p_exc = p_exc0 + dtk*c_dil
      p_exc = p_exc*exp(-max(plambda,0.d0)*dt)

      !recapture p, rho, h
      rho = m*(rho_s-rho_f)+rho_f
      h = rhoh/rho
      p = rho_f*gz*h + p_exc

      !call qfix_cmass(m,p,h,rho,u,v,rhoh,gz)
      
      dtk = dt
      !if (dtk<dt/1000.d0) then !integration stalled at boundary of physical space move on.
      !   dtk = dt
      !   write(*,*) 'exiting src integration:'
      !   write(*,*) 'h0,h,m0,m,rho0, rho,rhoh', h0,h,m_0,m,rho0,rho,rhoh
      !   write(*,*) 'dt,dtm,dtp', dt,dtm,dtp
      !   write(*,*) 'mkrate',mkrate
      !endif
      
      
      return
      end subroutine mp_update_FEexp

   !====================================================================
   ! subroutine mp_update_trapezoid: integrate dp/dt,dm/dt by a hybrid 
   ! trapezoid rule - some FE, some exponential, some implicit
   ! rationale: A-stability is only desired in relaxation cases I think.
   !====================================================================

      subroutine mp_update_trapezoid(dt,h,u,v,m,p,rhoh,gz)

      use digclaw_module, only: rho_f,rho_s,sigma_0,mu,alpha,setvars,qfix,qfix_cmass
      use geoclaw_module, only: grav,drytolerance

      implicit none

      !i/o
      real(kind=8), intent(inout) :: h,m,p
      real(kind=8), intent(in)  :: u,v,rhoh,dt
      real(kind=8), intent(in)  :: gz
      

      !local
      real(kind=8) :: h0,p0,m_0,p_eq0,p_exc0,sig_eff,sig_0,vnorm,m_eq
      real(kind=8) :: kappa,S,rho,tanpsi,D,tau,sigbed,kperm,phi
      real(kind=8) :: mkrate,plambda,dtk,alphainv,c_dil,p_exc1,m1
      real(kind=8) :: mstar,hstar,pstar,rhostar,h_n,m_n,p_n,p_eq,rho_n
      real(kind=8) :: p_exc_star,p_exc_n,p_exc,p_excstar,hu,hv,hm

      real(kind=8) :: convtol,deltaiter
      integer :: maxiter,iter,exitstatus

      vnorm = sqrt(u**2 + v**2)

      !explicit integration (hybrid FE and explicit exponential solution)---------------------------------------------------
      ! q* = q0 + 1/2dt*f(q0)
      call setvars(h,u,v,m,p,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)
      h0 = h
      m_0 = m
      p0 = p
      p_eq0 = rho_f*gz*h0
      p_exc0 = p0 - p_eq0
      !sig_0 = sigma_0  !
      !sig_0 = 0.5d0*alpha*rho_f*gz*rhoh*(rho_s-rho_f)/(rho**2)
      !sig_0 = 0.5d0*alpha*gz*rhoh*(rho_s-rho_f)/(rho)
       !later possibly needed to ensure stability of dp_exc/dt.
      !sig_eff = max(rhoh*gz-p0,0.d0) !max is fix small rounding error for double precision problems

      ! dm/dt = (2*k*rho^2/mu(rhoh)^2)*p_exc*m
      dtk = 0.5d0*dt
      mkrate = ((2.d0*kperm*rho**2)/(mu*rhoh**2))*p_exc0
      if (mkrate<=0.d0) then !integrate exponential
         mstar = m_0*exp(mkrate*dtk)
      else 
         mstar = min(m_0 + dtk*mkrate*m_0,1.d0)
      endif
      !dp_eq/dt = see George & Iverson 2014
      !alphainv = m_0*(sig_eff + sig_0)/alpha
      c_dil = -3.d0*vnorm*(alphainv*rho/(rhoh))*tanpsi
      p_excstar = p_exc0 + dtk*c_dil
      plambda = (2.d0*kperm/(h*mu))*(((6.d0*alphainv*rho)/(4.d0*rhoh)) &
               - ((3.d0*rho_f*gz*h*(rho-rho_f))/(4.d0*rhoh))) !should be always >= 0.d0, but fix if small rounding error below
      p_excstar = p_excstar*exp(-max(plambda,0.d0)*dtk)
      rhostar = mstar*(rho_s-rho_f)+rho_f
      hstar = rhoh/rhostar
      !recapture p for fix
      pstar = rho_f*gz*hstar + p_excstar
      call qfix_cmass(hstar,mstar,p,rhostar,p_excstar,hu,hv,hm,u,v,rhoh,gz)
      !implicit integration ---------------------------------------------------------------------------------------------
      ! q* = q0 + 1/2dt*f(q0)
      ! q^n+1 = q^* +  1/2dt*f(q^n+1)
      !note: h,m,p = h*,m*,p* ```````````````
      dtk = 0.5d0*dt
      maxiter = 1000
      convtol = 1.d-3
      exitstatus = 0
      do iter = 1,maxiter
         
         call setvars(h,u,v,m,p,gz,rho,kperm, &
                        alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)
         p_eq = rho_f*gz*h
         p_exc = p - p_eq
         !sig_0 = 0.5d0*alpha*rho_f*gz*rhoh*(rho_s-rho_f)/(rho**2)
         !sig_eff = max(rhoh*gz-p,0.d0) !max is fix small rounding error for double precision problems
         ! dm/dt = f_m(m,p) = (2*k*rho^2/mu(rhoh)^2)*p_exc*m
         mkrate = ((2.d0*kperm*rho**2)/(mu*rhoh**2))*p_exc
         !dp_eq/dt = f_p(m,p)  (see George & Iverson 2014 for f_p(m,p))
         !alphainv = m*(sig_eff + sig_0)/alpha
         c_dil = -3.d0*vnorm*(alphainv*rho/(rhoh))*tanpsi
         plambda = (2.d0*kperm/(h*mu))*(((6.d0*alphainv*rho)/(4.d0*rhoh)) &
               - ((3.d0*rho_f*gz*h*(rho-rho_f))/(4.d0*rhoh))) !should be always >= 0.d0, but fix if small rounding error below
         
         m_n = mstar + dtk*mkrate*m  
         p_exc_n = p_exc_star -dtk*max(plambda,0.d0)*p_exc + dtk*c_dil
         !check for convergence of fixed-point iteration
         deltaiter = (m-m_n)**2 + ((p_exc-p_exc_n)/(rho_f*gz*h))**2
         if (deltaiter<convtol) then
            m = min(m_n,1.d0)
            m = max(m_n,0.d0)
            rho = m*(rho_s-rho_f) + rho_f
            h = rhoh/rho
            p = min(rho_f*gz*h + p_exc_n,rhoh*gz)
            p = max(p,0.d0)
            exitstatus = 1
            exit
         else
            m = min(m_n,1.d0)
            m = max(m_n,0.d0)
            rho = m*(rho_s-rho_f) + rho_f
            h = rhoh/rho
            p = rho_f*gz*h + p_exc_n
         endif
      enddo

      !if fixed point fails
      if (exitstatus==0) then
         dtk = dt
         write(*,*) 'fixed-point fail'
         write(*,*) 'delta: ',deltaiter
         !write(*,*) 'h0,h,m:',h0,h,m
         !trapezoid failed use midpoint
         !solve q^{n+1} = q0 + dt*f(qstar)
         call setvars(hstar,u,v,mstar,pstar,gz,rhostar, &
                  kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)
         p_eq = rho_f*gz*hstar
         !p_exc = pstar - p_eq
         !sig_0 = 0.5d0*alpha*rho_f*gz*rhoh*(rho_s-rho_f)/(rho**2)
         !sig_eff = max(rhoh*gz-pstar,0.d0) !max is fix small rounding error for double precision problems
         ! dm/dt = f_m(m,p) = (2*k*rho^2/mu(rhoh)^2)*p_exc*m
         mkrate = ((2.d0*kperm*rhostar**2)/(mu*rhoh**2))*p_excstar
         m = m_0 + dtk*mkrate*mstar
         h = rhoh/(m*(rho_s-rho_f)+rho_f)
         !dp_eq/dt = f_p(m,p)  (see George & Iverson 2014 for f_p(m,p))
         !alphainv = mstar*(sig_eff + sig_0)/alpha
         c_dil = -3.d0*vnorm*(alphainv*rhostar/(rhoh))*tanpsi
         p_exc = p_exc0 + dtk*c_dil
         plambda = (2.d0*kperm/(hstar*mu))*(((6.d0*alphainv*rhostar)/(4.d0*rhoh)) &
               - ((3.d0*rho_f*gz*hstar*(rhostar-rho_f))/(4.d0*rhoh))) !should be always >= 0.d0, but fix if small rounding error below
         
         p_exc = p_exc0 -dtk*max(plambda,0.d0)*p_excstar + c_dil*dtk
         p = rho_f*gz*h + p_exc

         m = min(m_n,1.d0)
         m = max(m_n,0.d0)
         rho = m*(rho_s-rho_f) + rho_f
         h = rhoh/rho
         p = min(rho_f*gz*h + p_exc_n,rhoh*gz)
         p = max(p,0.d0)

      else
         !if(iter>1) write(*,*) 'iter ', iter
      endif

      return
      end subroutine mp_update_trapezoid


