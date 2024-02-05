

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
      double precision :: gacc,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi,dti,gz,gx
      double precision :: D,tau,sigbed,kperm,compress,pm,p_exc,sig_0
      double precision :: zeta,p_eq,mkrate,pkrate,gamma,dgamma,c_dil,alphainv
      double precision :: vnorm,hvnorm,theta,dtheta,w,hvnorm0
      double precision :: shear,sigebar,pmtanh01,rho_fp,seg
      double precision :: b_xx,b_yy,b_xy,chi,beta
      double precision :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      double precision :: vlow,m2,vreg,slopebound
      double precision :: b_eroded,b_remaining,dtcoeff
      integer :: i,j,ii,jj,jjend,icount
      logical :: ent

      !source fountain
      integer :: numCellsX,numCellsY,numCells,rhoh
      integer :: numCellsHalfX,numCellsHalfY,numCellsHalf
      integer :: srcI,srcJ,srcJlo,srcIlo,srcIhi,srcJhi
      double precision :: src_xloc, src_yloc
      double precision :: s_xloclo,s_yloclo
      double precision :: s_xlochi,s_ylochi
      double precision :: s_m0,s_q,s_Qp
      double precision :: s_tend,s_Vtot,s_xloc,s_yloc
      double precision :: srcXComp,srcYComp,x,y

      ! level awareness
      double precision :: dxmin, dymin
      common /comfine/ dxmin,dymin

      double precision, allocatable :: moll(:,:)

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

            !integrate momentum source term
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            rhoh = h*rho !this is invariant in src and always >0 below
            hvnorm0 = sqrt(hu**2 + hv**2)
            vnorm = hvnorm0/h

            if (hvnorm0>0.d0) then
               !integrate dynamic friction !DIG: TO DO - move dynamic friction to Riemann solver
               vnorm = dmax1(0.d0,vnorm - dt*tau/rhoh) !exact solution for Coulomb friction
               vnorm = vnorm*exp(-(1.d0-m)*2.0d0*mu*dt/(h*rhoh)) !exact solution (prior to h change) for effective viscous friction
               ! velocity determined, calculate directions etc. from vnorm
               hvnorm = h*vnorm
               hu = hvnorm*hu/hvnorm0 + gx*h*dt !gx=0 unless bed-normal !DIG: last term should ultimately be in Riemann solver
               hv = hvnorm*hv/hvnorm0
               u = hu/h
               v = hv/h
               vnorm = sqrt(u**2 + v**2)
               ! velocity now constant for remainder of src2
            endif

            if (p_initialized==0) cycle !# dlg: I can't remember why we don't just return in this case but don't want to change now

            !call admissibleq(h,hu,hv,hm,p,u,v,m,theta) ! dlg: don't see why this is needed again
            ! u,v change need new aux vals
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            !integrate m and p, keep rhoh ------------------------------------------------------------
            p_eq = rho_f*gz*h
            p_exc = p - p_eq

            !determine relaxation portion rhs of dp_exc/dt
            sig_0 = sigma_0  !0.5d0*alpha*p_eq*(rho_s-rho_f)/rho !latter possibly needed to ensure krate>0.
            alphainv = m*(rhoh*gz-p + sig_0)
            pkrate = (2.d0*kperm/(h*mu))*(((3.d0*alphainv*rho)/(2.d0*rhoh)) &
               - ((3.d0*rho_f*(gz*rhoh-p_eq))/(4.d0*rhoh)))
            !integrate shear induced dilatancy: dp_exc/dt = c_dil * (m-m_eq)
            c_dil = 3.d0*(alphainv*vnorm/h)*tanpsi
            p_exc = p_exc - dt*c_dil
            !relax p_exc
            p_exc = p_exc*exp(-dt*pkrate)

            !integrate m explicitly
            mkrate = (p-p_eq)*(2.d0*kperm)/(mu*h**2)
            m = m*exp(dt*mkrate)
            m = min(1.d0,m)
            !determine new rho then h
            rho = m*rho_s + (1.d0-m)*rho_f
            h = rhoh/rho
            !recapture p
            p = rho_f*gz*h + p_exc
            !recapture hu,hv,hm
            hu = h*u
            hv = h*v
            hm = h*m

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
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

            !enddo

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



      ! Manning friction------------------------------------------------
      if (ifriction==0) return
      if (coeffmanning>0.d0.and.frictiondepth>0.d0) then
         do i=1,mx
            do j=1,my

               h=q(i,j,1)
               if (h<=frictiondepth) then
                  !# apply friction source term only in shallower water
                  hu=q(i,j,2)
                  hv=q(i,j,3)
                  hm = q(i,j,4)
                  p =  q(i,j,5)
                  phi = aux(i,j,i_phi)
                  theta = aux(i,j,i_theta)
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  if (h<drytolerance) cycle
                  pm = q(i,j,6)/h
                  pm = max(0.0d0,pm)
                  pm = min(1.0d0,pm)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

                  if (h.lt.drytolerance) then
                     q(i,j,1)=0.d0
                     q(i,j,2)=0.d0
                     q(i,j,3)=0.d0
                  else
                     beta = 1.0d0-m
                     gamma= beta*dsqrt(hu**2 + hv**2)*(gz*coeffmanning**2)/(h**(7.d0/3.d0))
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

      return
      end
