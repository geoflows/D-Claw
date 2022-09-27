

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
      double precision :: zeta,p_hydro,p_litho,p_eq,krate,gamma,dgamma
      double precision :: vnorm,hvnorm,theta,dtheta,w,taucf,fsphi,hvnorm0
      double precision :: shear,sigebar,pmtanh01,rho_fp,seg
      double precision :: b_xx,b_yy,b_xy,chi,beta
      double precision :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      double precision :: vlow,m2,vreg,slopebound
      double precision :: b_eroded,b_remaining,dtcoeff
      integer :: i,j,ii,jj,jjend,icount
      logical :: ent

            !source fountain
      integer :: iBeg,iEnd,jBeg,jEnd,numCellsX,numCellsY,srcI,srcJ,numCells
      integer :: xCells,yCells
      double precision :: s_angle,s_h,s_m0,s_q,s_Qp,s_qx,s_qy
      double precision :: s_slope,s_slope_x,s_slope_y,s_tend
      double precision :: s_vel,s_velx,s_vely,s_Vtot,s_xloc,s_yloc
      double precision :: srcXComp,srcYComp,x,y

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

            jjend = 1
            dti = dt!/real(jjend,kind=8)
            !do jj=1,jjend


            !Ocean material removal
            !RPJ: 04/28/22
            ! If material reaches the ocean (topology = 0 elevation)
            ! and has moved away from the shoreline (all surrounding topo = 0)
            ! then remove the material to stop the slow creep into the flat ocean topo
            !topo = 0.0
            !do ii = -1,1
            !   do jj = -1,1
            !      topo = topo + aux(i+ii,j+jj,1)
            !   enddo
            !enddo
            !if (topo .lt. 1.0) then
            !   q(i,j,1) = 0.0
            !   q(i,j,2) = 0.0
            !   q(i,j,3) = 0.0
            !   q(i,j,4) = 0.0
            !   q(i,j,5) = 0.0
            !   cycle
            !endif


            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            !integrate momentum source term
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

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

            if (p_initialized==0) cycle

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)


            vnorm = sqrt(u**2.0 + v**2.0)

            !integrate shear-induced dilatancy
            sigebar = rho*gmod*h - p + sigma_0
            shear = 2.0*vnorm/h
            krate = 1.5*shear*m*tanpsi/alpha
            sigebar = sigebar*exp(krate*dti)
            !p = rho*gmod*h + sigma_0 - sigebar
            if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
               p = p - dti*3.0*vnorm*tanpsi/(h*compress)
            endif

            !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            !call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
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

            !integrate pressure relaxation
            !if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
            !   zeta = 3.d0/(compress*h*2.0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            !else
            !   zeta = (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            !endif
            zeta = ((m*(sigbed +  sigma_0))/alpha)*3.d0/(h*2.0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            krate=-zeta*2.0*kperm/(h*max(mu,1.d-16))
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
            !p_eq = max(p_eq,0.0)
            !p_eq = min(p_eq,p_litho)

            p = p_eq + (p-p_eq)*exp(krate*dti)


            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)


            krate = D*(rho-rho_fp)/rho
            hu = hu*exp(dti*krate/h)
            hv = hv*exp(dti*krate/h)
            hm = hm*exp(-dti*D*rho_fp/(h*rho))
            h = h + krate*dti

            !enddo

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

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
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  if (h<tol) cycle
                  pm = q(i,j,6)/h
                  pm = max(0.0,pm)
                  pm = min(1.0,pm)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

                  if (h.lt.tol) then
                     q(i,j,2)=0.d0
                     q(i,j,3)=0.d0
                  else
                     beta = 1.0-m !tan(1.5*p/(rho*gmod*h))/14.0
                     gamma= beta*dsqrt(hu**2 + hv**2)*(gmod*coeff**2)/(h**(7.0/3.0))
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
      ! Specifically for post-fire debris flow hazard assessment
      ! assumes triangle hydrograph as several locations
      ! Peak volume and concentration given
      !------------------------------
      !write(*,*) ' '
      !!write(*,*) ' SFA= ', src_fountain_active
      !write(*,*) ' SFNsh= ', src_ftn_num
      !write(*,*) ' SFNsr= ', src_ftn_num_sr
      !write(*,*) ' SFt= ', src_ftn_end_time
      !write(*,*) ' t= ', t

      if (src_fountain_active .eqv. .True.) then
         if (src_ftn_num .eq. 0) then
            src_fountain_active = .False.
         endif
         if (t .gt. src_ftn_end_time) then
            src_fountain_active = .False.
         endif


         do ii = 1,src_ftn_num

            ! inefficient testing of needing to be at the finest grid
            if (dx .gt. 21.0) then
               cycle
            endif
            s_Vtot = src_ftn_vtot(ii)
            s_Qp = 0.1 * (s_Vtot ** 0.833)
            if (s_Qp .lt. 0.0001) then
               cycle
            endif
            s_tend = 2.0 * s_Vtot / s_Qp

            if (t .gt. s_tend) then
               cycle
            endif

            s_xloc = src_ftn_xloc(ii)
            s_yloc = src_ftn_yloc(ii)

            srcI = floor((s_xloc-xlower)/dx)
            srcJ = floor((s_yloc-ylower)/dy)

            if(srcI.lt.1 .or. srcJ.lt.1 .or. srcI.gt.mx .or. srcJ.gt.my) then
               !source in a different grid
               cycle
            endif


            s_slope_x = (aux(srcI+1,srcJ,1)-aux(srcI-1,srcJ,1))/(2.d0*dx)
            s_slope_y = (aux(srcI,srcJ+1,1)-aux(srcI,srcJ-1,1))/(2.d0*dy)
            !s_slope = sqrt((s_slope_x * srcXComp)**2 + (s_slope_y * srcYComp)**2)
            s_slope = sqrt((s_slope_x)**2 + (s_slope_y)**2)
            if(s_slope .lt. 0.0001) then
               ! if no slope use a value for the component calculations (which will still be 0 in each direction)
               s_slope = 0.1
            endif
            srcXComp = s_slope_x/s_slope
            srcYComp = s_slope_y/s_slope



            !!!! Only adding volume
            s_q = 2.0*s_Qp*(t/s_tend-0.5)*sign(1.0d0,(0.5-t/s_tend))+s_Qp ! triangle function centered at 0.5

            if (s_q .lt. 1e-4) then
               cycle
            endif

            s_vel = abs(2.1 * s_Qp**0.33 * s_slope**0.33) !Rickenmann Eq 21
            s_velx = s_vel * srcXComp
            s_vely = s_vel * srcYComp
            s_m0 = src_ftn_m0(ii)


            if (dx .lt. 4.0) then
               ! 3.3 grid size
               !numCells = 9
               !check for grid edge, if needed adjust to keep all source on same grid
               if(srcI.gt.1 .and. srcI.lt.mx) then
                  xcells = -1
               else if(srcI.eq.1) then
                  xcells = 0
               else
                  xcells = -2
               endif
               if(srcJ.gt.1 .and. srcJ.lt.my) then
                  ycells = -1
               else if(srcI.eq.1) then
                  ycells = 0
               else
                  ycells = -2
               endif

               ! add volume equally across all grid cells
               do i=xcells,xcells+2
                  do j=ycells,ycells+2
                     q(srcI+i,srcJ+j,1) = q(srcI+i,srcJ+j,1) + s_q*dt/(9*dx*dy)
                     q(srcI+i,srcJ+j,2) = q(srcI+i,srcJ+j,2) + s_q*dt/(9*dx*dy)*s_velx
                     q(srcI+i,srcJ+j,3) = q(srcI+i,srcJ+j,3) + s_q*dt/(9*dx*dy)*s_vely
                     q(srcI+i,srcJ+j,4) = q(srcI+i,srcJ+j,4) + s_q*dt/(9*dx*dy)*s_m0
                  end do
               end do
            else if(dx .lt. 6.0) then
               ! 5.0 m cells
               !numcells = 4
               ! check for grid edge, if not find direction of cell closest to source loc
               if(srcI.gt.1 .and. srcI.lt.mx) then
                  if(((s_xloc-xlower)-(srcI*dx))/dx .lt. 0.5) then
                     xcells = -1
                  else
                     xcells = 0
                  endif
               else if(srcI.eq.1) then
                  xcells = 0
               else
                  xcells = -1
               end if

               if(srcJ.gt.1 .and. srcJ.lt.my) then
                  if(((s_yloc-ylower)-(srcJ*dy))/dy .lt. 0.5) then
                     ycells = -1
                  else
                     ycells = 0
                  endif
               else if(srcJ.eq.1) then
                  ycells = 0
               else
                  ycells = -1
               end if

               do i=xcells,xcells+1
                  do j=ycells,ycells+1
                     q(srcI+i,srcJ+j,1) = q(srcI+i,srcJ+j,1) + s_q*dt/(4*dx*dy)
                     q(srcI+i,srcJ+j,2) = q(srcI+i,srcJ+j,2) + s_q*dt/(4*dx*dy)*s_velx
                     q(srcI+i,srcJ+j,3) = q(srcI+i,srcJ+j,3) + s_q*dt/(4*dx*dy)*s_vely
                     q(srcI+i,srcJ+j,4) = q(srcI+i,srcJ+j,4) + s_q*dt/(4*dx*dy)*s_m0
                  end do
               end do

            else
               ! grid is 10 or 20m
               ! numcells = 1
               q(srcI,srcJ,1) = q(srcI,srcJ,1) + s_q*dt/(dx*dy)
               q(srcI,srcJ,2) = q(srcI,srcJ,2) + s_q*dt/(dx*dy)*s_velx
               q(srcI,srcJ,3) = q(srcI,srcJ,3) + s_q*dt/(dx*dy)*s_vely
               q(srcI,srcJ,4) = q(srcI,srcJ,4) + s_q*dt/(dx*dy)*s_m0
            end if



            !!!!
!            s_angle = src_ftn_angle(ii)
!            ! north is 0 degrees for aspect
!            srcXComp = sin(s_angle*3.14/180)
!            srcYComp = cos(s_angle*3.14/180)

!            s_slope_x = (aux(srcI+1,srcJ,1)-aux(srcI-1,srcJ,1))/(2.d0*dx)
!            s_slope_y = (aux(srcI,srcJ+1,1)-aux(srcI,srcJ-1,1))/(2.d0*dy)
!            !s_slope = sqrt((s_slope_x * srcXComp)**2 + (s_slope_y * srcYComp)**2)
!            s_slope = sqrt((s_slope_x)**2 + (s_slope_y)**2)
!            !!

!            if(s_slope .lt. 0.0001) then
!               !write(*,*) '\n slope =  ,ii= ', s_slope,ii
!               s_slope = 0.1
!            endif
!            srcXComp = s_slope_x/s_slope
!            srcYComp = s_slope_y/s_slope

!            s_q = 2.0*s_Qp*(t/s_tend-0.5)*sign(1.0d0,(0.5-t/s_tend))+s_Qp ! triangle function centered at 0.5
!            if (s_q .lt. 1e-6) then
!               cycle
!            endif
!            s_qx = s_q * srcXComp
!            s_qy = s_q * srcYComp
!            s_vel = abs(2.1 * s_q**0.33 * s_slope**0.33) !Rickenmann Eq 21
!            s_velx = s_vel * srcXComp
!            s_vely = s_vel * srcYComp

!            !write(*,*) '\n sq= ,ii= , t= ,dx =', s_q,ii, t, dx
!            !write(*,*) '\n sl= ', src_ftn_length
!            !write(*,*) '\n dx= ', dx
!            if (src_ftn_length .lt. dx) then
!               s_h = abs(s_q / (dx * s_vel))
!               !write(*,*) '\n sq= ,sh= ,ii= , t= ,dx =', s_q, s_h,ii, t, dx
!               if (s_h .lt. tol) then
!                  cycle
!               endif
!               q(srcI,srcJ,1) = s_h
!               !calculate velocity for correct flux along each axis
!               q(srcI,srcJ,2) = q(srcI,srcJ,1) * s_qx/(dx*s_h)
!               q(srcI,srcJ,3) = q(srcI,srcJ,1) * s_qy/(dy*s_h)
!               q(srcI,srcJ,4) = q(srcI,srcJ,1) * s_m0
!               q(srcI,srcJ,5) = (rho_F*(1-s_m0)+rho_S*s_m0)*grav*q(srcI,srcJ,1)
!            else
!               xCells = floor(src_ftn_length*srcXComp/dx)
!               yCells = floor(src_ftn_length*srcYComp/dy)
!               iBeg = max(srcI - xCells,1)
!               iEnd = min(srcI + xCells,mx)
!               jBeg = max(srcJ - yCells,1)
!               jEnd = min(srcJ + yCells,my)

!               !Loop through both directions to count faces
!               !Sum topography for volume flux calcs
!               numCellsX = 0
!               numCellsY = 0
!               do i=iBeg,iEnd
!                  x=xlower + (i-0.5d0)*dx
!                  if (abs(x-s_xloc).lt.(src_ftn_length*srcXComp)) then
!                     numCellsY = numCellsY + 1
!                  endif
!               enddo

!               !Loop other order to find other face count
!               do j=jBeg,jEnd
!                  y=ylower + (j-0.5d0)*dy
!                  if (abs(y-s_yloc)<(src_ftn_length*srcYComp)) then
!                     numCellsX = numCellsX + 1
!                  endif
!               enddo

!               ! find h for total flux across x and y faces to equal flow rate
!               ! gives correct total flux into domain
!               s_h = abs(s_q /(numCellsX*dx*s_velx + numCellsY*dy*s_vely))

!               if (s_h .lt. tol) then
!                  cycle
!               endif

!               do i=iBeg,iEnd
!                  x=xlower + (i-0.5d0)*dx
!                  if (abs(x-s_xloc).lt.(src_ftn_length*srcXComp)) then
!                     do j=jBeg,jEnd
!                        y=ylower + (j-0.5d0)*dy
!                        if (abs(y-s_yloc).lt.(src_ftn_length*srcYComp)) then
!                           q(i,j,1) = s_h
!                           q(i,j,2) = q(i,j,1) * s_velx
!                           q(i,j,3) = q(i,j,1) * s_vely
!                           q(i,j,4) = q(i,j,1) * s_m0
!                           q(i,j,5) = (rho_F*(1-s_m0)+rho_S*s_m0)*grav*q(i,j,1)
!                        endif
!                     enddo
!                  endif
!               enddo

!            endif  ! src length

         enddo ! source number

      endif ! source fountain active

      return
      end
