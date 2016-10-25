

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
      double precision :: shear,sigebar
      double precision :: b_xx,b_yy,b_xy,chi,beta
      integer :: i,j,ii,jj,jjend,icount,curvature

      double precision, allocatable :: moll(:,:)

      ! check for NANs in solution:
      call check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,2)

      gmod=grav
      coeff = coeffmanning
      tol = 1.e-30  !# to prevent divide by zero in gamma
      curvature = 0 !add friction due to curvature acceleration
      !write(*,*) 'src:init,value',p_initialized,init_pmin_ratio

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
            dti = dt/real(jjend,kind=8)
            do jj=1,jjend

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


            !integrate pressure relaxation
            if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
               zeta = 3.d0/(compress*h*2.0)  + (rho-rho_f)*rho_f*gmod/(4.d0*rho)
            else
               zeta = (rho-rho_f)*rho_f*gmod/(4.d0*rho)
            endif

            krate=-zeta*2.0*kperm/(h*max(mu,1.d-16))
            p_hydro = h*rho_f*gmod
            p_litho = (rho_s*m + (1.d0-m)*rho_f)*gmod*h

            !if (abs(compress*krate)>0.0) then
            !   p_eq = p_hydro + 3.0*vnorm*tanpsi/(compress*h*krate)
            !else
            !   p_eq = p_hydro
            !endif
            !if (abs(pm-.5)>.49) then

            p_eq = p_hydro
            !p_eq = max(p_eq,0.0)
            !p_eq = min(p_eq,p_litho)

            p = p_eq + (p-p_eq)*exp(krate*dti)


            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
            

            krate = D*(rho-rho_f)/rho
            hu = hu*exp(dti*krate/h)
            hv = hv*exp(dti*krate/h)
            hm = hm*exp(-dti*D*rho_f/(h*rho))
            h = h + krate*dti

            enddo

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

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
                     beta = 1.0 !tan(1.5*p/(rho*gmod*h))/14.0
                     gamma= dsqrt(hu**2 + hv**2)*(beta*gmod*coeff**2)/(h**(7.0/3.0))
                     dgamma=1.d0 + dt*gamma
                     q(i,j,2)= q(i,j,2)/dgamma
                     q(i,j,3)= q(i,j,3)/dgamma
                  endif
               endif
            enddo
         enddo
      endif
     ! ----------------------------------------------------------------

      return
      end
