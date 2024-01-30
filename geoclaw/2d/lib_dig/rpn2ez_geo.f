c======================================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  fwave,s,amdq,apdq)
c======================================================================
c
c Solves normal Riemann problems for debris flow equations

c On input, ql contains the state vector at the left edge of each cell
c     qr contains the state vector at the right edge of each cell
c
c This data is along a slice in the x-direction if ixy=1
c     or the y-direction if ixy=2.

c  Note that the i'th Riemann problem has left state qr(i-1,:)
c     and right state ql(i,:)
c  From the basic clawpack routines, this routine is called with
c     ql = qr
c
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      # This Riemann solver is for debris flow equations
!
!
!
!           for information contact
!           David George <dgeorge@uw.edu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use geoclaw_module
      use digclaw_module

      implicit none

      !i/o
      integer maxm,meqn,mwaves,mbc,mx,ixy

      double precision  fwave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision  s(1-mbc:maxm+mbc, mwaves)
      double precision  ql(1-mbc:maxm+mbc, meqn)
      double precision  qr(1-mbc:maxm+mbc, meqn)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  amdq(1-mbc:maxm+mbc, meqn)
      double precision  auxl(1-mbc:maxm+mbc, *)
      double precision  auxr(1-mbc:maxm+mbc, *)


      !local
      integer m,i,mw,maxiter,mhu,nhv,mcapa,icom,jcom,waves
      double precision dtcom,dxcom,dycom,tcom
      double precision wall(3),fw(6,3),sw(3),wave(6,3)
      double precision lamL(3),lamR(3),beta(3)
      logical entropy(5)
      logical rare1,rare2,wallprob,drystate
      !logical entropycorr1,entropycorr2

      double precision drytol,gmod,veltol
      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL
      double precision pR,pL,hmL,hmR,mL,mR,phi_bedL,phi_bedR
      double precision hstar,hstartest,s1m,s2m,bL,bR
      double precision dxdc,dx,taudirL,taudirR
      double precision theta,thetaL,thetaR
      double precision h1M,h2M,hu1M,hu2M,u1M,u2M,heR,heL
      double precision sE1,sE2
      double precision chiHL,chiHR,chiL,chiR,fsL,fsR



      common /cmcapa/  mcapa
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      gmod=grav
      drytol=drytolerance

      !set to true to use an entropy correction
      !entropycorr1 = .false.  !Harten (bound nonlinear waves away from s=0)
      !entropycorr2 = .false. !Harten-Hyman (split expansion shocks) not recommended

      !near-zero velocity tolerance
      veltol = 1.d-3

      waves = 3
      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

!-----------------------Initializing-----------------------------------
         !inform of a bad Riemann problem from the start
c         if((qr(i-1,1).lt.0.d0).or.(ql(i,1) .lt. 0.d0)) then
c            write(*,*) 'Negative input: hl,hr,i=',qr(i-1,1),ql(i,1),i
c         endif

         !Initialize Riemann problem for grid interface

         do mw=1,mwaves
            s(i,mw)=0.d0
            entropy(mw)=.false.
            do m=1,meqn
               fwave(i,m,mw)=0.d0
            enddo
         enddo
         do mw=1,waves
            sw(mw) = 0.d0
            do m=1,6
               wave(m,mw) = 0.d0
               fw(m,mw) = 0.d0
            enddo
         enddo

         !skip problem if in a completely dry area
         if (qr(i-1,1).le.drytol.and.ql(i,1).le.drytol) then
            go to 30
         endif

c        !set normal direction
         if (ixy.eq.1) then
            mhu=2
            nhv=3
            dx = dxcom
            taudirR = auxl(i,i_taudir_x)
            taudirL = auxr(i-1,i_taudir_x)
         else
            mhu=3
            nhv=2
            dx = dycom
            taudirR = auxl(i,i_taudir_y)
            taudirL = auxr(i-1,i_taudir_y)
         endif

         fsL = auxr(i-1,i_fsphi)
         fsR = auxl(i,i_fsphi)

         if (bed_normal.eq.1) then
            thetaL = auxr(i-1,i_theta)
            thetaR = auxl(i,i_theta)
            theta = 0.5d0*(thetaL+thetaR)
            gmod = grav*dcos(0.5d0*(thetaL+thetaR))
         else
            thetaL = 0.d0
            thetaR = 0.d0
            theta = 0.d0
         endif

         !zero (small) negative values if they exist and set velocities
         call admissibleq(ql(i,1),ql(i,mhu),ql(i,nhv),
     &            ql(i,4),ql(i,5),uR,vR,mR,thetaR)

         call admissibleq(qr(i-1,1),qr(i-1,mhu),qr(i-1,nhv),
     &            qr(i-1,4),qr(i-1,5),uL,vL,mL,thetaL)


         !Riemann problem variables
         hL = qr(i-1,1)
         hR = ql(i,1)
         huL = qr(i-1,mhu)
         huR = ql(i,mhu)
         hvL=qr(i-1,nhv)
         hvR=ql(i,nhv)
         hmL = qr(i-1,4)
         hmR = ql(i,4)
         pL = qr(i-1,5)
         pR = ql(i,5)
         bL = auxr(i-1,1)  - qr(i-1,7)
         bR = auxl(i,1) - ql(i,7)
         phi_bedL = auxr(i-1,i_phi)
         phi_bedR = auxl(i,i_phi)
         chiHL = qr(i-1,6)
         chiHR = ql(i,6)

         if (hL.ge.drytol) then
            chiL = chiHL/hL
         endif
         if (hR.ge.drytol) then
            chiR = chiHR/hR
         endif

         !test for wall problem vs. inundation problem
         do mw=1,waves
            wall(mw) = 1.d0
         enddo
         drystate=.false.
         wallprob = .false.
         if (hR.le.drytol) then
            hR = 0.d0
            pR = 0.d0
            mR = mL
            chiR = chiL
            drystate=.true.
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,
     &                                 rare1,rare2,1,drytol,gmod)
            hstartest=dmax1(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
c                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               wallprob=.true.
               hR=hL
               huR=-huL
               hvR=hvL
               hmR=hmL
               bR=bL
               uR=-uL
               vR=vL
               mR=mL
               pR=pL
               chiHR=chiHL
               chiR = chiL
               !thetaL = 0.d0
               !thetaR = 0.d0
            !elseif (hL+bL.lt.bR) then
               !bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            hL = 0.d0
            pL = 0.d0
            mL = mR
            chiL= chiR
            drystate=.true.
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,gmod)
            hstartest=dmax1(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
c               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               wallprob=.true.
               hL=hR
               huL=-huR
               hvL=hvR
               hmL=hmR
               mL = mR
               bL=bR
               uL=-uR
               vL=vR
               pL=pR
               chiHL=chiHR
               chiL = chiR
               !thetaL = 0.d0
               !thetaR = 0.d0
            !elseif (hR+bR.lt.bL) then
               !bL=hR+bR
            endif
         endif

c=================begin digclaw-auxset =================================

c         pm = dsign(1.d0,uR-uL)
c
c         if (hL.gt.drytol.and.hR.gt.drytol) then
c            !this is a completely wet 'normal' Riemann problem
c            call auxeval(hR,uR,vR,mR,pR,phi_bedR,theta,
c     &        kappa,SN,rhoR,tanpsi,DR,tauR,sigbed,kpermR,compressR,pm)
c            call auxeval(hL,uL,vL,mL,pL,phi_bedL,theta,
c     &        kappa,SN,rhoL,tanpsi,DL,tauL,sigbed,kpermL,compressL,pm)
c            D = 0.5d0*(DL + DR)
c            compress = 0.5d0*(compressL + compressR)
c            kperm = 0.5d0*(kpermR + kpermL)
c            tau = 0.5d0*(tauL + tauR)
c            rho = 0.5d0*(rhoR + rhoL)
c         elseif (hR.gt.drytol) then
c            !inundation problem
c            call auxeval(hR,uR,vR,mR,pR,phi_bedR,theta,
c     &        kappa,SN,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
c         elseif (hL.gt.drytol) then
c            !inundation problem
c            call auxeval(hL,uL,vL,mL,pL,phi_bedL,theta,
c     &        kappa,SN,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
c         endif



c================end digclaw-aux========================================

         !modify for inundation problem
c         if (hR.lt.drytol) then
c            sw(1) = sL
c            sw(3) = uL + 2.d0*sqrt(gmod*hL*eps)
c            sw(2) = 0.5d0*(sw(1)+sw(3))
c         endif
c         if (hL.lt.drytol) then
c            sw(1) = uR - 2.d0*sqrt(gmod*hR*eps)
c            sw(3) = sR
c            sw(2) = 0.5d0*(sw(1)+sw(3))
c         endif

         !-- solve Riemann problem
c         call riemann_dig2_aug_sswave(ixy,meqn,mwaves,hL,hR,huL,huR,
c     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
c     &         kappa,rho,kperm,compress,tanpsi,D,tau,
c     &         theta,gamma,eps,dx,sw,fw,wave)

         call riemann_dig2_aug_sswave_ez(ixy,6,3,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         thetaL,thetaR,phi_bedL,phi_bedR,dx,sw,fw,wave,wallprob,
     &         taudirL,taudirR,chiL,chiR,fsL,fsR)


c         call riemann_aug_JCP(1,3,3,hL,hR,huL,
c     &        huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,
c     &                                    drytolerance,gmod,sw,fw)



         !--------------------------------------------------------------

c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)
            do m=1,6
               fw(m,mw)=fw(m,mw)*wall(mw)
            enddo
         enddo

         !Entropy correction (Harten-Hyman) for nonlinear waves---------
c         if (entropycorr1.and.entropycorr2.and.(.not.drystate)) then
c            h1M =  max(hL  + wave(1,1),0.d0)
c            h2M =  max(hR  - wave(1,3),0.d0)
c            hu1M = huL + wave(2,1)
c            hu2M = huR - wave(2,3)
c            if (h1M.gt.drytol) then
c               u1M = hu1M/h1M
c            else
c               u1M = 0.d0
c            endif
c            if (h2M.gt.drytol) then
c               u2M = hu2M/h2M
c            else
c               u2M = 0.d0
c            endif
            !s1M = u1M - sqrt(gmod*h1M*eps)
            !s2M = u2M + sqrt(gmod*h2M*eps)

c            heL = max(hL + bL - max(bL,bR),0.d0)
c            heR = max(hR + bR - max(bL,bR),0.d0)

c            call riemanntype(heL,heR,uL,uR,hstar,s1M,s2M,
c     &                                 rare1,rare2,1,drytol,g)

c            if (sL.lt.-0.5d0*veltol.and.s1M.gt.0.5d0*veltol) then
c               entropy(1) = .true.
c               beta(1) = max(0.d0,s1M-sw(1))/(s1M - sL)
c               lamL(1) = sL
c               lamR(1) = s1M
c               do m=1,meqn
c                  fw(m,1) = fw(m,1)/sw(1)
c               enddo
c            endif
c            if (sR.gt.0.5d0*veltol.and.s2M.lt.-0.5d0*veltol) then
c               entropy(3) = .true.
c               beta(3) = max(0.d0,sR-sw(3))/(sR - s2M)
c               lamL(3) = s2M
c               lamR(3) = sR
c               do m=1,meqn
c                  fw(m,3) = fw(m,3)/sw(3)
c               enddo
c            endif
c         endif

c=======================================================================
c         do mw=1,mwaves
c            s(i,mw)=sw(mw)
c            fwave(i,1,mw)=fw(1,mw)
c            fwave(i,mhu,mw)=fw(2,mw)
c            fwave(i,nhv,mw)=fw(3,mw)
c            fwave(i,4,mw) = fw(4,mw)
c            fwave(i,5,mw) = fw(5,mw)
c         enddo

c============segregation================================================

         s(i,1) = sw(1)
         s(i,2) = sw(2)
         s(i,3) = sw(2)
         s(i,4) = sw(2)
         s(i,5) = sw(3)

         fwave(i,1,1) =   fw(1,1)
         fwave(i,mhu,1) = fw(2,1)
         fwave(i,nhv,1) = fw(3,1)
         fwave(i,4,1)   = fw(4,1)
         fwave(i,5,1) =   fw(5,1)
         fwave(i,6,1) =   fw(6,1)

         fwave(i,1,5) =   fw(1,3)
         fwave(i,mhu,5) = fw(2,3)
         fwave(i,nhv,5) = fw(3,3)
         fwave(i,4,5)   = fw(4,3)
         fwave(i,5,5) =   fw(5,3)
         fwave(i,6,5) =   fw(6,3)

         fwave(i,1,2) =   fw(1,2)
         fwave(i,mhu,2) = fw(2,2)
         fwave(i,nhv,2) = fw(3,2)
         fwave(i,4,2)   = 0.0
         fwave(i,5,2) =  0.0
         fwave(i,6,2) = fw(6,2)

         fwave(i,1,3) =   0.0
         fwave(i,mhu,3) = 0.0
         fwave(i,nhv,3) = 0.0
         fwave(i,4,3)   = fw(4,2)
         fwave(i,5,3) =  0.0
         fwave(i,6,3) =  0.0

         fwave(i,1,4) =   0.0
         fwave(i,mhu,4) = 0.0
         fwave(i,nhv,4) = 0.0
         fwave(i,4,4)   = 0.0
         fwave(i,5,4) =  fw(5,2)
         fwave(i,6,4) =  0.0


c==========Capacity for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
          if (ixy.eq.1) then
             dxdc=(Rearth*pi/180.d0)
          else
             dxdc=Rearth*pi*cos(auxl(i,3))/180.d0
          endif

          do mw=1,mwaves
             !if (s(i,mw) .gt. 316.d0) then
               ! shouldn't happen unless h > 10 km!
              !  write(6,*) 'speed > 316: i,mw,s(i,mw): ',i,mw,s(i,mw)
              !  endif
             s(i,mw)=dxdc*s(i,mw)
             do m=1,meqn
               fwave(i,m,mw)=dxdc*fwave(i,m,mw)
             enddo
          enddo
        endif

c============= compute fluctuations=============================================
 30      continue
         do m=1,meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do  mw=1,mwaves
               if (s(i,mw).lt.0.d0.and.(.not.entropy(mw))) then
                  amdq(i,m) = amdq(i,m) + fwave(i,m,mw)
               elseif (s(i,mw).gt.0.d0.and.(.not.entropy(mw))) then
                  apdq(i,m) = apdq(i,m) + fwave(i,m,mw)
               elseif (entropy(mw)) then !note fwave has already been divided by s
c                  amdq(i,m) = amdq(i,m) +
c     &               beta(mw)*lamL(mw)*fwave(i,m,mw)
c                  apdq(i,m) = apdq(i,m) +
c     &              (1.d0-beta(mw))*lamR(mw)*fwave(i,m,mw)
               else
                  !amdq(i,m) = amdq(i,m) + .5d0*fwave(i,m,mw)
                  !apdq(i,m) = apdq(i,m) + .5d0*fwave(i,m,mw)
               endif
            enddo
         enddo

      enddo !-- end loop on i

      return
      end subroutine









