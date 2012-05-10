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
      integer m,i,mw,maxiter,mhu,nhv,mcapa,icom,jcom
      double precision dtcom,dxcom,dycom,tcom
      double precision wall(3),fw(5,3),sw(3),wave(5,3)
      double precision lamL(3),lamR(3),beta(3)
      logical entropy(3)
      logical rare1,rare2,wallprob,drystate
      !logical entropycorr1,entropycorr2

      double precision drytol,g,gmod,veltol
      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL
      double precision pR,pL,hmL,hmR,mL,mR,phi_bedL,phi_bedR
      double precision hstar,hstartest,s1m,s2m,bL,bR
      double precision dxdc,dx,pm
      double precision sigbed,SN,rho,kappa,kperm,compress,tau,D,tanpsi
      double precision rhoL,kpermL,compressL,tauL,DL,tanpsiL
      double precision rhoR,kpermR,compressR,tauR,DR,tanpsiR
      double precision gamma,eps
      double precision h1M,h2M,hu1M,hu2M,u1M,u2M,heR,heL
      double precision phiL,phiR,sE1,sE2



      common /cmcapa/  mcapa
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      g=grav
      drytol=phys_tol

      !set to true to use an entropy correction
      !entropycorr1 = .false.  !Harten (bound nonlinear waves away from s=0)
      !entropycorr2 = .false. !Harten-Hyman (split expansion shocks) not recommended

      !near-zero velocity tolerance
      veltol = 1.d-3

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
            sw(mw) = 0.d0
            entropy(mw)=.false.
            do m=1,meqn
               fwave(i,m,mw)=0.d0
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
         else
            mhu=3
            nhv=2
            dx = dycom
         endif

         !zero (small) negative values if they exist and set velocities
         call admissibleq(ql(i,1),ql(i,mhu),ql(i,nhv),
     &            ql(i,4),ql(i,5),uR,vR,mR)

         call admissibleq(qr(i-1,1),qr(i-1,mhu),qr(i-1,nhv),
     &            qr(i-1,4),qr(i-1,5),uL,vL,mL)


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
         bL = auxr(i-1,1)
         bR = auxl(i,1)
         phi_bedL = auxr(i-1,i_phi)
         phi_bedR = auxl(i,i_phi)
         phiR = 0.5d0*g*hR**2 + hR*uR**2
         phiL = 0.5d0*g*hL**2 + hL*uL**2

         !test for wall problem vs. inundation problem
         do mw=1,mwaves
            wall(mw) = 1.d0
         enddo
         drystate=.false.
         wallprob = .false.
         if (hR.le.drytol) then
            drystate=.true.
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,
     &                                 rare1,rare2,1,drytol,g)
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
            !elseif (hL+bL.lt.bR) then
               !bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            drystate=.true.
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
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
            !elseif (hR+bR.lt.bL) then
               !bL=hR+bR
            endif
         endif

c=================begin digclaw-auxset =================================

         pm = dsign(1.d0,uR-uL)

         if (hL.gt.drytol.and.hR.gt.drytol) then
            !this is a completely wet 'normal' Riemann problem
            call auxeval(hR,uR,vR,mR,pR,phi_bedR,
     &        kappa,SN,rhoR,tanpsi,DR,tauR,sigbed,kpermR,compressR,pm)
            call auxeval(hL,uL,vL,mL,pL,phi_bedL,
     &        kappa,SN,rhoL,tanpsi,DL,tauL,sigbed,kpermL,compressL,pm)
            D = 0.5d0*(DL + DR)
            compress = 0.5d0*(compressL + compressR)
            kperm = 0.5d0*(kpermR + kpermL)
            tau = 0.5d0*(tauL + tauR)
            rho = 0.5d0*(rhoR + rhoL)
         elseif (hR.gt.drytol) then
            !inundation problem
            call auxeval(hR,uR,vR,mR,pR,phi_bedR,
     &        kappa,SN,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
         elseif (hL.gt.drytol) then
            !inundation problem
            call auxeval(hL,uL,vL,mL,pL,phi_bedL,
     &        kappa,SN,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
         endif

         gamma = 1.5d0*(rho_f/(6.d0*rho)+0.5d0)
         eps = kappa + (1.d0-kappa)*gamma
         gmod = grav*eps
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
         call riemann_dig2_aug_sswave(meqn,mwaves,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         kappa,rho,kperm,compress,tanpsi,D,tau,
     &         gamma,gmod,dx,sw,fw,wave)


c         call riemann_aug_JCP(1,3,3,hL,hR,huL,
c     &        huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,
c     &                                    drytolerance,gmod,sw,fw)



         !--------------------------------------------------------------

c        !eliminate ghost fluxes for wall
         do mw=1,mwaves
            sw(mw)=sw(mw)*wall(mw)
            do m=1,meqn
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
         do mw=1,mwaves
            s(i,mw)=sw(mw)
            fwave(i,1,mw)=fw(1,mw)
            fwave(i,mhu,mw)=fw(2,mw)
            fwave(i,nhv,mw)=fw(3,mw)
            fwave(i,4,mw) = fw(4,mw)
            fwave(i,5,mw) = fw(5,mw)
         enddo


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
                  amdq(i,m) = amdq(i,m) +
     &               beta(mw)*lamL(mw)*fwave(i,m,mw)
                  apdq(i,m) = apdq(i,m) +
     &              (1.d0-beta(mw))*lamR(mw)*fwave(i,m,mw)
               else
c                  amdq(i,m) = amdq(i,m) + .5d0*fwave(i,m,mw)
c                  apdq(i,m) = apdq(i,m) + .5d0*fwave(i,m,mw)
               endif
            enddo
         enddo

      enddo !-- end loop on i

      return
      end subroutine









