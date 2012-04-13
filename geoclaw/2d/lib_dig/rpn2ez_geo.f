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
      logical harten,rare1,rare2,wallprob

      double precision drytol,g,gmod
      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
      double precision pR,pL,hmL,hmR,mL,mR,phi_bedL,phi_bedR
      double precision bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision hstar,hstartest,hstarHLL,sLtest,sRtest,s1m,s2m
      double precision dxdc,dx,pm
      double precision sigbed,SN,rho,kappa,kperm,compress,tau,D,tanpsi
      double precision gamma,eps



      common /cmcapa/  mcapa
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      g=grav
      drytol=drytolerance

      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

!-----------------------Initializing-----------------------------------
         !inform of a bad riemann problem from the start
         if((qr(i-1,1).lt.0.d0).or.(ql(i,1) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(i-1,1),ql(i,1),i
         endif

         !Initialize Riemann problem for grid interface
         do mw=1,mwaves
            s(i,mw)=0.d0
            sw(mw) = 0.d0
            do m=1,meqn
               fwave(i,m,mw)=0.d0
               wave(m,mw) = 0.d0
               fw(m,mw) = 0.d0
            enddo
         enddo

         do m=1,3
            entropy(m)=.false.
         enddo

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

         !zero (small) negative values if they exist
         if (qr(i-1,1).lt.0.d0) then
            do m=1,meqn
               qr(i-1,m)=0.d0
            enddo
         endif

         if (ql(i,1).lt.0.d0) then
            do m=1,meqn
               ql(i,m)=0.d0
            enddo
         endif

         !skip problem if in a completely dry area
         if (qr(i-1,1).le.drytol.and.ql(i,1).le.drytol) then
            go to 30
         endif

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

         !check for wet/dry boundary
         if (hR.gt.drytol) then
            uR=huR/hR
            vR=hvR/hR
            mR=hmR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
         else
            uR = 0.d0
            vR = 0.d0
            mR = 0.d0
            phiR = 0.d0
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            vL=hvL/hL
            mL=hmL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
         else
            uL=0.d0
            vL=0.d0
            mL=0.d0
            phiL = 0.d0
         endif

         wall(1) = 1.d0
         wall(2) = 1.d0
         wall(3) = 1.d0
         wallprob = .false.
         if (hR.le.drytol) then
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
c                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               wallprob = .true.
               hR=hL
               huR=-huL
               bR=bL
               phiR=phiL
               uR=-uL
               vR=vL
            elseif (hL+bL.lt.bR) then
               bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,g)
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
c               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               wallprob = .true.
               hL=hR
               huL=-huR
               bL=bR
               phiL=phiR
               uL=-uR
               vL=vR
            elseif (hR+bR.lt.bL) then
               bL=hR+bR
            endif
         endif

c=================begin digclaw-auxset =================================

         pm = sign(1.d0,uR-uL)

         if (hL.gt.drytol.and.hR.gt.drytol) then
            !this is a completely wet 'normal' Riemann problem
            call auxeval(0.5d0*(hL+hR),0.5d0*(uL+uR),0.5d0*(vL+vR),
     &         0.5d0*(mL+mR),0.5d0*(pL+pR),0.5d0*(phi_bedL+phi_bedR),
     &         kappa,SN,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
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
         eps = gamma + (1.d0-gamma)*kappa
         gmod = grav
c================end digclaw-aux========================================


         !determine wave speeds
         sL=uL-sqrt(gmod*hL*eps) ! 1 wave speed of left state
         sR=uR+sqrt(gmod*hR*eps) ! 2 wave speed of right state
         uhat=(sqrt(hL)*uL + sqrt(hR)*uR)/(sqrt(hR)+sqrt(hL)) ! Roe average
         chat=sqrt(gmod*0.5d0*(hR+hL)*eps) ! Roe average
         sRoe1=uhat-chat ! Roe wave speed 1 wave
         sRoe2=uhat+chat ! Roe wave speed 2 wave

         sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
         sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

         !--------------------end initializing...finally----------
         !solve Riemann problem.

         sw(1) = sRoe1
         sw(3) = sRoe2
         sw(2) = uhat

         if (hR.lt.drytol) then
            sw(3) = uL + 2.d0*sqrt(gmod*hL)
         endif
         if (hL.lt.drytol) then
            sw(1) = uR - 2.d0*sqrt(gmod*hR)
         endif

         call riemann_dig2_conservative(meqn,mwaves,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uhat,uL,uR,vL,vR,mL,mR,
     &         kappa,rho,kperm,compress,tanpsi,D,tau,
     &         gamma,gmod,dx,sw,fw,wave)


c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)
            do m=1,meqn
               fw(m,mw)=fw(m,mw)*wall(mw)
            enddo
         enddo

         do mw=1,mwaves
            s(i,mw)=sw(mw)
            fwave(i,1,mw)=fw(1,mw)
            fwave(i,mhu,mw)=fw(2,mw)
            fwave(i,nhv,mw)=fw(3,mw)
            fwave(i,4,mw) = fw(4,mw)
            fwave(i,5,mw) = fw(5,mw)
         enddo

 30      continue
      enddo


c==========Capacity for mapping from latitude longitude to physical space====

        if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
          if (ixy.eq.1) then
             dxdc=(Rearth*pi/180.d0)
          else
             dxdc=Rearth*pi*cos(auxl(i,3))/180.d0
          endif

          do mw=1,mwaves
c             if (s(i,mw) .gt. 316.d0) then
c               # shouldn't happen unless h > 10 km!
c                write(6,*) 'speed > 316: i,mw,s(i,mw): ',i,mw,s(i,mw)
c                endif
             s(i,mw)=dxdc*s(i,mw)
             do m=1,meqn
               fwave(i,m,mw)=dxdc*fwave(i,m,mw)
             enddo
          enddo
         enddo
        endif

c===============================================================================


c============= compute fluctuations=============================================
         do i=1-mbc,mx+mbc
            do m=1,meqn
               amdq(i,m)=0.0d0
               apdq(i,m)=0.0d0
               do  mw=1,mwaves
                  if (s(i,mw).lt.0.d0) then
                     amdq(i,m)=amdq(i,m) + fwave(i,m,mw)
                  elseif (s(i,mw).gt.0.d0) then
                     apdq(i,m)=apdq(i,m) + fwave(i,m,mw)
                  else
                  amdq(i,m) = amdq(i,m) + .5d0*fwave(i,m,mw)
                  apdq(i,m) = apdq(i,m) + .5d0*fwave(i,m,mw)
                  endif
               enddo
            enddo
         enddo


      return
      end subroutine









