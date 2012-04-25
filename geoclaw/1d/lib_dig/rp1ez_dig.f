c =====================================================================
      subroutine rp1(maxmx,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &             fwave,s,amdq,apdq)
c =====================================================================
c
c
c
c     # On input,
c     # ql contains the state vector at the left edge of each cell
c     # qr contains the state vector at the right edge of each cell
c     #
c     # On output, wave contains the fwaves/s,
c     #            s the speeds,
c     #            amdq the  left-going flux difference  A**- \Delta q
c     #            apdq the right-going flux difference  A**+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp1 is called with ql=qr=q.
c
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for Iverson's debris flow equations         !
!                                                                           !
!                                                                           !
!           David George, Vancouver WA, Dec. 2009                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !input
      integer maxmx,meqn,mwaves,mbc,mx

      double precision ql(1-mbc:maxmx+mbc, meqn)
      double precision qr(1-mbc:maxmx+mbc, meqn)
      double precision s(1-mbc:maxmx+mbc, mwaves)
      double precision fwave(1-mbc:maxmx+mbc, meqn, mwaves)
      double precision amdq(1-mbc:maxmx+mbc, meqn)
      double precision apdq(1-mbc:maxmx+mbc, meqn)
      double precision auxl(1-mbc:maxmx+mbc, *)
      double precision auxr(1-mbc:maxmx+mbc, *)

      include "digparamsdec.i"
      include "digparamscommon.i"

      !local only
      integer i,j,m,mw
      double precision fw(4,3),sw(3),wave(4,3),iflag(4)
      double precision lamL(3),lamR(3),beta(3)
      logical entropy(3)
      logical harten

      double precision drytol,g,stol,em,u
      double precision hL,hR,huL,huR,hmL,hmR,pL,pR,gmod
      double precision theta,thetaL,thetaR,kL,kR,rhoL,rhoR,mR,mL,uR,uL
      double precision sL,sR,uhat,chat,sRoe1,sRoe2,sE1,sE2
      double precision s1M,s1L,s2M,s2R,h1M,hu1M,u1M,h2M,hu2M,u2M
      double precision speedtol,sRoe1mod,sRoe2mod
      double precision gamma,eps,rho,kappa,kperm,compress,tau,D,tanpsi

      g=grav
      drytol=2.d0*dry_tol

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
                 fw(m,mw) = 0.d0
                 wave(m,mw) = 0.d0
                 fwave(i,m,mw)=0.d0
              enddo
         enddo
         do m=1,3
            entropy(m) = .false.
         enddo

         !skip problem if in a completely dry area
         if (qr(i-1,1).le.drytol.and.ql(i,1).le.drytol) then
            go to 30
         endif

         !Riemann problem variables
         hL = qr(i-1,1)
         hR = ql(i,1)
         huL = qr(i-1,2)
         huR = ql(i,2)
         hmL = qr(i-1,3)
         hmR = ql(i,3)
         pL = qr(i-1,4)
         pR = ql(i,4)

         if (hR.gt.drytol) then
            uR=huR/hR
            mR = hmR/hR
         else
            uR = 0.d0
            mR = m0
            hR = 0.d0
            huR = 0.d0
            hmR = 0.d0
            pR = 0.d0
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            mL = hmL/hL
         else
            uL = 0.d0
            mL = m0
            hL=0.d0
            huL=0.d0
            hmL = 0.d0
            pL = 0.d0
         endif

         if (hL.gt.drytol.and.hR.gt.drytol) then
            kperm = 0.5d0*(auxr(i-1,i_kperm) + auxl(i,i_kperm))
            compress = 0.5d0*(auxr(i-1,i_alpha) + auxl(i,i_alpha))
            D = 0.5d0*(auxr(i-1,i_D) + auxl(i,i_D))
            tau = 0.5d0*(auxr(i-1,i_tau) + auxl(i,i_tau))
            tanpsi = 0.5d0*(auxr(i-1,i_tanpsi) + auxl(i,i_tanpsi))
            rho = 0.5d0*(auxr(i-1,i_rho) + auxl(i,i_rho))
         elseif (hR.gt.drytol) then
            kperm = auxl(i,i_kperm)
            compress = auxl(i,i_alpha)
            D = auxl(i,i_D)
            tau =auxl(i,i_tau)
            tanpsi = auxl(i,i_tanpsi)
            rho = auxl(i,i_rho)
         elseif (hL.gt.drytol) then
            kperm = auxr(i-1,i_kperm)
            compress = auxr(i-1,i_alpha)
            D = auxr(i-1,i_D)
            tau = auxr(i-1,i_tau)
            tanpsi = auxr(i-1,i_tanpsi)
            rho = auxr(i-1,i_rho)
         endif
         thetaL =auxr(i-1,i_theta)
         thetaR =auxl(i,i_theta)
         theta = 0.5d0*(auxr(i-1,i_theta) + auxl(i,i_theta))
         kL=auxr(i-1,i_kappa)
         kR=auxl(i,i_kappa)
         kappa = kR
         rhoL=auxr(i-1,i_rho)
         rhoR=auxl(i,i_rho)
         gamma = 1.5d0*(rho_f/(6.d0*rho)+0.5d0)
         eps = gamma + (1.d0-gamma)*kappa
         gmod = grav*cos(theta)

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

         call riemann_dig1_conservative(i,meqn,mwaves,hL,hR,huL,huR,
     &        hmL,hmR,pL,pR,uhat,uL,uR,mL,mR,kL,kR,rho_f,theta,
     &        thetaL,thetaR,rhoL,rhoR,gamma,kperm,compress,tanpsi,D,tau,
     &        g,gmod,drytol,sw,fw,wave)


         do mw=1,mwaves
            s(i,mw) = sw(mw)
            do m=1,meqn
               fwave(i,m,mw) = fw(m,mw)
               if (i.ge.47.and.i.le.184.and..false.) then
                  write(*,*) 'i,m,mw,fw',i,m,mw,fwave(i,m,mw)
               endif
            enddo
         enddo

         !Entropy correction -------------------------------------------

         s1L = sL
         s2R = sR

         h1M =  max(hL  + wave(1,1),0.d0)
         h2M =  max(hR  - wave(1,3),0.d0)
         hu1M = huL + wave(2,1)
         hu2M = huR - wave(2,3)
         if (h1M.gt.drytol) then
            u1M = hu1M/h1M
         else
            u1M = 0.d0
         endif
         if (h2M.gt.drytol) then
            u2M = hu2M/h2M
         else
            u2M = 0.d0
         endif
         s1M = u1M - sqrt(gmod*h1M)
         s2M = u2M + sqrt(gmod*h2M)

         if (s1L.le.0.d0.and.s1M.gt.0.d0) then
            entropy(1) = .true.
            beta(1) = (s1M - sw(1))/(s1M - s1L)
            lamL(1) = s1L
            lamR(1) = s1M
         endif
         if (s2R.ge.0.d0.and.s2M.lt.0.d0) then
            entropy(3) = .true.
            beta(3) = (s2R - sw(3))/(s2R - s2M)
            lamL(3) = s2M
            lamR(3) = s2R
         endif
         !do m=1,mwaves
         !   entropy(m)=.false.
         !enddo
c-----------------------------------------------------------------------

 30      continue

         do m=1,meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do  mw=1,mwaves
               if (s(i,mw).lt.0.d0.and.(.not.entropy(mw))) then
                  amdq(i,m) = amdq(i,m) + fwave(i,m,mw)
               elseif (s(i,mw).gt.0.d0.and.(.not.entropy(mw))) then
                  apdq(i,m) = apdq(i,m) + fwave(i,m,mw)
               elseif (s(i,mw).eq.0.d0) then
                  apdq(i,m) = apdq(i,m) + 0.5d0*fwave(i,m,mw)
                  amdq(i,m) = amdq(i,m) + 0.5d0*fwave(i,m,mw)
               elseif (entropy(mw).and.abs(s(i,mw)).gt.0.d0) then
                  amdq(i,m) = amdq(i,m) +
     &               beta(mw)*lamL(mw)*fwave(i,m,mw)/s(i,mw)
                  apdq(i,m) = apdq(i,m) +
     &              (1.d0-beta(mw))*lamR(mw)*fwave(i,m,mw)/s(i,mw)
               endif
            enddo
         enddo

      enddo




      return
      end subroutine
