

c-----------------------------------------------------------------------
      subroutine riemann_dig2_aug_sswave(meqn,mwaves,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         kappa,rho,kperm,compress,tanpsi,D,tau,
     &         gamma,gmod,dx,sw,fw,w)

      ! solve the dig Riemann problem for debris flow eqn
      ! this is for 2d version

      use geoclaw_module
      use digclaw_module

      implicit none

*     !i/o
      integer meqn,mwaves

      double precision hL,hR,huL,huR,hvL,hvR,hmL,hmR,pL,pR
      double precision bL,bR,uL,uR,vL,vR,mL,mR
      double precision kappa,rho,kperm,compress,tanpsi,D,tau
      double precision gamma,gmod,dx

      double precision fw(meqn,mwaves),w(meqn,mwaves)
      double precision sw(mwaves)
      double precision psi(4)

*     !local
      integer m,mw,ks,mp
      double precision h,u,v,mbar
      double precision det1,det2,det3,detR
      double precision R(0:4,0:4),A(5,5),del(0:4)
      double precision beta(0:4),betas(4)
      double precision sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision delb,s1m,s2m,hm,heL,heR,criticaltol
      double precision s1s2bar,s1s2tilde,hbar,source2dx,veltol1,veltol2
      double precision hstarHLL,deldelh,drytol
      logical sonic,rare1,rare2
      logical entropycorr1

      veltol1=1.d-6
      veltol2=0.d0
      criticaltol=1.d-6
      drytol=phys_tol

      entropycorr1 = .false.

      if (hL.ge.drytol.and.hR.ge.drytol) then
         h = 0.5d0*(hL + hR)
         u = uhat
         v = 0.5d0*(vL + vR)
         mbar = 0.5d0*(mL + mR)
      elseif (hL.ge.drytol) then
         h = hL
         u = uhat
         v = vL
         mbar = mL
      else
         h = hR
         u = uhat
         v = vR
         mbar = mR
      endif

      !determine wave speeds
      sL=uL-dsqrt(gmod*hL) ! 1 wave speed of left state
      sR=uR+dsqrt(gmod*hR) ! 2 wave speed of right state
      uhat=(dsqrt(hL)*uL + dsqrt(hR)*uR)/(dsqrt(hR)+dsqrt(hL)) ! Roe average
      chat=dsqrt(gmod*0.5d0*(hR+hL)) ! Roe average
      sRoe1=uhat-chat ! Roe wave speed 1 wave
      sRoe2=uhat+chat ! Roe wave speed 2 wave

      sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
      sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave
      sw(1) = sE1
      sw(3) = sE2
      !sw(2) = 0.5d0*(sw(1)+sw(3))!uhat


      call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,
     &                                          1,drytol,gmod)
      sw(1)= min(sw(1),s2m) !Modified Einfeldt speed
      sw(3)= max(sw(3),s1m) !Modified Einfeldt speed
      !sw(2) = 0.5d0*(sw(3)+sw(1))
      sw(2) = uhat

      !----Harten entropy fix for near zero-speed nonlinear waves
      ! Note: this might change near zero-speed shocks as well
      ! smoothed absolute value function between veltol and 2*veltol
      if (entropycorr1) then
         do mw=1,mwaves
            if (dabs(sw(mw)).le.2.d0*veltol1) then
              sw(mw) = dsign((sw(mw)**2 + (2.d0*veltol1)**2)/
     &            (4.d0*veltol1),sw(mw))
            endif
         enddo
      endif

      delb=bR-bL

      !determine ss-wave
      hbar =  0.5d0*(hL+hR)
      s1s2bar = 0.25d0*(uL+uR)**2- gmod*hbar
      s1s2tilde= max(0.d0,uL*uR) - gmod*hbar

c     !find if sonic problem
      sonic=.false.
      if (dabs(s1s2bar).le.criticaltol) sonic=.true.
      if (s1s2bar*s1s2tilde.le.criticaltol) sonic=.true.
      if (s1s2bar*sE1*sE2.le.criticaltol) sonic = .true.
      if (min(dabs(sE1),dabs(sE2)).lt.criticaltol) sonic=.true.
      if (sE1.lt.0.d0.and.s1m.gt.0.d0) sonic = .true.
      if (sE2.gt.0.d0.and.s2m.lt.0.d0) sonic = .true.
      if ((uL+dsqrt(gmod*hL))*(uR+dsqrt(gmod*hR)).lt.0.d0) sonic=.true.
      if ((uL-dsqrt(gmod*hL))*(uR-dsqrt(gmod*hR)).lt.0.d0) sonic=.true.

      if (sonic) then
         source2dx = -grav*hbar*delb
      else
         source2dx = -delb*grav*hbar*s1s2tilde/s1s2bar
      endif

      source2dx=min(source2dx,gmod*max(-hL*delb,-hR*delb))
      source2dx=max(source2dx,gmod*min(-hL*delb,-hR*delb))
      if (dabs(u).le.veltol2) then
         source2dx=-h*grav*delb
      endif

c     !find bounds in case of critical state resonance, or negative states
c     !find jump in h, deldelh
      hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve
      if (sonic) then
         deldelh =  -delb
      else
         deldelh = delb*gmod*hbar/s1s2bar
      endif
c     !find bounds in case of critical state resonance, or negative states
      if (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) then
         deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
         deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
      elseif (sE1.ge.criticaltol) then
         deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
         deldelh = max(deldelh,-hL)
      elseif (sE2.le.-criticaltol) then
         deldelh = min(deldelh,hR)
         deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
      endif


*     !determine R
      R(0,0) = 0.d0
      R(1,0) = 0.d0
      R(2,0) = 1.d0
      R(3,0) = 0.d0
      R(4,0) = 0.d0

      R(0,1) = 1.d0
      R(1,1) = sw(1)
      R(2,1) = sw(1)**2
      R(3,1) = gamma*rho*grav
      R(4,1) = gamma*rho*grav*sw(1)

      R(0,2) = kappa-1.d0
      R(1,2) = sw(2)*(kappa-1.d0)
      R(2,2) = (kappa-1.d0)*sw(2)**2
      R(3,2) = kappa*rho*grav
      R(4,2) = kappa*rho*grav*sw(2)

      R(0,3) = 1.d0
      R(1,3) = sw(3)
      R(2,3) = sw(3)**2
      R(3,3) = gamma*rho*grav
      R(4,3) = gamma*rho*grav*sw(3)

      R(0,4) = 0.d0
      R(1,4) = 0.d0
      R(2,4) = 0.d0
      R(3,4) = 0.d0
      R(4,4) = 1.d0

      !determine del
      del(0) = hR- hL - deldelh
      del(1) = huR - huL
      del(2) = hR*uR**2 + 0.5d0*kappa*grav*hR**2 -
     &      (hL*uL**2 + 0.5d0*kappa*grav*hL**2)
      del(2) = del(2) + (1.d0-kappa)*h*(pR-pL)/rho
      del(3) = pR - pL - gamma*rho*grav*deldelh
      del(4) = -gamma*rho*grav*u*(hR-hL) + gamma*rho*grav*del(1)
     &         + u*(pR-pL)

*     !determine the source term
      call psieval(tau,rho,D,tanpsi,kperm,compress,h,u,mbar,psi)

      if (dabs(u).gt.veltol2) then
         del(2) = del(2) -source2dx - dx*psi(2)
      else
         if (dabs(del(2)-source2dx).ge.dabs(dx*tau/rho)) then
            del(2)=dsign(dabs(dabs(del(2)-source2dx)
     &                   -dabs(dx*tau/rho)),del(2)-source2dx)
         else
            del(0)=0.d0
            del(2)=0.d0
            del(3)=0.d0
         endif
      endif
      del(1) = del(1) - 0.5d0*dx*psi(1)
      del(4) = del(4) - 0.5d0*dx*psi(4)

*     !R beta = del
*     !gauss routine replaces del with beta and R with it's inverse
      !want to keep R, so replacing with A
      do mw=0,4
         beta(mw) = del(mw)
         do m=0,5
            A(m+1,mw+1)=R(m,mw)
         enddo
      enddo
      call gaussj(A,5,5,beta,1,5,1)

      do mw=1,3
         do m=1,2
            fw(m,mw) = beta(mw)*R(m,mw)
         enddo
            !fw(5,mw)  = sw(mw)*beta(mw)*R(3,mw)
            fw(5,mw)  = beta(mw)*R(4,mw)
      enddo

      !corrector wave
      fw(2,2) = fw(2,2) + beta(0)
      fw(5,2) = fw(5,2) + beta(4)

      !waves and fwaves for delta hum
      fw(4,1) = fw(1,1)*mL
      fw(4,3) = fw(1,3)*mR
      fw(4,2) = hmR*uR-hmL*uL - fw(4,1)- fw(4,3)-dx*0.5d0*psi(3)

      !waves and fwaves for delta huv

      fw(3,1) = fw(1,1)*vL
      fw(3,3) = fw(1,3)*vR
      fw(3,2) = hvR*uR-hvL*uL -fw(3,1) -fw(3,3) !- dx*D*v*(rho-rho_f)/rho

      return
      end !subroutine riemann_dig2_aug_sswave

c-----------------------------------------------------------------------
      subroutine riemann_dig2_sswave(meqn,mwaves,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         kappa,rho,kperm,compress,tanpsi,D,tau,
     &         gamma,gmod,dx,veltol,sw,fw,w)

      ! solve the dig Riemann problem for debris flow eqn
      ! this is for 2d version

      use geoclaw_module
      use digclaw_module

      implicit none

*     !i/o
      integer meqn,mwaves

      double precision hL,hR,huL,huR,hvL,hvR,hmL,hmR,pL,pR
      double precision bL,bR,uL,uR,vL,vR,mL,mR
      double precision kappa,rho,kperm,compress,tanpsi,D,tau
      double precision gamma,gmod,dx,veltol

      double precision fw(meqn,mwaves),w(meqn,mwaves)
      double precision sw(mwaves)
      double precision psi(4)

*     !local
      integer m,mw,ks,mp
      double precision h,u,mbar
      double precision det1,det2,det3,detR
      double precision R(3,3),A(3,3),alph(3),delq(3),delf(3),del(3)
      double precision beta(3)
      double precision sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision delb,s1m,s2m,hm,heL,heR,criticaltol
      double precision hLstar,hRstar,uLstar,uRstar,huLstar,huRstar
      double precision s1s2bar,s1s2tilde,hbar,source2dx,veltol2,seL,seR
      double precision rare1st,rare2st,raremin,raremax,hstarHLL,sdelta
      logical rare1,rare2,rarecorrector,rarecorrectortest,sonic

      veltol2=0.d0
      criticaltol=1.d-6

      !determine wave speeds
      sL=uL-dsqrt(gmod*hL) ! 1 wave speed of left state
      sR=uR+dsqrt(gmod*hR) ! 2 wave speed of right state
      uhat=(dsqrt(hL)*uL + dsqrt(hR)*uR)/(dsqrt(hR)+dsqrt(hL)) ! Roe average
      chat=dsqrt(gmod*0.5d0*(hR+hL)) ! Roe average
      sRoe1=uhat-chat ! Roe wave speed 1 wave
      sRoe2=uhat+chat ! Roe wave speed 2 wave

      if (hL.ge.drytolerance.and.hR.ge.drytolerance) then
         h = 0.5d0*(hL + hR)
         u = uhat
         mbar = 0.5d0*(mL + mR)
      elseif (hL.ge.drytolerance) then
         h = hL
         u = uhat
         mbar = mL
      else
         h = hR
         u = uhat
         mbar = mR
      endif

      sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
      sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

      sw(1) = sE1
      sw(3) = sE2
      !sw(2) = 0.5d0*(sw(1)+sw(3))!uhat
      sw(2) = uhat

      !determine Riemann speeds
      heL = max(hL + bL - max(bL,bR),0.d0)
      heR = max(hR + bR - max(bL,bR),0.d0)
      seL = uL-dsqrt(gmod*heL)
      seR = uR+dsqrt(gmod*heR)
      call riemanntype(heL,heR,uL,uR,hm,s1m,s2m,rare1,rare2,
     &                                          1,drytolerance,gmod)
      sw(1)= min(sw(1),s2m) !Modified Einfeldt speed
      sw(3)= max(sw(3),s1m) !Modified Einfeldt speed

      sw(1) = min(sw(1),seL)
      sw(3) = max(sw(3),seR)
      sw(2) = 0.5d0*(sw(3)+sw(1))


      !determine delq
      delq(1) = hR-hL
      delq(2) = huR - huL
      delq(3) = pR - pL

      delf(1) = delq(2)
      delf(2) = hR*uR**2 + 0.5d0*kappa*grav*hR**2 -
     &      (hL*uL**2 + 0.5d0*kappa*grav*hL**2)
      delf(3) = 0.d0

      del(1) = delf(1)
      del(2) = delf(2) + (1.d0-kappa)*delq(3)/rho
      del(3) = -rho*u*grav*gamma*delq(1) + rho*grav*gamma*delq(2)
     &         + u*delq(3)

      delb=bR-bL

      !determine ss-wave
      hLstar=hL
      hRstar=hR
      uLstar=uL
      uRstar=uR
      huLstar=uLstar*hLstar
      huRstar=uRstar*hRstar
      hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
      s1s2bar = 0.25d0*(uLstar+uRstar)**2 - gmod*hbar
      s1s2tilde= max(0.d0,uLstar*uRstar) - gmod*hbar

c     !find if sonic problem
      sE1=min(sw(1),uL-dsqrt(gmod*hL))
      sE2=max(sw(3),uR+dsqrt(gmod*hr))
      sonic=.false.
      if (dabs(s1s2bar).le.criticaltol) sonic=.true.
      if (s1s2bar*s1s2tilde.le.criticaltol) sonic=.true.
      if (s1s2bar*sE1*sE2.le.criticaltol) sonic = .true.
      if (min(dabs(sE1),dabs(sE2)).lt.criticaltol) sonic=.true.
      if (sE1.lt.0.d0.and.s1m.gt.0.d0) sonic = .true.
      if (sE2.gt.0.d0.and.s2m.lt.0.d0) sonic = .true.
      if ((uL+dsqrt(gmod*hL))*(uR+dsqrt(gmod*hR)).lt.0.d0) sonic=.true.
      if ((uL-dsqrt(gmod*hL))*(uR-dsqrt(gmod*hR)).lt.0.d0) sonic=.true.

      if (sonic) then
         source2dx = -grav*hbar*delb
      else
         source2dx = -delb*grav*hbar*s1s2tilde/s1s2bar
      endif
c     !find bounds in case of critical state resonance, or negative states
      source2dx=min(source2dx,gmod*max(-hLstar*delb,-hRstar*delb))
      source2dx=max(source2dx,gmod*min(-hLstar*delb,-hRstar*delb))
      if (dabs(u).le.veltol2) then
         source2dx=-h*grav*delb
      endif

c     !determine the middle entropy corrector wave------------------------
      rarecorrectortest=.true.
      rarecorrector=.false.
      if (rarecorrectortest) then
         hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve
         sdelta=sw(3)-sw(1)
         raremin = 0.5d0
         raremax = 0.9d0
         if (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2d0
         if (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2d0
         if (rare1.or.rare2) then
            !see which rarefaction is larger
            rare1st=3.d0*(dsqrt(gmod*hL)-dsqrt(gmod*hm))
            rare2st=3.d0*(dsqrt(gmod*hR)-dsqrt(gmod*hm))
            if (max(rare1st,rare2st).gt.raremin*sdelta.and.
     &         max(rare1st,rare2st).lt.raremax*sdelta) then
                  rarecorrector=.true.
               if (rare1st.gt.rare2st) then
                  sw(2)=s1m
               elseif (rare2st.gt.rare1st) then
                  sw(2)=s2m
               else
                  sw(2)=0.5d0*(s1m+s2m)
               endif
            endif
         endif
         if (hstarHLL.lt.min(hL,hR)/5.d0) rarecorrector=.false.
      endif

*     !determine R
      R(1,1) = 1.d0
      R(2,1) = sw(1)
      R(3,1) = gamma*rho*grav

      R(1,2) = kappa-1.d0
      R(2,2) = uhat*(kappa-1.d0)
      R(3,2) = kappa*rho*grav

      R(1,3) = 1.d0
      R(2,3) = sw(3)
      R(3,3) = gamma*rho*grav

*     !determine the source term
      do m=1,4
            psi(m)=0.d0
      enddo
      call psieval(tau,rho,D,tanpsi,kperm,compress,h,u,mbar,psi)

      del(1) = del(1) - dx*psi(1)
      del(3) = del(3) - dx*psi(4)
      if (dabs(u).gt.veltol2) then
         del(2) = del(2) -source2dx - dx*psi(2)
      else
         if (dabs(del(2)-source2dx).ge.dabs(dx*tau/rho)) then
            del(2)=dsign(dabs(dabs(delf(2)-source2dx)
     &                   -dabs(dx*tau/rho)),del(2)-source2dx)
         else
            del(2)=0.d0
         endif
      endif

*     R alpha = delq
*     solve for alpha1 using Cramer's rule

      det1=R(1,1)*(R(2,2)*R(3,3)-R(2,3)*R(3,2))
      det2=R(1,2)*(R(2,1)*R(3,3)-R(2,3)*R(3,1))
      det3=R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1))
      detR = (det1-det2+det3)


      do ks = 1,3 !solve for alpha(k)
         !build the A matrix for Cramer's rule
         do mw = 1,3
            do m = 1,3
               A(m,mw)  = R(m,mw)
               A(m,ks)  = delq(m)
            enddo
         enddo
         det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
         det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
         det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

         alph(ks) = (det1-det2+det3)/detR
      enddo

      do ks = 1,3 !solve for beta(k)
         !build the A matrix for Cramer's rule
         do mw = 1,3
            do m = 1,3
               A(m,mw)  = R(m,mw)
               A(m,ks)  = del(m)
            enddo
         enddo
         det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
         det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
         det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

         beta(ks) = (det1-det2+det3)/detR
      enddo

      do m=1,3
         mp=m
         if (m.eq.3) mp = 5
         do mw=1,3
            w(mp,mw)  = alph(mw)*R(m,mw)
            fw(mp,mw) = beta(mw)*R(m,mw)
         enddo
      enddo

      !waves and fwaves for delta hum
      w(4,1) = alph(1)*mL
      w(4,3) = alph(3)*mR
      w(4,2) = hmR-hmL - mL*alph(1) - mR*alph(3)

      fw(4,1) = beta(1)*mL
      fw(4,3) = beta(3)*mR
      fw(4,2) = hmR*uR-hmL*uL - mL*beta(1)- mR*beta(3)-dx*psi(3)

      !waves and fwaves for delta huv
      w(3,1) = alph(1)*vL
      w(3,3) = alph(3)*vR
      w(3,2) = hvR-hvL - vL*alph(1) - vR*alph(3)

      fw(3,1) = beta(1)*vL
      fw(3,3) = beta(3)*vR
      fw(3,2) = hvR*uR - hvL*uL -beta(1)*vL - beta(3)*vR

      return
      end !subroutine riemann_dig2_sswave


c-----------------------------------------------------------------------
      subroutine riemann_dig2_conservative(meqn,mwaves,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         kappa,rho,kperm,compress,tanpsi,D,tau,
     &         gamma,gmod,dx,veltol,sw,fw,w)

      ! solve the dig Riemann problem for debris flow eqn
      ! this is for 2d version

      use geoclaw_module
      use digclaw_module

      implicit none

*     !i/o
      integer meqn,mwaves

      double precision hL,hR,huL,huR,hvL,hvR,hmL,hmR,pL,pR
      double precision bL,bR,uhat,uL,uR,vL,vR,mL,mR
      double precision kappa,rho,kperm,compress,tanpsi,D,tau
      double precision gamma,gmod,dx,veltol

      double precision fw(meqn,mwaves),w(meqn,mwaves)
      double precision sw(mwaves)
      double precision psi(4)

*     !local
      integer m,mw,ks,mp
      double precision h,u,mbar
      double precision det1,det2,det3,detR
      double precision R(3,3),A(3,3),alph(3),delq(3),delf(3),del(3)
      double precision beta(3)

      if (hL.ge.drytolerance.and.hR.ge.drytolerance) then
         h = 0.5d0*(hL + hR)
         u = uhat
         mbar = 0.5d0*(mL + mR)
      elseif (hL.ge.drytolerance) then
         h = hL
         u = uhat
         mbar = mL
      else
         h = hR
         u = uhat
         mbar = mR
      endif

*     !determine R
      R(1,1) = 1.d0
      R(2,1) = sw(1)
      R(3,1) = gamma*rho*gmod

      R(1,2) = kappa-1.d0
      R(2,2) = uhat*(kappa-1.d0)
      R(3,2) = kappa*rho*gmod

      R(1,3) = 1.d0
      R(2,3) = sw(3)
      R(3,3) = gamma*rho*gmod

      !determine delq
      delq(1) = hR-hL
      delq(2) = huR - huL
      delq(3) = pR - pL

      delf(1) = delq(2)
      delf(2) = hR*uR**2 + 0.5d0*kappa*gmod*hR**2 -
     &      (hL*uL**2 + 0.5d0*kappa*gmod*hL**2)
      delf(3) = 0.d0

      del(1) = delf(1)
      del(2) = delf(2) + (1.d0-kappa)*delq(3)/rho
      del(3) = -rho*u*gmod*gamma*delq(1) + rho*gmod*gamma*delq(2)
     &         + u*delq(3)


*     !determine the source term
      do m=1,4
            psi(m)=0.d0
      enddo
      call psieval(tau,rho,D,tanpsi,kperm,compress,h,u,mbar,psi)

      del(1) = del(1) - dx*psi(1)
      del(3) = del(3) - dx*psi(4)
      if (dabs(u).gt.0.d0*veltol) then
         del(2) = del(2) + gmod*h*(bR-bL) - dx*psi(2)
      else
         if (dabs(del(2)+gmod*h*(bR-bL)).ge.dabs(dx*tau/rho)) then
            del(2)=dsign(dabs(dabs(delf(2)+gmod*h*(bR-bL))
     &                   -dabs(dx*tau/rho)),del(2)+gmod*h*(bR-bL))
         else
            del(2)=0.d0
         endif
      endif

*     R alpha = delq
*     solve for alpha1 using Cramer's rule

      det1=R(1,1)*(R(2,2)*R(3,3)-R(2,3)*R(3,2))
      det2=R(1,2)*(R(2,1)*R(3,3)-R(2,3)*R(3,1))
      det3=R(1,3)*(R(2,1)*R(3,2)-R(2,2)*R(3,1))
      detR = (det1-det2+det3)


      do ks = 1,3 !solve for alpha(k)
         !build the A matrix for Cramer's rule
         do mw = 1,3
            do m = 1,3
               A(m,mw)  = R(m,mw)
               A(m,ks)  = delq(m)
            enddo
         enddo
         det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
         det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
         det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

         alph(ks) = (det1-det2+det3)/detR
      enddo

      do ks = 1,3 !solve for beta(k)
         !build the A matrix for Cramer's rule
         do mw = 1,3
            do m = 1,3
               A(m,mw)  = R(m,mw)
               A(m,ks)  = del(m)
            enddo
         enddo
         det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
         det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
         det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

         beta(ks) = (det1-det2+det3)/detR
      enddo

      do m=1,3
         mp=m
         if (m.eq.3) mp = 5
         do mw=1,3
            w(mp,mw)  = alph(mw)*R(m,mw)
            fw(mp,mw) = beta(mw)*R(m,mw)
         enddo
      enddo

      !waves and fwaves for delta hum
      w(4,1) = alph(1)*mL
      w(4,3) = alph(3)*mR
      w(4,2) = hmR-hmL - mL*alph(1) - mR*alph(3)

      fw(4,1) = beta(1)*mL
      fw(4,3) = beta(3)*mR
      fw(4,2) = hmR*uR-hmL*uL - mL*beta(1)- mR*beta(3)-dx*psi(3)

      !waves and fwaves for delta huv
      w(3,1) = alph(1)*vL
      w(3,3) = alph(3)*vR
      w(3,2) = hvR-hvL - vL*alph(1) - vR*alph(3)

      fw(3,1) = beta(1)*vL
      fw(3,3) = beta(3)*vR
      fw(3,2) = hvR*uR - hvL*uL -beta(1)*vL - beta(3)*vR

      return
      end !subroutine riemann_dig2_conservative



c=============================================================================
      subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,
     &             maxiter,drytol,g)

      !determine the Riemann structure (wave-type in each family)


      implicit none

      !input
      double precision hL,hR,uL,uR,drytol,g
      integer maxiter

      !output
      double precision s1m,s2m
      logical rare1,rare2

      !local
      double precision hm,u1m,u2m,um,delu
      double precision h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
      integer iter



c     !Test for Riemann structure

      h_min=min(hR,hL)
      h_max=max(hR,hL)
      delu=uR-uL

      if (h_min.le.drytol) then
         hm=0.d0
         um=0.d0
         s1m=uR+uL-2.d0*dsqrt(g*hR)+2.d0*dsqrt(g*hL)
         s2m=uR+uL-2.d0*dsqrt(g*hR)+2.d0*dsqrt(g*hL)
         if (hL.le.0.d0) then
            rare2=.true.
            rare1=.false.
         else
            rare1=.true.
            rare2=.false.
         endif

      else
         F_min= delu+2.d0*(dsqrt(g*h_min)-dsqrt(g*h_max))
         F_max= delu +
     &         (h_max-h_min)*(dsqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

         if (F_min.gt.0.d0) then !2-rarefactions

            hm=(1.d0/(16.d0*g))*
     &               max(0.d0,-delu+2.d0*(dsqrt(g*hL)+dsqrt(g*hR)))**2
            um=dsign(1.d0,hm)*(uL+2.d0*(dsqrt(g*hL)-dsqrt(g*hm)))

            s1m=uL+2.d0*dsqrt(g*hL)-3.d0*dsqrt(g*hm)
            s2m=uR-2.d0*dsqrt(g*hR)+3.d0*dsqrt(g*hm)

            rare1=.true.
            rare2=.true.

         elseif (F_max.le.0.d0) then !2 shocks

c           !root finding using a Newton iteration on dsqrt(h)===
            h0=h_max
            do iter=1,maxiter
               gL=dsqrt(.5d0*g*(1/h0 + 1/hL))
               gR=dsqrt(.5d0*g*(1/h0 + 1/hR))
               F0=delu+(h0-hL)*gL + (h0-hR)*gR
               dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+
     &                   gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
               slope=2.d0*dsqrt(h0)*dfdh
               h0=(dsqrt(h0)-F0/slope)**2
            enddo
               hm=h0
               u1m=uL-(hm-hL)*dsqrt((.5d0*g)*(1/hm + 1/hL))
               u2m=uR+(hm-hR)*dsqrt((.5d0*g)*(1/hm + 1/hR))
               um=.5d0*(u1m+u2m)

               s1m=u1m-dsqrt(g*hm)
               s2m=u2m+dsqrt(g*hm)
               rare1=.false.
               rare2=.false.

         else !one shock one rarefaction
            h0=h_min

            do iter=1,maxiter
               F0=delu + 2.d0*(dsqrt(g*h0)-dsqrt(g*h_max))
     &                  + (h0-h_min)*dsqrt(.5d0*g*(1/h0+1/h_min))
               slope=(F_max-F0)/(h_max-h_min)
               h0=h0-F0/slope
            enddo

            hm=h0
            if (hL.gt.hR) then
               um=uL+2.d0*dsqrt(g*hL)-2.d0*dsqrt(g*hm)
               s1m=uL+2.d0*dsqrt(g*hL)-3.d0*dsqrt(g*hm)
               s2m=uL+2.d0*dsqrt(g*hL)-dsqrt(g*hm)
               rare1=.true.
               rare2=.false.
            else
               s2m=uR-2.d0*dsqrt(g*hR)+3.d0*dsqrt(g*hm)
               s1m=uR-2.d0*dsqrt(g*hR)+dsqrt(g*hm)
               um=uR-2.d0*dsqrt(g*hR)+2.d0*dsqrt(g*hm)
               rare2=.true.
               rare1=.false.
            endif
         endif
      endif

      return

      end ! subroutine riemanntype----------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine riemann_aug_JCP(maxiter,meqn,mwaves,hL,hR,huL,huR,
     &   hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

      ! solve shallow water equations given single left and right states
      ! This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
      ! Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation

      ! To use the original solver call with maxiter=1.

      ! This solver allows iteration when maxiter > 1. The iteration seems to help with
      ! instabilities that arise (with any solver) as flow becomes transcritical over variable topo
      ! due to loss of hyperbolicity.



      implicit none

      !input
      double precision fw(5,3)
      double precision sw(3)
      double precision hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR
      double precision hvL,hvR,vL,vR
      double precision drytol,g
      integer meqn,mwaves,maxiter


      !local
      integer m,mw,k,iter
      double precision A(3,3)
      double precision r(3,3)
      double precision lambda(3)
      double precision del(3)
      double precision beta(3)

      double precision delh,delhu,delphi,delb,delnorm
      double precision rare1st,rare2st,sdelta,raremin,raremax
      double precision criticaltol,convergencetol,raretol
      double precision s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
      double precision huRstar,huLstar,uRstar,uLstar,hstarHLL
      double precision deldelh,deldelphi
      double precision s1m,s2m,hm
      double precision det1,det2,det3,determinant
      double precision sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat

      logical rare1,rare2,rarecorrector,rarecorrectortest,sonic

      !determine wave speeds
      sL=uL-dsqrt(g*hL) ! 1 wave speed of left state
      sR=uR+dsqrt(g*hR) ! 2 wave speed of right state
      uhat=(dsqrt(hL)*uL + dsqrt(hR)*uR)/(dsqrt(hR)+dsqrt(hL)) ! Roe average
      chat=dsqrt(g*0.5d0*(hR+hL)) ! Roe average
      sRoe1=uhat-chat ! Roe wave speed 1 wave
      sRoe2=uhat+chat ! Roe wave speed 2 wave
      sE1=min(sRoe1,sL)
      sE2=max(sRoe2,sR)

      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL
      delnorm = delh**2 + delphi**2

      call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,
     &                                          1,drytol,g)


      lambda(1)= min(sE1,s2m) !Modified Einfeldt speed
      lambda(3)= max(sE2,s1m) !Modified Eindfeldt speed
      sE1=lambda(1)
      sE2=lambda(3)
      hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

c     !determine the middle entropy corrector wave------------------------
      rarecorrectortest=.false.
      rarecorrector=.false.
      if (rarecorrectortest) then
         sdelta=lambda(3)-lambda(1)
         raremin = 0.5d0
         raremax = 0.9d0
         if (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2d0
         if (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2d0
         if (rare1.or.rare2) then
            !see which rarefaction is larger
            rare1st=3.d0*(dsqrt(g*hL)-dsqrt(g*hm))
            rare2st=3.d0*(dsqrt(g*hR)-dsqrt(g*hm))
            if (max(rare1st,rare2st).gt.raremin*sdelta.and.
     &         max(rare1st,rare2st).lt.raremax*sdelta) then
                  rarecorrector=.true.
               if (rare1st.gt.rare2st) then
                  lambda(2)=s1m
               elseif (rare2st.gt.rare1st) then
                  lambda(2)=s2m
               else
                  lambda(2)=0.5d0*(s1m+s2m)
               endif
            endif
         endif
         if (hstarHLL.lt.min(hL,hR)/5.d0) rarecorrector=.false.
      endif

      do mw=1,mwaves
         r(1,mw)=1.d0
         r(2,mw)=lambda(mw)
         r(3,mw)=(lambda(mw))**2
      enddo
      if (.not.rarecorrector) then
         lambda(2) = 0.5d0*(lambda(1)+lambda(3))
c         lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
         r(1,2)=0.d0
         r(2,2)=0.d0
         r(3,2)=1.d0
      endif
c     !---------------------------------------------------

c     !determine the steady state wave -------------------
      criticaltol = 1.d-6
      deldelh = -delb
      deldelphi = -g*0.5d0*(hR+hL)*delb

c     !determine a few quanitites needed for steady state wave if iterated
      hLstar=hL
      hRstar=hR
      uLstar=uL
      uRstar=uR
      huLstar=uLstar*hLstar
      huRstar=uRstar*hRstar

      !iterate to better determine the steady state wave
      convergencetol=1.d-6
      do iter=1,maxiter
         !determine steady state wave (this will be subtracted from the delta vectors)
         if (min(hLstar,hRstar).lt.drytol.and.rarecorrector) then
            rarecorrector=.false.
            hLstar=hL
            hRstar=hR
            uLstar=uL
            uRstar=uR
            huLstar=uLstar*hLstar
            huRstar=uRstar*hRstar
            lambda(2) = 0.5d0*(lambda(1)+lambda(3))
c           lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
            r(1,2)=0.d0
            r(2,2)=0.d0
            r(3,2)=1.d0
         endif

         hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
         s1s2bar = 0.25d0*(uLstar+uRstar)**2 - g*hbar
         s1s2tilde= max(0.d0,uLstar*uRstar) - g*hbar

c        !find if sonic problem
         sonic=.false.
         if (dabs(s1s2bar).le.criticaltol) sonic=.true.
         if (s1s2bar*s1s2tilde.le.criticaltol) sonic=.true.
         if (s1s2bar*sE1*sE2.le.criticaltol) sonic = .true.
         if (min(dabs(sE1),dabs(sE2)).lt.criticaltol) sonic=.true.
         if (sE1.lt.0.d0.and.s1m.gt.0.d0) sonic = .true.
         if (sE2.gt.0.d0.and.s2m.lt.0.d0) sonic = .true.
         if ((uL+dsqrt(g*hL))*(uR+dsqrt(g*hR)).lt.0.d0) sonic=.true.
         if ((uL-dsqrt(g*hL))*(uR-dsqrt(g*hR)).lt.0.d0) sonic=.true.

c        !find jump in h, deldelh
         if (sonic) then
            deldelh =  -delb
         else
            deldelh = delb*g*hbar/s1s2bar
         endif
c        !find bounds in case of critical state resonance, or negative states
         if (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) then
            deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
            deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
         elseif (sE1.ge.criticaltol) then
            deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
            deldelh = max(deldelh,-hL)
         elseif (sE2.le.-criticaltol) then
            deldelh = min(deldelh,hR)
            deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
         endif

c        !find jump in phi, deldelphi
         if (sonic) then
            deldelphi = -g*hbar*delb
         else
            deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
         endif
c        !find bounds in case of critical state resonance, or negative states
         deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
         deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))

         del(1)=delh-deldelh
         del(2)=delhu
         del(3)=delphi-deldelphi

c        !Determine determinant of eigenvector matrix========
         det1=r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2))
         det2=r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1))
         det3=r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
         determinant=det1-det2+det3

c        !solve for beta(k) using Cramers Rule=================
         do k=1,3
            do mw=1,3
               do m=1,3
                  A(m,mw)=r(m,mw)
                  A(m,k)=del(m)
               enddo
            enddo
            det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
            det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
            det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
            beta(k)=(det1-det2+det3)/determinant
         enddo

         !exit if things aren't changing
         if (dabs(del(1)**2+del(3)**2-delnorm).lt.convergencetol) exit
         delnorm = del(1)**2+del(3)**2
         !find new states qLstar and qRstar on either side of interface
         hLstar=hL
         hRstar=hR
         uLstar=uL
         uRstar=uR
         huLstar=uLstar*hLstar
         huRstar=uRstar*hRstar
         do mw=1,mwaves
            if (lambda(mw).lt.0.d0) then
               hLstar= hLstar + beta(mw)*r(1,mw)
               huLstar= huLstar + beta(mw)*r(2,mw)
            endif
         enddo
         do mw=mwaves,1,-1
            if (lambda(mw).gt.0.d0) then
               hRstar= hRstar - beta(mw)*r(1,mw)
               huRstar= huRstar - beta(mw)*r(2,mw)
            endif
         enddo

         if (hLstar.gt.drytol) then
            uLstar=huLstar/hLstar
         else
            hLstar=max(hLstar,0.d0)
            uLstar=0.d0
         endif
         if (hRstar.gt.drytol) then
            uRstar=huRstar/hRstar
         else
            hRstar=max(hRstar,0.d0)
            uRstar=0.d0
         endif

      enddo ! end iteration on Riemann problem

      do mw=1,mwaves
         sw(mw)=lambda(mw)
         fw(1,mw)=beta(mw)*r(2,mw)
         fw(2,mw)=beta(mw)*r(3,mw)
         fw(3,mw)=beta(mw)*r(2,mw)
      enddo
      !find transverse components (ie huv jumps).
      fw(3,1)=fw(3,1)*vL
      fw(3,3)=fw(3,3)*vR
      fw(3,2)= hR*uR*vR - hL*uL*vL - fw(3,1)- fw(3,3)

      return

      end !subroutine riemann_aug_JCP-------------------------------------------------

c=======================================================================
c  file: gaussj
c  routine: gaussj
c  solves linear system ax=b directly using Guassian Elimination
c  a(1:n,1:n) is the input matrix of size np by np
c  b(1:n,1:m) is the right hand side (can solve for m vectors)
c  On output a is replaced by it's inverse, and b replaced by x.

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (dabs(a(j,k)).ge.big)then
                  big=dabs(a(j,k))
                  irow=j
                  icol=k
                endif

              else if (ipiv(k).gt.1) then
                write(*,*) 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) write(*,*) 'singular matrix in gaussj'

        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then

          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END


