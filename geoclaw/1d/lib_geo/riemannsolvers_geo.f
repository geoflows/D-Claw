c-----------------------------------------------------------------------
      subroutine riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
     &       hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,fw)

      ! solve shallow water equations given single left and right states
      ! steady state wave is subtracted from delta [q,f]^T before decomposition 

      implicit none 

      !input
      integer meqn,mwaves,maxiter

      double precision hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
      double precision vL,vR,hvL,hvR
      double precision drytol,g

      !local
      integer iter

      logical subcritical,superpos,superneg,critical
      
      double precision delh,delhu,delphi,delb,delhdecomp,delphidecomp
      double precision s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
      double precision uRstar,uLstar,hstarHLL
      double precision bound,bound1,bound2,hLb,hRb
      double precision deldelh,deldelphi
      double precision alpha1,alpha2,beta1,beta2,delalpha1,delalpha2
      double precision criticaltol,convergencetol,rat
      double precision sL,sR
      double precision uhat,chat,sRoe1,sRoe2

      double precision fw(meqn,mwaves)

      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL

      convergencetol= abs(delb)/1.d2
      criticaltol = 1.d-6

      deldelh = -delb
      deldelphi = -g*0.5d0*(hR+hL)*delb

!     !if no source term, skip determining steady state wave
      if (abs(delb).gt.0.d0) then
!
         !determine a few quanitites needed for steady state wave if iterated
         hLstar=hL
         hRstar=hR
         uLstar=uL
         uRstar=uR

         alpha1=0.d0
         alpha2=0.d0
         
!        !iterate to better determine Riemann problem
         do iter=1,maxiter

            !determine steady state wave (this will be subtracted from the delta vectors)
            hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
            s1s2bar = 0.25d0*(uLstar+uRstar)**2 - g*hbar
            s1s2tilde= max(0.d0,uLstar*uRstar) - g*hbar
 
!           !bounds in case of critical state resonance, or negative states
            hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

            if (s1s2bar.eq.0.d0) then
               deldelh =  -delb
            else
               deldelh = delb*g*hbar/s1s2bar
            endif
               
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

            deldelphi = -g*hbar*delb*s1s2tilde/s1s2bar
            deldelphi = -deldelh*s1s2tilde

            deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
            deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))

!---------determine fwaves ------------------------------------------

!           !first decomposition
            delhdecomp = delh-deldelh
            delalpha1 = (sE2*delhdecomp - delhu)/(sE2-sE1)-alpha1
            alpha1 = alpha1 + delalpha1
            delalpha2 = (delhu - sE1*delhdecomp)/(sE2-sE1)-alpha2
            alpha2 = alpha2 + delalpha2

            !second decomposition
            delphidecomp = delphi - deldelphi
            beta1 = (sE2*delhu - delphidecomp)/(sE2-sE1)
            beta2 = (delphidecomp - sE1*delhu)/(sE2-sE1)
             
            if ((delalpha2**2+delalpha1**2).lt.convergencetol**2) then
               exit
            endif
!
            if (sE2.gt.criticaltol.and.sE1.lt.-criticaltol) then
               hLstar=hL+alpha1
               hRstar=hR-alpha2
c               hustar=huL+alpha1*sE1
               hustar = huL + beta1
            elseif (sE1.ge.criticaltol) then
               hLstar=hL
               hustar=huL
               hRstar=hR - alpha1 - alpha2
            elseif (sE2.le.-criticaltol) then
               hRstar=hR
               hustar=huR
               hLstar=hL + alpha1 + alpha2
            endif
!
            if (hLstar.gt.drytol) then
               uLstar=hustar/hLstar
            else
               hLstar=max(hLstar,0.d0)
               uLstar=0.d0
            endif
!
            if (hRstar.gt.drytol) then
               uRstar=hustar/hRstar
            else
               hRstar=max(hRstar,0.d0)
               uRstar=0.d0
            endif

         enddo
      endif

      delhdecomp = delh - deldelh
      delphidecomp = delphi - deldelphi

      !first decomposition
      alpha1 = (sE2*delhdecomp - delhu)/(sE2-sE1)
      alpha2 = (delhu - sE1*delhdecomp)/(sE2-sE1)

      !second decomposition
      beta1 = (sE2*delhu - delphidecomp)/(sE2-sE1)
      beta2 = (delphidecomp - sE1*delhu)/(sE2-sE1)

      ! 1st nonlinear wave
      fw(1,1) = alpha1*sE1
      fw(2,1) = beta1*sE1
      fw(3,1) = 0.d0
      ! advection of transverse wave
      fw(1,2) = 0.d0
      fw(2,2) = 0.d0
      fw(3,2) = hR*uR*vR - hL*uL*vL
      ! 2nd nonlinear wave
      fw(1,3) = alpha2*sE2
      fw(2,3) = beta2*sE2
      fw(3,3) = 0.d0


      return

      end subroutine !-------------------------------------------------


c-----------------------------------------------------------------------
      subroutine riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
     &            bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,fw)

      ! solve shallow water equations given single left and right states
      ! solution has two waves.
      ! flux - source is decomposed.

      implicit none 

      !input
      integer meqn,mwaves

      double precision hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
      double precision hvL,hvR,vL,vR
      double precision drytol,g

      double precision fw(meqn,mwaves)

      !local
      double precision delh,delhu,delphi,delb,delhdecomp,delphidecomp
      double precision deldelh,deldelphi
      double precision beta1,beta2


      !determine del vectors
      delh = hR-hL
      delhu = huR-huL
      delphi = phiR-phiL
      delb = bR-bL

      deldelphi = -g*0.5d0*(hR+hL)*delb
      delphidecomp = delphi - deldelphi

      !flux decomposition
      beta1 = (sE2*delhu - delphidecomp)/(sE2-sE1)
      beta2 = (delphidecomp - sE1*delhu)/(sE2-sE1)

      ! 1st nonlinear wave
      fw(1,1) = beta1
      fw(2,1) = beta1*sE1
      fw(3,1) = 0.d0
      ! advection of transverse wave
      fw(1,2) = 0.d0
      fw(2,2) = 0.d0
      fw(3,2) = hR*uR*vR - hL*uL*vL
      ! 2nd nonlinear wave
      fw(1,3) = beta2
      fw(2,3) = beta2*sE2
      fw(3,3) = 0.d0
      return

      end subroutine !-------------------------------------------------
