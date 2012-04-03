


c-----------------------------------------------------------------------
      subroutine riemann_dig1_conservative(i,meqn,mwaves,hL,hR,huL,huR,
     &      hmL,hmR,pL,pR,uhat,uL,uR,mL,mR,kL,kR,rho_f,theta,
     &      thetaL,thetaR,rhoL,rhoR,gamma,kperm,compress,tanpsi,D,tau,
     &      g,gmod,drytol,sw,fw,w)

      ! solve the dig Riemann problem for debris flow on flume

      implicit none

      common /comxt/ dtcom,dxcom,tcom

      common /startparams/ pimin

*     !input
      integer i,meqn,mwaves
      double precision pimin
      double precision dtcom,dxcom,tcom
      double precision hL,hR,huL,huR,
     &            hmL,hmR,pL,pR,uhat,uL,uR,mL,mR,kL,kR,rho_f,theta,
     &        thetaL,thetaR,rhoL,rhoR,gamma,kperm,compress,tanpsi,D,tau,
     &        g,gmod,drytol

      double precision fw(meqn,mwaves),w(meqn,mwaves)
      double precision sw(mwaves)
      double precision psi(4)

*     !local
      integer m,mw,ks,mp
      double precision rho,k,h,u,mbar,dx
      double precision det1,det2,det3,detR
      double precision R(3,3),A(3,3),alpha(3),delq(3),delf(3),del(3)
      double precision beta(3)


      dx = dxcom

      !write(*,*) 'riemann solver is called'
      if (hL.ge.drytol.and.hR.ge.drytol) then
         rho = 0.5d0*(rhoL + rhoR)
         k = kR
         h = 0.5d0*(hL + hR)
         u = uhat
         mbar = 0.5d0*(mL + mR)
      elseif (hL.ge.drytol) then
         rho = rhoL
         k = kR
         h = hL
         u = uL
         mbar = mL
      else
         rho = rhoR
         k = kR
         h = hR
         u = uR
         mbar = mR
      endif

*     !determine R
      R(1,1) = 1.d0
      R(2,1) = sw(1)
      R(3,1) = gamma*rho*gmod

      R(1,2) = k-1.d0
      R(2,2) = uhat*(k-1.d0)
      R(3,2) = k*rho*gmod

      R(1,3) = 1.d0
      R(2,3) = sw(3)
      R(3,3) = gamma*rho*gmod

      !determine delq
      delq(1) = hR-hL
      delq(2) = huR - huL
      delq(3) = pR - pL

      delf(1) = delq(2)
      delf(2) = hR*uR**2 + 0.5d0*k*g*cos(thetaR)*hR**2 -
     &      (hL*uL**2 + 0.5d0*k*g*cos(thetaL)*hL**2)
      delf(3) = 0.d0

      del(1) = delf(1)
      del(2) = delf(2) + (1.d0-k)*delq(3)/rho
      del(3) = -rho*u*gmod*gamma*delq(1) + rho*gmod*gamma*delq(2)
     &         + u*delq(3)

*     !determine the source term
      do m=1,4
            psi(m)=0.d0
      enddo
      call psieval(theta,tau,rho,D,tanpsi,k,kperm,
     &                       compress,h,u,mbar,psi)


      del(1) = del(1) - dx*psi(1)
      del(3) = del(3) - dx*psi(4)
      if (abs(u).gt.1.d-6) then
         del(2) = del(2) - dx*psi(2)
      else
         if (abs(delf(2)-g*h*dx*sin(theta)).ge.abs(dx*tau/rho)) then
            del(2)=sign(abs(abs(delf(2)-g*h*dx*sin(theta))
     &                   -abs(dx*tau/rho)),delf(2)-g*h*dx*sin(theta))
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

         alpha(ks) = (det1-det2+det3)/detR
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
         if (m.eq.3) mp = 4
         do mw=1,3
            w(mp,mw)  = alpha(mw)*R(m,mw)
            fw(mp,mw) = beta(mw)*R(m,mw)
         enddo
      enddo

      if (i.ge.48.and.i.le.183.and..false.) then
         do m= 1,3
            write(*,*) 'i,m,beta',i,m,beta(m)
         enddo
      endif


      w(3,1) = alpha(1)*mbar
      w(3,3) = alpha(3)*mbar
      w(3,2) = hmR-hmL - mbar*(alpha(1)+alpha(3))

      fw(3,1) = beta(1)*mbar
      fw(3,3) = beta(3)*mbar
      fw(3,2) = hmR*uR-hmL*uL - mbar*(beta(1)+ beta(3))-dx*psi(3)

      return
      end

c-----------------------------------------------------------------------
      subroutine riemann_dig1(meqn,mwaves,hL,hR,huL,huR,hmL,hmR,pL,pR,
     &         uhat,uL,uR,mL,mR,kL,kR,rho_f,rhoL,rhoR,g,gmod,drytol,
     &         sw,fw,w)

      ! solve the dig Riemann problem for debris flow on flume

      implicit none

*     !input
      integer meqn,mwaves

      double precision hL,hR,huL,huR,hmL,hmR,pL,pR,
     &            uhat,uL,uR,mL,mR,kL,kR,rho_f,rhoL,rhoR,g,gmod,drytol

      double precision fw(meqn,mwaves),w(meqn,mwaves)
      double precision sw(mwaves)
      double precision efsum(4)

*     !local
      integer m,mw,ks,mp
      double precision rho,k,h,u,mbar,gamma
      double precision det1,det2,det3,detR
      double precision R(3,3),A(3,3),alpha(3),delq(3),dq(4)

      if (hL.ge.drytol.and.hR.ge.drytol) then
         rho = 0.5d0*(rhoL + rhoR)
         k = kR
         h = 0.5d0*(hL + hR)
         u = uhat
         mbar = 0.5d0*(mL + mR)
      elseif (hL.ge.drytol) then
         rho = rhoL
         k = kR
         h = hL
         u = uhat
         mbar = mL
      else
         rho = rhoR
         k = kR
         h = hR
         u = uhat
         mbar = mR
      endif

      gamma = 1.5d0*(rho_f/(6.d0*rho) + (1.d0+k)/4.d0)

*     !determine R
      R(1,1) = 1.d0
      R(2,1) = sw(1)
      R(3,1) = rho*gmod

      R(1,2) = k-1.d0
      R(2,2) = uhat*(k-1.d0)
      R(3,2) = k*rho*gmod

      R(1,3) = 1.d0
      R(2,3) = sw(3)
      R(3,3) = gamma*rho*gmod

      !determine delq
      delq(1) = hR-hL
      delq(2) = huR - huL
      delq(3) = pR - pL

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

         alpha(ks) = (det1-det2+det3)/detR
      enddo

      do m=1,3
         mp=m
         if (m.eq.3) mp = 4
         do mw=1,3
            w(mp,mw) = alpha(mw)*R(m,mw)
         enddo
      enddo

      w(3,1) = alpha(1)*mbar
      w(3,3) = alpha(3)*mbar
      w(3,2) = hmR-hmL - mbar*(alpha(1)+alpha(3))

      do mw=1,mwaves
         fw(1,mw)=w(2,mw)
      enddo

      do m=2,meqn
         do mw=1,mwaves
            fw(m,mw) = sw(mw)*w(m,mw)
         enddo
      enddo

      return
      end
