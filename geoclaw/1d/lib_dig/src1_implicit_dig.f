
c ======================================================================
      subroutine one_step(convtol,tn,dt,dtmax,theta,q0,qn,istep,iflag)
c     !find maximum possible dt and take it.
c     !failure is based on error, or leaving space of
c     !admissible solutions
c ======================================================================

      implicit none

*     !i/o
      integer iflag(4)
      logical istep
      double precision convtol,tn,dt,dtmax,theta
      double precision q0(4),qn(4)

*     !local
      integer md,mdmax,maxiter,miniter,iter,meq
      double precision sgnu,tau,rho,D,tanpsi,m,u
      include "digparamsdec.i"
      include "digparamscommon.i"

*     !try implicit solve
      mdmax = 10
      miniter = 5
      maxiter = 50
      do md = 1,mdmax !find a stable time step and take it if possible
         call backwardeuler(iter,maxiter,convtol,dt,theta,q0,qn,istep,
     &                        iflag)
         if (istep) then
            tn = tn + dt
             if (qn(1).gt.dry_tol.and.abs(qn(2)).gt.0.d0) then
               call friction(theta,dt,qn(1),qn(2),qn(3),qn(4))
c               call fluidfriction(theta,dt,qn(1),qn(2),qn(3),qn(4))
               do meq = 1,4
                  q0(meq) = qn(meq)
               enddo
               call be_dil(iter,maxiter,convtol,dt,theta,q0,qn,istep,
     &                     iflag)
             endif
c            if (iter.le.miniter) dt = 2.d0*dt
            exit
         else
            dt = 0.5d0*dt
            endif

         enddo

      return
      end

c ======================================================================
      subroutine trapezoid(iter,maxiter,convtol,
     &                                       dt,theta,q0,qn,istep,iflag)

c ======================================================================

      implicit none

*     !i/o
      integer iflag(4)
      integer maxiter,iter
      double precision theta,convtol,dt
      double precision q0(4),qn(4)
      logical istep

*     !local
      integer mrow,mcol,meq,k,ip,ihm,ih
      double precision tau,rho,D,m,u,tanpsi,norm,xnorm
      double precision xk(4),psi0(4),psik(4),F(4)
      double precision jac(4,4),Fprime(4,4)

      include "digparamsdec.i"
      include "digparamscommon.i"

      call admissibleq(theta,q0(1),q0(2),q0(3),q0(4),u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,q0(1),u,m,q0(4))
      call psieval(theta,tau,rho,D,tanpsi,q0(1),u,m,q0(4),psi0)

      do meq = 1,4
         xk(meq) = q0(meq)
         enddo

      istep = .false.
      do k = 1,maxiter
         call admissibleq(theta,xk(1),xk(2),xk(3),xk(4),u,m,iflag)
         if (iflag(1).eq.1) exit
c         if (iflag(3).eq.1) exit
c         if (iflag(4).eq.1) exit
         call auxeval(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4))
         call psieval(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),psik)
         call jacobianA(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),jac)
         call jacobian_dil(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),jac)
         call jacobian_tau(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),jac)

         do mrow = 1,4
            do mcol = 1,4
                  Fprime = 0.5d0*dt*jac(mrow,mcol)
               enddo
            enddo
         do meq = 1,4
            Fprime(meq,meq) = Fprime(meq,meq) - 1.d0
            F(meq) = 0.5d0*dt*(psik(meq)+psi0(meq)) + (q0(meq)-xk(meq))
            enddo
*        !we are solving: Fprime * (x_k - x_k+1) = F(x_k)
*        !the answer, x, to Ax = b, is returned by gausssj in b
*        !=> x_k+1 = x_k - F(x_k)
         call gaussj(Fprime,4,4,F,1,1)

         do meq = 1,4
            xk(meq) = xk(meq) - F(meq)
            enddo
         norm = abs(F(1))
         do meq = 2,4
            norm = norm +abs(F(meq))*dble(1-iflag(meq))
            enddo
         call admissibleq(theta,xk(1),xk(2),xk(3),xk(4),u,m,iflag)

         if (norm.lt.convtol) then
            istep = .true.
            do meq=1,4
               qn(meq) = xk(meq)
               enddo
            exit
            endif
         enddo
      iter = k
      return
      end


c ======================================================================
      subroutine backwardeuler(iter,maxiter,convtol,
     &                                       dt,theta,q0,qn,istep,iflag)

c ======================================================================

      implicit none

*     !i/o
      integer iflag(4)
      integer maxiter,iter
      double precision theta,convtol,dt
      double precision q0(4),qn(4)
      logical istep

*     !local
      integer mrow,mcol,meq,k,ip,ihm,ih
      double precision tau,rho,D,m,u,tanpsi,norm,xnorm
      double precision xk(4),psi0(4),psik(4),F(4)
      double precision jac(4,4),Fprime(4,4)

      include "digparamsdec.i"
      include "digparamscommon.i"

      call admissibleq(theta,q0(1),q0(2),q0(3),q0(4),u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,q0(1),u,m,q0(4))

      do meq = 1,4
         xk(meq) = q0(meq)
         enddo

      istep = .false.
      do k = 1,maxiter
         call admissibleq(theta,xk(1),xk(2),xk(3),xk(4),u,m,iflag)
         if (iflag(1).eq.1) exit
c         if (iflag(3).eq.1) exit
c         if (iflag(4).eq.1) exit
         call auxeval(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4))
         call psievalA(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),psik)
         call jacobianA(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),jac)
c         call jacobian_dil(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),jac)
c         call jacobian_tau(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),jac)

         do mrow = 1,4
            do mcol = 1,4
                  Fprime = 0.5d0*dt*jac(mrow,mcol)
               enddo
            enddo
         do meq = 1,4
            Fprime(meq,meq) = Fprime(meq,meq) - 1.d0
            F(meq) = dt*psik(meq) + (q0(meq)-xk(meq))
            enddo
*        !we are solving: Fprime * (x_k - x_k+1) = F(x_k)
*        !the answer, x, to Ax = b, is returned by gausssj in b
*        !=> x_k+1 = x_k - F(x_k)
         call gaussj(Fprime,4,4,F,1,1)

         do meq = 1,4
            xk(meq) = xk(meq) - F(meq)
            enddo
         norm = abs(F(1))
         do meq = 2,4
            norm = norm +abs(F(meq))*dble(1-iflag(meq))
            enddo
         call admissibleq(theta,xk(1),xk(2),xk(3),xk(4),u,m,iflag)

         if (norm.lt.convtol) then
            istep = .true.
            do meq=1,4
               qn(meq) = xk(meq)
               enddo
            exit
            endif
         enddo
      iter = k
      return
      end

c ======================================================================
      subroutine be_dil(iter,maxiter,convtol,
     &                                       dt,theta,q0,qn,istep,iflag)

c ======================================================================

      implicit none

*     !i/o
      integer iflag(4)
      integer maxiter,iter
      double precision theta,convtol,dt
      double precision q0(4),qn(4)
      logical istep

*     !local
      integer mrow,mcol,meq,k,ip,ihm,ih
      double precision tau,rho,D,m,u,tanpsi,norm,xnorm
      double precision xk(4),psi0(4),psik(4),F(4)
      double precision jac(4,4),Fprime(4,4)

      include "digparamsdec.i"
      include "digparamscommon.i"

      call admissibleq(theta,q0(1),q0(2),q0(3),q0(4),u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,q0(1),u,m,q0(4))

      do meq = 1,4
         xk(meq) = q0(meq)
         enddo

      istep = .false.
      do k = 1,maxiter
         call admissibleq(theta,xk(1),xk(2),xk(3),xk(4),u,m,iflag)
         if (iflag(1).eq.1) exit
c         if (iflag(3).eq.1) exit
c         if (iflag(4).eq.1) exit
         call auxeval(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4))
         call psieval_dil(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),psik)
         call jacobian_dil(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),jac)
c         call jacobian_tau(theta,tau,rho,D,tanpsi,xk(1),u,m,xk(4),jac)

         do mrow = 1,4
            do mcol = 1,4
                  Fprime = 0.5d0*dt*jac(mrow,mcol)
               enddo
            enddo
         do meq = 1,4
            Fprime(meq,meq) = Fprime(meq,meq) - 1.d0
            F(meq) = dt*psik(meq) + (q0(meq)-xk(meq))
            enddo
*        !we are solving: Fprime * (x_k - x_k+1) = F(x_k)
*        !the answer, x, to Ax = b, is returned by gausssj in b
*        !=> x_k+1 = x_k - F(x_k)
         call gaussj(Fprime,4,4,F,1,1)

         do meq = 1,4
            xk(meq) = xk(meq) - F(meq)
            enddo
         norm = abs(F(1))
         do meq = 2,4
            norm = norm +abs(F(meq))!*dble(1-iflag(meq))
            enddo
         call admissibleq(theta,xk(1),xk(2),xk(3),xk(4),u,m,iflag)

         if (norm.lt.convtol) then
            istep = .true.
            do meq=1,4
               qn(meq) = xk(meq)
               enddo
            exit
            endif
         enddo
      iter = k
      return
      end

c ======================================================================
      subroutine psieval(theta,tau,rho,D,tanpsi,h,u,m,p,psi)
c     !evaluate the source term that is being integrated
c ======================================================================

      implicit none

*     !i/o
      double precision psi(4)
      double precision theta,tau,rho,D,tanpsi,h,u,m,p

*     !local
      double precision taushear

      include "digparamsdec.i"
      include "digparamscommon.i"

      taushear = (tau/rho)*tanh(u/smallu)
      psi(1) = D*(rho-rho_f)/rho
      psi(2) = u*D*(rho-rho_f)/rho + grav*h*sin(theta)- taushear
      psi(3) = -D*m*(rho_f/rho)
      if (h.gt.phys_tol) then
         psi(4) = (2.d0/(alpha*h))*(D - 0.*abs(u)*tanpsi)
      else
         psi(4) = 0.d0
      endif


      return
      end

c ======================================================================
      subroutine psievalA(theta,tau,rho,D,tanpsi,h,u,m,p,psi)
c     !evaluate the source term that is being integrated
c ======================================================================

      implicit none

*     !i/o
      double precision psi(4)
      double precision theta,tau,rho,D,tanpsi,h,u,m,p

*     !local
      double precision taushear

      include "digparamsdec.i"
      include "digparamscommon.i"

c      taushear = (tau/rho)*tanh(u/smallu)
      psi(1) = D*(rho-rho_f)/rho
      psi(2) = u*D*(rho-rho_f)/rho + grav*h*sin(theta)
      psi(3) = -D*m*(rho_f/rho)
      if (h.gt.phys_tol) then
         psi(4) = (2.d0/(alpha*h))*D
      else
         psi(4) = 0.d0
      endif


      return
      end

c ======================================================================
      subroutine psieval_dil(theta,tau,rho,D,tanpsi,h,u,m,p,psi)
c     !evaluate the source term that is being integrated
c ======================================================================

      implicit none

*     !i/o
      double precision psi(4)
      double precision theta,tau,rho,D,tanpsi,h,u,m,p

*     !local
      double precision taushear

      include "digparamsdec.i"
      include "digparamscommon.i"

c      taushear = (tau/rho)*tanh(u/smallu)
      psi(1) = 0.d0
      psi(2) = 0.d0
      psi(3) = 0.d0
      if (h.gt.phys_tol) then
         psi(4) = -(2.d0/(alpha*h))*abs(u)*tanpsi
      else
         psi(4) = 0.d0
      endif


      return
      end


c ======================================================================
      subroutine jacobianA(theta,tau,rho,D,tanpsi,h,u,m,p,jac)
c     !evaluate the source term Jacobian
c     !does not include basal shear stress or shear induced dilatency
c ======================================================================

      implicit none

*     !i/o
      double precision jac(4,4)
      double precision theta,tau,rho,D,tanpsi,h,u,m,p

*     !local
      double precision rho_rhofdrho,rhos_rhofdrho,rhofdrho,kdmuh

      include "digparamsdec.i"
      include "digparamscommon.i"

      rho_rhofdrho = (rho-rho_f)/rho
      rhos_rhofdrho = (rho_s-rho_f)/rho
      rhofdrho = rho_f/rho
      kdmuh = kappita/(mu*h)

      jac(1,1) = (rho_rhofdrho/h)*(kdmuh*p - rhofdrho*D)
      jac(1,2) = 0.d0
      jac(1,3) = (D/h)*rhofdrho*rhos_rhofdrho
      jac(1,4) = -kdmuh*rho_rhofdrho

      jac(2,1) = grav*sin(theta) +
     &                  (u/h)*rho_rhofdrho*(kdmuh*p + D*rho_rhofdrho)
      jac(2,2) = rho_rhofdrho*D/h
      jac(2,3) = (u*D/h)*rhofdrho*rhos_rhofdrho
      jac(2,4) = -kdmuh*u*rho_rhofdrho

      jac(3,1) = rhofdrho*(m/h)*(D - D*m*rhos_rhofdrho - kdmuh*p)
      jac(3,2) = 0.d0
      jac(3,3) = (D/h)*rhofdrho*(m*rhos_rhofdrho - 1.d0)
      jac(3,4) = m*kdmuh*rhofdrho

      jac(4,1) = (2.d0/(alpha*h*h))*(kdmuh*p - D)
      jac(4,2) = 0.d0
      jac(4,3) = 0.d0
      jac(4,4) = -2.d0*kdmuh/(alpha*h)


      return
      end

c ======================================================================
      subroutine jacobian_dil(theta,tau,rho,D,tanpsi,h,u,m,p,jac)
c     !add shear dilatency to the Jacobian
c ======================================================================

      implicit none

*     !i/o
      double precision jac(4,4)
      double precision theta,tau,rho,D,tanpsi,h,u,m,p

*     !local
      double precision g,signu,dtanpsidh,dtanpsidu,dtanpsidm,dtanpsidp
      double precision S,sig,delm,s2ah,us2ah,dtanhs,dsdh,dsdu,dsdm,dsdp

      include "digparamsdec.i"
      include "digparamscommon.i"

      g = grav

*     !shear induced dilatency
*     !source term is -(2*abs(u)/(alpha*h))*tan(psi)
      if (u.eq.0.d0) then
         return
      endif
      sig = tau/tan(phi_bed)
      delm = m_crit - m_min
      signu = sign(1.d0,u)
      s2ah = 2.d0*signu/(alpha*h)
      us2ah = u*s2ah

      jac(4,1) = jac(4,1) + (2.d0/h)*us2ah
      jac(4,2) = jac(4,2) - s2ah*tanpsi/h

      if (sig.gt.0.d0) then
         S = 2.d0*(rho_s/sig)*(u*delta/h)**2
      else
         S = 1.d16
      endif
      if (S.lt.100.d0) then
         dtanhs = 1.d0/(cosh(S)**2)
         dsdh = -S*(-4.d0/h + rho_f*g*cos(theta)/sig)
         dsdu = S*(2.d0/u)
         dsdm = -S*(rho_s-rho_f)*g*h*cos(theta)/sig
         dsdp = S/sig

         dtanpsidh = -c1*m/h + c1*delm*dtanhs*dsdh
         dtanpsidu = c1*delm*dtanhs*dsdu
         dtanpsidm = c1 + c1*delm*dtanhs*dsdm
         dtanpsidp = c1*delm*dtanhs*dsdp

      else
         dtanpsidh = -c1*m/h
         dtanpsidu = 0.d0
         dtanpsidm = c1
         dtanpsidp = 0.d0
      endif

      jac(4,1) = jac(4,1) - us2ah*dtanpsidh
      jac(4,2) = jac(4,2) - us2ah*dtanpsidu/h
      jac(4,3) = jac(4,3) - us2ah*dtanpsidm/h
      jac(4,4) = jac(4,4) - us2ah*dtanpsidp

      return
      end

c ======================================================================
      subroutine jacobian_tau(theta,tau,rho,D,tanpsi,h,u,m,p,jac)
c     !add basal shear stress
c ======================================================================

      implicit none

*     !i/o
      double precision jac(4,4)
      double precision theta,tau,rho,D,tanpsi,h,u,m,p

*     !local
      double precision rhos_rhofdrho,dtanhgudu
      double precision gamma,dtaudh,dtaudm,dtaudp,g

      include "digparamsdec.i"
      include "digparamscommon.i"

      g = grav
      rhos_rhofdrho = (rho_s-rho_f)/rho

*     !basal shear stress
      gamma = 1.d0/smallu
      dtaudh = rho_f*grav*cos(theta)*tan(phi_bed)
      dtaudm = (rho_s - rho_f)*g*h*cos(theta)*tan(phi_bed)
      dtaudp = -tan(phi_bed)
      if (abs(gamma*u).lt.100.d0) then
         dtanhgudu = gamma/(cosh(gamma*u)**2)
      else
         dtanhgudu = 0.d0
         endif

      jac(2,1) = jac(2,1) + (-tanh(gamma*u)/rho)*dtaudh
     &           - (tau*m*(rho_s-rho_f)*tanh(gamma*u)/(h*rho**2) )
     &           + dtanhgudu*tau*u/(rho*h)
      jac(2,2) = jac(2,2) - tau*dtanhgudu/(h*rho)
      jac(2,3) = jac(2,3) -
     &            (tanh(gamma*u)/(h*rho))*(dtaudm - tau*rhos_rhofdrho)
      jac(2,4) = jac(2,4) -tanh(gamma*u)*dtaudp/rho

      return
      end

c ======================================================================
      subroutine frictionshear(theta,dt,h,hu,hm,p)
c     !evaluate a potential solution and return the nearest admissible one
c ======================================================================


      implicit none

*     !i/o
      double precision theta,dt,h,hu,hm,p

*     !local
      integer iflag(4)
      double precision u,m,tau,rho,D,tanpsi

      include "digparamsdec.i"
      include "digparamscommon.i"

      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,h,u,m,p)

      if (hu.gt.0.d0) then
         hu = max(hu - dt*tau/rho, 0.d0)
      else
         hu = min(hu + dt*tau/rho, 0.d0)
      endif
      p = p - dt*(2.d0/alpha)*(abs(hu)/h**2)*tanpsi
      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      return
      end

c ======================================================================
      subroutine friction(theta,dt,h,hu,hm,p)
c     !evaluate a potential solution and return the nearest admissible one
c ======================================================================


      implicit none

*     !i/o
      double precision theta,dt,h,hu,hm,p

*     !local
      integer iflag(4)
      double precision u,m,tau,rho,D,tanpsi

      include "digparamsdec.i"
      include "digparamscommon.i"

      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,h,u,m,p)

      if (hu.gt.0.d0) then
         hu = max(hu - dt*tau/rho, 0.d0)
      else
         hu = min(hu + dt*tau/rho, 0.d0)
      endif
      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      return
      end

c ======================================================================
      subroutine fluidfriction(theta,dt,h,hu,hm,p)
c     !evaluate a potential solution and return the nearest admissible one
c ======================================================================


      implicit none

*     !i/o
      double precision theta,dt,h,hu,hm,p

*     !local
      integer iflag(4)
      double precision u,m,tau,rho,D,tanpsi,fterm

      include "digparamsdec.i"
      include "digparamscommon.i"

      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,h,u,m,p)

      fterm = 3.d0*(1.d0-m)*mu*u/(h*rho)
      if (hu.gt.0.d0) then
         hu = max(hu - dt*fterm, 0.d0)
      else
         hu = min(hu - dt*fterm, 0.d0)
      endif
      call admissibleq(theta,h,hu,hm,p,u,m,iflag)
      return
      end