
c
c
c =========================================================
      subroutine src1(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
c =========================================================
      implicit none

*     !i/o
      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, *)

      integer maxmx,meqn,mbc,mx,maux
      double precision xlower,dx,t,dt

*     !local
      double precision psi(4)
      integer iflag(4)
      integer i,k,maxsteps,meq,ifail
      double precision h,hu,hm,p,m,u,D,tau,rho,theta,tanpsi,g
      double precision hk,huk,hmk,pk,norm,tol,netforce
      double precision t0,tf

      include "digparamsdec.i"
      include "digparamscommon.i"

      g=grav
      maxsteps = 1
      tol = 10.d0
      t0 = t
      tf = t0 + dt
c      call calcaux(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
      do i=0,mx+1
         call admissibleq(theta,q(i,1),q(i,2),q(i,3),q(i,4),u,m,iflag)
         if (q(i,1).gt.phys_tol) then
            theta = aux(i,i_theta)
*           !note: hyperbolic step might have returned inadmissible q
            call auxeval(theta,tau,rho,D,tanpsi,q(i,1),u,m,q(i,4))

*           !begin integrating using adaptive timestepping
            ifail = 0

            call odeint_adapt(maxsteps,tol,t0,tf,theta,q(i,1:4),u,m,
     &                              ifail)

            if (ifail.gt.0) then
c               write(*,*) 'SRC1: Failure integrating source term'
c               write(*,*) 'SRC1: either increase tol or maxiter'
c               write(*,*) 'SRC1: i,ifail:', i,ifail
c               write(*,*) 'SRC1: q:', (q(i,k),k=1,4)
c               stop
            endif
            call admissibleq(theta,q(i,1),q(i,2),q(i,3),q(i,4),u,m,
     &                                                         iflag)
         elseif (q(i,1).gt.dry_tol) then
            call auxeval(theta,tau,rho,D,tanpsi,q(i,1),u,m,q(i,4))
            q(i,2) = q(i,2) + grav*q(i,1)*sin(theta)*dt
            call friction(theta,dt,q(i,1),q(i,2),q(i,3),q(i,4))
            call fluidfriction(theta,dt,q(i,1),q(i,2),q(i,3),q(i,4))
c            q(i,3) = m0*q(i,1)
c            q(i,4) = rho_f*grav*q(i,1)*cos(theta)
         endif

      enddo

      call calcaux(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)

      return
      end


c ======================================================================
      subroutine odeint_adapt(maxsteps,tol,t0,tf,theta,q,u,m,ifail)
c     !integrate stiff ODEs with adaptive timestep and error estimation
c     !proceed from t0 to tf. Call a seperate routine to take 1 step
c ======================================================================

      implicit none

*     !i/o
      integer maxsteps,ifail
      double precision tol,t0,tf,theta,u,m
      double precision q(4)

*     !local
      integer iflag(4)
      logical istep
      integer meq,n,mid,j
      double precision tn,dt,tau,rho,D,tanpsi,convtol,dtmax,dtmin
      double precision lambda,kprate
      double precision q0(4),qn(4)

      include "digparamsdec.i"
      include "digparamscommon.i"

      return

      do meq=1,4
         q0(meq) = q(meq)
         enddo

*     !set-up needed variables
      tn = t0
      dtmax = tf-t0
      dtmin = (tf-t0)/real(maxsteps)
      call admissibleq(theta,q0(1),q0(2),q0(3),q0(4),u,m,iflag)
      call auxeval(theta,tau,rho,D,tanpsi,q0(1),u,m,q0(4))
      call admissibledt(dtmax,theta,tau,rho,D,tanpsi,q0(1),u,m,q0(4))
      dt = dtmax
      dt = max(dtmax,dtmin)
      convtol = tol*dtmax**2

      do n=1,maxsteps + 1 !timestepping loop
         !begin next timestep, trying dt to completion
         call admissibleq(theta,q0(1),q0(2),q0(3),q0(4),u,m,iflag)
         if (iflag(1).eq.1) exit !cannot integrate undefined source

         call one_step(convtol,tn,dt,dtmax,theta,q0,qn,istep,iflag)

         if (tn.ge.tf) then !we are done
            do meq=1,4
               q(meq) = qn(meq)
               enddo
            ifail = 0
            exit
         elseif (istep) then !continue with a new step
            dtmax = tf - tn
            do meq = 1,4
               q0(meq) = qn(meq)
               enddo
            call auxeval(theta,tau,rho,D,tanpsi,q0(1),u,m,q0(4))
            call admissibledt(dtmax,theta,tau,rho,D,tanpsi,q0(1),
     &                                                   u,m,q0(4))
            dt = max(dtmax,dtmin)
         else !last step failed. Exit as is
            ifail = 1
            exit
            endif

         enddo

      if (tn.lt.tf) then
         do meq=1,4
            q(meq) = q0(meq)
            enddo
         ifail = ifail + 1
         endif

      return
      end