
c ===============================================================
       subroutine calcaux(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c ===============================================================
c
      implicit none

      double precision q(1-mbc:maxmx+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, maux)

      integer maxmx,meqn,mbc,mx,maux
      double precision xlower,dx

      integer i,ma,ibc
      logical iflag(4)
      double precision h,u,m,pbed,rho,S,m_eqn,tanpsi,phi,compress,pms,
     &              kperm,tau,up,um,dudx,pm,sqrtarg,kappa,D,theta,g
      double precision sigbed


      include "digparamsdec.i"
      include "digparamscommon.i"

      g=grav

      do i=1,mx
         theta = aux(i,i_theta)
         phi = aux(i,i_phi)
         call admissibleq(theta,q(i,1),q(i,2),q(i,3),q(i,4),u,m,iflag)
         h=q(i,1)
         if (h.le.dry_tol) then
            aux(i,i_S) = 0.d0
            aux(i,i_rho) = 0.d0
            aux(i,i_tanpsi) = 0.d0
            aux(i,i_D) = 0.d0
            aux(i,i_tau) = 0.d0
            aux(i,i_kappa) = 1.d0
            aux(i,i_sigbed) = 0.d0
            aux(i,i_kperm) = 0.d0
            aux(i,i_alpha) = 0.d0
            cycle
         endif
         pbed = q(i,4)
         !write(*,*) 'pbed:',pbed
         if (q(i,1).gt.0.d0.and.q(i-1,1).gt.0.d0) then
            dudx = (q(i,2)/q(i,1) - q(i-1,2)/q(i-1,1))/dx
         else
            dudx = 0.d0
         endif
         pms = abs(tanh(dudx/dudx_eps))
         pm = sign(1.d0,dudx)

         call auxeval(h,u,m,pbed,kappa,theta,
     &              S,rho,tanpsi,D,tau,phi,sigbed,kperm,compress,pm)

         aux(i,i_S) = S
         aux(i,i_rho) = rho
         aux(i,i_tanpsi) = tanpsi
         aux(i,i_D) = D
         aux(i,i_tau) = tau
         aux(i,i_kappa) = kappa
         aux(i,i_sigbed) = sigbed
         aux(i,i_kperm) = kperm
         aux(i,i_alpha) = compress

      enddo

      do ibc = 1,mbc
         do ma = 1,maux
            aux(1-ibc,ma) = aux(1,ma)
            aux(mx+ibc,ma) = aux(mx,ma)
            enddo
         enddo

      return
      end