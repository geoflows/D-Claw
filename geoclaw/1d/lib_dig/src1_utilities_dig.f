c ======================================================================
      subroutine admissibleq(theta,h,hu,hm,p,u,m,iflag)
c     !evaluate a potential solution and return the nearest admissible one
c ======================================================================

      implicit none

*     !i/o
      integer iflag(4)
      double precision theta,h,hu,hm,p,u,m

*     !local
      integer i
      double precision rho
      double precision pmax,mlo,mhi,phi,plo,hlo

      include "digparamsdec.i"
      include "digparamscommon.i"


      do i = 1,4
         iflag(i) = 0
      enddo

      if (h.le.dry_tol) then
         iflag(1) = 1
         h = max(0.d0,h)
         hu = 0.d0
         hm = h*m0
         p  = h*rho_f*grav*cos(theta)
         u = 0.d0
         m = m0
         return
      endif

      hlo = 2.d0*dry_tol
      if (h.lt.hlo) then
         h = (h**2 + hlo**2)/(2.d0*hlo)
      endif

      u  = hu/h
      m = hm/h
      m = min(m,1.d0)
      m = max(m,0.d0)
      !mlo =  phys_tol
      mlo = m_min
      mhi = 1.d0 - 0.5d0*mlo
      if (m.le.mlo) then
         iflag(3) = 1
         !write(*,*) 'mbefore',m
         m = (m**2 + mlo**2)/(2.d0*mlo)
         hm = h*m
         !write(*,*) 'm,mlo:',m,mlo
      elseif (m.ge.mhi) then
         iflag(3) = 1
         m = 1.d0 - ((1.d0-mhi)**2 + (1.d0-m)**2)/(2.d0*(1.d0-mhi))
         hm = h*m
      endif

      rho = (rho_s*m + (1.d0-m)*rho_f)
      pmax = rho*grav*h*cos(theta)
      p = min(pmax,p)
      p = max(0.d0,p)
      plo = phys_tol*rho*h*grav*cos(theta)
      phi = pmax - plo
      !write(*,*) 'p,pmax,phi,plo',p,pmax,phi,plo
      if (p.lt.plo) then
         iflag(4) = 1
         p = (p**2 + plo**2)/(2.d0*plo)
      elseif (p.gt.phi) then
         iflag(4) = 1
         p = pmax - ((pmax-p)**2+ (pmax-phi)**2)/(2.d0*(pmax-phi))
      endif
      if (p.gt.pmax) write(*,*) 'pnew,pmax',p
      return
      end

c ======================================================================
      subroutine auxeval(h,u,m,p,kappa,theta,
     &              S,rho,tanpsi,D,tau,phi,sigbed,kperm,compress,pm)
c     !evaluate the aux variables for this single grid cell
c ======================================================================
      implicit none

*     !i/o
      double precision, intent(in)  :: h,u,m,p,theta,pm
      double precision, intent(out) :: S,rho,tanpsi,D,tau,kappa
      double precision, intent(out) :: phi,sigbed,kperm,compress

*     !local
      double precision g,m_eqn,hdel,hpos,sqrtarg
      double precision sigmin,SN,SNN,S0,s02,m_eqnS,m_eqnN,m_eqnC,sigpos
      double precision m_eqncrit,m_eqnboyer,mpos

      include "digparamsdec.i"
      include "digparamscommon.i"

      g=grav

      rho = rho_s*m + rho_f*(1.d0-m)
      sigbed = max(0.d0,rho*g*h*cos(theta) - p)
      SN = mu*abs(u)/((h**2)*(rho_s-rho_f)*m*g*cos(theta))
      SNN = mu*abs(u)/(h*sigbed)
      S = SNN

      m_eqncrit = m_crit*(1.d0 - tanh(c2*SN))
      !m_eqnboyer = m_crit/(1.d0+sqrt(S))
      if (h*sigbed.gt.0.d0) then
         m_eqnboyer = sqrt(h*sigbed)*m_crit/
     &       (sqrt(h*sigbed)+ sqrt(mu*abs(u)))
      else
         m_eqnboyer = 0.d0
      endif
      m_eqn = m_eqnboyer
      tanpsi = c1*(m-m_eqn)
      tau = max(0.d0,sigbed*tan(phi))!+atan(tanpsi))
      !kperm = (kappita**2*(1.d0-m)**3)/(180.d0*m**2)
      kperm = kappita**2*exp(max(0.d0,m-m_min)/(-0.03))/40.0
      compress = alpha/((m)*(sigbed + 1.d5))
      if (p_initialized.gt.0) then
         D = (2.d0*kperm/(h*mu))*(rho_f*g*h*cos(theta)-p)
      else
         D = 0.d0
      endif

c----------------------------------------------------------------
      kappa = 1.d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c-----------------------------------------------------------------

      return
      end

c ======================================================================
      subroutine psieval(theta,tau,rho,D,tanpsi,kappa,kperm,
     &                       compress,h,u,m,psi)
c     !evaluate the source term that is being integrated
c ======================================================================

      implicit none

*     !i/o
      double precision psi(4)
      double precision theta,tau,rho,D,tanpsi,kappa,kperm,compress
      double precision h,u,m,drytol

*     !local
      double precision taushear,zeta

      include "digparamsdec.i"
      include "digparamscommon.i"
      drytol = 2.d0*dry_tol
      taushear = (tau/rho)*sign(1.d0,u)
      if (h.gt.drytol) then
         taushear = taushear + (1.d0-m)*mu*u/(h*rho)
      endif

      !psi(2) = grav*h*sin(theta) - taushear + u*D*(rho-rho_f)/rho
      if (h.gt.drytol) then
         psi(1) = D*(rho-rho_f)/rho
         psi(2) = grav*h*sin(theta) - taushear + u*D*(rho-rho_f)/rho
         psi(3) = -D*m*(rho_f/rho)
         zeta = 6.d0/(compress*h*(1.d0+kappa))  +
     &        1.5d0*(rho-rho_f)*rho_f*grav*sin(theta)/(6.d0*rho)
         psi(4) = - 3.d0*abs(u)*tanpsi/(h*compress*(1.d0+kappa))
      else
         psi(1) = 0.d0
         psi(2) = 0.d0
         psi(3) = 0.d0
         psi(4) = 0.d0
      endif

      return
      end

c ======================================================================
      subroutine auxeval_old(h,u,m,p,kappa,theta,
     &              S,rho,tanpsi,D,tau,phi,sigbed,kperm,compress,pm)
c     !evaluate the aux variables for this single grid cell
c ======================================================================
      implicit none

*     !i/o
      double precision, intent(in)  :: h,u,m,p,theta,pm
      double precision, intent(out) :: S,rho,tanpsi,D,tau,kappa
      double precision, intent(out) :: phi,sigbed,kperm,compress

*     !local
      double precision g,m_eqn,hdel,hpos,sqrtarg
      double precision sigmin,SN,SNN,S0,s02,m_eqnS,m_eqnN,m_eqnC,sigpos
      double precision m_eqncrit,mpos

      include "digparamsdec.i"
      include "digparamscommon.i"

      g=grav

      rho = rho_s*m + rho_f*(1.d0-m)
      sigbed = max(0.d0,rho*g*h*cos(theta) - p)
c      if (h.ge.phys_tol.and.sigbed.gt.0.d0) then
c         S = 2.d0*(rho_s/sigbed)*(u*delta/h)**2
c      elseif (h.lt.phys_tol.and.sigbed.gt.0.d0) then
c         S = 2.d0*(rho_s/sigbed)*u**2
c      else
c         S = 1.d16
c      endif
c      S = 2.d0*(rho_s/(sigbed+delta))*(u*delta/(h+delta))**2
c      S = 2.d0*(rho_s/(sigbed+p+delta))*(u*delta/(h+delta))**2 !!!
c      S = 2.d0*(rho_s/(sigbed + p))*(u*delta/(h+delta))**2 !!!!
      hdel = 2.d0*delta
      !hdel = dry_tol
      mpos = max(m,1.d-6)
      if (h.gt.hdel) then
         hpos= h
      else
         hpos = (h**2 + hdel**2)/(2.d0*hdel)
      endif
c            hpos = 1.d0
      sigmin = rho_f*hdel*grav*cos(theta)
      if (sigbed.gt.sigmin) then
         sigpos = sigbed
      else
         sigpos = (sigbed**2 + sigmin**2)/(2.d0*sigmin)
      endif

      sigpos = sigbed + delta
      hpos = h + delta
      SN = mu*abs(u)/((hpos**2)*(rho_s-rho_f)*mpos*g*cos(theta))
c      SNN = mu*abs(u)/((hpos**2)*(sigpos))
c      S0 = 2.d0*(rho_s/(sigpos))*(u*delta/(hpos))**2
c      S02 = 2.d0*(rho_s/(sigpos+p))*(u*delta)**2
      S = SN

      m_eqnS = m_crit - (m_crit - m_min)*tanh(S)
      m_eqnN = m_crit - c2*S
      m_eqnC = ((rho - rho_f)/(rho_s - rho_f))*(1.d0 - c2*tanh(S))
      m_eqncrit = m_crit*(1.d0 - tanh(c2*S))
      m_eqn = m_eqncrit

      tanpsi = c1*(m-m_eqn)
      !phi = max(1.d0,phi + atan(tanpsi))
      tau = max(0.d0,sigbed*tan(phi))!+atan(tanpsi))
      kperm = (kappita**2*max(1.d0-m,1.d-6)**3)/(180.d0*max(m,1.d-6)**2)
      compress = alpha/(max(1.d-16,m)*
     &               (max(rho*grav*cos(theta)*h - p,0.d0) + 1.d5))
      if (h.gt.dry_tol.and.p_initialized.gt.0) then
         D = (2.d0*kperm/(h*mu))*(rho_f*g*h*cos(theta)-p)
      else
         D = 0.d0
      endif
      if (phi_int.eq.phi_bed) then
         sqrtarg = 0.d0
      else
         sqrtarg = 1.d0-(cos(phi_int)**2)*(1.d0 + tan(phi_bed)**2)
      endif

      kappa = (2.d0 - pm*2.d0*sqrt(sqrtarg))/(cos(phi_int)**2)
      kappa = kappa - 1.d0
c----------------------------------------------------------------
      kappa = 1.d0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !D = 0.d0
      !kperm = 1.d-8
      !compress = 5.d-5
c-----------------------------------------------------------------

      return
      end

c ======================================================================
      subroutine admissibleq_old(theta,h,hu,hm,p,u,m,iflag)
c     !evaluate a potential solution and return the nearest admissible one
c ======================================================================

      implicit none

*     !i/o
      integer iflag(4)
      double precision theta,h,hu,hm,p,u,m

*     !local
      integer i
      double precision pmax,rho

      include "digparamsdec.i"
      include "digparamscommon.i"

      do i = 1,4
         iflag(i) = 0
         enddo

      if (h.le.dry_tol) then
         iflag(1) = 1
         !if (h.lt.0.d0) write(*,*) 'h',h
         h = max(0.d0,h)
         hu = 0.d0
         hm = h*m0 !max(0.d0,min(hm,h))
         p  = h*rho_f*grav*cos(theta) !max(0.d0,min(p,rho_f*grav*h*cos(theta)))
         u = 0.d0
         m = m0
      else
         if (hm.gt.h.or.hm.lt.0.d0) then
            iflag(3) = 1
            hm = max(hm,0.d0)
            hm = min(hm,h)
            endif
         u  = hu/h
         m  = hm/h
         pmax = (rho_s*m + (1.d0-m)*rho_f)*grav*h*cos(theta)
         if (p.lt.0.d0.or.p.gt.pmax) then
             iflag(4) = 1
             p = max(p,0.d0)
             p  = min(p,pmax)
            endif
         endif
      return
      end

