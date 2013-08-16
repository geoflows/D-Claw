! ============================================================================
!  Program:     digclaw_mod
!  File:        geoclaw_mod.f90
!  Created:     2012-04-10
!  Author:      David George
! ============================================================================
!      Copyright (C)  2012-04-10 David George <dgeorge@uw.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module digclaw_module

   use geoclaw_module

   implicit none

    ! ========================================================================
    ! General digclaw parameters
    ! ========================================================================
    double precision :: rho_s,rho_f,phi_bed,phi_int,delta,kappita
    double precision :: mu,alpha,m_crit,c1,m0,dudx_eps,phys_tol

    integer :: init_ptype,p_initialized,bed_normal
    double precision :: init_pmax_ratio,init_ptf2,init_ptf,init_pmin_ratio

    integer, parameter ::  i_dig    = 4 !Start of digclaw aux variables
    integer, parameter ::  i_phi    = i_dig
    integer, parameter ::  i_theta  = i_dig + 1
    integer, parameter ::  DIG_PARM_UNIT = 78


contains

    ! ========================================================================
    !  set_dig(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_dig(fname)

         implicit none

         ! Input
         character*25, intent(in), optional :: fname

         ! Locals
         double precision :: deg2rad
         integer, parameter :: iunit = 127
         character*25 :: file_name
         logical :: found_file


         deg2rad = pi/180.d0

         ! Read user parameters from setgeo.data
         if (present(fname)) then
            file_name = fname
         else
            file_name = 'setdig.data'
         endif
         inquire(file=file_name,exist=found_file)
         if (.not. found_file) then
            print *, 'You must provide a file ', file_name
            stop
         endif

         call opendatafile(iunit, file_name)

         read(iunit,*) rho_s
         read(iunit,*) rho_f
         read(iunit,*) phi_bed
         phi_bed = deg2rad*phi_bed
         read(iunit,*) phi_int
         phi_int = deg2rad*phi_int
         read(iunit,*) delta
         read(iunit,*) kappita
         read(iunit,*) mu
         read(iunit,*) alpha
         read(iunit,*) m_crit
         read(iunit,*) c1
         read(iunit,*) m0
         read(iunit,*) dudx_eps
         read(iunit,*) phys_tol
         read(iunit,*) bed_normal
         close(iunit)

         open(unit=DIG_PARM_UNIT,file='fort.dig',status="unknown",action="write")

         write(DIG_PARM_UNIT,*) ' '
         write(DIG_PARM_UNIT,*) '--------------------------------------------'
         write(DIG_PARM_UNIT,*) 'SETDIG:'
         write(DIG_PARM_UNIT,*) '---------'
         write(DIG_PARM_UNIT,*) '    rho_s:',rho_s
         write(DIG_PARM_UNIT,*) '    rho_f:',rho_f
         write(DIG_PARM_UNIT,*) '    phi_bed:', phi_bed/deg2rad
         write(DIG_PARM_UNIT,*) '    phi_int:', phi_int/deg2rad
         write(DIG_PARM_UNIT,*) '    delta:', delta
         write(DIG_PARM_UNIT,*) '    kappita:', kappita
         write(DIG_PARM_UNIT,*) '    mu:', mu
         write(DIG_PARM_UNIT,*) '    alpha:', alpha
         write(DIG_PARM_UNIT,*) '    m_crit:', m_crit
         write(DIG_PARM_UNIT,*) '    c1:', c1
         write(DIG_PARM_UNIT,*) '    m0:', m0
         write(DIG_PARM_UNIT,*) '    dudx_eps:', dudx_eps
         write(DIG_PARM_UNIT,*) '    phys_tol:', phys_tol


   end subroutine set_dig

    ! ========================================================================
    !  set_pinit(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
   subroutine set_pinit(fname)

        implicit none

        ! Input
        character*25, intent(in), optional :: fname

        ! Locals
        integer, parameter :: iunit = 127
        character*25 :: file_name
        logical :: found_file


         ! Read user parameters from setgeo.data
         if (present(fname)) then
            file_name = fname
         else
            file_name = 'setpinit.data'
         endif
         inquire(file=file_name,exist=found_file)
         if (.not. found_file) then
            print *, 'You must provide a file ', file_name
            stop
         endif

         call opendatafile(iunit, file_name)
         read(iunit,*) init_ptype
         read(iunit,*) init_pmax_ratio
         read(iunit,*) init_ptf
         read(iunit,*) init_ptf2
         close(unit=iunit)

         p_initialized = 0
         init_pmin_ratio = 1.d16


         write(DIG_PARM_UNIT,*) ' '
         write(DIG_PARM_UNIT,*) '--------------------------------------------'
         write(DIG_PARM_UNIT,*) 'SETPINIT:'
         write(DIG_PARM_UNIT,*) '---------'
         write(DIG_PARM_UNIT,*) '    init_ptype:',init_ptype
         write(DIG_PARM_UNIT,*) '    init_pmax_ratio:',init_pmax_ratio
         write(DIG_PARM_UNIT,*) '    init_ptf:',init_ptf
         close(DIG_PARM_UNIT)



   end subroutine set_pinit


   !====================================================================
   !subroutine admissibleq
   !accept solution q, return q in admissible space
   !====================================================================

   subroutine admissibleq(h,hu,hv,hm,p,u,v,m,theta)

      implicit none

      !Input
      double precision :: h,hu,hv,hm,p,u,v,m,theta

      !Locals
      double precision :: mlo,mhi,hlo,pmax,phi,plo,rho,dry_tol,m_min,gmod

      gmod = grav
      dry_tol = drytolerance
      if (bed_normal.eq.1) gmod = grav*dcos(theta)

      if (h.le.dry_tol) then
         h =  0.0*max(h,0.d0)
         hu = 0.d0
         hv = 0.d0
         hm = h*m0
         p  = h*gmod*rho_f
         u = 0.d0
         v = 0.d0
         m = m0
         return
      endif

      u = hu/h
      v = hv/h
      m = hm/h

      !mlo = 1.d-3
      mlo = 1.d-16
      mhi = 1.d0 - mlo

      if (m.lt.mlo) then
         m = dmax1(m,mlo)
         !m = (m**2 + mlo**2)/(2.d0*mlo)
         hm = h*m
      elseif (m.gt.mhi) then
         m = dmin1(m,1.d0)
         !m = 1.d0 - ((1.d0-mhi)**2 + (1.d0-m)**2)/(2.d0*(1.d0-mhi))
         hm = h*m
      endif

      rho = rho_s*m + (1.d0-m)*rho_f
      pmax = rho*gmod*h
      plo = rho_f*dry_tol*gmod*dry_tol
      phi = pmax - plo
      if (p.lt.plo) then
         p = dmax1(0.d0,p)
         !p = (p**2 + plo**2)/(2.d0*plo)
      elseif (p.gt.phi) then
         p = dmin1(pmax,p)
         !p = pmax - ((pmax-p)**2+ (pmax-phi)**2)/(2.d0*(pmax-phi))
      endif

      return

   end subroutine admissibleq

   !====================================================================
   ! subroutine auxeval: evaluates the auxiliary variables as functions
   !                     of the solution vector q
   !====================================================================

   subroutine auxeval(h,u,v,m,p,phi_bed,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

      implicit none

      !i/o
      double precision, intent(in)  :: h,u,v,m,p,pm,phi_bed,theta
      double precision, intent(out) :: S,rho,tanpsi,D,tau,kappa
      double precision, intent(out) :: sigbed,kperm,compress

      !local
      double precision m_eqn,vnorm,gmod,sigbedc,sqrtarg


      gmod=grav
      if (bed_normal.eq.1) gmod=grav*dcos(theta)
      vnorm = dsqrt(u**2 + v**2)
      rho = rho_s*m + rho_f*(1.d0-m)
      sigbed = dmax1(0.d0,rho*gmod*h - p)
      sigbedc = (rho-rho_f)*gmod*h**2
      !sigbedc = sigbed
      if (sigbedc.gt.0.d0) then
         !Note: m_eqn = m_crit/(1+sqrt(N))
         !From Boyer et. al.
         !and N = mu*venorm/(h*sigbed)
         !S = mu*abs(vnorm)/(h*sigbed)
         S = mu*dabs(vnorm)/(h*sigbedc)
         m_eqn = dsqrt(h*sigbedc)*m_crit/(dsqrt(h*sigbedc)+ dsqrt(mu*dabs(vnorm)))
         m_eqn = dmax1(m_eqn,m_crit/(1.d0 + dsqrt(2.d0)))
      else
         S = 0.d0
         m_eqn= 0.d0
      endif
      tanpsi = c1*(m-m_eqn)
      kperm = (kappita**2*(1.d0-m)**3)/(180.d0*m**2)

      !kperm = kappita*exp(-(m-0.6)/(0.04))
      compress = alpha/(sigbed + 1.d5)
      if (p_initialized.eq.0.and.vnorm.le.0.d0) then
         D = 0.d0
         tanpsi = 0.d0
      elseif (h*mu.gt.0.d0) then
         D = (kperm/(mu*h))*(rho_f*gmod*h - p)
      else
         D = 0.d0
      endif

      tau = dmax1(0.d0,sigbed*dtan(phi_bed + datan(tanpsi)))

      !kappa: earth pressure coefficient
      !if (phi_int.eq.phi_bed) then
      !   sqrtarg = 0.d0
      !else
      !   sqrtarg = 1.d0-(dcos(phi_int)**2)*(1.d0 + dtan(phi_bed)**2)
      !endif

      !kappa = (2.d0 - pm*2.d0*dsqrt(sqrtarg))/(dcos(phi_int)**2)
      !kappa = kappa - 1.d0
      kappa = 1.d0

   end subroutine auxeval


   !====================================================================
   !subroutine psieval: evaluate the source term
   !====================================================================

   subroutine psieval(tau,rho,D,tanpsi,kperm,compress,h,u,m,psi)

      implicit none

      !i/o
      double precision, intent(out) :: psi(4)
      double precision, intent(in)  :: tau,rho,D,tanpsi,kperm,compress
      double precision, intent(in)  :: h,u,m

      !local
      double precision :: taushear,drytol,vnorm

      drytol = drytolerance

      taushear = (tau/rho)*dsign(1.d0,u)
      vnorm = dabs(u)
      if (h.lt.drytol.or..true.) then
         psi(1) = 0.d0
         psi(2) = 0.d0
         psi(3) = 0.d0
         psi(4) = 0.d0
      else
         psi(1) =  D*(rho-rho_f)/rho
         psi(2) =  u*D*(rho-rho_f)/rho
         psi(3) = -D*m*(rho_f/rho)
         psi(4) = 0.d0
      endif

   end subroutine psieval

   ! ========================================================================
   !  calc_pmin
   ! ========================================================================
   !  Determines minimum pore pressure for mobilization
   ! ========================================================================

   subroutine calc_pmin(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      !Input
      double precision :: dx,dy,xlower,ylower
      double precision :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      integer :: maxmx,maxmy,mx,my,mbc,meqn,maux

      !Locals
      double precision :: forcemag,pcrit,rho,h,h_r,h_l,b_r,b_l,dry_tol
      double precision :: tanpsi,gmod,thetaL,thetaR,phiL,phiR
      double precision :: forcemagL,forcemagR,pcritL,pcritR
      integer :: i,j

      gmod = grav

      !if (init_ptype.eq.2) then
      !   init_pmin_ratio = 1.0d0
      !   return
      !endif
      dry_tol = drytolerance
      rho = m0*rho_s + (1.d0-m0)*rho_f
      tanpsi = max(c1*(m0 - m_crit),0.d0)
      !tanpsi = c1*(m0-m_crit)

      do i=1,mx
         do j=1,my
            h_r = q(i,j,1)
            h_l = q(i-1,j,1)
            if (h_l.le.dry_tol.or.h_r.le.dry_tol) then
               cycle
            endif
            b_r = aux(i,j,1)
            b_l = aux(i-1,j,1)
            phiL = aux(i-1,j,i_phi)
            phiR = aux(i,j,i_phi)
            phiL = dmax1(0.d0,phiL + datan(tanpsi))
            phiR = dmax1(0.d0,phiR + datan(tanpsi))
            if (bed_normal.eq.1) then
               thetaR = aux(i,j,i_theta)
               thetaL = aux(i-1,j,i_theta)
               gmod = grav*dcos(0.5d0*(thetaL+thetaR))
            else
               thetaL = 0.d0
               thetaR = 0.d0
            endif

            if (h_l.gt.dry_tol.and.h_r.gt.dry_tol) then
               h = 0.5d0*(h_r + h_l)
            elseif (h_r.gt.dry_tol) then
               h = h_r
            else
               h = h_l
            endif

            !determine pressure min ratio
            forcemagL = abs(-grav*h*dsin(thetaL)*dx + gmod*(b_r-b_l)*h + 0.5d0*gmod*(h_r**2 - h_l**2))
            forcemagR = abs(-grav*h*dsin(thetaR)*dx + gmod*(b_r-b_l)*h + 0.5d0*gmod*(h_r**2 - h_l**2))
            pcritR = (rho*h_r*gmod - rho*forcemagR/(dx*tan(phiR)))/(rho_f*gmod*h_r)
            pcritL = (rho*h_l*gmod - rho*forcemagL/(dx*tan(phiR)))/(rho_f*gmod*h_l)
            pcrit = max(pcritR,pcritL)
            init_pmin_ratio = min(init_pmin_ratio,pcrit)
            init_pmin_ratio = max(init_pmin_ratio,0.d0)

            !repeat for y-Riemann problems
            h_r = q(i,j,1)
            h_l = q(i,j-1,1)
            if (h_l.le.dry_tol.or.h_r.le.dry_tol) then
               cycle
            endif
            b_r = aux(i,j,1)
            b_l = aux(i,j-1,1)
            phiL = aux(i,j-1,i_phi)
            phiR = aux(i,j,i_phi)
            phiL = dmax1(0.d0,phiL + datan(tanpsi))
            phiR = dmax1(0.d0,phiR + datan(tanpsi))

            if (bed_normal.eq.1) then
               thetaR = aux(i,j,i_theta)
               thetaL = aux(i,j-1,i_theta)
               gmod = grav*dcos(0.5d0*(thetaL+thetaR))
            endif

            if (h_l.gt.dry_tol.and.h_r.gt.dry_tol) then
               h = 0.5d0*(h_r + h_l)
            elseif (h_r.gt.dry_tol) then
               h = h_r
            else
               h = h_l
            endif

           !determine pressure min ratio
            forcemag = abs(gmod*(b_r-b_l)*h + 0.5d0*gmod*(h_r**2 - h_l**2))
            pcritR = (rho*h_r*gmod - rho*forcemag/(dy*tan(phiR)))/(rho_f*gmod*h_r)
            pcritL = (rho*h_l*gmod - rho*forcemag/(dy*tan(phiR)))/(rho_f*gmod*h_l)
            pcrit = max(pcritR,pcritL)
            init_pmin_ratio = min(init_pmin_ratio,pcrit)
            init_pmin_ratio = max(init_pmin_ratio,0.d0)


         enddo
      enddo
      write(*,*) 'init_pmin_ratio:',init_pmin_ratio
   end subroutine calc_pmin



end module digclaw_module
