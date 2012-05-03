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

    integer init_ptype,p_initialized
    double precision init_pmax_ratio,init_ptf,init_pmin_ratio

    integer, parameter ::  i_dig    = 4 !Start of digclaw aux variables
    integer, parameter ::  i_phi    = i_dig
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
         close(unit=iunit)

         p_initialized = 0
         init_pmin_ratio = 1.d0


         write(DIG_PARM_UNIT,*) ' '
         write(DIG_PARM_UNIT,*) '--------------------------------------------'
         write(DIG_PARM_UNIT,*) 'SETPINIT:'
         write(DIG_PARM_UNIT,*) '---------'
         write(DIG_PARM_UNIT,*) '    init_ptype:',init_ptype
         write(DIG_PARM_UNIT,*) '    init_pmax_ratio:',init_pmax_ratio
         write(DIG_PARM_UNIT,*) '    init_ptf:',init_ptf
         close(DIG_PARM_UNIT)



   end subroutine set_pinit

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
      double precision forcemag,pcrit,rho,h,h_r,h_l,b_r,b_l,dry_tol,phi
      integer :: i,j

      dry_tol = drytolerance

      do i=2-mbc,mx+mbc
         do j=2-mbc,my+mbc
            h_r = q(i,j,1)
            h_l = q(i-1,j,1)
            b_r = aux(i,j,1)
            b_l = aux(i-1,j,1)
            phi = 0.5d0*(aux(i-1,j,i_phi)+ aux(i,j,i_phi))
            rho = m0*rho_s + (1.d0-m0)*rho_f

            if (h_l.le.dry_tol.and.h_r.le.dry_tol) then
               cycle
            endif

            if (h_l.gt.dry_tol.and.h_r.gt.dry_tol) then
               h = 0.5d0*(h_r + h_l)
            elseif (h_r.gt.dry_tol) then
               h = h_r
            else
               h = h_l
            endif

            !determine pressure min ratio
            forcemag = abs(grav*(b_r-b_l)*h -0.5d0*grav*(h_r**2 - h_l**2))

            pcrit = rho*h*grav - rho*forcemag/(dx*tan(phi))
            pcrit = max(pcrit,0.0)
            init_pmin_ratio = min(init_pmin_ratio,pcrit/(rho_f*grav*h))
            init_pmin_ratio = max(init_pmin_ratio,0.d0)

            !repeat for y-Riemann problems

            h_r = q(i,j,1)
            h_l = q(i,j-1,1)
            b_r = aux(i,j,1)
            b_l = aux(i,j-1,1)
            phi = 0.5d0*(aux(i,j-1,i_phi)+ aux(i,j,i_phi))
            rho = m0*rho_s + (1.d0-m0)*rho_f

            if (h_l.le.dry_tol.and.h_r.le.dry_tol) then
               cycle
            endif

            if (h_l.gt.dry_tol.and.h_r.gt.dry_tol) then
               h = 0.5d0*(h_r + h_l)
            elseif (h_r.gt.dry_tol) then
               h = h_r
            else
               h = h_l
            endif

           !determine pressure min ratio
            forcemag = abs(grav*(b_r-b_l)*h -0.5d0*grav*(h_r**2 - h_l**2))

            pcrit = rho*h*grav - rho*forcemag/(dy*tan(phi))
            pcrit = max(pcrit,0.0)
            init_pmin_ratio = min(init_pmin_ratio,pcrit/(rho_f*grav*h))
            init_pmin_ratio = max(init_pmin_ratio,0.d0)

         enddo
      enddo

   end subroutine calc_pmin


   !====================================================================
   !subroutine admissibleq
   !accept solution q, return q in admissible space
   !====================================================================

   subroutine admissibleq(h,hu,hv,hm,p,u,v,m)

      implicit none

      !Input
      double precision :: h,hu,hv,hm,p,u,v,m

      !Locals
      double precision :: mlo,mhi,hlo,pmax,phi,plo,rho,dry_tol,m_min

      dry_tol = drytolerance

      if (h.le.dry_tol) then
         h = 0.d0
         hu = 0.d0
         hv = 0.d0
         hm = 0.d0
         p  = 0.d0
         u = 0.d0
         v = 0.d0
         m = m0
         return
      endif

      hlo = dry_tol

      if (h.lt.hlo) then
         u = 0.d0!hu/h
         v = 0.d0!hv/h
         m = m0!hm/m
         h = (h**2 + hlo**2)/(2.d0*hlo)
         p  = h*rho_f*grav
         hu = 0.d0!h*u
         hv = 0.d0!h*v
         hm = 0.d0!h*m
         return
      endif

      u = hu/h
      v = hv/h
      m = hm/h
      m = min(m,1.d0)
      m = max(m,0.d0)

      mlo = 1.d-3
      mhi = 1.d0 - mlo

      if (m.le.mlo) then
         m = (m**2 + mlo**2)/(2.d0*mlo)
         hm = h*m
      elseif (m.ge.mhi) then
         m = 1.d0 - ((1.d0-mhi)**2 + (1.d0-m)**2)/(2.d0*(1.d0-mhi))
         hm = h*m
      endif

      rho = (rho_s*m + (1.d0-m)*rho_f)
      pmax = rho*grav*h
      p = min(pmax,p)
      p = max(0.d0,p)
      plo = rho_f*dry_tol*grav
      phi = pmax - plo
      if (p.lt.plo) then
         p = (p**2 + plo**2)/(2.d0*plo)
      elseif (p.gt.phi) then
         p = pmax - ((pmax-p)**2+ (pmax-phi)**2)/(2.d0*(pmax-phi))
      endif

      return

   end subroutine admissibleq

   !====================================================================
   ! subroutine auxeval: evaluates the auxiliary variables as functions
   !                     of the solution vector q
   !====================================================================

   subroutine auxeval(h,u,v,m,p,phi_bed,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

      implicit none

      !i/o
      double precision, intent(in)  :: h,u,v,m,p,pm,phi_bed
      double precision, intent(out) :: S,rho,tanpsi,D,tau,kappa
      double precision, intent(out) :: sigbed,kperm,compress

      !local
      double precision m_eqn,vnorm,g,sigbedc,sqrtarg


      g=grav
      vnorm = dsqrt(u**2 + v**2)
      rho = rho_s*m + rho_f*(1.d0-m)
      sigbed = max(0.d0,rho*g*h - p)
      sigbedc = (rho-rho_f)*grav*h**2
      !sigbedc = sigbed
      if (sigbedc.gt.0.d0) then
         !Note: m_eqn = m_crit/(1+sqrt(N))
         !From Boyer et. al.
         !and N = mu*venorm/(h*sigbed)
         !S = mu*abs(vnorm)/(h*sigbed)
         S = mu*dabs(vnorm)/(h*sigbedc)
         m_eqn = dsqrt(h*sigbedc)*m_crit/(dsqrt(h*sigbedc)+ dsqrt(mu*dabs(vnorm)))
         m_eqn = dmax1(m_eqn,0.5*m_crit)
      else
         S = 0.d0
         m_eqn=0.d0
      endif
      tanpsi = c1*(m-m_eqn)
      tau = max(0.d0,sigbed*tan(phi_bed + atan(tanpsi)))
      kperm = (kappita**2*(1.d0-m)**3)/(180.d0*m**2)
      !kperm = kappita**2*exp(max(0.d0,m-m_crit)/(-0.03))/40.0
      compress = alpha/((m)*(sigbed + 1.d5))
      if (p_initialized.gt.0.and.h*mu.gt.0.d0) then
         D = (kperm/(h*mu))*(rho_f*g*h - p)
      else
         D = 0.d0
      endif
      !kappa: earth pressure coefficient
      !if (phi_int.eq.phi_bed) then
      !   sqrtarg = 0.d0
      !else
      !   sqrtarg = 1.d0-(cos(phi_int)**2)*(1.d0 + tan(phi_bed)**2)
      !endif

      !kappa = (2.d0 - pm*2.d0*sqrt(sqrtarg))/(cos(phi_int)**2)
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
      double precision :: taushear,drytol

      drytol = drytolerance

      taushear = (tau/rho)*dsign(1.d0,u)

      if (h.gt.drytol.and.p_initialized.eq.1) then
         psi(1) =  D*(rho-rho_f)/rho
         psi(2) =  u*D*(rho-rho_f)/rho
         psi(3) = -D*m*(rho_f/rho)
         psi(4) = 0.d0
      else
         psi(1) = 0.d0
         psi(2) = 0.d0
         psi(3) = 0.d0
         psi(4) = 0.d0
      endif

   end subroutine psieval

   !====================================================================
   ! subroutine calcaux: calculate all auxiliary variables
   !====================================================================

   subroutine calcaux(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      integer, parameter ::  i_rho    = i_dig+1
      integer, parameter ::  i_tanpsi = i_dig+2
      integer, parameter ::  i_D      = i_dig+3
      integer, parameter ::  i_tau    = i_dig+4
      integer, parameter ::  i_kappa  = i_dig+5
      integer, parameter ::  i_S      = i_dig+6
      integer, parameter ::  i_sigbed = i_dig+7
      integer, parameter ::  i_theta  = i_dig+8
      integer, parameter ::  i_kperm  = i_dig+9
      integer, parameter ::  i_alpha  = i_dig+10

      !i/o
      double precision :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      double precision, intent(in) :: xlower,ylower,dx,dy
      integer, intent(in) :: maxmx,maxmy,mx,my,meqn,mbc,maux

      !local
      double precision :: h,u,v,m,pbed,rho,S,m_eqn,tanpsi,phi,compress,pms
      double precision :: kperm,tau,up,um,dudx,dudy,div,pm,kappa,D,theta,g
      double precision :: sigbed,dry_tol
      integer :: i,j,ma,ibc

      g= grav
      dry_tol = drytolerance

      do i=1,mx
         do j=1,my
            phi = aux(i,j,i_phi)
            call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),q(i,j,5),u,v,m)
            h=q(i,j,1)
            if (h.le.dry_tol) then
               aux(i,j,i_S) = 0.d0
               aux(i,j,i_rho) = 0.d0
               aux(i,j,i_tanpsi) = 0.d0
               aux(i,j,i_D) = 0.d0
               aux(i,j,i_tau) = 0.d0
               aux(i,j,i_kappa) = 1.d0
               aux(i,j,i_sigbed) = 0.d0
               aux(i,j,i_kperm) = 0.d0
               aux(i,j,i_alpha) = 0.d0
               cycle
            endif
            pbed = q(i,j,5)

            if (q(i,j,1).gt.0.d0.and.q(i-1,j,1).gt.0.d0) then
               dudx = (q(i,j,2)/q(i,j,1) - q(i-1,j,2)/q(i-1,j,1))/dx
            else
               dudx = 0.d0
            endif
            if (q(i,j,1).gt.0.d0.and.q(i-1,j-1,1).gt.0.d0) then
               dudy = (q(i,j,2)/q(i,j,1) - q(i-1,j,2)/q(i-1,j,1))/dy
            else
               dudy = 0.d0
            endif
            div = dudx + dudy
            pm = sign(1.d0,div)

            call auxeval(h,u,v,m,pbed,kappa,S,rho,tanpsi,D,tau,phi,sigbed,kperm,compress,pm)

            aux(i,j,i_S) = S
            aux(i,j,i_rho) = rho
            aux(i,j,i_tanpsi) = tanpsi
            aux(i,j,i_D) = D
            aux(i,j,i_tau) = tau
            aux(i,j,i_kappa) = kappa
            aux(i,j,i_sigbed) = sigbed
            aux(i,j,i_kperm) = kperm
            aux(i,j,i_alpha) = compress

         enddo
      enddo

      !right and left boundary
      do j=1,my
         do ma = i_dig,maux
            do ibc=1,mbc
               aux(1-ibc,j,ma) = aux(1,j,ma)
               aux(mx+ibc,j,ma) = aux(mx,j,ma)
            enddo
         enddo
      enddo
      !top and bottom boundary
      do i=1,mx
         do ma = i_dig,maux
            do ibc=1,mbc
               aux(i,1-ibc,ma) = aux(i,1,ma)
               aux(i,my+ibc,ma) = aux(i,my,ma)
            enddo
         enddo
      enddo

   end subroutine calcaux


end module digclaw_module
