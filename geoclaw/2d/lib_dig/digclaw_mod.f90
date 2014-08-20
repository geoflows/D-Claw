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
    double precision :: rho_s,rho_f,phi_bed,theta_input,delta,kappita
    double precision :: mu,alpha,m_crit,c1,m0,phys_tol,sigma_0

    integer :: init_ptype,p_initialized,bed_normal
    double precision :: init_pmax_ratio,init_ptf2,init_ptf,init_pmin_ratio
    double precision :: grad_eta_max,cohesion_max,grad_eta_ave,eta_cell_count

    integer, parameter ::  i_dig    = 4 !Start of digclaw aux variables
    integer, parameter ::  i_phi    = i_dig
    integer, parameter ::  i_theta  = i_dig + 1
    integer, parameter ::  i_fs     = i_dig + 2
    integer, parameter ::  i_cohesion  = i_dig + 3
    integer, parameter ::  i_taudir_x = i_dig + 4
    integer, parameter ::  i_taudir_y = i_dig + 5
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
         read(iunit,*) theta_input
         theta_input = deg2rad*theta_input
         read(iunit,*) delta
         read(iunit,*) kappita
         read(iunit,*) mu
         read(iunit,*) alpha
         read(iunit,*) m_crit
         read(iunit,*) c1
         read(iunit,*) m0
         read(iunit,*) sigma_0
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
         write(DIG_PARM_UNIT,*) '    theta_input:', theta_input/deg2rad
         write(DIG_PARM_UNIT,*) '    delta:', delta
         write(DIG_PARM_UNIT,*) '    kappita:', kappita
         write(DIG_PARM_UNIT,*) '    mu:', mu
         write(DIG_PARM_UNIT,*) '    alpha:', alpha
         write(DIG_PARM_UNIT,*) '    m_crit:', m_crit
         write(DIG_PARM_UNIT,*) '    c1:', c1
         write(DIG_PARM_UNIT,*) '    m0:', m0
         write(DIG_PARM_UNIT,*) '    sigma_0:', sigma_0
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
         grad_eta_max = 0.0
         cohesion_max = 0.0
         grad_eta_ave = 0.0
         eta_cell_count = 1.e-6


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

      !i/o
      double precision, intent(in) :: theta
      double precision, intent(inout) :: h,hu,hv,hm,p
      double precision, intent(out) :: u,v,m

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
         !p = dmax1(-5.0*pmax,p)
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
      double precision :: m_eqn,vnorm,gmod,sigbedc,hbounded,shear,tanphi

      if (h.lt.drytolerance) return

      hbounded = h !max(h,delta)
      gmod=grav
      if (bed_normal.eq.1) gmod=grav*dcos(theta)
      vnorm = dsqrt(u**2 + v**2)
      rho = rho_s*m + rho_f*(1.d0-m)
      shear = 2.0*vnorm/hbounded
      sigbed = dmax1(0.d0,rho*gmod*h - p)
      sigbedc = rho_s*(shear*delta)**2 + sigbed

      if (sigbedc.gt.0.0) then
         S = mu*shear/(sigbedc)
      else
         S = 0.d0
      endif
      !Note: m_eqn = m_crit/(1+sqrt(S))
      !From Boyer et. al
      m_eqn = m_crit/(1.d0 + sqrt(S))
      tanpsi = c1*(m-m_eqn)*tanh(shear/0.1)
      kperm = (1.0 + 0.0*vnorm/sqrt(gmod*h))*kappita*exp(-(m-0.60)/(0.04))
      !compress = alpha/(sigbed + 1.d5)
      compress = alpha/(m*(sigbed +  sigma_0))

      if (m.le.1.d-99) then
         kperm = 0.0
         tanpsi = 0.0
      endif

      if (p_initialized.eq.0.and.vnorm.le.0.d0) then
      !if (vnorm.le.0.d0) then
         tanpsi = 0.d0
         D = 0.d0
      elseif (h*mu.gt.0.d0) then
         D = 2.0*(kperm/(mu*h))*(rho_f*gmod*h - p)
      else
         D = 0.d0
      endif
      tanphi = dtan(phi_bed + datan(tanpsi))
      !if (S.gt.0.0) then
      !   tanphi = tanphi + 0.38*mu*shear/(shear + 0.005*sigbedc)
      !endif

      tau = dmax1(0.d0,sigbed*dtan(phi_bed + datan(tanpsi)))

      if (p_initialized.eq.0) then
         !tau = tau + rho*gmod*h*max(0.0,grad_eta_max*tan(phi_bed)-tan(phi_bed))
      !   tau = tau + pm
      endif
      !tau = (grav/gmod)*dmax1(0.d0,sigbed*tanphi)
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
   !  calc_taudir
   ! ========================================================================
   !  Determines the resistive force vector for static cells
   !  outputs direction cosines at each interface
   ! ========================================================================

subroutine calc_taudir(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      !Input
      double precision :: dx,dy,xlower,ylower
      double precision :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      integer :: maxmx,maxmy,mx,my,mbc,meqn,maux

      !Locals
      double precision :: hL,hR,huL,huR,hvL,hvR,bL,bR
      double precision :: phiL,phiR,phi,thetaL,thetaR,theta
      double precision :: gmod,dry_tol
      double precision :: EtaL,EtaR,EtaTL,EtaTR,EtaBL,EtaBR
      double precision :: detadx,detadxL,detadxR,detadxTR,detadxTL,detadxBR,detadxBL
      double precision :: detady,detadyL,detadyR,detadyTR,detadyTL,detadyBR,detadyBL
      double precision :: hTL,hTR,hBL,hBR,bTL,bTR,bBL,bBR

      integer :: i,j

      dry_tol = drytolerance
      gmod = grav


      do i=2-mbc,mx+mbc
         do j=2-mbc,my+mbc-1
            !note: for edge valued aux, aux(i,. _x) is at i-1/2.

            hR = q(i,j,1)
            hL = q(i-1,j,1)
            if ((hL<=dry_tol).and.(hR<=dry_tol)) then
               aux(i,j,i_taudir_x) = 1.0
               cycle
            endif

            huL = q(i-1,j,2)
            huR = q(i,j,2)
            hvL = q(i-1,j,3)
            hvR = q(i,j,3)

            if ((huL**2+huR**2)>0.0) then
               aux(i,j,i_taudir_x) = 1.0
               cycle
            endif

            bR = aux(i,j,1)
            bL = aux(i-1,j,1)
            phiL = aux(i-1,j,i_phi)
            phiR = aux(i,j,i_phi)

            hTL = q(i-1,j+1,1)
            bTL = aux(i-1,j+1,1)
            hBL = q(i-1,j-1,1)
            bBL = aux(i-1,j-1,1)

            hTR = q(i,j+1,1)
            bTR = aux(i,j+1,1)
            hBR = q(i,j-1,1)
            bBR = aux(i,j-1,1)

            if (bed_normal.eq.1) then
               thetaR = aux(i,j,i_theta)
               thetaL = aux(i-1,j,i_theta)
               gmod = grav*cos(0.5d0*(thetaL+thetaR))
            else
               thetaL = 0.d0
               thetaR = 0.d0
            endif
            phi   = 0.5*(phiL + phiR)
            theta = 0.5*(thetaL + thetaR)

            if ((hR>dry_tol).and.(hL>dry_tol))then
               EtaR = hR+bR
               EtaL = hL+bL
            elseif (hR>dry_tol) then
               EtaR = hR+bR
               EtaL = min(EtaR,hL+bL)
            else
               EtaL = hL+bL
               EtaR = min(EtaL,hR+bR)
            endif

            detadx = (EtaR-EtaL)/dx - sin(theta)

            !---------left minmod deta/dy-------------------
            if (hTL>dry_tol) then
               EtaTL = hTL+bTL
            else
               EtaTL = min(EtaL,hTL+bTL)
            endif
            detadyTL = (EtaTL-EtaL)/dy

            if (hBL>dry_tol) then
               EtaBL = hBL+bBL
            else
               EtaBL = min(EtaL,hBL+bBL)
            endif
            detadyBL = (EtaL-EtaBL)/dy

            if (detadyTL*detadyBL.gt.0.0) then
               detadyL = min(abs(detadyTL),abs(detadyBL))*sign(1.0,etaTL-etaBL)
            else
               detadyL = 0.0
            endif

            !-----------right minmod deta/dy------------------
            if (hTR>dry_tol) then
               EtaTR = hTR+bTR
            else
               EtaTR = min(EtaR,hTR+bTR)
            endif
            detadyTR = (EtaTR-EtaR)/dy

            if (hBR>dry_tol) then
               EtaBR = hBR+bBR
            else
               EtaBR = min(EtaR,hBR+bBR)
            endif
            detadyBR = (EtaR-EtaBR)/dy

            if (detadyTR*detadyBR.gt.0.0) then
               detadyR = min(abs(detadyTR),abs(detadyBR))*sign(1.0,etaTR-etaBR)
            else
               detadyR = 0.0
            endif

            !---------minmod deta/dy--------------------------
            if (detadyR*detadyL.gt.0.0) then
               detady = min(abs(detadyR),abs(detadyL))*sign(1.0,detadyR)
            else
               detady = 0.0
            endif

            !---------direction cosine for resistive static friction--
            !---------factor of safety for dry static material (p=0)--
            if (detady>0.0) then
               aux(i,j,i_taudir_x) = abs(detadx)/sqrt(detadx**2 + detady**2)
            elseif (detadx>0.0) then
               aux(i,j,i_taudir_x) = 1.0
            else
               aux(i,j,i_taudir_x) = 1.0
            endif
         enddo
         aux(i,my+mbc,i_taudir_x) = aux(i,my+mbc-1,i_taudir_x)
      enddo

      do j=2-mbc,my+mbc
         do i=2-mbc,mx+mbc-1
            !Note: edge value aux(.,j._y) is at j-1/2

            hR = q(i,j,1)
            hL = q(i,j-1,1)
            if ((hL<=dry_tol).and.(hR<=dry_tol)) then
               aux(i,j,i_taudir_y) = 1.0
               cycle
            endif

            huL = q(i,j-1,2)
            huR = q(i,j,2)
            hvL = q(i,j-1,3)
            hvR = q(i,j,3)

            if ((huL**2+huR**2)>0.0) then
               aux(i,j,i_taudir_y) = 1.0
               cycle
            endif

            bR = aux(i,j,1)
            bL = aux(i,j-1,1)
            phiL = aux(i,j-1,i_phi)
            phiR = aux(i,j,i_phi)

            hTL = q(i+1,j-1,1)
            bTL = aux(i+1,j-1,1)
            hBL = q(i-1,j-1,1)
            bBL = aux(i-1,j-1,1)

            hTR = q(i+1,j,1)
            bTR = aux(i+1,j,1)
            hBR = q(i-1,j,1)
            bBR = aux(i-1,j,1)

            if (bed_normal.eq.1) then
               thetaR = aux(i,j,i_theta)
               thetaL = aux(i,j-1,i_theta)
               gmod = grav*cos(0.5d0*(thetaL+thetaR))
            else
               thetaL = 0.d0
               thetaR = 0.d0
            endif
            phi   = 0.5*(phiL + phiR)
            theta = 0.5*(thetaL + thetaR)

            if ((hR>dry_tol).and.(hL>dry_tol))then
               EtaR = hR+bR
               EtaL = hL+bL
            elseif (hR>dry_tol) then
               EtaR = hR+bR
               EtaL = min(EtaR,hL+bL)
            else
               EtaL = hL+bL
               EtaR = min(EtaL,hR+bR)
            endif

            detady = (EtaR-EtaL)/dy

            !---------left minmod deta/dx-------------------
            if (hTL>dry_tol) then
               EtaTL = hTL+bTL
            else
               EtaTL = min(EtaL,hTL+bTL)
            endif
            detadxTL = (EtaTL-EtaL)/dx -sin(theta)

            if (hBL>dry_tol) then
               EtaBL = hBL+bBL
            else
               EtaBL = min(EtaL,hBL+bBL)
            endif
            detadxBL = (EtaL-EtaBL)/dx -sin(theta)

            if (detadxTL*detadxBL.gt.0.0) then
               detadxL = min(abs(detadxTL),abs(detadxBL))*sign(1.0,etaTL-etaBL)
            else
               detadxL = 0.0
            endif

            !-----------right minmod deta/dy------------------
            if (hTR>dry_tol) then
               EtaTR = hTR+bTR
            else
               EtaTR = min(EtaR,hTR+bTR)
            endif
            detadxTR = (EtaTR-EtaR)/dx -sin(theta)

            if (hBR>dry_tol) then
               EtaBR = hBR+bBR
            else
               EtaBR = min(EtaR,hBR+bBR)
            endif
            detadxBR = (EtaR-EtaBR)/dx - sin(theta)

            if (detadxTR*detadxBR.gt.0.0) then
               detadxR = min(abs(detadxTR),abs(detadxBR))*sign(1.0,etaTR-etaBR)
            else
               detadxR = 0.0
            endif

            !---------minmod deta/dy--------------------------
            if (detadxR*detadxL.gt.0.0) then
               detadx = min(abs(detadxR),abs(detadxL))*sign(1.0,detadxR)
            else
               detadx = 0.0
            endif

            !---------direction cosine for resistive static friction--
            !---------factor of safety for dry static material (p=0)--
            if (detadx>0) then
               aux(i,j,i_taudir_y) = abs(detady)/sqrt(detady**2 + detadx**2)
            endif

         enddo
         aux(mx+mbc,j,i_taudir_y) = aux(mx+mbc-1,j,i_taudir_y)
      enddo


end subroutine calc_taudir

   ! ========================================================================
   !  calc_pmin
   ! ========================================================================
   !  Determines minimum pore pressure for mobilization
   !  Determines factor of safety and cohesion for static states
   ! ========================================================================

subroutine calc_pmin(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)


      implicit none

      !Input
      double precision :: dx,dy,xlower,ylower
      double precision :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision :: aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      integer :: maxmx,maxmy,mx,my,mbc,meqn,maux

      !Locals
      double precision :: h,hL,hR,hu,hv,b,bL,bR,bT,bB,hT,hB
      double precision :: phi,theta,rho
      double precision :: gmod,dry_tol
      double precision :: EtaL,EtaR,EtaT,EtaB,Eta
      double precision :: detadx,detadxL,detadxR,detady,detadyT,detadyB
      double precision :: grad_eta,init_pmin_ratio_noc


      integer :: i,j

      dry_tol = drytolerance
      gmod = grav
      rho = m0*rho_s + (1.0-m0)*rho_f

      do i=1,mx
         do j=1,my

            h = q(i,j,1)
            hL = q(i-1,j,1)
            hR = q(i+1,j,1)
            if (h<dry_tol) then
               aux(i,j,i_fs) = 10.0
               aux(i,j,i_cohesion) = 0.0
               cycle
            endif

            hu = q(i,j,2)
            hv = q(i,j,3)

            if ((hu**2+hv**2)>0.0) then
               aux(i,j,i_fs) = 0.0
               aux(i,j,i_cohesion) = 0.0
               cycle
            endif

            b = aux(i,j,1)
            bR = aux(i+1,j,1)
            bL = aux(i-1,j,1)
            phi = aux(i,j,i_phi)

            if ((phi)==0.0) then
               aux(i,j,i_fs) = 0.0
               aux(i,j,i_cohesion) = 0.0
               init_pmin_ratio = 0.0
               cycle
            endif

            hT = q(i,j+1,1)
            bT = aux(i,j+1,1)
            hB = q(i,j-1,1)
            bB = aux(i,j-1,1)

            if (bed_normal.eq.1) then
               theta = aux(i,j,i_theta)
               gmod = grav*cos(theta)
            else
               theta = 0.d0
            endif

            Eta  = h+b
            !---------max deta/dx-------------------
            EtaR = hR+bR
            EtaL = hL+bL
            if (hR<=dry_tol) then
               EtaR = min(Eta,bR)
            endif
            if (hL<=dry_tol) then
               EtaL = min(Eta,bL)
            endif
            detadxR = (EtaR-Eta)/dx -sin(theta)
            detadxL = (Eta-EtaL)/dx -sin(theta)
            if (detadxR*detadxL<=0.0) then
               detadx = 0.0
            elseif (abs(detadxR)>abs(detadxL)) then
               detadx = detadxL
            else
               detadx = detadxR
            endif


            !---------max deta/dy-------------------
            EtaT = hT+bT
            EtaB = hB+bB
            if (hT<=dry_tol) then
               EtaT = min(Eta,bT)
            endif
            if (hB<=dry_tol) then
               EtaB = min(Eta,bB)
            endif
            detadyT = (EtaT-Eta)/dy
            detadyB = (Eta-EtaB)/dy
            if (detadyT*detadyB<=0.0) then
               detady = 0.0
            elseif (abs(detadyT)>abs(detadyB)) then
               detady = detadyB
            else
               detady = detadyT
            endif

            grad_eta = sqrt(detadx**2 + detady**2)
            grad_eta_ave = grad_eta_ave + grad_eta/tan(phi)
            eta_cell_count = eta_cell_count + 1.0

            grad_eta_max = max(grad_eta_max,grad_eta/tan(phi))


            !aux(i,j,i_cohesion) = max(0.0,rho*gmod*h*(grad_eta-tan(phi))+0.0*rho_f*gmod*h*tan(phi))

            aux(i,j,i_cohesion) = 1.0-grad_eta/tan(phi)
            init_pmin_ratio = min(init_pmin_ratio, 1.0-grad_eta/tan(phi))

            if (grad_eta>0.0) then
               aux(i,j,i_fs) = tan(phi)/grad_eta
            else
               aux(i,j,i_fs) = 10.0
            endif
         enddo
      enddo

      if (init_ptype==2.or.init_ptype==4) then
         init_pmin_ratio = 1.0-(grad_eta_ave/eta_cell_count)
      endif
      if (init_ptype>0) then
         write(*,*) '--------------------------------------------'
         write(*,*) 'hydrostatic liquefaction ratio:', rho_f/rho
         write(*,*) 'initiation liquefaction  ratio:',init_pmin_ratio, grad_eta_ave
         write(*,*) 'maximum surface slope angle:',180.*atan(tan(phi)*grad_eta_max)/3.14, grad_eta_max
         write(*,*) 'average failure liquefaction ratio:', 1.0-(grad_eta_ave/eta_cell_count) , eta_cell_count
         write(*,*) '--------------------------------------------'
      endif
   end subroutine calc_pmin


end module digclaw_module
