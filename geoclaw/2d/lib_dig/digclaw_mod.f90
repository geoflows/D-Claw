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
    double precision :: rho_s,rho_f,phi_bed,phi_int,delta,kappita,mu,alpha
    double precision :: m_crit,m_min,m_eqn0,c1,m0,dudx_eps,phys_tol

    integer init_htype,init_ptype,p_initialized
    double precision init_pmax_ratio,init_ptf,init_pmin_ratio

    integer, parameter ::  i_S =4
    integer, parameter ::  i_rho = 5
    integer, parameter ::  i_tanpsi = 6
    integer, parameter ::  i_D = 7
    integer, parameter ::  i_tau = 8
    integer, parameter ::  i_kappa = 9
    integer, parameter ::  i_phi = 10
    integer, parameter ::  i_sigbed = 11
    integer, parameter ::  i_theta = 12
    integer, parameter ::  i_kperm = 13
    integer, parameter ::  i_alpha = 14
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


         open(unit=DIG_PARM_UNIT,file='fort.dig',status="unknown",action="write")

         write(DIG_PARM_UNIT,*) ' '
         write(DIG_PARM_UNIT,*) '--------------------------------------------'
         write(DIG_PARM_UNIT,*) 'SETDIG:'
         write(DIG_PARM_UNIT,*) '---------'

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
         read(iunit,*) m_min
         read(iunit,*) m_eqn0
         read(iunit,*) c1
         read(iunit,*) m0
         read(iunit,*) dudx_eps
         read(iunit,*) phys_tol
         close(iunit)

         write(DIG_PARM_UNIT,*) '    rho_s:',rho_s
         write(DIG_PARM_UNIT,*) '    rho_f:',rho_f
         write(DIG_PARM_UNIT,*) '    phi_bed:', phi_bed
         write(DIG_PARM_UNIT,*) '    phi_int:', phi_int
         write(DIG_PARM_UNIT,*) '    delta:', delta
         write(DIG_PARM_UNIT,*) '    kappita:', kappita
         write(DIG_PARM_UNIT,*) '    mu:', mu
         write(DIG_PARM_UNIT,*) '    alpha:', alpha
         write(DIG_PARM_UNIT,*) '    m_crit:', m_crit
         write(DIG_PARM_UNIT,*) '    m_min:', m_min
         write(DIG_PARM_UNIT,*) '    m_eqn0:', m_eqn0
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


         open(unit=DIG_PARM_UNIT,file='fort.geo',status="unknown",action="write")

         write(DIG_PARM_UNIT,*) ' '
         write(DIG_PARM_UNIT,*) '--------------------------------------------'
         write(DIG_PARM_UNIT,*) 'SETPINIT:'
         write(DIG_PARM_UNIT,*) '---------'


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
         read(iunit,*) init_htype
         read(iunit,*) init_ptype
         read(iunit,*) init_pmax_ratio
         read(iunit,*) init_ptf
         if (init_ptype.gt.0) p_initialized = 0
         close(unit=iunit)
         init_pmin_ratio = 1.d0


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

            pcrit = rho*h*grav - rho*forcemag/(dx*tan(phi))
            pcrit = max(pcrit,0.0)
            init_pmin_ratio = min(init_pmin_ratio,pcrit/(rho_f*grav*h))
            init_pmin_ratio = max(init_pmin_ratio,0.d0)

         enddo
      enddo



   end subroutine calc_pmin




end module digclaw_module
