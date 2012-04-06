! ============================================================================
!  Program:     qinit_module
!  File:        qinit_mod.f90
!  Created:     2012-04-05
!  Author:      David L George
! ============================================================================
!      Copyright (C) 2012-04-05 David George <dgeorge@uw.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================
!  Module for initialization data that might come from files
! ============================================================================
module qinit_module

      implicit none

      ! Work array
      double precision, allocatable :: qinitwork(:)

      ! Topography file data
      character*150, allocatable :: qinitfname(:)
      integer :: mqinitfiles,mqinitsize
      double precision, allocatable :: xlowqinit(:), ylowqinit(:)
      double precision, allocatable :: xhiqinit(:), yhiqinit(:)
      double precision, allocatable :: dxqinit(:), dyqinit(:)

      integer, allocatable ::  mxqinit(:), myqinit(:)
      integer, allocatable :: i0qinit(:), mqinit(:)
      integer, allocatable :: iqinit(:), qinitftype(:)
      integer, allocatable ::  minlevelqinit(:), maxlevelqinit(:)




contains


    ! ========================================================================
    ! Read qinit files as specified in setqinit.data
    !
    ! Each file has a type stored in qinitftype(i).
    !   qinittype = 1:  standard GIS format: 3 columns: x,y,z(m)
    !   qinittype = 2:  Header as in DEM file, height(m) one value per line
    !   qinittype = 3:  Header as in DEM file, height(m) one row per line
    ! For other formats modify readqinit routine.
    !
    ! advancing northwest to northeast then from north to south. Values should
    ! be uniformly spaced.
    !
    ! Associated with each file is a initialization type, iqinit(file):
    !     as follows:
    !     1,2,3: perturbation to q(1,2,3)
    !     4:     file defines eta, and h =q(i,j,1)= max(eta-b,0)
    ! ========================================================================

   subroutine set_qinit(fname)

      use geoclaw_module

      implicit none

      ! Input arguments
      character*25, intent(in), optional :: fname

      ! Locals
      integer, parameter :: iunit = 7
      integer :: i,j,iqinitfile
      character*25 :: file_name
      logical :: found_file


      ! Open and begin parameter file output
      write(GEO_PARM_UNIT,*) ' '
      write(GEO_PARM_UNIT,*) '--------------------------------------------'
      write(GEO_PARM_UNIT,*) 'SETQINIT:'
      write(GEO_PARM_UNIT,*) '---------'

      if (present(fname)) then
         file_name = fname
      else
         file_name  = 'setqinit.data'
      endif
      inquire(file=file_name,exist=found_file)
      if (.not. found_file) then
         print *, 'You must provide a file ', file_name
         stop
      endif

      call opendatafile(iunit, file_name)

      read(iunit,*) mqinitfiles

      if (mqinitfiles==0) then
         write(GEO_PARM_UNIT,*) '   mqinitfiles = 0'
         write(GEO_PARM_UNIT,*) '   no initial perturbation = 0'
         write(GEO_PARM_UNIT,*) '   h will be set max(0-b,0)   '
         return
      endif

      write(GEO_PARM_UNIT,*) '   mqinitfiles = ',mqinitfiles

      ! Read and allocate data parameters for each file
      allocate(mxqinit(mqinitfiles),myqinit(mqinitfiles))
      allocate(xlowqinit(mqinitfiles),ylowqinit(mqinitfiles))
      allocate(xhiqinit(mqinitfiles),yhiqinit(mqinitfiles))
      allocate(dxqinit(mqinitfiles),dyqinit(mqinitfiles))
      allocate(minlevelqinit(mqinitfiles),maxlevelqinit(mqinitfiles))
      allocate(qinitfname(mqinitfiles),qinitftype(mqinitfiles))
      allocate(iqinit(mqinitfiles))
      allocate(i0qinit(mqinitfiles),mqinit(mqinitfiles))

      do i=1,mqinitfiles
         read(iunit,*) qinitfname(i)
         read(iunit,*) qinitftype(i),iqinit(i),minlevelqinit(i), maxlevelqinit(i)

         write(GEO_PARM_UNIT,*) '   '
         write(GEO_PARM_UNIT,*) '   ',qinitfname(i)
         write(GEO_PARM_UNIT,*) '  qinitftype = ', qinitftype(i)
         write(GEO_PARM_UNIT,*) '  iqinit = ', iqinit(i)
         write(GEO_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                                  minlevelqinit(i), maxlevelqinit(i)

         call read_qinit_header(qinitfname(i),qinitftype(i),mxqinit(i), &
                myqinit(i),xlowqinit(i),ylowqinit(i),xhiqinit(i),yhiqinit(i), &
                dxqinit(i),dyqinit(i))
            mqinit(i) = mxqinit(i)*myqinit(i)
      enddo

      ! Indexing into work array
      i0qinit(1)=1
      if (mqinitfiles > 1) then
         do i=2,mqinitfiles
            i0qinit(i)=i0qinit(i-1) + mqinit(i-1)
         enddo
      endif

      ! Read and allocate space in work array for each file
      mqinitsize = sum(mqinit)
      allocate(qinitwork(mqinitsize))

      do i=1,mqinitfiles
            call read_qinit(mxqinit(i),myqinit(i),qinitftype(i),qinitfname(i), &
                qinitwork(i0qinit(i):i0qinit(i)+mqinit(i)-1))
      enddo

   end subroutine set_qinit


    ! ========================================================================
    !
    !  Read qinit file.
    !
    ! ========================================================================
    subroutine read_qinit(mx,my,filetype,fname,qinit)

        use geoclaw_module

        implicit none

        ! Arguments
        integer, intent(in) :: mx,my,filetype
        character*150, intent(in) :: fname
        double precision, intent(inout) :: qinit(1:mx*my)

        ! Locals
        integer, parameter :: iunit = 19, miss_unit = 17
        double precision, parameter :: qinit_missing = -150.d0
        logical, parameter :: maketype2 = .false.
        integer :: i,j,num_points,missing,status,qinit_start
        double precision :: no_data_value,x,y,z

        print *, ' '
        print *, 'Reading qinit file  ', fname

        open(unit=iunit, file=fname, status='unknown',form='formatted')

        select case(abs(filetype))
            ! ASCII file with x,y,z values on each line.
            ! (progressing from upper left corner across rows, then down)
            ! Assumes a uniform rectangular grid of data values.
            case(1)
                i = 0
                status = 0
                do while (status == 0)
                    i = i + 1
                    read(iunit,fmt=*,iostat=status) x,y,qinit(i)
                enddo

            ! ================================================================
            ! ASCII file with header followed by z data
            ! (progressing from upper left corner across rows, then down)
            ! one value per line if filetype=2 or
            ! mx values per line if filetype=3
            ! ================================================================
            case(2:3)
                ! Read header
                do i=1,5
                    read(iunit,*)
                enddo
                read(iunit,*) no_data_value

                ! Read in data
                missing = 0
                select case(abs(filetype))
                    case(2)
                        do i=1,mx*my
                            read(iunit,*) qinit(i)
                            if (qinit(i) == no_data_value) then
                                missing = missing + 1
                                qinit(i) = qinit_missing
                            endif
                        enddo
                    case(3)
                        do j=1,my
                            read(iunit,*) (qinit((j-1)*mx + i),i=1,mx)
                            do i=1,mx
                                if (qinit((j-1)*mx + i) == no_data_value) then
                                    missing = missing + 1
                                    qinit((j-1)*mx + i) = qinit_missing
                                endif
                            enddo
                        enddo
                end select

                ! Write a warning if we found and missing values
                if (missing > 0)  then
                    print *, '   WARNING...some missing data values this file'
                    print *, '       ',missing,' missing data values'
                    print *, '              (see fort.missing)'
                    print *, '   These values have arbitrarily been set to ',&
                        qinit_missing
                endif
        end select

        close(unit=iunit)

   end subroutine read_qinit


    ! ========================================================================
    ! subroutine read_qinit_header(fname,qinit_type,mx,my,xll,yll,xhi,yhi,dx,dy)
    ! ========================================================================
    !  Read qinit file header to determine space needed in allocatable array
    !
    !  :Input:
    !   - fname - (char) Name of file
    !   - qinit_type - (int) Type of file format (1 < qinit_type < 3)
    !
    !  :Output:
    !   - mx,my - (int) Number of grid points
    !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
    !   - dx,dy - (float) Spatial resolution of grid
    ! ========================================================================
    subroutine read_qinit_header(fname,qinit_type,mx,my,xll,yll,xhi,yhi,dx,dy)

        use geoclaw_module

        implicit none

        ! Input and Output
        character*150, intent(in) :: fname
        integer, intent(in) :: qinit_type
        integer, intent(out) :: mx,my
        double precision, intent(out) :: xll,yll,xhi,yhi,dx,dy

        ! Local
        integer, parameter :: iunit = 19
        integer :: qinit_size, status
        double precision :: x,y,z,nodata_value
        logical :: found_file

        inquire(file=fname,exist=found_file)
        if (.not. found_file) then
            print *, 'Missing qinitgraphy file:'
            print *, '   ', fname
            stop
        endif

        open(unit=iunit, file=fname, status='unknown',form='formatted')

        select case(abs(qinit_type))
            ! ASCII file with 3 columns
            ! determine data size
            case(1)
                ! Initial size variables
                qinit_size = 0
                mx = 0

                ! Read in first values, determines xlow and yhi
                read(iunit,*) xll,yhi
                qinit_size = qinit_size + 1
                mx = mx + 1

                ! Go through first row figuring out mx, continue to count
                y = yhi
                do while (yhi == y)
                    read(iunit,*) x,y,z
                    qinit_size = qinit_size + 1
                    mx = mx + 1
                enddo
                mx = mx - 1
                ! Continue to count the rest of the lines
                status = 0
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) x,y,z
                    qinit_size = qinit_size + 1
                enddo
                if (status > 0) then
                    print *, "IO error occured in ",fname,", aborting!"
                    stop
                endif

                ! Calculate remaining values
                my = qinit_size / mx
                xhi = x
                yll = y
                dx = (xhi-xll) / (mx-1)
                dy = (yhi-yll) / (my-1)

            ! ASCII file with header followed by z data
            case(2:3)
                read(iunit,*) mx
                read(iunit,*) my
                read(iunit,*) xll
                read(iunit,*) yll
                read(iunit,*) dx
                read(iunit,*) nodata_value
                dy = dx
                xhi = xll + (mx-1)*dx
                yhi = yll + (my-1)*dy

            case default
                print *, 'ERROR:  Unrecognized qinit_type'
                print *, '    qinit_file_type = ',qinit_type
                print *, '  for qinit file:'
                print *, '   ', fname
                stop
        end select

        close(iunit)

        write(GEO_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
        write(GEO_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
        write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

    end subroutine read_qinit_header

end module qinit_module
