! ============================================================================
!  Program:     auxinit_module
!  File:        auxinit_mod.f90
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
module auxinit_module

      implicit none

      ! Work array
      double precision, allocatable :: auxinitwork(:)

      ! file data
      character*150, allocatable :: auxinitfname(:)
      integer :: mauxinitfiles,mauxinitsize
      double precision, allocatable :: xlowauxinit(:), ylowauxinit(:)
      double precision, allocatable :: xhiauxinit(:), yhiauxinit(:)
      double precision, allocatable :: dxauxinit(:), dyauxinit(:)

      integer, allocatable ::  mxauxinit(:), myauxinit(:)
      integer, allocatable :: i0auxinit(:), mauxinit(:)
      integer, allocatable :: iauxinit(:), auxinitftype(:)
      integer, allocatable ::  minlevelauxinit(:), maxlevelauxinit(:)




contains


    ! ========================================================================
    ! Read auxinit files as specified in auxinit.data
    !
    ! Each file has a type stored in auxinitftype(i).
    !   auxinitftype = 1:  standard GIS format: 3 columns: x,y,z(m)
    !   auxinitftype = 2:  Header as in DEM file, height(m) one value per line
    !   auxinitftype = 3:  Header as in DEM file, height(m) one row per line
    ! For other formats modify readauxinit routine.
    !
    ! advancing northwest to northeast then from north to south. Values should
    ! be uniformly spaced.
    !
    ! Associated with each file is a initialization type, iauxinit(file):
    !     as follows:
    !     defines a perturbation (or definition of) aux(i,j,iauxinit)
    ! ========================================================================

   subroutine set_auxinit(fname)

      use geoclaw_module

      implicit none

      ! Input arguments
      character*25, intent(in), optional :: fname

      ! Locals
      integer, parameter :: iunit = 7
      integer :: i,j,iauxinitfile
      character*25 :: file_name
      logical :: found_file


      ! Open and begin parameter file output
      write(GEO_PARM_UNIT,*) ' '
      write(GEO_PARM_UNIT,*) '--------------------------------------------'
      write(GEO_PARM_UNIT,*) 'SETAUXINIT:'
      write(GEO_PARM_UNIT,*) '---------'

      if (present(fname)) then
         file_name = fname
      else
         file_name  = 'setauxinit.data'
      endif
      inquire(file=file_name,exist=found_file)
      if (.not. found_file) then
         print *, 'You must provide a file ', file_name
         stop
      endif

      call opendatafile(iunit, file_name)

      read(iunit,*) mauxinitfiles

      if (mauxinitfiles==0) then
         write(GEO_PARM_UNIT,*) '   mauxinitfiles = 0'
         write(GEO_PARM_UNIT,*) '   no initial perturbation = 0'
         write(GEO_PARM_UNIT,*) '   h will be set max(0-b,0)   '
         return
      endif

      write(GEO_PARM_UNIT,*) '   mauxinitfiles = ',mauxinitfiles

      ! Read and allocate data parameters for each file
      allocate(mxauxinit(mauxinitfiles),myauxinit(mauxinitfiles))
      allocate(xlowauxinit(mauxinitfiles),ylowauxinit(mauxinitfiles))
      allocate(xhiauxinit(mauxinitfiles),yhiauxinit(mauxinitfiles))
      allocate(dxauxinit(mauxinitfiles),dyauxinit(mauxinitfiles))
      allocate(minlevelauxinit(mauxinitfiles),maxlevelauxinit(mauxinitfiles))
      allocate(auxinitfname(mauxinitfiles),auxinitftype(mauxinitfiles))
      allocate(iauxinit(mauxinitfiles))
      allocate(i0auxinit(mauxinitfiles),mauxinit(mauxinitfiles))

      do i=1,mauxinitfiles
         read(iunit,*) auxinitfname(i)
         read(iunit,*) auxinitftype(i),iauxinit(i),minlevelauxinit(i), maxlevelauxinit(i)

         write(GEO_PARM_UNIT,*) '   '
         write(GEO_PARM_UNIT,*) '   ',auxinitfname(i)
         write(GEO_PARM_UNIT,*) '  auxinitftype = ', auxinitftype(i)
         write(GEO_PARM_UNIT,*) '  iauxinit = ', iauxinit(i)
         write(GEO_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                                  minlevelauxinit(i), maxlevelauxinit(i)

         call read_auxinit_header(auxinitfname(i),auxinitftype(i),mxauxinit(i), &
                myauxinit(i),xlowauxinit(i),ylowauxinit(i),xhiauxinit(i),yhiauxinit(i), &
                dxauxinit(i),dyauxinit(i))
            mauxinit(i) = mxauxinit(i)*myauxinit(i)
      enddo

      ! Indexing into work array
      i0auxinit(1)=1
      if (mauxinitfiles > 1) then
         do i=2,mauxinitfiles
            i0auxinit(i)=i0auxinit(i-1) + mauxinit(i-1)
         enddo
      endif

      ! Read and allocate space in work array for each file
      mauxinitsize = sum(mauxinit)
      allocate(auxinitwork(mauxinitsize))

      do i=1,mauxinitfiles
            call read_auxinit(mxauxinit(i),myauxinit(i),auxinitftype(i),auxinitfname(i), &
                auxinitwork(i0auxinit(i):i0auxinit(i)+mauxinit(i)-1))
      enddo

   end subroutine set_auxinit


    ! ========================================================================
    !
    !  Read auxinit file.
    !
    ! ========================================================================
    subroutine read_auxinit(mx,my,filetype,fname,auxinit)

        use geoclaw_module

        implicit none

        ! Arguments
        integer, intent(in) :: mx,my,filetype
        character*150, intent(in) :: fname
        double precision, intent(inout) :: auxinit(1:mx*my)

        ! Locals
        integer, parameter :: iunit = 19, miss_unit = 17
        double precision, parameter :: auxinit_missing = -150.d0
        logical, parameter :: maketype2 = .false.
        integer :: i,j,num_points,missing,status,auxinit_start
        double precision :: no_data_value,x,y,z

        print *, ' '
        print *, 'Reading auxinit file  ', fname

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
                    read(iunit,fmt=*,iostat=status) x,y,auxinit(i)
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
                            read(iunit,*) auxinit(i)
                            if (auxinit(i) == no_data_value) then
                                missing = missing + 1
                                auxinit(i) = auxinit_missing
                            endif
                        enddo
                    case(3)
                        do j=1,my
                            read(iunit,*) (auxinit((j-1)*mx + i),i=1,mx)
                            do i=1,mx
                                if (auxinit((j-1)*mx + i) == no_data_value) then
                                    missing = missing + 1
                                    auxinit((j-1)*mx + i) = auxinit_missing
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
                        auxinit_missing
                endif
        end select

        close(unit=iunit)

   end subroutine read_auxinit


    ! ========================================================================
    ! subroutine read_auxinit_header(fname,auxinit_type,mx,my,xll,yll,xhi,yhi,dx,dy)
    ! ========================================================================
    !  Read auxinit file header to determine space needed in allocatable array
    !
    !  :Input:
    !   - fname - (char) Name of file
    !   - auxinitftype - (int) Type of file format (1 < auxinitftype < 3)
    !
    !  :Output:
    !   - mx,my - (int) Number of grid points
    !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
    !   - dx,dy - (float) Spatial resolution of grid
    ! ========================================================================
    subroutine read_auxinit_header(fname,auxinit_type,mx,my,xll,yll,xhi,yhi,dx,dy)

        use geoclaw_module

        implicit none

        ! Input and Output
        character*150, intent(in) :: fname
        integer, intent(in) :: auxinit_type
        integer, intent(out) :: mx,my
        double precision, intent(out) :: xll,yll,xhi,yhi,dx,dy

        ! Local
        integer, parameter :: iunit = 19
        integer :: auxinit_size, status
        double precision :: x,y,z,nodata_value
        logical :: found_file

        inquire(file=fname,exist=found_file)
        if (.not. found_file) then
            print *, 'Missing auxinit file:'
            print *, '   ', fname
            stop
        endif

        open(unit=iunit, file=fname, status='unknown',form='formatted')

        select case(abs(auxinit_type))
            ! ASCII file with 3 columns
            ! determine data size
            case(1)
                ! Initial size variables
                auxinit_size = 0
                mx = 0

                ! Read in first values, determines xlow and yhi
                read(iunit,*) xll,yhi
                auxinit_size = auxinit_size + 1
                mx = mx + 1

                ! Go through first row figuring out mx, continue to count
                y = yhi
                do while (yhi == y)
                    read(iunit,*) x,y,z
                    auxinit_size = auxinit_size + 1
                    mx = mx + 1
                enddo
                mx = mx - 1
                ! Continue to count the rest of the lines
                status = 0
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) x,y,z
                    auxinit_size = auxinit_size + 1
                enddo
                if (status > 0) then
                    print *, "IO error occured in ",fname,", aborting!"
                    stop
                endif

                ! Calculate remaining values
                my = auxinit_size / mx
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
                print *, 'ERROR:  Unrecognized auxinit_type'
                print *, '    auxinit_file_type = ',auxinit_type
                print *, '  for auxinit file:'
                print *, '   ', fname
                stop
        end select

        close(iunit)

        write(GEO_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
        write(GEO_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
        write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

    end subroutine read_auxinit_header

end module auxinit_module
