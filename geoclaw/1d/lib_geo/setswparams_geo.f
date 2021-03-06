      subroutine setswparams
      implicit double precision (a-h,o-z)
      character*25 fname
      logical foundFile

      common /geo/ grav,abs_tol
      common /sourceterms/ frictioncoeff,ifrictiontype

      fname = 'setswparams.data'
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname
        write(*,*) 'file can be generated by setrun.py'
        write(*,*) '>> python setrun.py'
        stop
      endif

      iunit = 8
      call opendatafile(iunit, fname)

      read(8,*) grav
      read(8,*) abs_tol
      read(8,*) ifrictiontype
      read(8,*) frictioncoeff

      close(unit=8)


      return
      end
