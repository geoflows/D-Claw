
      subroutine opendatafile(iunit, fname)
c     #
c     # Open the file fname and determine how many leading lines are
c     # comments.  Then rewind and skip over the comment lines so that the
c     # file is ready for reading data from.
c     #
c     # All comment lines must start with # in the first column.

      integer iunit, line, commentlines
      character*25 fname
      character*12 fname12
      character*1 firstchar
      logical foundFile
      
c     # Debug statement to test fname:
c     write(6,*) '+++ fname = XXX',fname,'XXX'

      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
c         # for backward compatability: maybe fname is declared character*12
c         # in calling routine...
          fname12 = fname(1:12)
c         #write(6,*) 'truncated fname12 = XXX',fname12,'XXX'
          inquire(file=fname12,exist=foundFile)
          if (.not. foundFile) then
            write(*,*) '*** in opendatafile, file not found:', fname 
            stop
          endif
          open(unit=iunit,file=fname12,status='old',form='formatted')
          write(6,*) 'Reading data file: ', fname12
      else
          open(unit=iunit,file=fname,status='old',form='formatted')
          write(6,*) 'Reading data file: ', fname
      endif

c     # this version may not work in f77
c     firstchar = '#'
c     commentlines = -1
c     do while (firstchar == '#')
c        read(iunit,*) firstchar
c        commentlines = commentlines + 1
c        enddo

      firstchar = '#'
      commentlines = -1
 10   continue
      if (firstchar .eq. '#') then
	 read(iunit,100) firstchar
 100     format(a1)
         commentlines = commentlines + 1
	 go to 10
      endif
      
      write(6,602) commentlines
  602 format('         first',i2,
     &       ' lines are comments and will be skipped')
      rewind(iunit)
      do line=1,commentlines
         read(iunit,*) firstchar
         enddo
      
      return
      end
