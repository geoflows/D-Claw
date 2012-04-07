c=========================================================================
      subroutine setflood
c=========================================================================


      implicit double precision (a-h,o-z)
      character*20 fname
      logical foundFile

      include "call.i"
      include "geo.i"
      include "flowgrades.i"

      write(parmunit,*) ' '
      write(parmunit,*) '--------------------------------------------'
      write(parmunit,*) 'SETFLOOD:'
      write(parmunit,*) '------------'

c       # read user parameters from setflood.data

      fname  = 'setflood.data'
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname
        write(*,*) 'Or comment out call setflood in setprob'
        stop
      endif

      open(unit=7,file=fname,status='old',form='formatted')

      read(7,*) mflowgrades

      if (mflowgrades .eq. 0) then
         write(parmunit,*) '  No flow grades specified'
         return
         endif

      if (mflowgrades.gt.maxflowgrades) then
           write(*,*) 'SETFLOOD: ERROR mflowgrades > maxflowgrades'
           write(*,*) 'Decrease the number of flowgrades or'
           write(*,*) 'Increase maxflowgrades in flowgrades.i'
           stop
           endif

      do i=1,mflowgrades
         read(7,*) flowgradevalue(i),iflowgradevariable(i),
     &         iflowgradetype(i),iflowgrademinlevel(i)
         enddo

      close(7)

      write(parmunit,*) '   mflowgrades:',  mflowgrades

      do i=1,mflowgrades
         write(parmunit,701) flowgradevalue(i),iflowgradevariable(i),
     &         iflowgradetype(i),iflowgrademinlevel(i)
  701    format(d12.3,3i4)
         enddo

      return
      end
