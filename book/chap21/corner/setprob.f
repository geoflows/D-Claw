      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comaux/ zl,cl,zr,cr
c
c     # Set the material parameters for the acoustic equations
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                

c
c     # Piecewise constant medium with single interface at x=0
c     # Impedance and sound speed to left and right:

      read(7,*) zl
      read(7,*) cl
      read(7,*) zr
      read(7,*) cr

      return
      end
