      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /cparam/ rho,bulk,cc,zz

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
c     # Density and sound speed:

      read(7,*) rho
      read(7,*) cc
c
c     # Compute bulk modulus and impedance:

      bulk = cc*cc*rho
      zz = rho*cc

      return
      end
