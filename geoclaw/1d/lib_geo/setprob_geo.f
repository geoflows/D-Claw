      subroutine setprob
      implicit double precision (a-h,o-z)

      logical relimit1,relimit2

      common /geo/ grav,abs_tol 
      common /fluxlimiting/ relimit1,relimit2
      common /sourceterms/ frictioncoeff,ifrictiontype

      open(unit=7,file='setprob.data',status='old',form='formatted')

c     # Graviational constant g:
      read(7,*) grav
      read(7,*) abs_tol
      read(7,*) relimit1
      read(7,*) relimit2
      read(7,*) ifrictiontype
      read(7,*) frictioncoeff

      close(unit=7)
      return
      end
