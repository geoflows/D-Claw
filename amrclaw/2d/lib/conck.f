c
c -----------------------------------------------------------
c
      subroutine conck(level, nvar, time, rest)
c
      implicit double precision (a-h,o-z)

      include  "call.i"
      logical  rest

      iadd(i,j,ivar)  = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      iaddaux(i,j) = locaux + i - 1 + mitot*(j-1) +
     .                        mitot*mjtot*(mcapa-1)
c
c ******************************************************************
c conck - conservation check  for specified level
c         mostly a debugging tool
c         this assumes grids don't overlap
c ******************************************************************
c
c
c  grid loop for given level
c
      hx      = hxposs(level)
      hy      = hyposs(level)
      dt      = possk(level)
      totmass = 0.d0
      totvol  = 0.d0

      mptr = lstart(level)
 20   if (mptr .eq. 0) go to 85
         loc    = node(store1,mptr)
         locaux = node(storeaux,mptr)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot  = nx + 2*nghost
         mjtot  = ny + 2*nghost
c
         if (mcapa .eq. 0) then
           do 50 j  = nghost+1, mjtot-nghost
           do 50 i  = nghost+1, mitot-nghost
              if (dabs(alloc(iadd(i,j,1))).gt.1.d-16) then
                  sv = alloc(iadd(i,j,4))/alloc(iadd(i,j,1))
                  rho = 1000.d0*(1.d0-sv) + 2700.d0*sv
                  totmass = totmass + rho*alloc(iadd(i,j,1))
              endif
              totvol = totvol + alloc(iadd(i,j,1))
 50           continue
          else
c          # with capa array:
           do 60 j  = nghost+1, mjtot-nghost
           do 60 i  = nghost+1, mitot-nghost
              if (dabs(alloc(iadd(i,j,1))).gt.1.d-16) then
                  sv = alloc(iadd(i,j,4))/alloc(iadd(i,j,1))
                  rho = 1000.d0*(1.d0-sv) + 2700.d0*sv
                  totmass = totmass
     &               + rho*alloc(iadd(i,j,1))*alloc(iaddaux(i,j))
              endif
              totvol = totvol + alloc(iadd(i,j,1))*alloc(iaddaux(i,j))
 60           continue
          endif
c
       mptr = node(levelptr,mptr)
       go to 20
c
 85    totmass = totmass * hx * hy
 86    totvol = totvol * hx * hy
       if (time.eq.tstart .and. (level.eq.1) .and. .not. rest) then
           tmass0 = totmass
           tvol0 = totvol
           write(6,*) 'Total mass at initial time: ',tmass0
            write(6,*) 'Total volume at initial time: ',tvol0
           endif
       !write(outunit,777) time, totmass, 100.0*(totmass-tmass0)/tmass0
       !write(outunit,778) time, totvol, 100.0*(totvol-tvol0)/tvol0
 777   format('time t = ',e12.5,',  total mass = ',e22.15, ' diff % = ',
     &         e11.4)
 778   format('time t = ',e12.5,',  total vol = ',e22.15, '  diff % = ',
     &         e11.4)
c
 99   return
      end
