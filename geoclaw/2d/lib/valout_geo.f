c
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c
      use digclaw_module, only : rho_f,rho_s,mom_autostop,mom_perc
      use digclaw_module, only: amidoneyet,globmaxmom

      implicit double precision (a-h,o-z)
      character*10  matname1, matname2, matname3
      double precision :: locmaxmom
      double precision momh, momvel, momm, momrho, momprop
      include  "call.i"

      logical outaux

      iadd(i,j,ivar) = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      iaddaux(i,j,iaux) = locaux + i - 1 + mitot*((iaux-1)*mjtot+j-1)

c ::::::::::::::::::::::::::: VALOUT ::::::::::::::::::::::::::::::::::;
c graphics output of soln values for contour or surface plots.
c modified for GeoClaw to output the surface level along with q.
c    surface = q(i,j,1) + aux(i,j,1)
c :::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::::::;

c
c     ### MATLAB graphics output
c
      outaux = .false.
      if (matlabout) then
c        ###  make the file names and open output files
         matname1 = 'fort.qxxxx'
         matname2 = 'fort.txxxx'
         matname3 = 'fort.axxxx'
         matunit1 = 50
         matunit2 = 60
         matunit3 = 70
         nstp     = matlabu
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            matname1(ipos:ipos) = char(ichar('0') + idigit)
            matname2(ipos:ipos) = char(ichar('0') + idigit)
            matname3(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue
         open(unit=matunit1,file=matname1,status='unknown',
     .       form='formatted')

         level = lst
         if (mom_autostop) then
           locmaxmom = 0. ! initialize local max momentum as zero. if
           ! using mom_autostop
         endif

         write(6,47) locmaxmom, globmaxmom
47       format('GeoClaw: Starting momentum calc ', d12.6, " ", d12.6)
         ngrids = 0
 65      if (level .gt. lfine) go to 90
            mptr = lstart(level)
 70         if (mptr .eq. 0) go to 80
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              write(matunit1,1001) mptr, level, nx, ny
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')


c  old        xcorn = rnode(cornxlo,mptr) - .5d0*hxposs(level)
c  old        ycorn = rnode(cornylo,mptr) - .5d0*hyposs(level)
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              write(matunit1,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy',/)

      do j = nghost+1, mjtot-nghost
         do i = nghost+1, mitot-nghost
            do ivar=1,nvar
               if (dabs(alloc(iadd(i,j,ivar))) .lt. 1d-90) then
                  alloc(iadd(i,j,ivar)) = 0.d0
               endif
            enddo
            if (mom_autostop) then
              if (level .eq. lst ) then
                ! calculate and add to momentum to get a level one momentum sum.
                momh = alloc(iadd(i,j,1))
                if (momh .gt. 0.0001) then ! if substantial thickness.

                  momvel = (  (alloc(iadd(i,j,2))/momh)**2.
     &                      + (alloc(iadd(i,j,3))/momh)**2.)**0.5
                  momm = alloc(iadd(i,j,4)) / momh
                  momrho = (rho_s * momm) + ((1.-momm) * rho_f)
                  ! hard code values for sediment and fluid density, but
                  ! better than nothing and prob approx right.
                  locmaxmom = locmaxmom ! momentum = (mass * velocity) = density * volume * velocity.
     &                        + (momh * hxposs(level) * hyposs(level)
     &                        * momrho * momvel)
                endif
              endif
            endif
            surface = alloc(iadd(i,j,1)) + alloc(iaddaux(i,j,1))
     &              - alloc(iadd(i,j,7))
            write(matunit1,109) (alloc(iadd(i,j,ivar)), ivar=1,nvar),
     &         surface
         enddo
         write(matunit1,*) ' '
      enddo
  109       format(9e26.16)


            mptr = node(levelptr, mptr)
            go to 70
 80      level = level + 1
         go to 65

 90     continue
        ! If using mom_autostop then make calculations.
        ! local step_max. if in excess of mom_perc of total max change amidoneyet.
        ! set max to new max.
        if (mom_autostop) then
          ! if current max is largest, increase the global max
          if (locmaxmom .gt. globmaxmom) then
            globmaxmom = locmaxmom
          endif

          ! calculate new momentum proportion.
          momprop = locmaxmom / globmaxmom

          ! test if done condition is met.
          if (momprop .le. mom_perc) then
            amidoneyet = .True. ! assign amidoneyet as true if done.
          endif

          ! write current status to the log.
          write(6,112) momprop, locmaxmom, globmaxmom, time
 112      format('GeoClaw: Current momentum proportion ', d12.6,
     &           ' ( ', d12.6 ' / ', d12.6,
     &           ') at t = ', d12.6,/)

      endif
      ! end if mom_autostop

        if (outaux) then
c        # output aux array to fort.aXXXX
         open(unit=matunit3,file=matname3,status='unknown',
     .       form='formatted')
         level = lst
         ngrids = 0
 165     if (level .gt. lfine) go to 190
            mptr = lstart(level)
 170        if (mptr .eq. 0) go to 180
              ngrids  = ngrids + 1
              nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
              ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
              loc     = node(store1, mptr)
              locaux  = node(storeaux,mptr)
              mitot   = nx + 2*nghost
              mjtot   = ny + 2*nghost
              if (ny.gt.1) then
                  write(matunit3,1001) mptr, level, nx, ny
                else
c                 # output in 1d format if ny=1:
                  write(matunit3,1003) mptr, level, nx
                endif
              xlow = rnode(cornxlo,mptr)
              ylow = rnode(cornylo,mptr)
              if (ny.gt.1) then
                  write(matunit3,1002)
     &              xlow,ylow,hxposs(level),hyposs(level)
                else
                  write(matunit3,1004)
     &              xlow,hxposs(level)
                endif

         do j = nghost+1, mjtot-nghost
            do i = nghost+1, mitot-nghost
               do ivar=1,naux
                  if (dabs(alloc(iaddaux(i,j,ivar))) .lt. 1d-90) then
                     alloc(iaddaux(i,j,ivar)) = 0.d0
                  endif
               enddo
               write(matunit3,109) (alloc(iaddaux(i,j,ivar)),
     &                              ivar=1,naux)
            enddo
            write(matunit3,*) ' '
         enddo

            mptr = node(levelptr, mptr)
            go to 170
 180     level = level + 1
         go to 165

 190    continue
        close(unit=matunit3)
        endif !# end outputting aux array

        open(unit=matunit2,file=matname2,status='unknown',
     .       form='formatted')

c     # nvar+1 variable printed since surface also printed

      !write(matunit2,1000) time,nvar+1,ngrids,3,2
      write(matunit2,1000) time,nvar+1,ngrids,naux,2
 1000 format(e18.8,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 naux'/,
     &       i5,'                 ndim'/,/)
c
 1003 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx')

 1004 format(e18.8,'    xlow', /,
     &       e18.8,'    dx', /)

      write(6,601) matlabu,time
  601 format('GeoClaw: Frame ',i4,
     &       ' output files done at time t = ', d12.6,/)

      matlabu = matlabu + 1

      close(unit=matunit1)
      close(unit=matunit2)
      endif

      return
      end
