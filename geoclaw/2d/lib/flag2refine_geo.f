c
c -------------------------------------------------------------------
      subroutine flag2refine(mx,my,mbc,meqn,maux,xlower,ylower,dx,dy,
     &                 t,level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
c -------------------------------------------------------------------

c
c ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
c
c User routine to control flagging of points for refinement.
c
c Specific for GeoClaw
c
c
c The logical function allowflag(x,y,t) is called to
c check whether further refinement at this level is allowed in this cell
c at this time.
c
c    q   = grid values including ghost cells (bndry vals at specified
c          time have already been set, so can use ghost cell values too)
c
c  aux   = aux array on this grid patch
c
c amrflags  = array to be flagged with either the value
c             DONTFLAG (no refinement needed)  or
c             DOFLAG   (refinement desired)
c

c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      use geoclaw_module
      use topo_module
      use dtopo_module
      use qinit_module

      implicit double precision (a-h, o-z)

      dimension   q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      dimension   aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
      dimension   amrflags(1-mbc:mx+mbc,1-mbc:my+mbc)
      logical     allowflag
      external  allowflag
      logical shoreregion,wave

      include 'regions.i'

c     # loop over interior points on this grid:

      do 200 j = 1,my
        y = ylower +  (j-0.5d0)*dy
        y1 = ylower + (j-1)*dy
        y2 = ylower + j*dy
        do 100 i = 1,mx
          x = xlower +  (i-0.5d0)*dx
          x1 = xlower +  (i-1)*dx
          x2 = xlower +  i*dx
c         # (i,j) grid cell is [x1,x2] x [y1,y2].

c         # default for each point is not to flag unless some condition
c         # below is satisfied:

          amrflags(i,j) = DONTFLAG


          do 30 m=1,mtopofiles
c           # check to see if refinement is forced in any topo file region:
            if (level .lt. minleveltopo(m) .and.
     &          t.ge.tlowtopo(m) .and. t.le.thitopo(m)) then
              xlow = xlowtopo(m)
              xhi = xhitopo(m)
              ylow = ylowtopo(m)
              yhi = yhitopo(m)
                 if (x2.gt.xlow.and.x1.lt.xhi.and.
     &               y2.gt.ylow.and.y1.lt.yhi) then
                    amrflags(i,j) = DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
   30       continue

           do 40 m=1,mregions
c           # check to see if refinement is forced in any other region:
            if (level .lt. minlevelregion(m) .and.
     &          t.ge.tlowregion(m) .and. t.le.thiregion(m)) then
              xlow = xlowregion(m)
              xhi = xhiregion(m)
              ylow = ylowregion(m)
              yhi = yhiregion(m)
                 if (x2.gt.xlow.and.x1.lt.xhi.and.
     &               y2.gt.ylow.and.y1.lt.yhi) then
                    amrflags(i,j) = DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
   40       continue

         do m = 1,num_dtopo
c           # check if we're in the dtopo region and need to refine:
c           # force refinement to level minleveldtopo
            t0dt = t0dtopo(m)
            tfdt = tfdtopo(m)
            minlevldt = minleveldtopo(m)
            if (level.lt.minleveldtopo(m).and.
     &              t.le.tfdtopo(m).and. !t.ge.t0dtopo(m).and.
     &              x2.gt.xlowdtopo(m).and.x1.lt.xhidtopo(m).and.
     &              y2.gt.ylowdtopo(m).and.y1.lt.yhidtopo(m)) then
                amrflags(i,j)=DOFLAG
                go to 100 !# flagged, so no need to check anything else
                endif
         enddo

         do m=1,mqinitfiles
            if (abs(t).lt.1.d0) then
c              # check if we're in the region where initial perturbation is
c              # specified and need to force refinement:
               if (level.lt.minlevelqinit(m).and.
     &              x2.gt.xlowqinit(m).and.x1.lt.xhiqinit(m).and.
     &              y2.gt.ylowqinit(m).and.y1.lt.yhiqinit(m)) then
                  amrflags(i,j)=DOFLAG
                go to 100 !# flagged, so no need to check anything else
                endif
             endif
         enddo

c        -----------------------------------------------------------------

c        # refinement not forced, so check if it is allowed, and if so,
c        # check if there is a reason to flag this point:

         if (allowflag(x,y,t,level)) then

             depth= q(i,j,1)
             momentum = sqrt(q(i,j,2)**2 + q(i,j,3)**2)
             surface = q(i,j,1) + aux(i,j,1)

c            # determine if flowgrades are used
            do iflow=1,mflowgrades
                  if (iflowgradevariable(iflow).eq.1) then
                    flowgradenorm=depth
                    flowgradegrad=depth

                  elseif (iflowgradevariable(iflow).eq.2) then
                    flowgradenorm=momentum
                    flowgradegrad=momentum

                  elseif (iflowgradevariable(iflow).eq.3) then

                    if (depth.gt.drytolerance) then
                      flowgradenorm=dabs(surface)
                      flowgradegrad=dabs(surface)
                    else
                      flowgradenorm=0.0
                      flowgradegrad=0.0
                    endif
                  endif

                  if (iflowgradetype(iflow).eq.1) then
                    flowgrademeasure=flowgradenorm
                  else
                    write(*,*) 'only flowgradetype = 1 supported'
                    stop
                    flowgrademeasure=flowgradegrad
                  endif

                  if (flowgrademeasure.gt.flowgradevalue(iflow)
     &                  .and.level.lt.iflowgrademinlevel(iflow)) then
                    amrflags(i,j)=DOFLAG
                    go to 100 !# flagged, so no need to check anything else
                  endif
            enddo

            if (mflowgrades.eq.0) then !tsunami-type refinement
               shoreregion = dabs(aux(i,j,1)) .lt. depthdeep
               wave = (dabs(surface-sealevel).gt.wavetolerance.and.
     &                q(i,j,1).gt.drytolerance)

c
               if (wave) then
c               # the surface is not at sea level here

                  if (level.lt.maxleveldeep) then
c                    # in deep water we can refine to this level
                     amrflags(i,j)=DOFLAG
                     go to 100 !# flagged, so no need to check anything else
                  endif

                  if (shoreregion.and.q(i,j,1).gt.drytolerance) then
c                    # following comment is regarding commented nested do loop above.
c                    # this cell is wet and a neighbor is dry ==> shoreline
                     amrflags(i,j)=DOFLAG
                     go to 100 !# flagged, so no need to check anything else
                  endif
               endif
            endif
         endif

 100    continue  !# end loop on i
 200    continue  !# end loop on j

      return
      end
