c     =========================================
      logical function allowflag(x,y,t,level)
c     =========================================

c     # Indicate whether the grid point at (x,y,t) at this refinement level
c     # is allowed to be flagged for further refinement.
c
c     # Modified for GeoClaw to check whether the point lies in any of
c     # the various regions specified in the data files.
c
c     # KRB modify 12/30/2023 to use overlap region istead of contain.
c
c     # This routine is called from routine flag2refine.
c
c     # If Richardson error estimates are used (if tol>0) then this routine
c     # is also called from errf1.

      use geoclaw_module
      use topo_module
      use dtopo_module
      use qinit_module

      implicit double precision (a-h,o-z)

      include 'regions.i' ! needed for hxposs
      include 'call.i'

c========================================================================


      allowflag=.false.

c      following commented by dlg on 10/9/08.
c      my understanding of maxleveldeep might be differnet
c      still shouldn't be allowed if maxlevel allowed in a region is less
c      than maxleveldeep
c      might want to allow high levels of refinement in some deep regions
c      but not others.
c
c      if (level .lt. maxleveldeep) then
c         # allow refinement to next level in deep water
c          allowflag = .true.
c          go to 900  !# no need to check anything else
c          endif


      x1 = x - 0.5*hxposs(level)
      x2 = x + 0.5*hxposs(level)
      y1 = y - 0.5*hxposs(level)
      y2 = y + 0.5*hxposs(level)

      do m=1,mtopofiles
       	if (level.lt.maxleveltopo(m)) then

          xlow = xlowtopo(m)
          xhi = xhitopo(m)
          ylow = ylowtopo(m)
          yhi = yhitopo(m)

          if (x2.gt.xlow.and.x1.lt.xhi.and.
     &        y2.gt.ylow.and.y1.lt.yhi.and.
     &  	    t.gt.tlowtopo(m).and.t.lt.thitopo(m)) then
            allowflag=.true.
            go to 900  !# no need to check anything else
      	  endif
      	endif
      enddo

      do m=1,mregions
       	if (level.lt.maxlevelregion(m)) then

          xlow = xlowregion(m)
          xhi = xhiregion(m)
          ylow = ylowregion(m)
          yhi = yhiregion(m)

          if (x2.gt.xlow.and.x1.lt.xhi.and.
     &        y2.gt.ylow.and.y1.lt.yhi.and.
     &  	    t.ge.tlowregion(m).and.t.le.thiregion(m)) then
            allowflag=.true.
            go to 900  !# no need to check anything else
          endif
  	    endif
      enddo

      do m=1,num_dtopo
        if (level.lt.maxleveldtopo(m)) then
          xlow = xlowdtopo(m)
          xhi = xhidtopo(m)
          ylow = ylowdtopo(m)
          yhi = yhidtopo(m)

          if (x2.gt.xlow.and.x1.lt.xhi.and.
     &      y2.gt.ylow.and.y1.lt.yhi.and.
     &      t.ge.t0dtopo(m).and.t.le.tfdtopo(m)) then

              allowflag=.true.
              go to 900  !# no need to check anything else
          endif
        endif
      enddo

      do m=1,mqinitfiles
        if (level.lt.maxlevelqinit(m)) then
          if (abs(t).lt.1.d0) then
            xlow = xlowqinit(m)
            xhi = xhiqinit(m)
            ylow = ylowqinit(m)
            yhi = yhiqinit(m)

            if (x2.gt.xlow.and.x1.lt.xhi.and.
     &          y2.gt.ylow.and.y1.lt.yhi.and.
     &          t.ge.t0dtopo(m).and.t.le.tfdtopo(m)) then

              allowflag=.true.
              go to 900  !# no need to check anything else
            endif
          endif
         endif
      enddo

  900 continue
!      if (level.gt.1) then
!        write(*,*) "overlap", allowflag, level
!      endif
      return
      end
