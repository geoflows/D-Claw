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

      double precision :: xlo,ylo,xhi,yhi

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


      xlo = x - 0.5*hxposs(level)
      xhi = x + 0.5*hxposs(level)
      ylo = y - 0.5*hxposs(level)
      yhi = y + 0.5*hxposs(level)

      do m=1,mtopofiles
     	if (level.lt.maxleveltopo(m)) then
      	  if (xlo.lt.xhitopo(m).and.xhi.gt.xlowtopo(m) .and.
     &	      ylo.lt.yhitopo(m).and.yhi.gt.ylowtopo(m) .and.
     &	      t.gt.tlowtopo(m).and.t.lt.thitopo(m)) then
                  allowflag=.true.
                  go to 900  !# no need to check anything else
	  endif
	endif
      enddo

      do m=1,mregions
     	if (level.lt.maxlevelregion(m)) then
      	  if (xlo.lt.xhiregion(m).and.xhi.gt.xlowregion(m).and.
     &	      ylo.lt.yhiregion(m).and.yhi.gt.ylowregion(m).and.
     &	      t.ge.tlowregion(m).and.t.le.thiregion(m)) then
                  allowflag=.true.
                  go to 900  !# no need to check anything else
	  endif
	endif
       enddo

      do m=1,num_dtopo
        if (x.gt.xlowdtopo(m).and.x.lt.xhidtopo(m).and.
     &      y.gt.ylowdtopo(m).and.y.lt.yhidtopo(m).and.
     &      t.ge.t0dtopo(m).and.t.le.tfdtopo(m)) then
            if (level.lt.maxleveldtopo(m)) then
                  allowflag=.true.
                  go to 900  !# no need to check anything else
                  endif
        endif
      enddo

      do m=1,mqinitfiles
         if (abs(t).lt.1.d0) then
            if (xlo.lt.xhiqinit(m).and.xhi.gt.xlowqinit(m).and.
     &	        ylo.lt.yhiqinit(m).and.yhi.gt.ylowqinit(m)) then
                if (level.lt.maxlevelqinit(m)) then
                     allowflag=.true.
                     go to 900  !# no need to check anything else
                 endif
            endif
         endif
      enddo

  900 continue
      return
      end
