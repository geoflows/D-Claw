
c
c ------------------------------------------------------------------
c
      subroutine bc2amr(val,aux,nrow,ncol,meqn,naux,
     1                  hx, hy, level, time,
     2                  xleft,  xright,  ybot, ytop,
     3                  xlower, ylower,xupper,yupper,
     4                  xperiodic, yperiodic,spheredom)

c
c    Specific to geoclaw:  extrapolates aux(i,j,1) at boundaries
c    to constant.
c
c :::::::::: bc2amr ::::::::::::::::::::::::::::::::::::::::::::::;
c
c     Take a grid patch with mesh widths hx,hy, of dimensions nrow by
c     ncol,  and set the values of any piece of
c     of the patch which extends outside the physical domain
c     using the boundary conditions.
c
c     ------------------------------------------------
c     # Standard boundary condition choices for amr2ez in clawpack
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     #            =  4  sphere bcs (left half maps to right half of same
c     #                  side, and vice versa), as if domain folded in half
c     ------------------------------------------------
c
c     The corners of the grid patch are at
c        (xleft,ybot)  --  lower left corner
c        (xright,ytop) --  upper right corner
c
c     The physical domain itself is a rectangle bounded by
c        (xlower,ylower)  -- lower left corner
c        (xupper,yupper)  -- upper right corner
c
c     the picture is the following:
c
c               _____________________ (xupper,yupper)
c              |                     |
c          _________ (xright,ytop)   |
c          |   |    |                |
c          |   |    |                |
c          |   |    |                |
c          |___|____|                |
c (xleft,ybot) |                     |
c              |                     |
c              |_____________________|
c   (xlower,ylower)
c
c
c     Any cells that lie outside the physical domain are ghost cells whose
c     values should be set in this routine.  This is tested for by comparing
c     xleft with xlower to see if values need to be set at the left, as in
c     the figure above, and similarly at the other boundaries.
c
c     Patches are guaranteed to have at least 1 row of cells filled
c     with interior values so it is possible to  extrapolate.
c     Fix trimbd if you want more than 1 row pre-set.
c
c     Make sure the order the boundaries are specified is correct
c     so that diagonal corner cells are also properly taken care of.
c
c     Periodic boundaries are set before calling this routine, so if the
c     domain is periodic in one direction only you
c     can safely extrapolate in the other direction.
c
c     Don't overwrite ghost cells in periodic directions!

c     This particular routine bc2amr_noslopesets auxillary values so
c     that no slope in topography occurs at the physical boundary.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      implicit double precision (a-h,o-z)

      common /combc2/ mthbc(4)

      dimension val(nrow,ncol,meqn), aux(nrow,ncol,naux)
      logical xperiodic, yperiodic, spheredom

      hxmarg = hx*.01
      hymarg = hy*.01

      if (xperiodic .and. (yperiodic .or. spheredom)) go to 499


c      Either by passing data arrays or looking at a file called something
c      like "setdbc.data" at this point learn about whether there are any
c      input-file specified boundary conditions specified by an external
c      file.
c
c      Follow the pattern of dtopo files in specifying these inputs.
c      That is, dtopo input file type 3
c      (http://www.clawpack.org/topo.html#topo-dtopo)
c
c      End goal, these files will be described in setrun.py with something like:
c
c      digdata.dbcfiles.append([
c              dtopotype,
c              minlevel,
c              maxlevel,
c              bc_rest,
c              iq,
c              fname])
c
c      where
c        dtopotype is the dtopo file type
c        minlevel and max level are the min and max AMR levels
c        bc_rest is the boundary condition type for the rest of the boundary.
c            (alternatively we could just have the bcs for the rest specified
c            by clawdata.mthbc_XXX in setrun.py) specify everything EXCEPT
c            what is by an overriding file. Then rather than complex
c            indexing selection based on what is in/out of the file, you
c            could fill the whole domain boundary strip based on mthbc_XXX
c            (using current routines) and then overwrite based on the file
c            at the very end.
c               This makes the most sense to me, so I'm going to stub out
c               below.
c      So if you wanted to specify n elements of q, then you would need n
c      elements appended to dbcfiles.
c
c      Some statements here about reasonable default values?
c
c      either h, hu, or hv must be set?
c      If not specified, default values are (the following are KRB thoughts based on 10 seconds of thought)
c        h = zero order extrapolation
c        hv = 0
c        hu = 0
c        hm = zero order extrapolation of h times m0
c        pb = hydrostatic

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      if (xleft .ge. xlower-hxmarg) then
c        # not a physical boundary -- no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 199
         endif
c
c     # number of grid cells from this patch lying outside physical domain:
      nxl = (xlower+hxmarg-xleft)/hx
c
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2amr'
      stop
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 m=1,meqn
         do 115 i=1,nxl
            do 115 j = 1,ncol
               aux(i,j,1) = aux(nxl+1,j,1)  !inserted for bc2amr_noslope
               val(i,j,m) = val(nxl+1,j,m)
  115       continue
      go to 199

  120 continue
c     # periodic:   handled elsewhere in amr
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 m=1,meqn
         do 135 i=1,nxl
            do 135 j = 1,ncol
               aux(i,j,1) = aux(2*nxl+1-i,j,1)  !inserted for bc2amr_noslope
               val(i,j,m) = val(2*nxl+1-i,j,m)
  135       continue
c     # negate the normal velocity:
      do 136 i=1,nxl
         do 136 j = 1,ncol
            val(i,j,2) = -val(i,j,2)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      if (xright .le. xupper+hxmarg) then
c        # not a physical boundary --  no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 299
         endif
c
c     # number of grid cells lying outside physical domain:
      nxr = (xright - xupper + hxmarg)/hx
      ibeg = max0(nrow-nxr+1, 1)
c
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2amr'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do 215 m=1,meqn
         do 215 i=ibeg,nrow
            do 215 j = 1,ncol
               aux(i,j,1) = aux(ibeg-1,j,1) !inserted for bc2amr_noslope
               val(i,j,m) = val(ibeg-1,j,m)
  215       continue
      go to 299

  220 continue
c     # periodic:   handled elsewhere in amr
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 m=1,meqn
         do 235 i=ibeg,nrow
            do 235 j = 1,ncol
               aux(i,j,1) = aux(2*ibeg-1-i,j,1) !inserted for bc2amr_noslope
               val(i,j,m) = val(2*ibeg-1-i,j,m)
  235       continue
c     # negate the normal velocity:
      do 236 i=ibeg,nrow
         do 236 j = 1,ncol
            val(i,j,2) = -val(i,j,2)
  236    continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      if (ybot .ge. ylower-hymarg) then
c        # not a physical boundary -- no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 399
         endif
c
c     # number of grid cells lying outside physical domain:
      nyb = (ylower+hymarg-ybot)/hy
c
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2amr'
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 m=1,meqn
         do 315 j=1,nyb
            do 315 i=1,nrow
                aux(i,j,1) = aux(i,nyb+1,1) !inserted for bc2amr_noslope
                val(i,j,m) = val(i,nyb+1,m)
  315       continue
      go to 399

  320 continue
c     # periodic:   handled elsewhere in amr
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 m=1,meqn
         do 335 j=1,nyb
            do 335 i=1,nrow
                aux(i,j,1) =  aux(i,2*nyb+1-j,1) !inserted for bc2amr_noslope
                val(i,j,m) =  val(i,2*nyb+1-j,m)
  335       continue
c     # negate the normal velocity:
      do 336 j=1,nyb
         do 336 i=1,nrow
            val(i,j,3) = -val(i,j,3)
  336    continue
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      if (ytop .le. yupper+hymarg) then
c        # not a physical boundary --  no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 499
         endif
c
c     # number of grid cells lying outside physical domain:
      nyt = (ytop - yupper + hymarg)/hy
      jbeg = max0(ncol-nyt+1, 1)
c
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2amr'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 m=1,meqn
         do 415 j=jbeg,ncol
            do 415 i=1,nrow
               aux(i,j,1) = aux(i,jbeg-1,1)  !inserted for bc2amr_noslope
               val(i,j,m) =  val(i,jbeg-1,m)
  415       continue
      go to 499

  420 continue
c     # periodic:   handled elsewhere in amr
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 m=1,meqn
         do 435 j=jbeg,ncol
            do 435 i=1,nrow
               aux(i,j,1) =  aux(i,2*jbeg-1-j,1)  !inserted for bc2amr_noslope
               val(i,j,m) =  val(i,2*jbeg-1-j,m)
  435       continue
c     # negate the normal velocity:
      do 436 j=jbeg,ncol
         do 436 i=1,nrow
            val(i,j,3) = -val(i,j,3)
  436    continue
      go to 499

  499 continue

      return
      end
! KRB PSEUDOCODE.
!
! this block would start with if setdbc.data exists or if setdbc_flag = True.
!
! Then files would be read and/or arrays already loaded would be used.
!
! for now...
!
! hard code some variables here (not declared above yet) that will
! eventually be what would be read in in setdbc.data based on the current
! value of time.
!
! ymin, ymax, xmin, xmax, q1, q2, q3, q4, q5 are all floats.
! jstart, jend, istart, iend are all ints.

ymin = XXX ! these four define the bounding box of the dbc file
ymax = XXX
xmin = XXX
xmax = XXX

! (note that this means the extent needs to cover the extent of the ghost cells.
! there should be a check/warning for this. Also careful precision leq/geq/lt/gt
! checking for overlap)

q1 = XXX # values for each of element of q within the dbc bounding box.
q2 = XXX # eventually these may be spatially variable. floats for now.
q3 = XXX
q4 = XXX
q5 = XXX

! Ensure it is possible for multiple bcs to be specified (e.g., two on
! left side, one on left and one on right so that multiple inflow
! hydrographs can be specified).
!
! Do some sort of check about sub/supercritical. If subcritical print a
! warning to output.
!
! For each bc domain for this timestep (including time interpolation
! already supported by dtopo capabilities?)
!
!       use domain bounds and dbc bounds to determine the i and j indice
!       values for val that are modified by this part of setdbc.data.

jstart = XXX
istart = XXX
jend = XXX
iend = XXX

! such that we are filling val[istart:iend+1, jstart:jend+1, :] = [q1, q2, q3, q4, q5] # interpret this as python

!       At present, fill all those with values for q1--q5.
!       Eventually, use existing interpolation schemes to interpolate
!       values correctly given the grid discretization and the extent of
!       the dbc domain.
!
!       use default values where values are not specified.
