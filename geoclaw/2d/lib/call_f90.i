
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::   data structure info.
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       integer    cornxlo,cornylo,cornxhi,cornyhi,timemult
       integer    store1,store2,storeaux
       integer    tempptr,errptr,ffluxptr,cfluxptr
       integer    rsize

       parameter (rsize =  5)
       parameter (nsize = 13)

!  :::::::   integer part of node descriptor
       parameter (levelptr  = 1)
       parameter (tempptr   = 2)
       parameter (errptr    = 3)
       parameter (nestlevel = 4)
       parameter (cfluxptr  = 5)
       parameter (ffluxptr  = 6)
       parameter (store1    = 7)
       parameter (store2    = 8)
       parameter (ndilo     = 9)
       parameter (ndihi     = 10)
       parameter (ndjlo     = 11)
       parameter (ndjhi     = 12)
       parameter (storeaux  = 13)

! :::::::  real part of node descriptor
       parameter (cornxlo  = 1)
       parameter (cornylo  = 2)
       parameter (cornxhi  = 3)
       parameter (cornyhi  = 4)
       parameter (timemult = 5)

! :::::::   for linking nodes
       parameter (nextfree = 2)
       parameter (null = 0)
       parameter (nil  = 0)

! :::::::  for flagging points
       parameter (goodpt = 0.0)
       parameter (badpt  = 2.0)
       parameter (badpro = 3.0)

       parameter (rinfinity = 10.e32)
       parameter (iinfinity = 999999)
       parameter (horizontal = 1)
       parameter (vertical = 2)
       parameter  (maxgr = 5000)
       parameter  (maxlv = 10)
       parameter  (maxcl = 500)
       parameter  (max1d = 60)
       parameter  (maxvar = 10)
       parameter  (maxaux = 20)
       parameter  (maxout = 5000)

       logical    printout,matlabout,ncarout


       common  /nodal/ &
              hxposs(maxlv), hyposs(maxlv),possk(maxlv), &
              rnode(rsize, maxgr), node(nsize, maxgr),   &
              icheck(maxlv),lstart(maxlv),newstl(maxlv),  &
              listsp(maxlv),intratx(maxlv),intraty(maxlv), &
              kratio(maxlv), iregsz(maxlv),jregsz(maxlv), &
              iregst(maxlv),jregst(maxlv), &
              iregend(maxlv),jregend(maxlv), &
              numgrids(maxlv),numcells(maxlv), &
              tol, tolsp, ibuff,  mstart, ndfree, lfine, &
              iorder,mxnest,kcheck,nghost,printout,matlabout,ncarout

       common /cmatlab/ matlabu,ngrids
!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      ::::  for alloc array/memory
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!       Static memory implementation
!        parameter  (memsize = 10000000)
!        common  /calloc/   alloc(memsize)

!      Dynamic memmory
       double precision, pointer, dimension(:) :: alloc
       common /calloc/ alloc, memsize


!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::   for space management of alloc array
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       parameter (lfdim=5000)

       common /space/ lfree(lfdim,2),lenf

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  domain description variables
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       logical xperdom, yperdom, spheredom

       common /cdom/ xupper,yupper,xlower,ylower, &
                      xperdom, yperdom, spheredom


!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  collect stats
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       common   /stats/  rvoll(10),evol,rvol,lentot,lenmax,lendim

!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      :::::  method parameters
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
       parameter (maxwave = 10)
       character * 10 auxtype(maxaux)
       common /cmeth/ method(7), mthlim(maxwave), mwaves
       common /cmcapa/  mcapa
       common /auxstuff/ auxtype
       common /ccfl/ cfl,cflmax,cflv1,cfl_level
!
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!      ::::  for i/o assignments
!      ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

       integer chkunit,dbugunit,inunit,outunit,pltunit1,rstunit
       integer matunit,parmunit

       parameter (parmunit = 12)
       parameter (chkunit = 10)
       parameter (inunit  = 5)
       parameter (outunit = 66)
       parameter (pltunit1 = 3)
       parameter (rstunit = 9)
       parameter (dbugunit = 11)
       parameter (matunit = 70)


!      ::::  common for debugging flags (verbose output)

       logical &
              dprint,    & !  domain flags output 
              eprint,    & !  error estimation output
              edebug,    & !  even more error estimation output
              gprint,    & !  verbose grid generation (clustering,colating...)
              nprint,    & !  nestck reporting
              pprint,    & !  projec tagged pts.
              rprint,    & !  regridding -  summary of new grids
              sprint,    & !  space (memory) output
              tprint,    & !  tick (time stepping) reporting
              uprint       !  updating/upbnding reporting


       common /bugflags/ dprint, eprint, edebug, gprint, nprint, pprint, &
                        rprint, sprint, tprint, uprint

!
!      ::::  common for conservation check
      common /ctstart/ tstart
      common /comconck/ tmass0,tvol0
