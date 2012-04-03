c
c
c =========================================================
      subroutine out2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,t,iframe,aux,maux)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 2 dimensions
c
c     # Write the results to the file fort.q<iframe>.nc
c     # Use format required by matlab script  plotclaw2.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c     # set outaux = .true. to also output the aux arrays to fort.a<iframe>
c

c
c -----------------------------------------------------
c     Routine to write netcdf files in the classic format
!        #jj-2011.03.29
!        # Each file written by the fortran code has 
!        # Dimensions:
!        #           timedimension : UNLIMITED
!        #           meqn          : The number of equations
!        #           dimx_<gridno> : X dimension for grid number <gridno>
!        #           dimy_<gridno> : Y dimension for grid number <gridno>
!        # Variables:
!        #           timedimension : Stores the time of the frame
!        #           ngrids        : Number of grids in this frame
!        #           naux          : Number of Auxilary Variables
!        #           ndim          : Number of Dimensions in the frame
!        #           grid_<gridno> : A grid of (dimx,dimy,meqn)
!        # Attributes:
!        # (grid_<no>) gridno      : The number of this grid <grid_no>
!        #           level         : The AMR level
!        #           dim_names     : a list of dimensions [dimx,dimy]
!        #           dim<x,y>.low  : The lowest dimension value 
!        #           dim<x,y>.d    : The distance between grid points 
c -----------------------------------------------------

c
      implicit double precision (a-h,o-z)
      include 'netcdf.inc'
      
      dimension   q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
      character*10 fname1, fname2, fname3
      logical outaux,do_ascii
            
      real(kind=8) time
      integer ncid,rcode
      integer timeid,tVarID,meqnID,ngridsVarID,nauxVarID,ndimVarID
      integer dimxid,dimyid,xlowid,ylowid,dxid,dyid
      integer gridid
      integer ntimes
      character*2 gridstr
      character*40 dim_names
      REAL(kind=8), ALLOCATABLE  ::grid(:,:,:)
      integer nx, ny
      real xlow,ylow
      

      nx=mx
      ny=my
      outaux = .false.
      do_ascii = .false.
c
c     # first create the file name and open file
c
         fname1 = 'fort.qxxxx'
         fname2 = 'fort.txxxx'
         fname3 = 'fort.axxxx'
         nstp = iframe
         do 55 ipos = 10, 7, -1
            idigit = mod(nstp,10)
            fname1(ipos:ipos) = char(ichar('0') + idigit)
            fname2(ipos:ipos) = char(ichar('0') + idigit)
            fname3(ipos:ipos) = char(ichar('0') + idigit)
            nstp = nstp / 10
 55      continue

         if (do_ascii) then  open(unit=50,file=fname1,status='unknown',form='formatted')
         if (do_ascii) then  open(unit=60,file=fname2,status='unknown',form='formatted')

c
c     # the following parameters are used in amrclaw where there are
c     # multiple grids.  Here they are all set to 1:
      ngrids = 1
      mptr = 1
      level = 1

      ntimes=1
c        write(50,1001) mptr,level,mx,my
c        write(50,1002) xlower,ylower,dx,dy
        !!!!Define netcdf file
         rcode=NF_CREATE(fname1//'.nc',NF_NOCLOBBER,ncid)
         if(rcode.ne.NF_NOERR) print *,'ERROR OPENING NETCDF FILE'
         rcode=NF_DEF_DIM(ncid,'timedimension',NF_UNLIMITED,timeid)
         rcode=NF_DEF_VAR(ncid,'timedimension',NF_DOUBLE,1,timeid,
     &   tVarID)
         rcode=NF_DEF_DIM(ncid,'meqn',meqn,meqnid)
         rcode=NF_DEF_VAR(ncid,'ngrids',NF_INT,0,0,ngridsVarID)
         rcode=NF_DEF_VAR(ncid,'naux',NF_INT,0,0,nauxVarID)
         rcode=NF_DEF_VAR(ncid,'ndim',NF_INT,0,0,ndimVarID)
        write(gridstr,67) mptr
              
67            format(I2.2)              
              rcode=NF_DEF_DIM(ncid,'dimx_'//trim(gridstr),nx,dimxid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE DIMS'
              rcode=NF_DEF_DIM(ncid,'dimy_'//trim(gridstr),ny,dimyid)
              if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE DIMS'
              

              rcode=NF_DEF_Var(ncid,'grid_'//trim(gridstr),NF_DOUBLE,4,
     &              (/dimxid,dimyid,meqnid,timeid/),gridid)
               if(rcode.ne.NF_NOERR) print *,'ERROR  DEFINE VAR'
               
              rcode=NF_PUT_ATT_INT(ncid,gridid,'gridno',NF_INT,1,
     &              mptr)
     
              rcode=NF_PUT_ATT_INT(ncid,gridid,'level',NF_INT,1,level)
              
              dim_names="['dimx','dimy']"
              rcode=NF_PUT_ATT_TEXT(ncid,gridid,'dim_names',
     &         LEN_TRIM(dim_names),TRIM(dim_names))
     
              rcode=NF_PUT_ATT_DOUBLE(ncid,gridid,'dimx.lower',NF_FLOAT,
     &              1,xlower)     
              rcode=NF_PUT_ATT_DOUBLE(ncid,gridid,'dimy.lower',NF_FLOAT,
     &              1,ylower)
     
              rcode=NF_PUT_ATT_DOUBLE(ncid,gridid,'dimx.d',NF_FLOAT,1,
     &          dx)
              rcode=NF_PUT_ATT_DOUBLE(ncid,gridid,'dimy.d',NF_FLOAT,1,
     &          dy) 
     
              rcode=NF_ENDDEF(ncid)
            write(0,*) "xlower, ylower:",xlower,ylower,sizeof(xlower)        
      allocate(grid(nx,ny,meqn))
      time=t      
c
      do j=1,my
        do i=1,mx
          do m=1,meqn
c            # exponents with more than 2 digits cause problems reading
c            # into matlab... reset tiny values to zero:
             if (dabs(q(i,j,m)) .lt. 1d-99) q(i,j,m) = 0.d0
             grid(i,j,m)=q(i,j,m)
             enddo
              if (do_ascii) then write(50,1005) (grid(i,j,m), m=1,meqn)
         enddo
         if (do_ascii) then write(50,*) ' '
      enddo
      if (do_ascii) then write(50,*) ' '
      
      rcode=NF_PUT_VARA_DOUBLE(ncid,gridid,(/1,1,1,1/),
     & (/nx,ny,meqn,1/),grid)
     
      deallocate(grid)
      
      rcode=NF_PUT_VAR_DOUBLE(ncid,tVarID,t)
      if(rcode.ne.NF_NOERR) print *,'ERROR  Write Time'
      rcode=NF_PUT_VAR_INT(ncid,ngridsVarID,int(ngrids))
      if(rcode.ne.NF_NOERR) print *,'ERROR  Write GridNo'
      rcode=NF_PUT_VAR_INT(ncid,nauxVarID,maux)
      rcode=NF_PUT_VAR_INT(ncid,ndimVarID,2)
      rcode=NF_CLOSE(ncid)

 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my')

 1002 format(e26.16,'    xlow', /,
     &       e26.16,'    ylow', /,
     &       e26.16,'    dx', /,
     &       e26.16,'    dy',/)
 1005     format(4e26.16)

            
      if (outaux) then 
c     # also output the aux arrays:
      open(unit=70,file=fname3,status='unknown',form='formatted')
      write(70,1001) mptr,level,mx,my
      write(70,1002) xlower,ylower,dx,dy
      do 120 j=1,my
         do 110 i=1,mx
            do m=1,maux
c              # exponents with more than 2 digits cause problems reading
c              # into matlab... reset tiny values to zero:
               if (dabs(aux(i,j,m)) .lt. 1d-99) aux(i,j,m) = 0.d0
            enddo
c
            write(70,1005) (aux(i,j,m), m=1,maux)
c
  110       continue
         write(70,*) ' '
  120    continue
      write(70,*) ' '
      close(unit=70)
      endif

      if (do_ascii) then write(60,1000) t,meqn,ngrids,maux,2

 1000 format(e26.16,'    time', /,
     &       i5,'                 meqn'/,
     &       i5,'                 ngrids'/,
     &       i5,'                 maux'/,
     &       i5,'                 ndim'/,/)
!c
      if (do_ascii) then close(unit=50)
      if (do_ascii) then close(unit=60)

      return
      end
