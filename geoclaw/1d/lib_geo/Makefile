#
#	Makefile for the geoclaw 1d library code:
#
#	To compile the claw1 library files type:  make lib
#       (from the this sub-directory)
#
#	To compile a single file.f type:	  make file.o
#
#
# Set your environment compiler $FC and compiler flags $FFLAGS
LFLAGS = $(FFLAGS) -O
LINK   = $(FC)

# Or comment above, uncomment below to use some defaults
#LFLAGS =  -funroll-loops -finline-functions -fdefault-real-8
#LINK = gfortran

# Build rules
%.o : %.f ; $(LINK) $(LFLAGS) -c $< -o $@

LIBOBJECTS = \
  $(GEOCLAW)/1d/lib/bc1_geo.o \
  $(GEOCLAW)/1d/lib/b4step1_geo.o \
  $(GEOCLAW)/1d/lib/novaccum1_geo.o \
  $(GEOCLAW)/1d/lib/src1_geo.o \
  $(GEOCLAW)/1d/lib/out1_geo.o \
  $(GEOCLAW)/1d/lib/rp1_geo.o \
  $(GEOCLAW)/1d/lib/riemannsolvers_geo.o \
  $(GEOCLAW)/1d/lib/rp1ez_geo.o \
  $(GEOCLAW)/1d/lib/step1fw_geo.o \
  $(GEOCLAW)/1d/lib/claw1.o \


# make list of .o files required from the sources above:
LIBSOURCES = $(subst .o,.f, $(LIBOBJECTS))

lib: $(LIBOBJECTS)

backup : $(LIBSOURCES)
	cat $(LIBSOURCES) > library.`date +%y%m%d`

clean:
	-rm -f $(LIBOBJECTS)

new:
				make clean
				make lib


### DO NOT remove this line - make depends on it ###
