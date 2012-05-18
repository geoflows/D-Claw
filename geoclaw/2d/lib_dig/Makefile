
# Library version of Makefile
# ---------------------------

# Only "make .objs" and "make .htmls" are useful in this directory.

# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# Note: This version uses dynamic memory allocation from f90 and so requires
# an f90 compiler such as gfortran.

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran44   # default if not set as environment variable

# Add any desired compiler flags such as -g here:
FFLAGS ?= -O2

DIGCLAW_LIB = $(shell pwd)
GEOCLAW_LIB = $(CLAW)/geoclaw/2d/lib
CLAW_SOURCES = \
  $(GEOCLAW_LIB)/geoclaw_mod.f90 \
  $(DIGCLAW_LIB)/digclaw_mod.f90 \
  $(GEOCLAW_LIB)/topo_mod.f90 \
  $(GEOCLAW_LIB)/dtopo_mod.f90 \
  $(GEOCLAW_LIB)/qinit_mod.f90 \
  $(GEOCLAW_LIB)/auxinit_mod.f90 \
  $(DIGCLAW_LIB)/setprob_geo.f \
  $(DIGCLAW_LIB)/b4step2_geo.f \
  $(DIGCLAW_LIB)/setaux_geo.f \
  $(DIGCLAW_LIB)/qinit_geo.f \
  $(DIGCLAW_LIB)/rpn2ez_geo.f \
  $(DIGCLAW_LIB)/riemannsolvers_geo.f \
  $(DIGCLAW_LIB)/rpt2_geo.f \
  $(DIGCLAW_LIB)/src2_geo.f

default: .objs

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###

