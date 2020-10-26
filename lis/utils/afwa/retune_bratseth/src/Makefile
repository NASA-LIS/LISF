# Makefile for Fortran utilities
#
# To use, first run lvt/configure script to generate configure.lis in lis/make.
#
# Eric Kemp, SSAI and NASA/GSFC, 26 Oct 2020
#------------------------------------------------------------------------------

# Use the configure.lis found in LISF/lis/make
include ../../../../make/configure.lis

# Optionally overwrite compiler and compiler flags
#FC = mpif90
##FFLAGS = -O0 -g -traceback -check all -warn
#FFLAGS = -O2 -fpe0

all: procOBA_NWP procOBA_Sat

procOBA_NWP: procOBA_NWP.o EMK_GridHashMod.o EMK_ReportsMod.o \
             EMK_StationsMod.o sharedMod.o
	$(FC) $^ -o procOBA_NWP -L$(LIB_ESMF) -lesmf \
                                -L$(LIB_NETCDF) -lnetcdff -lnetcdf

procOBA_Sat: procOBA_Sat.o EMK_GridHashMod.o EMK_ReportsMod.o \
             EMK_StationsMod.o sharedMod.o
	$(FC) $^ -o procOBA_Sat -L$(LIB_ESMF) -lesmf \
                                -L$(LIB_NETCDF) -lnetcdff -lnetcdf

procOBA_NWP.o: procOBA_NWP.F90 EMK_GridHashMod.o EMK_ReportsMod.o \
                     EMK_StationsMod.o sharedMod.o
	$(FC) -c $(FFLAGS) $< -I$(MOD_ESMF)

procOBA_Sat.o: procOBA_Sat.F90 EMK_GridHashMod.o EMK_ReportsMod.o \
                     EMK_StationsMod.o sharedMod.o
	$(FC) -c $(FFLAGS) $< -I$(MOD_ESMF)

EMK_GridHashMod.o : EMK_GridHashMod.F90
	$(FC) -c $(FFLAGS) $<

EMK_ReportsMod.o: EMK_ReportsMod.F90
	$(FC) -c $(FFLAGS) $<

EMK_StationsMod.o: EMK_StationsMod.F90
	$(FC) -c $(FFLAGS) $<

sharedMod.o: sharedMod.F90
	$(FC) -c $(FFLAGS) $< -I$(MOD_ESMF)


.PHONY: clean
clean:
	rm -f *.o *.mod procOBA_NWP procOBA_Sat
