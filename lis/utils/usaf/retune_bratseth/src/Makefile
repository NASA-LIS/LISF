# Makefile for Fortran utilities
#
# To use, first run lis/configure script to generate configure.lis in lis/make.
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

procOBA_NWP: procOBA_NWP.o USAF_GridHashMod.o USAF_ReportsMod.o \
             USAF_StationsMod.o USAF_sharedMod.o
	$(FC) $^ -o procOBA_NWP -L$(LIB_ESMF) -lesmf \
                                -L$(LIB_NETCDF) -lnetcdff -lnetcdf

procOBA_Sat: procOBA_Sat.o USAF_GridHashMod.o USAF_ReportsMod.o \
             USAF_StationsMod.o USAF_sharedMod.o
	$(FC) $^ -o procOBA_Sat -L$(LIB_ESMF) -lesmf \
                                -L$(LIB_NETCDF) -lnetcdff -lnetcdf

procOBA_NWP.o: procOBA_NWP.F90 USAF_GridHashMod.o USAF_ReportsMod.o \
                     USAF_StationsMod.o USAF_sharedMod.o
	$(FC) $(FFLAGS) $< -I$(MOD_ESMF)

procOBA_Sat.o: procOBA_Sat.F90 USAF_GridHashMod.o USAF_ReportsMod.o \
                     USAF_StationsMod.o USAF_sharedMod.o
	$(FC) $(FFLAGS) $< -I$(MOD_ESMF)

USAF_GridHashMod.o : USAF_GridHashMod.F90
	$(FC) $(FFLAGS) $<

USAF_ReportsMod.o: USAF_ReportsMod.F90
	$(FC) $(FFLAGS) $<

USAF_StationsMod.o: USAF_StationsMod.F90
	$(FC) $(FFLAGS) $<

USAF_sharedMod.o: USAF_sharedMod.F90
	$(FC) $(FFLAGS) $< -I$(MOD_ESMF)

.PHONY: clean
clean:
	rm -f *.o *.mod procOBA_NWP procOBA_Sat

.PHONY: realclean
realclean:
	rm -f *.o *.mod procOBA_NWP procOBA_Sat
