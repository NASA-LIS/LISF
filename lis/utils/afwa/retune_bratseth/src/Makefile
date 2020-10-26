#ESMF = /discover/nobackup/projects/lis/libs/esmf/7.1.0r_intel-18.0.3.222_impi-18.0.3.222/
ESMF = /discover/nobackup/projects/lis/libs/esmf/7.1.0r_intel-19.1.0.166_impi-20.0.0.166_sles12.3
ESMF_LIB = $(ESMF)/lib/libO/Linux.intel.64.intelmpi.default
ESMF_MOD = $(ESMF)/mod/modO/Linux.intel.64.intelmpi.default

NETCDF = /discover/nobackup/projects/lis/libs/netcdf/4.5.0_intel-19.1.0.166_sles12.3
NETCDF_LIB = $(NETCDF)/lib

FC = mpif90
#FFLAGS = -O0 -g -traceback -check all -warn
FFLAGS = -O2 -fpe0

all: procOBA_NWP procOBA_Sat


procOBA_NWP: procOBA_NWP.o EMK_GridHashMod.o EMK_ReportsMod.o \
             EMK_StationsMod.o sharedMod.o
	$(FC) $^ -o procOBA_NWP -L$(ESMF_LIB) -lesmf \
                                -L$(NETCDF_LIB) -lnetcdff -lnetcdf

procOBA_Sat: procOBA_Sat.o EMK_GridHashMod.o EMK_ReportsMod.o \
             EMK_StationsMod.o sharedMod.o
	$(FC) $^ -o procOBA_Sat -L$(ESMF_LIB) -lesmf \
                                -L$(NETCDF_LIB) -lnetcdff -lnetcdf

procOBA_NWP.o: procOBA_NWP.F90 EMK_GridHashMod.o EMK_ReportsMod.o \
                     EMK_StationsMod.o sharedMod.o
	$(FC) -c $(FFLAGS) $< -I$(ESMF_MOD)

procOBA_Sat.o: procOBA_Sat.F90 EMK_GridHashMod.o EMK_ReportsMod.o \
                     EMK_StationsMod.o sharedMod.o
	$(FC) -c $(FFLAGS) $< -I$(ESMF_MOD)

EMK_GridHashMod.o : EMK_GridHashMod.F90
	$(FC) -c $(FFLAGS) $<

EMK_ReportsMod.o: EMK_ReportsMod.F90
	$(FC) -c $(FFLAGS) $<

EMK_StationsMod.o: EMK_StationsMod.F90
	$(FC) -c $(FFLAGS) $<

sharedMod.o: sharedMod.F90
	$(FC) -c $(FFLAGS) $< -I$(ESMF_MOD)


.PHONY: clean
clean:
	rm -f *.o *.mod procOBA_NWP procOBA_Sat
