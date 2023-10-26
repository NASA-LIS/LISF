# Set default goal.
.DEFAULT_GOAL := LIS

# Define shell commands
RM := rm -f

# Define functions
make-clean := $(RM) *.o *.mod LIS
make-real-clean := $(RM) *.o *.mod LIS *.d Filepath
COUPLED:=$(shell grep -E '^\#[[:space:]]*define[[:space:]]+COUPLED' LIS_misc.h)
append_offline := $(if $(COUPLED),,vim -X -N -u NONE -i NONE -e -s -c '$$s;$$; ../offline;|wq' Filepath)

# Ensure that Filepath exists and/or is up to date.
#
# Do not include Filepath when the goal is either
# clean, realclean, or help.
require-filepath := $(if $(filter clean realclean help,$(MAKECMDGOALS)),false,true)
user.cfg:
	touch $@
Filepath: default.cfg user.cfg LIS_misc.h
	@echo Generating Filepath
	$(call make-real-clean)
	./plugins.py
	$(call append_offline)
ifeq "$(require-filepath)" "true"
   include ./Filepath
endif

# Expand any tildes in directory names. Change spaces to colons.
space := $(subst ,, )
VPATH := $(foreach dir,$(dirs),$(wildcard $(dir)))
VPATH := $(subst $(space),:,$(VPATH))

# Find header files and build an include directory list.
FIND_HEADERS = $(wildcard $(dir)/*.h)
HEADER_FILES = $(foreach dir, $(dirs),$(FIND_HEADERS))
HEADER_DIRS1 := $(dir $(HEADER_FILES))
uniq = $(if $1,$(firstword $1) $(call uniq,$(filter-out $(firstword $1),$1)))
HEADER_DIRS2 = $(call uniq,$(HEADER_DIRS1))
HEADER_DIRS = $(foreach dir,$(HEADER_DIRS2),-I$(dir))

# Get list of files and determine objects and dependency files
FIND_FILES = $(wildcard $(dir)/*.F $(dir)/*.f $(dir)/*.F90 $(dir)/*.f90 $(dir)/*.c $(dir)/*.cc)
FILES      = $(foreach dir, $(dirs),$(FIND_FILES))
SOURCES   := $(sort $(notdir $(FILES)))
DEPS      := $(addsuffix .d, $(basename $(SOURCES)))
OBJS      := $(addsuffix .o, $(basename $(SOURCES)))
DOCS      := $(addsuffix .tex, $(basename $(SOURCES)))
FILES     := $(subst $(space),:,$(strip $(FILES)))

# Define the dependency generator.
#
# To override LIS_MAKEDEP_FLAGS below, set a LIS_MAKEDEP_FLAGS environment
# variable.
LIS_MAKEDEP_FLAGS ?= --verbose status
DEPGEN := ./makedep.py
DEPGEN_FLAGS := $(LIS_MAKEDEP_FLAGS)
DEPGEN_FLAGS += --exclude-fmod-files esmf,netcdf,grib_api,hdf5,mpi
DEPGEN_FLAGS += --exclude-finc-files mpif.h,netcdf.inc,CXMLDEF.FOR,mkl_vml.f90

include ./configure.lis

# Define goals.
LIS: $(OBJS)
	 $(FC) -o $@ $(OBJS) $(FOPTS) $(LDFLAGS)

LIBTARGET    =  lislib
TARGETDIR    =  ./
$(LIBTARGET) : 
		gmake -i -r lis_contrib ; \
		$(AR) lislib.a $(OBJS)  ; \

LIBTARGET    =  explis
TARGETDIR    =  ./
$(LIBTARGET) : 
		gmake -i -r lis_contrib                  ; \
		$(AR) -ru ../../main/libwrflib.a $(OBJS) ; \

LIBTARGET    =  gcelis
TARGETDIR    =  ./
$(LIBTARGET) : 
		gmake -i -r lis_contrib ; \
		$(AR) lislib.a $(OBJS)  ; \

lis_contrib: $(OBJS)

.PHONY: debug
debug:
	@echo "FILES: |$(FILES)|"
	@echo "FFLAGS: $(FFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "OBJS: $(OBJS)"
	@echo "VPATH: |$(VPATH)|"

.PHONY: depend
depend:
	@$(DEPGEN) $(DEPGEN_FLAGS) $(FILES) $(VPATH)

.PHONY: clean
clean:
	$(call make-clean)

.PHONY: realclean
realclean:
	$(call make-real-clean)


# Define pattern rules.
.SUFFIXES:

%.d: %.F
	@$(DEPGEN) $(DEPGEN_FLAGS) $< $(VPATH)

%.d: %.f
	@$(DEPGEN) $(DEPGEN_FLAGS) $< $(VPATH)

%.d: %.F90
	@$(DEPGEN) $(DEPGEN_FLAGS) $< $(VPATH)

%.d: %.f90
	@$(DEPGEN) $(DEPGEN_FLAGS) $< $(VPATH)

%.d: %.c
	@$(DEPGEN) $(DEPGEN_FLAGS) $< $(VPATH)

%.d: %.cc
	@$(DEPGEN) $(DEPGEN_FLAGS) $< $(VPATH)

# Additional prerequisite rules not found by makedep.py.
# NONE

ifeq ($(LIS_ARCH),linux_ifc)
include ./Crocus.make
endif

%.o: %.F
	$(FC77) -cpp $(FFLAGS77) $(HEADER_DIRS) $<

%.o: %.f
	$(FC77) $(FFLAGS77) $<

%.o: %.F90
	$(FC) $(FFLAGS) $(HEADER_DIRS) $<

%.o: %.f90
	$(FC) -cpp $(FFLAGS) $(HEADER_DIRS) $<

%.o: %.c
	$(CC) $(CFLAGS) $(HEADER_DIRS) $<

%.o: %.cc
	$(CC) $(CFLAGS) $(HEADER_DIRS) $<

version.o: version.F90
	$(FC) $(FFLAGS) $(HEADER_DIRS) $<
	$(RM) $<

version.F90:
	@echo "Building version file $@"
	sh build_version.sh

ifeq ($(LIS_ARCH),cray_cray)
model_interface_mod.o: model_interface_mod.F90
	$(FC) $(FFLAGS) -O0 $(HEADER_DIRS) $<
soil_htc_jls_mod.o: soil_htc_jls_mod.F90
	$(FC) $(FFLAGS) -O0 $(HEADER_DIRS) $<
endif

ifeq ($(LIS_ARCH),linux_gfortran)
# For MPI
LIS_historyMod.o: LIS_historyMod.F90
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(HEADER_DIRS) $<
LIS_DAobservationsMod.o: LIS_DAobservationsMod.F90
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(HEADER_DIRS) $<
USAF_bratsethMod.o: USAF_bratsethMod.F90
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(HEADER_DIRS) $<
gmaopert_Mod.o: gmaopert_Mod.F90
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(HEADER_DIRS) $<
HYMAP2_routingMod.o: HYMAP2_routingMod.F90
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(HEADER_DIRS) $<
read_SMOSNRTNNL2sm.o: read_SMOSNRTNNL2sm.F90
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(HEADER_DIRS) $<
HYMAP2_model.o: HYMAP2_model.F90
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(HEADER_DIRS) $<
# For HDF-EOS2
read_GLASSlai.o: read_GLASSlai.F90
	$(FC) $(FFLAGS) -fallow-argument-mismatch $(HEADER_DIRS) $<
endif

.PHONY: help
help:
	@echo ""
	@echo "usage: make <target>"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "clean"
	@echo "   $(make-clean)"
	@echo ""
	@echo "realclean"
	@echo "   $(make-real-clean)"
	@echo ""
	@echo "depend"
	@echo "   generate all dependencies at once (instead of file by file)"
	@echo ""
	@echo "LIS"
	@echo "   build LIS -- DEFAULT target"
	@echo ""
	@echo "lislib"
	@echo "   build lislib for ?"
	@echo ""
	@echo "explis"
	@echo "   build explis for coupling LIS with NU-WRF"
	@echo ""
	@echo "gcelis"
	@echo "   build gcelis for coupling LIS with GCE"
	@echo ""
	@echo "debug"
	@echo "   print some debugging information"
	@echo ""


# Do not include DEPS when the goal is either
# debug, depend, clean, realclean, or help.
require-deps := $(if $(filter debug depend clean realclean help,$(MAKECMDGOALS)),false,true)
ifeq "$(require-deps)" "true"
-include $(DEPS)
endif
