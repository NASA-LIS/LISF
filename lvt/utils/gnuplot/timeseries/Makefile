# Set up special characters
include $(LVT_SRC)/make/configure.lvt

LDFLAGS := -L$(LIB_MPI) -lmpi -L$(LIB_ESMF) -lesmf -L$(LIB_NETCDF) -lnetcdf -lstdc++ -lrt
#AIX: LDFLAGS := -q64 -L$(LIB_ESMF) -lesmf -lrt -lC_r -lmpiuni
#LDFLAGS := -q64 -L$(LIB_ESMF) -lesmf -lnetcdf_stubs -lC_r

null  :=
space := $(null) $(null)

# Check for directory in which to put executable
ifeq ($(MODEL_EXEDIR),$(null))
MODEL_EXEDIR := .
endif

# Check for name of executable
ifeq ($(EXENAME),$(null))
EXENAME := plotTS
endif

# Determine platform
UNAMES := $(shell uname -s)
UMACHINE := $(shell uname -m)

# Load dependency search path.
dirs := . 

#------------------------------------------------------------------------
# Primary target: build the model
#------------------------------------------------------------------------
all: $(MODEL_EXEDIR)/$(EXENAME)

# Get list of files and determine objects and dependency files
FIND_FILES = $(wildcard $(dir)/*.F $(dir)/*.f $(dir)/*.F90 $(dir)/*.c $(dir)/*.f)
FILES      = $(foreach dir, $(dirs),$(FIND_FILES))
SOURCES   := $(sort $(notdir $(FILES)))
DEPS      := $(addsuffix .d, $(basename $(SOURCES)))
OBJS      := $(addsuffix .o, $(basename $(SOURCES)))
DOCS      := $(addsuffix .tex, $(basename $(SOURCES)))

$(MODEL_EXEDIR)/$(EXENAME): $(OBJS)
	 $(FC) -o $@ $(OBJS) $(FOPTS) $(LDFLAGS)
debug: $(OBJS)
	echo "FFLAGS: $(FFLAGS)"
	echo "LDFLAGS: $(LDFLAGS)"
	echo "OBJS: $(OBJS)"

#***********************************************************************
#********** Architecture-specific flags and rules***********************
#***********************************************************************

.SUFFIXES:
.SUFFIXES: .f .f90 .F90 .c .o

.f.o:
	$(FC77) $(FFLAGS77) $<
.F90.o:
	$(FC) $(FFLAGS) $<
.c.o:
	cc -c $(cpp_path) $(CFLAGS) $<
.f90.o:
	$(FC) $(FFLAGS) $<


# Add user defined compiler flags if set, and replace FC if USER option set.
FFLAGS  += $(USER_FFLAGS)
ifneq ($(USER_FC),$(null))
FC := $(USER_FC)
endif

clean:
	$(RM) -f *.o *.mod *.stb  $(MODEL_EXEDIR)/$(EXENAME)
realclean:
	$(RM) -f *.o *.d *.mod *.stb  $(MODEL_EXEDIR)/$(EXENAME)

#------------------------------------------------------------------------
#!!!!!!!!!!!!!!!!DO NOT EDIT BELOW THIS LINE.!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#------------------------------------------------------------------------
# These rules cause a dependency file to be generated for each source
# file.  It is assumed that the tool "makdep" (provided with this
# distribution in clm2/tools/makdep) has been built and is available in
# the user's $PATH.  Files contained in the clm2 distribution are the
# only files which are considered in generating each dependency.  The
# following filters are applied to exclude any files which are not in
# the distribution (e.g. system header files like stdio.h).
#
#  1) Remove full paths from dependencies. This means gnumake will not break
#     if new versions of files are created in the directory hierarchy
#     specified by VPATH.
#
#  2) Because of 1) above, remove any file dependencies for files not in the
#     clm2 source distribution.
#
# Finally, add the dependency file as a target of the dependency rules.  This
# is done so that the dependency file will automatically be regenerated
# when necessary.
#
#     i.e. change rule
#       make.o : make.c make.h
#       to:
#       make.o make.d : make.c make.h
#------------------------------------------------------------------------
DEPGEN := $(LIS_SRC)/make/MAKDEP/makdep -s F
%.d : %.c
	@echo "Building dependency file $@"
	@$(DEPGEN) -f $(cpp_path) $< > $@
%.d : %.f
	@echo "Building dependency file $@"
	@$(DEPGEN) -f $(cpp_path) $< > $@
%.d : %.F90
	@echo "Building dependency file $@"
	@$(DEPGEN) -f $(cpp_path) $<  > $@
%.d : %.F
	@echo "Building dependency file $@"
	@$(DEPGEN) -f $(cpp_path) $< > $@
#
# if goal is clean or realclean then don't include .d files
# without this is a hack, missing dependency files will be created
# and then deleted as part of the cleaning process
#
INCLUDE_DEPS=TRUE
ifeq ($(MAKECMDGOALS), realclean)
 INCLUDE_DEPS=FALSE
endif
ifeq ($(MAKECMDGOALS), clean)
 INCLUDE_DEPS=FALSE
endif

ifeq ($(INCLUDE_DEPS), TRUE)
-include $(DEPS)
endif







