#
# --- LIS makefile 
# --- Requires ESMVv7+
# --- LIS ESMF component.
#

# ###########################
# Include ESMFMKFILE fragment
# ###########################

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif
include $(ESMFMKFILE)

# ##############
# Local Variable
# ##############
HR="================================================================================"

ifndef LISDIR
LISDIR=../../.
endif

LIS_DIR=$(LISDIR)
LIS_TDIR=$(notdir $(LIS_DIR))
LIS_MODDIR=make
LIS_LIBDIR=make
LIS_OBJDIR=make
LIS_MKDIR=make

include $(LIS_DIR)/$(LIS_MKDIR)/configure.lis

LIS_LIB=libLIS.a
LIS_MODS=lis_coremod.mod \
	lis_daobservationsmod.mod \
	lis_dataassimmod.mod \
	lis_domainmod.mod \
	lis_forc_attributesmod.mod \
	lis_logmod.mod \
	lis_metforcingmod.mod \
	lis_paramsmod.mod \
	lis_surfacemodelmod.mod \
	lis_tbotadjustmod.mod \
	lis_timemgrmod.mod \
	liswrfgridcompmod.mod \
	map_utils.mod
LIS_FILES=$(addprefix $(LIS_DIR)/$(LIS_MODDIR)/,$(LIS_MODS)) \
	$(LIS_LIB)

LISNUOPC_DIR=$(shell pwd)
LISNUOPC_TDIR=$(notdir $(LISNUOPC_DIR))

LISNUOPC_LIB=liblis_nuopc.a
LISNUOPC_MK=lis.mk
LISNUOPC_DEP_FRONT=lis_nuopc
LISNUOPC_MODS=lis_nuopc.mod
LISNUOPC_FILES= $(LISNUOPC_MK) \
	$(LISNUOPC_LIB) \
	$(LISNUOPC_MODS)

# ########
# Make all
# ########
all: build_lis_nuopc

# ###########
# LIS Library
# ###########
#$(LIS_DIR)/$(LIS_LIBDIR)/$(LIS_LIB): 
#	ar rc $@ $(LIS_DIR)/$(LIS_OBJDIR)/*.o

$(LIS_LIB): 
	ar rc $@ $(LIS_DIR)/$(LIS_OBJDIR)/*.o

# ##############
# Build Settings
# ##############

.SUFFIXES: 
.SUFFIXES: .c .C .f90 .F90 .F .f

#LIS_F90COMPILEPATHS =-I$(LISHYDRO_DIR)/$(LISHYDRO_MODDIR)
LIS_F90COMPILEPATHS=-I$(LIS_DIR)/$(LIS_MODDIR)

.F:
	@echo "Must have an explicit rule for" $*
.f:
	@echo "Must have an explicit rule for" $*
.C:
	@echo "Must have an explicit rule for" $*
.c: 
	@echo "Must have an explicit rule for" $*

%.o : %.f90
	@echo $(HR)
	@echo "Compiling $@..."
	@echo
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(LIS_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) $<

%.o : %.F90
	@echo $(HR)
	@echo "Compiling $@..."
	@echo
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(LIS_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

# #####################################
# Generate Local ESMF Makefile fragment
# #####################################
%.mk :
	@echo $(HR)	
	@echo "Generating LIS Makefile fragment $@"
	@echo
	@echo "# ESMF self-describing build dependency makefile fragment" > $@
	@echo >> $@
	@echo "ESMF_DEP_FRONT     = $(LISNUOPC_DEP_FRONT)"                  >> $@
	@echo "ESMF_DEP_INCPATH   = $(LISNUOPC_DIR)"                  >> $@
	@echo "ESMF_DEP_CMPL_OBJS = $(LISNUOPC_DIR)/$(LISNUOPC_LIB)" >> $@
	@echo "ESMF_DEP_LINK_OBJS = $(LISNUOPC_DIR)/$(LISNUOPC_LIB)" >> $@
	@echo "ESMF_DEP_SHRD_PATH = $(LIB_GRIBAPI) $(LIB_JASPER)" >> $@
	@echo "ESMF_DEP_SHRD_LIBS = grib_api_f90 grib_api jasper" >> $@

# #####################
# Build NUOPC Component
# #####################

.PHONY: lis_nuopc clean nuopcinstall

build_lis_nuopc: $(LISNUOPC_FILES)

$(LISNUOPC_LIB): $(LIS_LIB) lis_nuopc_gluecode.o \
	NUOPC_FileReadUtility.o NUOPC_FileWriteUtility.o \
	NUOPC_FillUtility.o NUOPC_CopyUtility.o \
	NUOPC_LogUtility.o lis_nuopc_cap.o 
	@echo $(HR)
	@echo "Copying static library $@..."
	@echo
	mv -f $(LIS_LIB) $@
	ar rc $@ *.o
	
# Dependencies
lis_nuopc.mod: lis_nuopc_cap.o
lis_nuopc_gluecode.mod: lis_nuopc_gluecode.o
nuopc_filereadutility.mod: NUOPC_FileReadUtility.o
nuopc_filewriteutility.mod: NUOPC_FileWriteUtility.o
nuopc_fillutility.mod: NUOPC_FillUtility.o
nuopc_copyutility.mod: NUOPC_CopyUtility.o
nuopc_logutility.mod: NUOPC_LogUtility.o

lis_nuopc_cap.o: lis_nuopc_gluecode.mod \
	nuopc_filereadutility.mod nuopc_filewriteutility.mod \
	nuopc_fillutility.mod nuopc_copyutility.mod nuopc_logutility.mod
lis_nuopc_gluecode.o: nuopc_filereadutility.mod nuopc_filewriteutility.mod \
	nuopc_fillutility.mod  nuopc_copyutility.mod \
	nuopc_logutility.mod $(LIS_FILES)

# ###########
# Check Build
# ###########
define checkfile
	@if [ ! -e $(1) ]; then \
	echo "File is missing:$(1)"; \
	exit 1; fi;

endef # blank line in checkfile is required

define checkdir
	@if [ ! -d $(1) ]; then \
	echo "Directory is missing:$(1)"; \
	exit 1; fi;
endef # blank line in checkdir is required

check: check_esmf check_lis_nuopc

# ##################
# Check ESMF Version
# ##################

check_esmf:
	@echo $(HR)
	@echo "Checking ESMFMKFILE file..."
	@echo
	@echo "ESMFMKFILE=$(ESMFMKFILE)"
	@if [ "$(ESMF_VERSION_MAJOR)" -lt 7 ]; then \
	echo "Please use ESMF version 7+"; \
	exit 1; fi;
	@echo "ESMF Version=$(ESMF_VERSION_STRING)"

# #####################
# Check LIS NUOPC
# #####################
check_lis_nuopc: 
	@echo $(HR)
	@echo "Checking for LIS NUOPC files..."
	@echo
	$(foreach FILENAME, $(LISNUOPC_FILES), $(call checkfile, $(FILENAME)))

# #########
# Clean all
# #########
clean: clean_lis_nuopc

# #########################
# Clean LIS NUOPC Component
# #########################
clean_lis_nuopc:
	@echo $(HR)
	@echo "Cleaning LIS NUOPC build..."
	@echo
	rm -f *.o *.mod *.a *.mk VERSION

# #########
# Clean LIS
# #########
clean_lis:
	@echo $(HR)
	@echo "Cleaning LIS build..."
	@echo ""
	$(call checkdir, $(LIS_DIR))
	(cd $(LIS_DIR)/$(LIS_MKDIR); make clean)

# ###########################
# Install LIS NUOPC Component
# ###########################

ifndef DESTDIR
DESTDIR := $(LISNUOPC_DIR)
endif

INSTDATE := $(shell date '+%Y-%m-%d-%H-%M-%S')
ifndef INSTDIR
INSTDIR  := LIS_$(INSTDATE)
endif

install:
	@mkdir -p $(DESTDIR)/$(INSTDIR)
	@cp $(LISNUOPC_FILES) $(DESTDIR)/$(INSTDIR)
	@echo "# ESMF self-describing build dependency makefile fragment" > $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_MK)
	@echo >> $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_MK)
	@echo "ESMF_DEP_FRONT     = $(LISNUOPC_DEP_FRONT)" >> $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_MK)
	@echo "ESMF_DEP_INCPATH   = $(DESTDIR)/$(INSTDIR)" >> $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_MK)
	@echo "ESMF_DEP_CMPL_OBJS = $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_LIB)" >> $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_MK)
	@echo "ESMF_DEP_LINK_OBJS = $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_LIB)" >> $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_MK)
	@echo "ESMF_DEP_SHRD_PATH = $(LIB_GRIBAPI) $(LIB_JASPER)" >> $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_MK)
	@echo "ESMF_DEP_SHRD_LIBS = grib_api_f90 grib_api jasper" >> $(DESTDIR)/$(INSTDIR)/$(LISNUOPC_MK)

	@touch VERSION
	@if [ -d .svn ]; then \
	  echo "SVN Repository" > VERSION; \
	  svn info . | grep URL >> VERSION; \
	  svn info . | grep "Last Changed Rev" >> VERSION; \
	fi
	@if [ -d .git ]; then \
	  echo "Git Repository" > VERSION; \
	  git show . | grep -m 1 "commit " >> VERSION; \
	  git show . | grep -m 1 "Author: " >> VERSION; \
	  git show . | grep -m 1 "Date: " >> VERSION; \
	fi
	@cp VERSION $(DESTDIR)/$(INSTDIR)/
	@echo Installation into \"$(DESTDIR)/$(INSTDIR)\" complete!
	@echo
# ------------------------------------------------------------------------------
