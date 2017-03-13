#
# --- LIS makefile 
# --- Requires ESMVv7+
# --- LIS ESMF component.
#

# ###############
# Local Variables
# ###############

HR    := ========================================
HR    := $(HR)$(HR)
COMMA := ,
DIR   := $(CURDIR)

# ###########################
# Include ESMFMKFILE fragment
# ###########################

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif
include $(ESMFMKFILE)

# #########################
# Determine Build Precision
# #########################

ifeq ($(BUILD_PREC),r4)
override ESMF_F90COMPILECPPFLAGS += -DREAL4
else
override ESMF_F90COMPILECPPFLAGS += -DREAL8
endif

# ###########################
# Determine Installation Path
# ###########################

ifndef DESTDIR
DESTDIR  := $(DIR)
endif

ifndef INSTDIR
INSTDIR  := LIS_$(shell date '+%Y-%m-%d-%H-%M-%S')
endif

INSTPATH := $(abspath $(DESTDIR)/$(INSTDIR))

# ###############
# Model Variables
# ###############

MODEL_DIR    := $(abspath $(DIR)/../..)
MODEL_OBJDIR := $(abspath $(MODEL_DIR)/make)
MODEL_MODDIR := $(abspath $(MODEL_DIR)/make)
MODEL_LIBDIR := $(abspath $(MODEL_DIR))
MODEL_LIB    := $(abspath $(MODEL_LIBDIR)/liblis.a)
MODEL_MKDIR  := $(abspath $(MODEL_DIR)/make)
MODEL_MK     := $(abspath $(MODEL_MKDIR)/Makefile)
MODEL_MKINC  := $(abspath $(MODEL_MKDIR)/configure.lis)
MODEL_DEPDIR := $(abspath $(MODEL_DIR)/make/MAKDEP)
MODEL_DEPMK  := $(abspath $(MODEL_DEPDIR)/Makefile)

MODEL_MODS   := $(abspath $(MODEL_MODDIR)/lis_coremod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_daobservationsmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_dataassimmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_domainmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_forc_attributesmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_logmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_metforcingmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_paramsmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_surfacemodelmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_tbotadjustmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/lis_timemgrmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/liswrfgridcompmod.mod)
MODEL_MODS   += $(abspath $(MODEL_MODDIR)/map_utils.mod)

MODEL_FILES  := $(MODEL_LIB) $(MODEL_MODS)

# #############
# Cap Variables
# #############

CAP_DIR       := $(abspath $(DIR))
CAP_LIB       := liblis_nuopc.a
CAP_MK        := lis.mk
CAP_DEP_FRONT := lis_nuopc
CAP_VERS      := version

CAP_OBJS      := LIS_NUOPC_Cap.o
CAP_OBJS      += LIS_NUOPC_Gluecode.o
CAP_OBJS      += LIS_NUOPC_DataCopy.o
CAP_OBJS      += beta_NUOPC_FileRead.o
CAP_OBJS      += beta_NUOPC_Auxiliary.o
CAP_OBJS      += beta_NUOPC_Fill.o
CAP_OBJS      += beta_NUOPC_Log.o
CAP_OBJS      += beta_NUOPC_Base.o

CAP_MODS      := lis_nuopc.mod
CAP_MODS      += lis_nuopc_gluecode.mod
CAP_MODS      += lis_nuopc_datacopy.mod
CAP_MODS      += beta_nuopc_fileread.mod
CAP_MODS      += beta_nuopc_auxiliary.mod
CAP_MODS      += beta_nuopc_fill.mod
CAP_MODS      += beta_nuopc_log.mod
CAP_MODS      += beta_nuopc_base.mod

CAP_FILES     := $(CAP_LIB) $(CAP_OBJS) $(CAP_MODS) $(CAP_MK) $(CAP_VERS)

# ###############################
# Include Model Makefile Fragment
# ###############################

include $(MODEL_MKINC)
ESMF_F90COMPILEPATHS    += -I$(MODEL_MODDIR)

ifeq ($(DEBUG),on)
ESMF_F90COMPILEOPTS     += -g -traceback
ESMF_CXXCOMPILEOPTS     += -g -traceback
ESMF_F90COMPILECPPFLAGS += -DDEBUG
ESMF_CXXCOMPILECPPFLAGS += -DDEBUG
endif

# #######################
# Primary Makefile Target
# #######################
.PHONY: nuopc nuopcinstall nuopcclean clean_cap clean_model

nuopc: $(CAP_FILES) nuopcinstall

# ############
# Dependencies
# ############

LIS_NUOPC_Cap.o: LIS_NUOPC_Macros.h LIS_NUOPC_Gluecode.o \
        beta_NUOPC_FileRead.o beta_NUOPC_Auxiliary.o \
        beta_NUOPC_Fill.o beta_NUOPC_Log.o \
        beta_NUOPC_Base.o
LIS_NUOPC_Gluecode.o: LIS_NUOPC_Macros.h LIS_NUOPC_DataCopy.o \
        beta_NUOPC_FileRead.o beta_NUOPC_Auxiliary.o \
        beta_NUOPC_Fill.o beta_NUOPC_Log.o \
        $(MODEL_FILES)
LIS_NUOPC_DataCopy.o: LIS_NUOPC_Macros.h

# ###########
# LIS Library
# ###########

$(MODEL_LIB): $(MODEL_MODS)
	@echo $(HR)
	@echo "Building Model Library..."
	@echo
	ar rc $@ $(MODEL_OBJDIR)/*.o

$(MODEL_MODS):
	@echo "Building Model..."
	@echo
	(cd $(MODEL_MKDIR); make -f $(MODEL_DEPMK))
	(cd $(MODEL_MKDIR); make -f $(MODEL_MK))

# ##############
# Build Settings
# ##############

.SUFFIXES: 
.SUFFIXES: .c .C .f90 .F90 .F .f

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
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) $<

%.o : %.F90
	@echo $(HR)
	@echo "Compiling $@..."
	@echo
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

# #####################
# Build NUOPC Component
# #####################

$(CAP_LIB): $(MODEL_LIB) $(CAP_OBJS)
	@echo $(HR)
	@echo "Copying static library $@..."
	@echo
	ar cruT $@ $^

$(CAP_VERS):
	@echo $(HR)
	@echo "Generating Version Information"
	@echo
	@echo "# NUOPC Cap Version" > $@
	@if [ -d .svn ]; then \
	  echo "SVN Repository" > $@; \
	  svn info . | grep URL >> $@; \
	  svn info . | grep "Last Changed Rev" >> $@; \
	fi
	@if [ -d `git rev-parse --git-dir` ]; then \
	  echo "Git Repository" > $@; \
	  git show . | grep -m 1 "commit " >> $@; \
	  git show . | grep -m 1 "Author: " >> $@; \
	  git show . | grep -m 1 "Date: " >> $@; \
	fi

$(CAP_MK):
	@echo $(HR)
	@echo "Generating NUOPC Makefile Fragment"
	@echo
	@echo "# ESMF self-describing build dependency makefile fragment" > $@
	@echo "" >> $@
	@echo "ESMF_DEP_FRONT     = $(CAP_DEP_FRONT)" >> $@
	@echo "ESMF_DEP_INCPATH   = $(CAP_DIR)" >> $@
	@echo "ESMF_DEP_CMPL_OBJS = " >> $@
	@echo "ESMF_DEP_LINK_OBJS = $(CAP_DIR)/$(CAP_LIB)" >> $@
	@echo "ESMF_DEP_SHRD_PATH = " >> $@
	@echo "ESMF_DEP_SHRD_LIBS = " >> $@

# -----------------------------------------------------------------------------
# Install Library, Modules, and Makefile Fragment
# -----------------------------------------------------------------------------

nuopcinstall: $(CAP_LIB) $(CAP_MODS) $(CAP_VERS) \
 $(addprefix $(INSTPATH)/,$(CAP_MODS)) \
 $(addprefix $(INSTPATH)/,$(CAP_LIB)) \
 $(addprefix $(INSTPATH)/,$(CAP_VERS)) \
 $(addprefix $(INSTPATH)/,$(CAP_MK))

$(INSTPATH)/%:
	@echo $(HR)
	@echo "Installing $(notdir $@)"
	@echo
	@mkdir -p $(INSTPATH)
	@cp $(notdir $@) $@

$(INSTPATH)/$(CAP_MK):
	@echo $(HR)
	@echo "Installing NUOPC Makefile Fragment"
	@echo
	@mkdir -p $(INSTPATH)
	@echo "# ESMF self-describing build dependency makefile fragment" > $@
	@echo "" >> $@
	@echo "ESMF_DEP_FRONT     = $(CAP_DEP_FRONT)" >> $@
	@echo "ESMF_DEP_INCPATH   = $(INSTPATH)" >> $@
	@echo "ESMF_DEP_CMPL_OBJS = " >> $@
	@echo "ESMF_DEP_LINK_OBJS = $(INSTPATH)/$(CAP_LIB)" >> $@
	@echo "ESMF_DEP_SHRD_PATH = " >> $@
	@echo "ESMF_DEP_SHRD_LIBS = " >> $@

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

check: check_esmf check_model check_cap

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

# ###########
# Check Model
# ###########

check_model:
	@echo $(HR)
	@echo "Checking for Model files..."
	@echo
	$(foreach FILENAME, $(MODEL_FILES), $(call checkfile, $(FILENAME)))

# #########
# Check Cap
# #########

check_cap:
	@echo $(HR)
	@echo "Checking for WRF-Hydro NUOPC files..."
	@echo
	$(foreach FILENAME, $(CAP_FILES), $(call checkfile, $(FILENAME)))

# #########
# Clean all
# #########

nuopcclean: clean_cap clean_model

# ##########
# Clean  Cap
# ##########

clean_cap:
	@echo $(HR)
	@echo "Cleaning Cap build..."
	@echo
	rm -f $(CAP_FILES)

# ###########
# Clean Model
# ###########

clean_model:
	@echo $(HR)
	@echo "Cleaning Model build..."
	@echo ""
	$(call checkdir, $(MODEL_DIR))
	(cd $(MODEL_DIR); make -f $(MODEL_MK) clean)

# ------------------------------------------------------------------------------
