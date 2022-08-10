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
else ifeq ($(BUILD_PREC),r8)
override ESMF_F90COMPILECPPFLAGS += -DREAL8
#else
#override ESMF_F90COMPILECPPFLAGS += -DREAL8
endif

# #################################
# Compile with Debugging Directives
# #################################

ifneq ($(filter $(BUILD_TYPE),DEBUG Debug debug),)
override ESMF_F90COMPILECPPFLAGS += -DDEBUG
override ESMF_CXXCOMPILECPPFLAGS += -DDEBUG
ifneq ($(filter $(COMPILER),INTEL Intel intel),)
ESMF_F90COMPILEOPTS += -O0 -traceback -check arg_temp_created,bounds,format,output_conversion,stack,uninit
endif
ifneq ($(filter $(COMPILER),GNU Gnu gnu),)
ESMF_F90COMPILEOPTS += -O0 -Wall -Wextra -Wconversion -Wno-unused -Wno-unused-dummy-argument -fbacktrace -fimplicit-none -fcheck=all,no-pointer
endif
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

ifndef INSTPATH
INSTPATH := $(abspath $(DESTDIR)/$(INSTDIR))
endif

# ###############
# Model Variables
# ###############

ifndef LIS_DIR
MODEL_RDIR   := $(abspath $(DIR)/../../..)
else
MODEL_RDIR   := $(abspath $(LIS_DIR))
endif
MODEL_DIR    := $(abspath $(MODEL_RDIR)/lis)
MODEL_OBJDIR := $(abspath $(MODEL_DIR)/make)
MODEL_MODDIR := $(abspath $(MODEL_DIR)/make)
MODEL_LIBDIR := $(abspath $(MODEL_DIR))
MODEL_LIB    := $(abspath $(MODEL_LIBDIR)/liblis.a)
MODEL_MKDIR  := $(abspath $(MODEL_DIR)/make)
MODEL_MK     := $(abspath $(MODEL_MKDIR)/Makefile)
MODEL_MKINC  := $(abspath $(MODEL_MKDIR)/configure.lis)

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
CAP_VERS      := VERSION

CAP_OBJS      := LIS_NUOPC_Cap.o
CAP_OBJS      += LIS_NUOPC_Gluecode.o
CAP_OBJS      += LIS_NUOPC_DataCopy.o
CAP_OBJS      += LIS_ESMF_Extensions.o

CAP_MODS      := lis_nuopc.mod
CAP_MODS      += lis_nuopc_gluecode.mod
CAP_MODS      += lis_nuopc_datacopy.mod
CAP_MODS      += lis_esmf_extensions.mod

CAP_FILES     := $(CAP_OBJS) $(CAP_MODS) $(CAP_LIB) $(CAP_VERS) $(CAP_MK)

# ###############################
# Include Model Makefile Fragment
# ###############################

include $(MODEL_MKINC)
override ESMF_F90COMPILEPATHS += -I$(MODEL_MODDIR)
override DEP_SHRD_PATH          = $(patsubst -L%,%,$(LIS_LIB_PATHS))
override DEP_SHRD_LIBS          = $(patsubst -l%,%,$(LIS_LIB_FLAGS))

# #################################
# Compile with WRF_HYDRO directives
# #################################

ifneq ($(WRF_HYDRO),)
ifneq ($(WRF_HYDRO),0)
override ESMF_F90COMPILECPPFLAGS += -DWRF_HYDRO
override ESMF_CXXCOMPILECPPFLAGS += -DWRF_HYDRO
endif
endif

# ##################################
# Compile with GSM_EXTLND directives
# ##################################
ifneq ($(GSM_EXTLND),)
ifneq ($(GSM_EXTLND),0)
override ESMF_F90COMPILECPPFLAGS += -DGSM_EXTLND
override ESMF_CXXCOMPILECPPFLAGS += -DGSM_EXTLND
endif
endif

# ##################################
# Compile with ESMF_TRACE directives
# ##################################
ifneq ($(ESMF_TRACE),)
ifneq ($(ESMF_TRACE),0)
override ESMF_F90COMPILECPPFLAGS += -DESMF_TRACE
override ESMF_CXXCOMPILECPPFLAGS += -DESMF_TRACE
endif
endif

# #######################
# Primary Makefile Target
# #######################
.PHONY: nuopc nuopcinstall nuopcdistclean nuopcclean install_mk

nuopc: $(CAP_FILES)

nuopcinstall: $(CAP_LIB) $(CAP_MODS) $(CAP_VERS) \
 $(addprefix $(INSTPATH)/,$(CAP_MODS)) \
 $(addprefix $(INSTPATH)/,$(CAP_LIB)) \
 $(addprefix $(INSTPATH)/,$(CAP_VERS)) \
 install_mk

# ############
# Dependencies
# ############

LIS_NUOPC_Cap.o: LIS_NUOPC_Macros.h LIS_NUOPC_Gluecode.o \
        LIS_ESMF_Extensions.o
LIS_NUOPC_Gluecode.o: LIS_NUOPC_Macros.h LIS_NUOPC_DataCopy.o \
        LIS_ESMF_Extensions.o $(MODEL_MODS)
LIS_NUOPC_DataCopy.o: LIS_NUOPC_Macros.h $(MODEL_MODS)

lis_nuopc.mod: LIS_NUOPC_Cap.o
lis_nuopc_gluecode.mod: LIS_NUOPC_Gluecode.o
lis_nuopc_datacopy.mod: LIS_NUOPC_DataCopy.o
lis_esmf_extensions.mod: LIS_ESMF_Extensions.o

# ###########
# Build model
# ###########

build_model:
	@echo $(HR)
	@echo "Building Model..."
	@echo ""
	$(call checkdir, $(MODEL_DIR))
	( cd $(MODEL_DIR); ./compile )

$(MODEL_MODS): build_model

$(MODEL_LIB): build_model
	@echo $(HR)
	@echo "Building Model Library..."
	@echo
	ar cr $@ $(MODEL_OBJDIR)/*.o

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
	$(call checkfile, $(MODEL_LIB))
	cp $(MODEL_LIB) $@
	ar cr $@ $(CAP_OBJS)

$(CAP_VERS): $(CAP_LIB) $(CAP_MODS)
	@echo $(HR)
	@echo "Generating Version Information"
	@echo
ifneq (,$(wildcard $(MODEL_RDIR)/.git))
	@echo "# NUOPC Cap Revision #" > $(CAP_VERS)
	@git log . | grep -m 1 "commit " >> $(CAP_VERS)
	@git log . | grep -m 1 "Author: " >> $(CAP_VERS)
	@git log . | grep -m 1 "Date: " >> $(CAP_VERS)
	@echo >> $(CAP_VERS)
	@echo "# Model Revision     #" >> $(CAP_VERS)
	@git log $(MODEL_RDIR) | grep -m 1 "commit " >> $(CAP_VERS)
	@git log $(MODEL_RDIR) | grep -m 1 "Author: " >> $(CAP_VERS)
	@git log $(MODEL_RDIR) | grep -m 1 "Date: " >> $(CAP_VERS)
	@echo >> $(CAP_VERS)
else ifneq (,$(wildcard $(MODEL_RDIR)/.svn))
	@echo "# NUOPC Cap Revision #" > $(CAP_VERS)
	@svn info . | grep URL >> $(CAP_VERS)
	@svn info . | grep "Last Changed Rev" >> $(CAP_VERS)
	@echo >> $(CAP_VERS)
	@echo "# Model Revision     #" >> $(CAP_VERS)
	@svn info $(MODEL_RDIR) | grep URL >> $(CAP_VERS)
	@svn info $(MODEL_RDIR) | grep "Last Changed Rev" >> $(CAP_VERS)
	@echo >> $(CAP_VERS)
else
	@echo "# Version Information Not Found" > $(CAP_VERS)
endif

$(CAP_MK): $(CAP_LIB) $(CAP_MODS)
	@echo $(HR)
	@echo "Generating NUOPC Makefile Fragment"
	@echo
	@echo "# ESMF self-describing build dependency makefile fragment" > $(CAP_MK)
	@echo "" >> $(CAP_MK)
	@echo "ESMF_DEP_FRONT     = $(CAP_DEP_FRONT)" >> $(CAP_MK)
	@echo "ESMF_DEP_INCPATH   = $(CAP_DIR)" >> $(CAP_MK)
	@echo "ESMF_DEP_CMPL_OBJS = " >> $(CAP_MK)
	@echo "ESMF_DEP_LINK_OBJS = $(CAP_DIR)/$(CAP_LIB)" >> $(CAP_MK)
	@echo "ESMF_DEP_SHRD_PATH = $(DEP_SHRD_PATH)" >> $(CAP_MK)
	@echo "ESMF_DEP_SHRD_LIBS = $(DEP_SHRD_LIBS)" >> $(CAP_MK)

# -----------------------------------------------------------------------------
# Install Library, Modules, and Makefile Fragment
# -----------------------------------------------------------------------------

$(INSTPATH)/%: %
	@echo $(HR)
	@echo "Installing $(notdir $@)"
	@echo
	@mkdir -p $(INSTPATH)
	@cp $(notdir $@) $@

install_mk:
	@echo $(HR)
	@echo "Installing NUOPC Makefile Fragment"
	@echo
	@mkdir -p $(INSTPATH)
	@echo "# ESMF self-describing build dependency makefile fragment" > $(INSTPATH)/$(CAP_MK)
	@echo "" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_FRONT     = $(CAP_DEP_FRONT)" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_INCPATH   = $(INSTPATH)" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_CMPL_OBJS = " >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_LINK_OBJS = $(INSTPATH)/$(CAP_LIB) $(DEP_LINK_OBJS)" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_SHRD_PATH = $(DEP_SHRD_PATH)" >> $(INSTPATH)/$(CAP_MK)
	@echo "ESMF_DEP_SHRD_LIBS = $(DEP_SHRD_LIBS)" >> $(INSTPATH)/$(CAP_MK)

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

# ###################
# Clean Cap and Model
# ###################

nuopcdistclean: nuopcclean
	@echo $(HR)
	@echo "Cleaning Model build..."
	@echo ""
	$(call checkdir, $(MODEL_MKDIR))
	make -C $(MODEL_MKDIR) -f $(MODEL_MK) clean
	rm -f $(MODEL_LIB)

# ##########
# Clean  Cap
# ##########

nuopcclean:
	@echo $(HR)
	@echo "Cleaning Cap build..."
	@echo
	rm -f $(CAP_FILES)

# ------------------------------------------------------------------------------
