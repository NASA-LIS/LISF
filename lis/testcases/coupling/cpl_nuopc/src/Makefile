# GNU Makefile template for user ESMF application

################################################################################
################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below. Following the ESMF User's Guide, a complete ESMF     ##
## installation should ensure that a single environment variable "ESMFMKFILE" ##
## is made available on the system. This variable should point to the         ##
## "esmf.mk" file.                                                            ##
##                                                                            ##
## This example Makefile uses the "ESMFMKFILE" environment variable.          ##
##                                                                            ##
## If you notice that this Makefile cannot find variable ESMFMKFILE then      ##
## please contact the person responsible for the ESMF installation on your    ##
## system.                                                                    ##
## As a work-around you can simply hardcode the path to "esmf.mk" in the      ##
## include line below. However, doing so will render this Makefile a lot less ##
## flexible and non-portable.                                                 ##
################################################################################

HR    := ========================================
HR    := $(HR)$(HR)
COMMA := ,

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)
ifneq ($(filter $(BUILD_TYPE),DEBUG Debug debug),)
ifneq ($(filter $(COMPILER),INTEL Intel intel),)
ESMF_F90COMPILEOPTS += -O0 -g -traceback -check arg_temp_created,bounds,format,output_conversion,stack,uninit
endif
ifneq ($(filter $(COMPILER),GNU Gnu gnu),)
ESMF_F90COMPILEOPTS += -O0 -g -Wall -Wextra -Wconversion -Wno-unused -Wno-unused-dummy-argument -fbacktrace -fimplicit-none -fcheck=all,no-pointer
endif
endif

# -----------------------------------------------------------------------------
# Files
# -----------------------------------------------------------------------------

TST_SRCS := $(wildcard *.F90)
TST_OBJS := $(TST_SRCS:.F90=.o)
TST_MODS := lis_nuopc_test_atm.mod
TST_MODS += lis_nuopc_test_drv.mod
TST_MODS += lis_nuopc_test_hyd.mod
TST_EXE  := lis_nuopc_test.exe

# -----------------------------------------------------------------------------
# Makefile fragments
# -----------------------------------------------------------------------------

DEP_FRONTS    :=
DEP_INCS      :=
DEP_CMPL_OBJS :=
DEP_LINK_OBJS :=
DEP_SHRD_PATH :=
DEP_SHRD_LIBS :=

ifdef COMP1_MK
  include $(COMP1_MK)
  DEP_FRONTS    += -DCOMP1_MOD=$(ESMF_DEP_FRONT)
  DEP_FRONTS    += -DCOMP1_STR=\"$(ESMF_DEP_FRONT)\"
  DEP_INCS      += $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS += $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS += $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH += $(addprefix -L, $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_PATH += $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS += $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifdef COMP2_MK
  include $(COMP2_MK)
  DEP_FRONTS    += -DCOMP2_MOD=$(ESMF_DEP_FRONT)
  DEP_FRONTS    += -DCOMP2_STR=\"$(ESMF_DEP_FRONT)\"
  DEP_INCS      += $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS += $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS += $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH += $(addprefix -L, $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_PATH += $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS += $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifdef COMP3_MK
  include $(COMP3_MK)
  DEP_FRONTS    += -DCOMP3_MOD=$(ESMF_DEP_FRONT)
  DEP_FRONTS    += -DCOMP3_STR=\"$(ESMF_DEP_FRONT)\"
  DEP_INCS      += $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS += $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS += $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH += $(addprefix -L, $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_PATH += $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS += $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifdef COMP4_MK
  include $(COMP4_MK)
  DEP_FRONTS    += -DCOMP4_MOD=$(ESMF_DEP_FRONT)
  DEP_FRONTS    += -DCOMP4_STR=\"$(ESMF_DEP_FRONT)\"
  DEP_INCS      += $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS += $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS += $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH += $(addprefix -L, $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_PATH += $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS += $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif
ifdef COMP5_MK
  include $(COMP5_MK)
  DEP_FRONTS    += -DCOMP5_MOD=$(ESMF_DEP_FRONT)
  DEP_FRONTS    += -DCOMP5_STR=\"$(ESMF_DEP_FRONT)\"
  DEP_INCS      += $(addprefix -I, $(ESMF_DEP_INCPATH))
  DEP_CMPL_OBJS += $(ESMF_DEP_CMPL_OBJS)
  DEP_LINK_OBJS += $(ESMF_DEP_LINK_OBJS)
  DEP_SHRD_PATH += $(addprefix -L, $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_PATH += $(addprefix -Wl$(COMMA)-rpath$(COMMA), $(ESMF_DEP_SHRD_PATH))
  DEP_SHRD_LIBS += $(addprefix -l, $(ESMF_DEP_SHRD_LIBS))
endif

################################################################################
################################################################################

# -----------------------------------------------------------------------------
# Compiler Arguments
# -----------------------------------------------------------------------------

.SUFFIXES: .f90 .F90 .c .C

%.o : %.f90
	@echo $(HR)
	@echo "Compiling $@"
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(DEP_FRONTS) $(DEP_INCS) $(ESMF_F90COMPILEFREENOCPP) $<

%.o : %.F90
	@echo $(HR)
	@echo "Compiling $@"
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(DEP_FRONTS) $(DEP_INCS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR) $<

%.o : %.c
	@echo $(HR)
	@echo "Compiling $@"
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(DEP_FRONTS) $(DEP_INCS) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

%.o : %.C
	@echo $(HR)
	@echo "Compiling $@"
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(DEP_FRONTS) $(DEP_INCS) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

# -----------------------------------------------------------------------------
# Build Targets
# -----------------------------------------------------------------------------

all: $(TST_EXE) install

$(TST_EXE): $(TST_OBJS) $(DEP_LINK_OBJS)
	@echo $(HR)
	@echo "Building Executable"
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) -o $@ $^ $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(DEP_SHRD_PATH) $(DEP_SHRD_LIBS) $(ESMF_F90ESMFLINKLIBS)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

lis_nuopc_test_app.o: lis_nuopc_test_drv.o
lis_nuopc_test_drv.o: lis_nuopc_test_atm.o lis_nuopc_test_hyd.o

# -----------------------------------------------------------------------------
# Install
# -----------------------------------------------------------------------------

install: $(TST_EXE) $(TST_MODS) $(TST_OBJS) \
 $(addprefix $(BUILD_DIR)/,$(TST_MODS)) \
 $(addprefix $(BUILD_DIR)/,$(TST_OBJS)) \
 $(addprefix $(RUN_DIR)/,$(TST_EXE))

$(RUN_DIR)/% $(BUILD_DIR)/%: %
	@echo $(HR)
	@echo "Installing $(notdir $@)"
	@echo
	@mkdir -p $(dir $@)
	@cp $(notdir $@) $@

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: clean install
clean:
	@echo $(HR)
	@echo "Deleting Build Files"
	rm -f $(TST_EXE) $(TST_OBJS) $(TST_MODS)
