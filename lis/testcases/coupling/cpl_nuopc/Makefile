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

HR      := ========================================
HR      := $(HR)$(HR)
COMMA   := ,
DIR     := $(CURDIR)
CAP_DIR ?= $(abspath $(DIR)/../../../runmodes/nuopc_cpl_mode/)
SRC_DIR ?= $(abspath $(DIR)/src)
EXP_DIR ?= $(abspath $(DIR)/example)
BLD_DIR ?= $(abspath $(DIR)/build)
RUN_DIR ?= $(abspath $(DIR)/run)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: all build distclean clean dust

all: build setup

build:
	@echo $(HR)
	@echo "Building"
	make -C $(CAP_DIR) nuopcinstall INSTPATH=$(BLD_DIR) BUILD_TYPE=Debug COMPILER=$(COMPILER)
	make -C $(SRC_DIR) COMP1_MK=$(BLD_DIR)/lis.mk BUILD_DIR=$(BLD_DIR) RUN_DIR=$(RUN_DIR) COMPILER=$(COMPILER)

setup:
	@echo $(HR)
	@echo "Run Setup"
	@echo
	@mkdir -p $(RUN_DIR)
	@cp -R $(EXP_DIR)/* $(RUN_DIR)

distclean:
	@echo $(HR)
	@echo "Dist Cleaning"
	make -C $(SRC_DIR) clean
	make -C $(CAP_DIR) nuopcdistclean
	rm -rf $(BLD_DIR) $(RUN_DIR)

clean:
	@echo $(HR)
	@echo "Cleaning"
	make -C $(SRC_DIR) clean
	make -C $(CAP_DIR) nuopcclean
	rm -rf $(BLD_DIR)

dust:
	@echo $(HR)
	@echo "Deleting Output Files"
	rm -f  $(RUN_DIR)/*.nc
	rm -f  $(RUN_DIR)/test.log
	rm -f  $(RUN_DIR)/core*
	rm -f  $(RUN_DIR)/PET*.ESMF_LogFile
	rm -rf $(RUN_DIR)/LIS_LOGS
	rm -rf $(RUN_DIR)/LIS_OUTPUT
	rm -rf $(RUN_DIR)/LND_CAP_OUTPUT

