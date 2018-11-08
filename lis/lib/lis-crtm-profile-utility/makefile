#================================================================================
#
# Makefile for Profile Utility library build
#
# The default build target will work for most workstation platforms. For linux
# systems where there are several compiler choices, specific compiler targets
# are also provided. To modify the default Linux compiler, edit the make.macros
# file definition of LINUX_FLAGS.
#
# Note there are production and debug targets. Default build is production.
#================================================================================

# -------------
# Define macros
# -------------
include make.macros

# The library name
LIBRARY_TAG = Profile_Utility
LIBRARY_NAME = lib$(LIBRARY_TAG).a

# The library directories
SOURCE_DIR = src
TEST_DIR = test
LIBRARY_DIR = lib
INCLUDE_DIR = include


# Default compiler build
# ----------------------
all:
	cd $(SOURCE_DIR); make; cd ..

# Install the library and include files in their respective directories
# ---------------------------------------------------------------------
install: install_lib install_include

install_lib:
	@if [ ! -d $(LIBRARY_DIR) ]; then \
	  mkdir $(LIBRARY_DIR); \
	fi; \
	cd $(SOURCE_DIR); \
	if [ -f $(LIBRARY_NAME) ]; then \
	  $(MOVE) $(LIBRARY_NAME) ../$(LIBRARY_DIR); \
	fi ; \
	cd ..

install_include:
	@if [ ! -d $(INCLUDE_DIR) ]; then \
	  mkdir $(INCLUDE_DIR); \
	fi; \
	cd $(SOURCE_DIR); \
	for FILE in `ls -1 | egrep "\.(mod|MOD|m|M)$$"`; do \
	  if [ -f $$FILE ]; then \
	    $(COPY) -f $$FILE ../$(INCLUDE_DIR); \
	  fi ; \
	done ; \
	cd ..


# Test code build
TEST_CASES = Units_Conversion Level_Layer_Conversion Geopotential Atmospheric_Properties
test: all install
	@cd $(TEST_DIR) ; \
	for test_case in ${TEST_CASES}; do \
	  echo ; \
	  echo "  Building test case: $$test_case" ; \
	  echo "  -------------------" ; \
	  cd $$test_case ; \
	  make run ; \
	  cd .. ; \
	done ; \
	cd ..


# --------
# Clean up
# --------
# All the build products
clean: clean_test
	@cd $(SOURCE_DIR); make clean; cd ..

clean_test:
	@cd $(TEST_DIR) ; \
	for test_case in ${TEST_CASES}; do \
	  cd $$test_case ; \
	  make clean ; \
	  cd .. ; \
	done ; \
	cd ..

# Everything
distclean: clean clean_lib clean_include

clean_lib:
	$(REMOVE)r $(LIBRARY_DIR)

clean_include:
	$(REMOVE)r $(INCLUDE_DIR)

