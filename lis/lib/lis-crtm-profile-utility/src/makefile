#==============================================================================
#
# Makefile for Profile Utility Build
#
#==============================================================================

# -------------
# Define macros
# -------------
include ../make.macros

# This makefile
MAKE_FILE = makefile

# The library name
LIBRARY_TAG = Profile_Utility
LIBRARY_NAME = lib$(LIBRARY_TAG).a

# The filelist
include make.filelist


# ------------------
# Define dependecies
# ------------------
# Default build based on OS type
all:
	@echo "OS type detected: "`uname -s`
	@case `uname -s` in \
	  "SunOS")   make -f $(MAKE_FILE) library $(SUNOS_FLAGS) ;; \
	  "AIX")     make -f $(MAKE_FILE) library $(AIX_FLAGS) ;; \
	  "IRIX64" ) make -f $(MAKE_FILE) library $(IRIX64_FLAGS) ;; \
	  "Linux"|"Darwin" )  make -f $(MAKE_FILE) library $(LINUX_FLAGS) ;; \
	  *) echo "This system is not supported" ;; \
       esac

# Build the library and include files
library: $(OBJ_FILES)
	$(ARCHIVER) $(ARCHIVER_FLAGS) $(LIBRARY_NAME) $(OBJ_FILES)

# Clean up the chaff
clean:
	$(REMOVE) $(OBJ_FILES) *.mod *.MOD *.stb *.a

# Define build dependencies
include make.dependencies

# Define suffix rules
include ../make.rules
