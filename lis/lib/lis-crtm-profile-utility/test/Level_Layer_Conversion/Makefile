#==============================================================================
#
# Makefile for Level to Layer conversion test code
#
# $Id: Makefile 1958 2008-03-18 16:35:36Z paul.vandelst@noaa.gov $
#==============================================================================

# Define the executable file
EXE_FILE = Level_Layer_Conversion_Test

# include the common test make components
include ../make.common

# The run target
OUTFILE = Profile_Conversion.asc
run: all
	$(EXE_FILE)
	@echo; \
        echo "     Differencing $(OUTFILE) with baseline:" ; \
	diff $(OUTFILE).Baseline $(OUTFILE)

# include the dependency list
include make.dependencies
