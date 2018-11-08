#==============================================================================
#
# Makefile for geopotential height test code
#
# $Id: Makefile 1960 2008-03-18 17:26:10Z paul.vandelst@noaa.gov $
#==============================================================================

# Define the executable file
EXE_FILE = Geopotential_Test

# include the common test make components
include ../make.common

# The run target
OUTFILE = Geopotential_Test.asc
run: all
	$(EXE_FILE)
	@echo; \
        echo "     Differencing $(OUTFILE) with baseline:" ; \
	diff $(OUTFILE).Baseline $(OUTFILE)

# include the dependency list
include make.dependencies
