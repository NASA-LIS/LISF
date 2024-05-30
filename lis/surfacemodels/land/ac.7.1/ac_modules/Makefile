# Compiler options

DEBUG = 0
STATIC = 0

FC = gfortran
FCFLAGS = $(CPPFLAGS)
LDFLAGS_EXE =
LDFLAGS_LIB =

ifeq ($(STATIC), 1)
    ifneq ($(DEBUG), 0)
        $(error DEBUG=0 is required for STATIC=1 build)
    endif
endif

ifeq ($(findstring gfortran, $(FC)), gfortran)
    ifeq ($(STATIC), 1)
        # For statically compiled (release) executables
        FCFLAGS += -O2 -mtune=generic -funroll-loops
        LDFLAGS_EXE += -static
        LDFLAGS_LIB += -static
    else
        FCFLAGS += -fPIC -fall-intrinsics
        LDFLAGS_LIB += -shared

        ifeq ($(DEBUG), 1)
            # For debugging and development
            PCFLAGS += -gl
            FCFLAGS += -g -Wall -Wextra -Wimplicit-interface -fmax-errors=1 \
                       -fcheck=all -fbacktrace -std=f2008
        else
            # For production
            FCFLAGS += -O2 -march=native -funroll-loops
        endif
    endif
else ifeq ($(FC), ifort)
    ifneq ($(STATIC), 0)
        $(error STATIC=0 is required for FC=ifort)
    endif

    LDFLAGS_LIB += -shared
    FCFLAGS += -fPIC

    ifeq ($(DEBUG), 1)
        # For debugging and development
        FCFLAGS += -g -fpe0 -warn -traceback -debug extended -stand f08 \
                   -check all -check noarg_temp_created
    else
        # For production
        FCFLAGS += -O2 -fp-model strict -xHost
    endif
endif


# Main targets

EXECUTABLE = aquacrop
STATIC_LIBRARY = libaquacrop.a
SHARED_LIBRARY = libaquacrop.so
ifeq ($(STATIC), 1)
    LIBRARY = $(STATIC_LIBRARY)
else
    LIBRARY = $(SHARED_LIBRARY)
endif

bin: $(EXECUTABLE)
lib: $(LIBRARY)
all: bin lib
.DEFAULT_GOAL = all


# Fortran-related targets

OBJ = climprocessing.o \
    defaultcropsoil.o \
    global.o  \
    inforesults.o \
    initialsettings.o \
    kinds.o \
    project_input.o \
    rootunit.o \
    run.o \
    simul.o \
    startunit.o \
    tempprocessing.o \
    utils.o

$(STATIC_LIBRARY): $(OBJ) \
    $(OBJ_INTERFACE)
	$(AR) cr $(STATIC_LIBRARY) $(OBJ) $(OBJ_INTERFACE)

$(SHARED_LIBRARY): $(OBJ) \
    $(OBJ_INTERFACE)
	$(FC) $(LDFLAGS_LIB) $(FCFLAGS) $(wildcard *.f90) startunit.F90 -o $(LIBRARY)

fortranclean:
	$(RM) $(wildcard *.mod) $(wildcard *.o) $(STATIC_LIBRARY) $(SHARED_LIBRARY)


aquacrop.o: aquacrop.F90 \
    startunit.o
	$(FC) $(FCFLAGS) -c $< -o $@

climprocessing.o: climprocessing.f90 \
    global.o \
    kinds.o
	$(FC) $(FCFLAGS) -c $< -o $@

defaultcropsoil.o: defaultcropsoil.f90 \
    global.o \
    kinds.o
	$(FC) $(FCFLAGS) -c $< -o $@

global.o: global.f90 \
    kinds.o \
    project_input.o \
    utils.o
	$(FC) $(FCFLAGS) -c $< -o $@

inforesults.o: inforesults.f90 \
    global.o \
    kinds.o \
    utils.o
	$(FC) $(FCFLAGS) -c $< -o $@

initialsettings.o: initialsettings.f90 \
    global.o \
    defaultcropsoil.o \
    kinds.o
	$(FC) $(FCFLAGS) -c $< -o $@

kinds.o: kinds.f90
	$(FC) $(FCFLAGS) -c $< -o $@

project_input.o: project_input.f90 \
    kinds.o \
    utils.o
	$(FC) $(FCFLAGS) -c $< -o $@

rootunit.o: rootunit.f90 \
    global.o \
    kinds.o \
    utils.o
	$(FC) $(FCFLAGS) -c $< -o $@

run.o: run.f90 \
    global.o \
    inforesults.o \
    kinds.o \
    rootunit.o \
    simul.o \
    tempprocessing.o \
    climprocessing.o \
    utils.o
	$(FC) $(FCFLAGS) -c $< -o $@

simul.o: simul.f90 \
    global.o \
    tempprocessing.o \
    kinds.o
	$(FC) $(FCFLAGS) -c $< -o $@

startunit.o: startunit.F90 \
    global.o \
    initialsettings.o \
    kinds.o \
    run.o \
    utils.o
	$(FC) $(FCFLAGS) -c $< -o $@

tempprocessing.o: tempprocessing.f90 \
    global.o \
    kinds.o \
    project_input.o \
    utils.o
	$(FC) $(FCFLAGS) -c $< -o $@

utils.o: utils.f90 \
    kinds.o
	$(FC) $(FCFLAGS) -c $< -o $@


# Interface Fortran code: always adding the parent module and the kinds and
# utils modules as dependencies for convenience, even if not actually used.

# Finally, the actual 'executable' and 'clean' targets

# Fortran-based executable
$(EXECUTABLE): aquacrop.o
	$(FC) $(LDFLAGS_EXE) aquacrop.o $(OBJ) -o $(EXECUTABLE)

clean: fortranclean
	$(RM) $(EXECUTABLE)
