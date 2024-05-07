
# Compiler/Linker
ifeq ($(LIS_ARCH),linux_ifc)
FC = ifort
LD = ifort 
else ifeq ($(LIS_ARCH),linux_gfortran)
FC = gfortran
LD = gfortran
else
$(error "Please specify the compiler via the LIS_ARCH env variable.  Only linux_ifc and linux_gfortran are supported. See the LISF Installation Guide for more information.")
endif

# compiler/linker flags
# Compile with DEBUG messages
#FCFLAGS = -DDEBUG -fpp -g -u -traceback -fpe0 -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl -I/home/ytian/proj-disk/libs/CRTM_Profile_Utility_intel_11_1_038/include
# Compile without DEBUG messages
ifeq ($(LIS_ARCH),linux_ifc)
FCFLAGS = -fpp -g -u -traceback -fpe0 -nomixed_str_len_arg -names lowercase -convert big_endian -assume byterecl 
else ifeq ($(LIS_ARCH),linux_gfortran)
FCFLAGS = -cpp -g -u -fbacktrace -ffpe-trap=invalid,zero,overflow -fconvert=big-endian
endif
#-I/home/ytian/proj-disk/libs/CRTM_Profile_Utility_intel_11_1_038/include

#LIBS = /home/ytian/proj-disk/libs/CRTM_Profile_Utility_intel_11_1_038/lib/libProfile_Utility.a

# files removal
RM = /bin/rm -f

#SRCS = atm_sub.F90 dielwat_sub.F90 lis_mem.F90 soil.F90 water.F90 dielice_sub.F90 fresnel_sub.F90 rghref_sub.F90 veg.F90 dielsoil_sub.F90 liebe89_sub.F90 snow.F90 veg_sub.F90

SRCS := $(wildcard *.F90) 

# program's object files
OBJS = $(SRCS:.F90=.o) 

# program's executable
PROG = liblis_mem.a 
#PROG = lis_mem_debug.a 

$(PROG): $(OBJS)
	ar crvs $(PROG) $(OBJS)

# compile F90 source files into object files.
%.o: %.F90 
	$(FC) $(FCFLAGS) -c $<

# clean everything
clean:
	$(RM) $(OBJS) $(PROG) *__genmod.???


