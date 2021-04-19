# system choice
 CONFIG = "intel"

ifeq ($(CONFIG),"intel")
FC      = ifort       # compiler
LD      = ifort       # linker
LDFLAGS =             # linker options
SFLAGS  = -c -r8 -i-dynamic -mcmodel=medium    # standard compiler flags
OFLAGS  = -O2  # optimization compiler flags
DFLAGS  =  # debug compiler flags
endif

# executable name
EXEC = ../NODALEP

# put compiler flags together
#CFLAGS = $(SFLAGS) $(OFLAGS) $(OPT) $(MPIINCL) $(HDF5INCL) # optimize
 CFLAGS = $(SFLAGS) $(OFLAGS) $(OPT) # optimize (basic)
# CFLAGS = $(SFLAGS) $(DFLAGS) $(OPT) $(MPIINCL) $(HDF5INCL) # debug


# put libraries together
LDLIBS = $(MPILIB) $(MPI) $(HDF5LIB) $(HDF5)

# source directory
source   = ../source/

OBJS =\
  NODALEP.o\
  units_module.o\
  fermi_integral_module.o\
  readin_nuprox_module.o\
  readin_stp_module.o\
  readin_lambda_module.o\
  readin_LNk_module.o\
  readin_distf_module.o\
  eos_module.o\
  egroup_module.o\
  heating_cooling_module.o\
  luminosity_estimates.o\
  nluminosity_estimates.o\
  lambda_theory.o\
  emission_testing.o
$(EXEC):\
  $(OBJS); $(LD) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(EXEC)

clean: 
	rm -v *.o *.mod ../NODALEP  

NODALEP.o:\
  units_module.o\
  fermi_integral_module.o\
  readin_nuprox_module.o\
  readin_stp_module.o\
  readin_lambda_module.o\
  readin_LNk_module.o\
  readin_distf_module.o\
  eos_module.o\
  egroup_module.o\
  heating_cooling_module.o\
  luminosity_estimates.o\
  nluminosity_estimates.o\
  lambda_theory.o\
  emission_testing.o\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)

units_module.o:\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)

fermi_integral_module.o:\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)

readin_nuprox_module.o:\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)
  
readin_stp_module.o:\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)
 
readin_distf_module.o:\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)

readin_lambda_module.o:\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)

readin_LNk_module.o:\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)
 
eos_module.o:\
  egroup_module.o\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)  

egroup_module.o:\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)    

heating_cooling_module.o:\
  lambda_theory.o\
  eos_module.o\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)
  
luminosity_estimates.o:\
  lambda_theory.o\
  fermi_integral_module.o\
  units_module.o\
  emission_testing.o\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)

nluminosity_estimates.o:\
  lambda_theory.o\
  fermi_integral_module.o\
  units_module.o\
  emission_testing.o\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)
  
lambda_theory.o:\
  units_module.o\
  fermi_integral_module.o\
  emission_testing.o\
  egroup_module.o\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)
  
emission_testing.o:\
  units_module.o\
  egroup_module.o\
  $(source)$(@:.o=.f90); $(FC) $(CFLAGS) $(source)$(@:.o=.f90)    
