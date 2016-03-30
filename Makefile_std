############################
# SOFUN master Makefile
# adopted from LPX Makefile
############################

##################################
## Select configuration profile ##
##################################
# pgf     - PGF95 compiler
# gfor    - gfortran compiler
# intel   - ifort compiler

PROFILE=pgi
# PROFILE=intel

##################
## pgf profile ##
##################
ifeq ($(PROFILE),pgi)
# Compiler and options
FCOM=pgf95
CPPFLAGS=-E
COMPFLAGS=-Mextend -Mfreeform -Mdalign -Kieee -Ktrap=fp -O2
DPCOMPFLAGS=-r8 -Mextend -Mfreeform  -Mdalign -Kieee -Ktrap=fp -O2
#COMPFLAGS= -Mextend -Mdalign -Kieee -Ktrap=fp -O2 -Mprof=lines # to analyze computation time by subroutines
DEBUGFLAGS=-g -O0 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform
DPDEBUGFLAGS=-g -O0 -r8 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform

# System libraries
#LIBS = -L $(NETCDF_LIB) -lnetcdf -lnetcdff
endif

#####################
## intel profile ##
#####################
ifeq ($(PROFILE),intel)
# Compiler and options
FCOM=ifort
CPPFLAGS=-e -fpp -preprocess_only -E
COMPFLAGS=-O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -extend_source -free -g -traceback ##-r8 -i4 -align -pc64 -fp-model strict 
DEBUGFLAGS=-O3 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -extend_source -free -warn all -implicitnone -g -traceback -fpe0 -fpstkchk -CU

# System libraries
# NETCDF_INC = /usr/local/include
# NETCDF_LIB = /usr/local/lib
# LIBS = -L $(NETCDF_LIB) -lnetcdf
endif

####################
## general config ##
####################

# Check if FCOM is set
ifeq ($(strip $(FCOM)),)
$(error 'ERROR. Select a valid configuration profile in the Makefile (e.g. PROFILE=gfor).')
endif

# Add flags for MPI parallelization (enable the following lines when the parallel_mpi feature is turned on)
#LIBS += $(shell mpif90 --showme:link)
#COMPFLAGS += $(shell mpif90 --showme:compile)
#DEBUGFLAGS += $(shell mpif90 --showme:compile)

# Add library include files to compiler flags
#COMPFLAGS += -I$(NETCDF_INC)
#DEBUGFLAGS += -I$(NETCDF_INC)

# name of executable
EXE        = runsofun
SPLASH_EXE = runsplash
PMODEL_EXE = runpmodel
CMODEL_EXE = runcmodel
CNMODEL_EXE = runcnmodel

ARCHIVES= ./src/sofun.a
# ARLPJ= ./lpj/lpj.a (archive names when compiling with different option)

# Export variables that are needed by Makefiles in the subdirectories (called below)
export FCOM CPPFLAGS COMPFLAGS DEBUGFLAGS #LIBS

# Targets
# -------
standard: 
	 $(MAKE) -C src
	 $(FCOM) -o $(EXE) $(COMPFLAGS) $(ARCHIVES)
#  include libraries when necessary
#	 $(FCOM) -o $(EXE) $(COMPFLAGS) $(ARCHIVES) $(LIBS)

# code for debugging:
debug: 
	$(MAKE) debug -C src
	$(FCOM) -o $(EXE) $(DEBUGFLAGS) $(ARCHIVES) #$(LIBS)

# including double precision flag
dp: 
	 $(MAKE) dp -C src
	 $(FCOM) -o $(EXE) $(DPCOMPFLAGS) $(ARCHIVES)

# code for debugging, including double precision flags:
dpdebug: 
	$(MAKE) dpdebug -C src
	$(FCOM) -o $(EXE) $(DPDEBUGFLAGS) $(ARCHIVES) #$(LIBS)

# reduced model setup: only SPLASH
splash: 
	 $(MAKE) splash -C src
	 $(FCOM) -o $(SPLASH_EXE) $(COMPFLAGS) $(ARCHIVES)

# reduced model setup: only SPLASH and PMODEL
pmodel: 
	 $(MAKE) pmodel -C src
	 $(FCOM) -o $(PMODEL_EXE) $(COMPFLAGS) $(ARCHIVES)

# reduced model setup: fixed allocation, no litter, soil and inorganic C and N dynamics
cmodel: 
	 $(MAKE) cmodel -C src
	 $(FCOM) -o $(CMODEL_EXE) $(COMPFLAGS) $(ARCHIVES)

# full model setup
cnmodel: 
	 $(MAKE) cnmodel -C src
	 $(FCOM) -o $(CNMODEL_EXE) $(COMPFLAGS) $(ARCHIVES)

# clean: remove exe and .o and .do files
.PHONY: clean
clean:
	-rm $(EXE) $(SPLASH_EXE) $(PMODEL_EXE) $(CMODEL_EXE) $(CNMODEL_EXE)
	$(MAKE) clean -C src
# include libraries when necessary
#	$(MAKE) clean -C lpj/cdfcode

#EOF
