############################
# SOFUN master Makefile
# adopted from LPX Makefile
############################

##################################
## Select configuration profile ##
##################################
# pgf     - PGF95 compiler
# gfor    - gfortran compiler

PROFILE=pgf

##################
## pgf profile ##
##################
ifeq ($(PROFILE),pgf)
# Compiler and options
FCOM=pgf95
CPPFLAGS=-E
COMPFLAGS=-Mextend -Mfreeform
#COMPFLAGS       = -Mextend -Mdalign -Kieee -Ktrap=fp -O2 -Mprof=lines # to analyze computation time by subroutines
DEBUGFLAGS=-g -O0 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform

# System libraries
#LIBS = -L $(NETCDF_LIB) -lnetcdf -lnetcdff
endif

#####################
## gfor profile ##
#####################
ifeq ($(PROFILE),gfor)
# Compiler and options
FCOM=gfortran
CPPFLAGS=-e
COMPFLAGS=-ffixed-line-length-0 -fdefault-real-8 -ffree-form

DEBUGFLAGS=-ffixed-line-length-0 -fdefault-real-8 -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ffree-form


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
$(error 'ERROR: Please select a configuration profile in the Makefile (e.g. PROFILE=gfor).')
endif

# Add flags for MPI parallelization (enable the following lines when the parallel_mpi feature is turned on)
#LIBS += $(shell mpif90 --showme:link)
#COMPFLAGS += $(shell mpif90 --showme:compile)
#DEBUGFLAGS += $(shell mpif90 --showme:compile)

# Add library include files to compiler flags
#COMPFLAGS += -I$(NETCDF_INC)
#DEBUGFLAGS += -I$(NETCDF_INC)

# name of executable
EXE = runsofun

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

# clean: remove exe and .o and .do files
.PHONY: clean
clean:
	-rm $(EXE)
	$(MAKE) clean -C src
# include libraries when necessary
#	$(MAKE) clean -C lpj/cdfcode

#EOF
