
# standard
FCOM=pgf95 -Mextend

# for debugging
# FCOM=pgf95 -Mextend -g -O0 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp 


# standard
# FCOM=gfortran -ffixed-line-length-0 -fdefault-real-8

# for debugging
# FCOM=gfortran -ffixed-line-length-0 -fdefault-real-8 -g -fbounds-check -Wall -fbacktrace -finit-real=nan

export FCOM

# name of executable
EXE = runsofun

standard: 
	 $(FCOM) -c params_core.mod.F
	 $(FCOM) -c classdefs.mod.F
	 $(FCOM) -c params_siml.mod.F
	 $(FCOM) -c params_site.mod.F
	 $(FCOM) -c params_modl.mod.F
	 $(FCOM) -c gridvars.mod.F
	 $(FCOM) -c vars.mod.F
	 $(FCOM) -c forcing.mod.F
	 $(FCOM) -c rates.mod.F
	 $(FCOM) -c pools.mod.F
	 $(FCOM) -c init.F
	 $(FCOM) -c update.F
	 $(FCOM) -c kill.F
	 $(FCOM) -c gpp.F
	 $(FCOM) -c npp.F
	 $(FCOM) -c littersom.F
	 $(FCOM) -c ntransform.F
	 $(FCOM) -c turnover.F
	 $(FCOM) -c allocation.F
	 $(FCOM) -c establishment.F
	 $(FCOM) -c getout.F
	 $(FCOM) -c nuptake.F
	 $(FCOM) params_core.mod.o classdefs.mod.o params_siml.mod.o params_site.mod.o params_modl.mod.o gridvars.mod.o vars.mod.o forcing.mod.o rates.mod.o pools.mod.o init.o update.o kill.o gpp.o npp.o littersom.o ntransform.o turnover.o allocation.o establishment.o getout.o nuptake.o sofun.F -o $(EXE)

.PHONY: clean
clean:
	rm $(EXE) *.o *.mod

