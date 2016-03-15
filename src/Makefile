# The following variables are imported from Makefile in ../
# and don't need to be set here explicitly:
#FCOM
#CPPFLAGS
#COMPFLAGS
#DEBUGFLAGS

# List of source code files for standard (full) setup:
MODS=params_core.mod.f90 classdefs.mod.f90 sofunutils.mod.f90 params_siml.mod.f90 params_soil.mod.f90 params_site.mod.f90 gridvars_siterun.mod.f90 forcing_siterun.mod.f90 params_modl.mod.f90 rates.mod.f90 vars_core.mod.f90 outvars.mod.f90 phenology_lpx.mod.f90 waterbal_stash.mod.f90 vegdynamics_fix.mod.f90 gpp_pmodel.mod.f90 nuptake_constexu.mod.f90 npp_lpj.mod.f90 findroot_fzeroin.mod.f90 allocation.mod.f90 turnover.mod.f90 littersom_lpj.mod.f90 ntransform_xuri.mod.f90 soiltemp_sitch.mod.f90 
SOURCES=sofun.f90 biosphere.f90 kill.f90 writeout_ascii.f90 

# List of source code files for reduced setup, executing only SPLASH:
SPLASH_MODS=params_core.mod.f90 classdefs.mod.f90 sofunutils.mod.f90 params_siml.mod.f90 params_soil.mod.f90 params_site.mod.f90 gridvars_siterun.mod.f90 forcing_siterun.mod.f90 params_modl.mod.f90 vars_core.mod.f90 outvars.mod.f90 waterbal_stash.mod.f90 soiltemp_sitch.mod.f90 
SPLASH_SOURCES=sofun.f90 biosphere_splash.f90

# List of source code files for reduced setup, executing SPLASH and P-MODEL:
PMODEL_MODS=params_core.mod.f90 classdefs.mod.f90 sofunutils.mod.f90 params_siml.mod.f90 params_soil.mod.f90 params_site.mod.f90 gridvars_siterun.mod.f90 forcing_siterun_pmodel.mod.f90 params_modl.mod.f90 vars_core.mod.f90 outvars.mod.f90 waterbal_stash.mod.f90 soiltemp_sitch.mod.f90 gpp_pmodel.mod.f90
PMODEL_SOURCES=sofun.f90 biosphere_pmodel.f90

# List of build dependencies (Makefile and include files):
DEPEND=Makefile ../Makefile sofun_module_control.inc

# Generate list of .o, .do files from SOURCES
OBJS=$(SOURCES:.f90=.o)
MODOBJS=$(MODS:.f90=.o)
DEBUGOBJS=$(SOURCES:.f90=.do)
DEBUGMODOBJS=$(MODS:.f90=.do)
SPLASH_OBJS=$(SPLASH_SOURCES:.f90=.o)
SPLASH_MODOBJS=$(SPLASH_MODS:.f90=.o)
PMODEL_OBJS=$(PMODEL_SOURCES:.f90=.o)
PMODEL_MODOBJS=$(PMODEL_MODS:.f90=.o)

# Archive file names
ARCHIVE=sofun.a
# CDFARCHIVE=./cdfcode/cdf.a

# Export variables that are needed by Makefiles in the subdirectories (called below)
export FCOM CPPFLAGS COMPFLAGS DEBUGFLAGS #LIBS

# Targets
# -------
# standard source code: (depends on object files)
all:	$(MODOBJS) $(OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(MODOBJS) $(OBJS) 

# code for debugging:
debug:	$(DEBUGMODOBJS) $(DEBUGOBJS)
	# $(MAKE) debug -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(DEBUGMODOBJS) $(DEBUGOBJS)

# reduced model setup: only SPLASH
splash: $(SPLASH_MODOBJS) $(SPLASH_OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(SPLASH_MODOBJS) $(SPLASH_OBJS) 

# reduced model setup: only SPLASH and PMODEL
pmodel: $(PMODEL_MODOBJS) $(PMODEL_OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(PMODEL_MODOBJS) $(PMODEL_OBJS) 

# default rules (.f90 -> .f -> .o/.do)
# ------------
%.f: %.f90 $(DEPEND)
	rm -f $*.f
	$(FCOM) $(CPPFLAGS) $*.f90 > $*.f 

$(MODOBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(OBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(DEBUGOBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(DEBUGMODOBJS): %.do: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(SPLASH_MODOBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(SPLASH_OBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(PMODEL_MODOBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(PMODEL_OBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

# don't delete the intermediate *.f files
.SECONDARY: $(SOURCES:.f90=.f) $(SPLASH_SOURCES:.f90=.f) $(PMODEL_SOURCES:.f90=.f)

# clean: remove .f, .o, .do, and .stb files
.PHONY: clean
clean:
	-rm -f *.f *.o *.do *.stb *.mod
#	rm $(ARCHIVE)

#EOF
