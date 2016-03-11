# The following variables are imported from Makefile in ../
# and don't need to be set here explicitly:
#FCOM
#CPPFLAGS
#COMPFLAGS
#DEBUGFLAGS

# List of source code files for standard (full) setup:
MODS=params_core.mod.F classdefs.mod.F sofunutils.mod.F params_siml.mod.F params_soil.mod.F params_site.mod.F gridvars_siterun.mod.F forcing_siterun.mod.F params_modl.mod.F rates.mod.F vars_core.mod.F outvars.mod.F phenology_lpx.mod.F waterbal_stash.mod.F vegdynamics_fix.mod.F gpp_pmodel.mod.F nuptake_constexu.mod.F npp_lpj.mod.F findroot_fzeroin.mod.F allocation.mod.F turnover.mod.F littersom_lpj.mod.F ntransform_xuri.mod.F soiltemp_sitch.mod.F 
SOURCES=sofun.F biosphere.F kill.F writeout_ascii.F 

# List of source code files for reduced setup, executing only SPLASH:
SPLASH_MODS=params_core.mod.F classdefs.mod.F sofunutils.mod.F params_siml.mod.F params_soil.mod.F params_site.mod.F gridvars_siterun.mod.F forcing_siterun.mod.F params_modl.mod.F vars_core.mod.F outvars.mod.F waterbal_stash.mod.F soiltemp_sitch.mod.F 
SPLASH_SOURCES=sofun.F biosphere_splash.F

# List of source code files for reduced setup, executing SPLASH and P-MODEL:
PMODEL_MODS=params_core.mod.F classdefs.mod.F sofunutils.mod.F params_siml.mod.F params_soil.mod.F params_site.mod.F gridvars_siterun.mod.F forcing_siterun_pmodel.mod.F params_modl.mod.F vars_core.mod.F outvars.mod.F waterbal_stash.mod.F soiltemp_sitch.mod.F gpp_pmodel.mod.F
PMODEL_SOURCES=sofun.F biosphere_pmodel.F

# List of build dependencies (Makefile and include files):
DEPEND=Makefile ../Makefile sofun_module_control.inc

# Generate list of .o, .do files from SOURCES
OBJS=$(SOURCES:.F=.o)
MODOBJS=$(MODS:.F=.o)
DEBUGOBJS=$(SOURCES:.F=.do)
DEBUGMODOBJS=$(MODS:.F=.do)
SPLASH_OBJS=$(SPLASH_SOURCES:.F=.o)
SPLASH_MODOBJS=$(SPLASH_MODS:.F=.o)
PMODEL_OBJS=$(PMODEL_SOURCES:.F=.o)
PMODEL_MODOBJS=$(PMODEL_MODS:.F=.o)

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

# default rules (.F -> .for -> .o/.do)
# ------------
%.for: %.F $(DEPEND)
	rm -f $*.for
	$(FCOM) $(CPPFLAGS) $*.F > $*.for 

$(MODOBJS): %.o: %.for
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.for

$(OBJS): %.o: %.for
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.for

$(DEBUGOBJS): %.do: %.for
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.for

$(DEBUGMODOBJS): %.do: %.for
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.for

$(SPLASH_MODOBJS): %.o: %.for
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.for

$(SPLASH_OBJS): %.o: %.for
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.for

$(PMODEL_MODOBJS): %.o: %.for
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.for

$(PMODEL_OBJS): %.o: %.for
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.for

# don't delete the intermediate *.for files
.SECONDARY: $(SOURCES:.F=.for) $(SPLASH_SOURCES:.F=.for) $(PMODEL_SOURCES:.F=.for)

# clean: remove .for, .o, .do, and .stb files
.PHONY: clean
clean:
	-rm -f *.for *.o *.do *.stb *.mod
#	rm $(ARCHIVE)

#EOF
