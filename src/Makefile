# The following variables are imported from Makefile in ../
# and don't need to be set here explicitly:
#FCOM
#CPPFLAGS
#COMPFLAGS
#DEBUGFLAGS

# List of source code files for standard (full) setup (equal to cnmodel):
MODS=params_core.mod.f90 classdefs.mod.f90 sofunutils.mod.f90 params_siml_cnmodel.mod.f90 params_soil.mod.f90 params_site.mod.f90 gridvars_siterun.mod.f90 forcing_siterun.mod.f90 interface_biosphere.mod.f90 rates.mod.f90 waterbal_stash.mod.f90 soiltemp_sitch.mod.f90 plant.mod.f90 gpp_pmodel.mod.f90 phenology_lpx.mod.f90 lambertw.mod.f90 findroot_fzeroin.mod.f90 npp.mod.f90 turnover.mod.f90 vegdynamics_dummy.mod.f90 ntransform_xuri.mod.f90 littersom_lpj.mod.f90 nuptake.mod.f90 allocation_cnmodel.mod.f90 landuse.mod.f90
SOURCES=sofun.f90 biosphere_cnmodel.f90

# List of source code files for reduced setup, executing only SPLASH: XXX UNDER CONSTRUCTION XXX
SPLASH_MODS=params_core.mod.f90 classdefs.mod.f90 sofunutils.mod.f90 params_siml.mod.f90 params_soil.mod.f90 params_site.mod.f90 gridvars_siterun.mod.f90 forcing_siterun.mod.f90 params_modl.mod.f90 interface_biosphere.mod.f90 vars_core.mod.f90 outvars.mod.f90 waterbal_stash.mod.f90 soiltemp_sitch.mod.f90 
SPLASH_SOURCES=sofun.f90 biosphere_splash.f90

# List of source code files for reduced setup, executing SPLASH and P-MODEL:
PMODEL_MODS=params_core.mod.f90 classdefs.mod.f90 sofunutils.mod.f90 params_siml_pmodel.mod.f90 params_soil.mod.f90 params_site.mod.f90 gridvars_siterun.mod.f90 forcing_siterun_pmodel.mod.f90 interface_biosphere.mod.f90 waterbal_stash.mod.f90 soiltemp_sitch.mod.f90 plant_pmodel.mod.f90 gpp_pmodel.mod.f90 
PMODEL_SOURCES=sofun.f90 biosphere_pmodel.f90

# List of source code files for reduced setup, executing the C-only model: XXX UNDER CONSTRUCTION XXX
CMODEL_MODS=params_core.mod.f90 classdefs.mod.f90 sofunutils.mod.f90 params_siml.mod.f90 params_soil.mod.f90 params_site.mod.f90 gridvars_siterun.mod.f90 forcing_siterun.mod.f90 interface_biosphere.mod.f90 rates.mod.f90 waterbal_stash.mod.f90 soiltemp_sitch.mod.f90 plant.mod.f90 gpp_pmodel.mod.f90 phenology_lpx.mod.f90 lambertw.mod.f90 findroot_fzeroin.mod.f90 npp.mod.f90 turnover.mod.f90 allocation_cmodel.mod.f90 vegdynamics_dummy.mod.f90 landuse.mod.f90
CMODEL_SOURCES=sofun.f90 biosphere_cmodel.f90

# List of source code files for reduced setup, executing the full CN-only model:
CNMODEL_MODS=params_core.mod.f90 classdefs.mod.f90 sofunutils.mod.f90 params_siml_cnmodel.mod.f90 params_soil.mod.f90 params_site.mod.f90 gridvars_siterun.mod.f90 forcing_siterun.mod.f90 interface_biosphere.mod.f90 rates.mod.f90 waterbal_stash.mod.f90 soiltemp_sitch.mod.f90 lambertw.mod.f90 plant.mod.f90 gpp_pmodel.mod.f90 phenology_lpx.mod.f90 findroot_fzeroin.mod.f90 turnover.mod.f90 npp.mod.f90 vegdynamics_dummy.mod.f90 ntransform_xuri.mod.f90 littersom_lpj.mod.f90 nuptake.mod.f90 allocation_cnmodel.mod.f90 landuse.mod.f90
CNMODEL_SOURCES=sofun.f90 biosphere_cnmodel.f90

# List of build dependencies (Makefile and include files):
DEPEND=Makefile ../Makefile sofun_module_control.inc

# Generate list of .o, .do files from SOURCES
OBJS=$(SOURCES:.f90=.o)
MODOBJS=$(MODS:.f90=.o)
DEBUGOBJS=$(SOURCES:.f90=.do)
DEBUGMODOBJS=$(MODS:.f90=.do)

SPLASH_OBJS=$(SPLASH_SOURCES:.f90=.o)
SPLASH_MODOBJS=$(SPLASH_MODS:.f90=.o)
DBGSPLASH_OBJS=$(SPLASH_SOURCES:.f90=.do)
DBGSPLASH_MODOBJS=$(SPLASH_MODS:.f90=.do)

PMODEL_OBJS=$(PMODEL_SOURCES:.f90=.o)
PMODEL_MODOBJS=$(PMODEL_MODS:.f90=.o)
DBGPMODEL_OBJS=$(PMODEL_SOURCES:.f90=.do)
DBGPMODEL_MODOBJS=$(PMODEL_MODS:.f90=.do)

CMODEL_OBJS=$(CMODEL_SOURCES:.f90=.o)
CMODEL_MODOBJS=$(CMODEL_MODS:.f90=.o)
DBGCMODEL_OBJS=$(CMODEL_SOURCES:.f90=.do)
DBGCMODEL_MODOBJS=$(CMODEL_MODS:.f90=.do)

CNMODEL_OBJS=$(CNMODEL_SOURCES:.f90=.o)
CNMODEL_MODOBJS=$(CNMODEL_MODS:.f90=.o)
DBGCNMODEL_OBJS=$(CNMODEL_SOURCES:.f90=.do)
DBGCNMODEL_MODOBJS=$(CNMODEL_MODS:.f90=.do)

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

# reduced model setup: only SPLASH
dbgsplash: $(DBGSPLASH_MODOBJS) $(DBGSPLASH_OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(DBGSPLASH_MODOBJS) $(DBGSPLASH_OBJS) 

# reduced model setup: only SPLASH and PMODEL
pmodel: $(PMODEL_MODOBJS) $(PMODEL_OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(PMODEL_MODOBJS) $(PMODEL_OBJS) 

# reduced model setup: only SPLASH and PMODEL
dbgpmodel: $(DBGPMODEL_MODOBJS) $DBG(PMODEL_OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(DBGPMODEL_MODOBJS) $DBG(PMODEL_OBJS) 

# reduced model setup: C-model only
cmodel: $(CMODEL_MODOBJS) $(CMODEL_OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(CMODEL_MODOBJS) $(CMODEL_OBJS) 

# reduced model setup: C-model only
dbgcmodel: $(DBGCMODEL_MODOBJS) $(DBGCMODEL_OBJS) 
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(DBGCMODEL_MODOBJS) $(DBGCMODEL_OBJS) 

# full model setup: CN-model
cnmodel: $(CNMODEL_MODOBJS) $(CNMODEL_OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(CNMODEL_MODOBJS) $(CNMODEL_OBJS) 

# full model setup: CN-model
dbgcnmodel: $(DBGCNMODEL_MODOBJS) $(DBGCNMODEL_OBJS)
	# $(MAKE) -C cdfcode
	-rm $(ARCHIVE)
	# cp $(CDFARCHIVE) $(ARCHIVE)
	ar r $(ARCHIVE) $(DBGCNMODEL_MODOBJS) $(DBGCNMODEL_OBJS) 

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
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(SPLASH_MODOBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(SPLASH_OBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(DBGSPLASH_MODOBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(DBGSPLASH_OBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(PMODEL_MODOBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(PMODEL_OBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(DBGPMODEL_MODOBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(DBGPMODEL_OBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(CMODEL_MODOBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(CMODEL_OBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(DBGCMODEL_MODOBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(DBGCMODEL_OBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(CNMODEL_MODOBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(CNMODEL_OBJS): %.o: %.f
	$(FCOM) -c -o $@ $(COMPFLAGS) $*.f

$(DBGCNMODEL_MODOBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

$(DBGCNMODEL_OBJS): %.do: %.f
	$(FCOM) -c -o $@ $(DEBUGFLAGS) $*.f

# don't delete the intermediate *.f files
.SECONDARY: $(SOURCES:.f90=.f) $(SPLASH_SOURCES:.f90=.f) $(PMODEL_SOURCES:.f90=.f) $(CMODEL_SOURCES:.f90=.f) $(CNMODEL_SOURCES:.f90=.f)

# clean: remove .f, .o, .do, and .stb files
.PHONY: clean
clean:
	-rm -f *.f *.o *.do *.stb *.mod
#	rm $(ARCHIVE)

#EOF
