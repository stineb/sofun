# The following variables are imported from Makefile in ../
# and don't need to be set here explicitly:
#FCOM
#CPPFLAGS
#COMPFLAGS
#DEBUGFLAGS

# List of source code files:  xxx where to put sofun (main) ??? xxx
MODS=params_core.mod.F classdefs.mod.F params_siml.mod.F params_site.mod.F gridvars_siterun.mod.F forcing_siterun.mod.F params_modl.mod.F rates.mod.F pools.mod.F input_site.mod.F fluxes.mod.F waterbal_stash.mod.F treegeometry.mod.F nuptake_fun.mod.F fixveg.mod.F outvars.mod.F 
SOURCES=sofun.F biosphere.F getpar_modl.F init.F initio.F kill.F gpp.F npp.F littersom.F ntransform.F turnover.F allocation.F getout.F writeout_ascii.F

# List of build dependencies (Makefile and include files):
DEPEND=Makefile ../Makefile sofun_module_control.inc

# Generate list of .o, .do files from SOURCES
OBJS=$(SOURCES:.F=.o)
MODOBJS=$(MODS:.F=.o)
DEBUGOBJS=$(SOURCES:.F=.do)
DEBUGMODOBJS=$(MODS:.F=.do)

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

# don't delete the intermediate *.for files
.SECONDARY:     $(SOURCES:.F=.for)

# clean: remove .for, .o, .do, and .stb files
.PHONY: clean
clean:
	-rm -f *.for *.o *.do *.stb *.mod
#	rm $(ARCHIVE)

#EOF
