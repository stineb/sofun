c //////////////////////////////////////////////////////////////////////
c SELECT MODULES FOR DIFFERENT PARAMETRISATIONS
c ----------------------------------------------------------------------

c ----------------------------------------------------------------------
c MAINTENANCE RESPIRATION
c Default: like in LPX
c Option:   After Mäkelä et al., 2008, New Phyt.
#define _resp_maint_makela 1
c ----------------------------------------------------------------------

c ----------------------------------------------------------------------
c ALLOCATION
c Default: annually called, like in LPX for grass PFT
c Option:  Pipe model as described in Mäkelä et al., 2008, New Phyt.
#define _allocation_pipe 0
c ----------------------------------------------------------------------

c //////////////////////////////////////////////////////////////////////
c SELECT MODEL INTEGRATION
c defining which model parts are to be calculated online and which to 
c be prescribed. Default is always to calculate it online (not prescr.)
c In general, if selected, daily values for the course of the whole 
c transient simulation are read in when initialising the model. This 
c requires input files formatted like ''
c This is only useful for site scale simulations. 
c ----------------------------------------------------------------------

c Prescribe GPP per day
#define _prescr_gpp_site 1

c //////////////////////////////////////////////////////////////////////
c SANITY
c Compile additional code for online sanity checks. This is recommended
c for simulations during development stage. De-activate later to reduce
c computational costs.
c ----------------------------------------------------------------------
#define _check_sanity 1
