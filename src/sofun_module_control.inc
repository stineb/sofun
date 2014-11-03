c //////////////////////////////////////////////////////////////////////
c SELECT MODULES FOR DIFFERENT PARAMETRISATIONS
c ----------------------------------------------------------------------

c ----------------------------------------------------------------------
c MAINTENANCE RESPIRATION
c Default: like in LPX
c Option:   After Mäkelä et al., 2008, New Phyt.
#define _resp_maint_makela 1
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

