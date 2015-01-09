! //////////////////////////////////////////////////////////////////////
! SELECT MODULES FOR DIFFERENT PARAMETRISATIONS
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! FIXED/DYNAMIC VEGETATION
! Default: fixed, reading in FPC for each PFT from site parameter file)
! Option:  dynamic, like in LPX
#define _fix_veg 1
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! MAINTENANCE RESPIRATION
! Default: like in LPX
! Option:   After Mäkelä et al., 2008, New Phyt.
#define _resp_maint_makela 1
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! PHOTOSYNTHESIS
! Default: xxx
! Option:  After Mäkelä et al., 2008, New Phyt.
#define _gpp_makela 1
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! N UPTAKE
! Default: FUN type after Fisher et al., 2011
! Option:  After Mäkelä et al., 2008, New Phyt.
#define _nuptake_makela 1
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ALLOCATION
! Default: annually called, like in LPX for grass PFT
! Option:  Pipe model as described in Mäkelä et al., 2008, New Phyt.
! Note that _allocation_pipe does not account for the separation between
! sapwood and heartwood. The pool 'psapw' is taken to represent 
! 'live wood'. This affects how turnover is formulated (see turnover.F)
#define _allocation_makela 1
! ----------------------------------------------------------------------

! //////////////////////////////////////////////////////////////////////
! SELECT MODEL INTEGRATION
! defining which model parts are to be calculated online and which to 
! be prescribed. Default is always to calculate it online (not prescr.)
! In general, if selected, daily values for the course of the whole 
! transient simulation are read in when initialising the model. This 
! requires input files formatted like ''
! This is only useful for site scale simulations. 
! ----------------------------------------------------------------------

! Prescribe GPP per day
#define _prescr_gpp_site 0

! //////////////////////////////////////////////////////////////////////
! SANITY
! Compile additional code for online sanity checks. This is recommended
! for simulations during development stage. De-activate later to reduce
! computational costs.
! ----------------------------------------------------------------------
#define _check_sanity 1
