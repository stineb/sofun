from subprocess import call
import os
import os.path

##--------------------------------------------------------------------
## Simulation suite for site-scale simulations. Chose any of ...
## - "swissface"
## - "fluxnet"
## - "fluxnet2015"
## - "fluxnet_cnmodel"
## - "gcme"
## - "campi"
## - "campi_cmodel"
## - "fluxnet_fixalloc"
## - "atkin"
## - "atkinfull"
## - "olson"
## - "olson_cmodel"
## - "swbm"
## - "ameriwue"
##--------------------------------------------------------------------
## For global simulations, set simsuite to 'global'.
## - "global"
##--------------------------------------------------------------------
## This links NetCDF input files from directories mirrored locally from
## /work/bstocker/labprentice/data on Imperial's HPC CX1 server into the 
## input directory structure required for SOFUN.
##--------------------------------------------------------------------
simsuite = 'global'

##--------------------------------------------------------------------
## For an example simulation (simulation name 'EXAMPLE_global'), set 
## this to true 
##--------------------------------------------------------------------
example = False

##--------------------------------------------------------------------
## Manually et the root directory for the local mirror of 
## /work/bstocker/labprentice/data
##--------------------------------------------------------------------
dataroot = '/Users/benjaminstocker/data/'
# dataroot = '/Users/benjaminstocker/alphadata01/bstocker/data/'

##--------------------------------------------------------------------
## Link directories
##--------------------------------------------------------------------
## link output direcories
## ASCII output
os.system( 'unlink output' )
os.system( 'ln -svf ../output_' + simsuite + '_sofun output'  )

## NetCDF output
os.system( 'unlink output_nc' )
os.system( 'ln -svf ../output_nc_' + simsuite + '_sofun output_nc'  )


## link NetCDF input files for global simulations
if simsuite == 'global':
	##--------------------------------------
	## GLOBAL SIMULATIONS
	##--------------------------------------

	## Grid information
	##--------------------------------------
	dirn = 'input/global/grid'
	os.system( 'mkdir -p ' + dirn )

	## elevation
	call(['ln', '-svf', dataroot + 'watch_wfdei/WFDEI-elevation.nc', dirn ])

	## land masks at 1 deg and 0.5 deg resolution
	call(['ln', '-svf', dataroot + 'landmasks/gicew_1x1deg.cdf', dirn ])
	call(['ln', '-svf', dataroot + 'landmasks/gicew_halfdeg.cdf', dirn ])

	## CO2
	##--------------------------------------
	dirn = 'input/global/co2'
	call(['ln', '-svf', dataroot + 'co2/cCO2_rcp85_const850-1765.dat', dirn ])

	## fapar (fapar3g)
	##--------------------------------------
	dirn = 'input/global/fapar'
	os.system( 'mkdir -p ' + dirn )
	call(['ln', '-svf', dataroot + 'fAPAR/fAPAR3g_v2/fAPAR3g_v2_1982_2016_FILLED.nc', dirn ])


	## fapar (modis evi)
	##--------------------------------------
	dirn = 'input/global/fapar'
	os.system( 'mkdir -p ' + dirn )
	call(['ln', '-svf', dataroot + 'modis_monthly-evi/zmaw_data/modis_vegetation__LPDAAC__v5__0.5deg_FILLED.nc', dirn ])


	## WATCH-WFDEI climate input data
	##--------------------------------------
	dirn = './input/global/climate'
	if not os.path.isdir( dirn ):
		os.system( 'mkdir -p ' + dirn )

	## temperature
	src = dataroot + 'watch_wfdei/Tair_daily/*'
	dst = 'input/global/climate/temp'
	if not os.path.isdir( dst ):
		os.system( 'mkdir -p ' + dst )
	os.system( 'ln -svf ' + src + ' ' + dst )

	## precipitation (rain and snow)
	dst = 'input/global/climate/prec'
	if not os.path.isdir( dst ):
		os.system( 'mkdir -p ' + dst )

	src = dataroot + 'watch_wfdei/Rainf_daily/*'
	os.system( 'ln -svf ' + src + ' ' + dst )

	src = dataroot + 'watch_wfdei/Snowf_daily/*'
	os.system( 'ln -svf ' + src + ' ' + dst )

	## humidity (specific humidity in the case of WATCH-WFDEI)
	src = dataroot + 'watch_wfdei/Qair_daily/*'
	dst = 'input/global/climate/humd'
	if not os.path.isdir( dst ):
		os.system( 'mkdir -p ' + dst )
	os.system( 'ln -svf ' + src + ' ' + dst )

	## solar (shortwave) radiation
	src = dataroot + 'watch_wfdei/SWdown_daily/*'
	dst = 'input/global/climate/srad'
	if not os.path.isdir( dst ):
		os.system( 'mkdir -p ' + dst )
	os.system( 'ln -svf ' + src + ' ' + dst )


	## CRU climate input data (only ccov)
	##--------------------------------------
	## cloud cover
	src = dataroot + 'cru/ts_3.23/cru_ts3.23.1901.2014.cld.dat.nc'
	dst = 'input/global/climate/ccov'
	if not os.path.isdir( dst ):
		os.system( 'mkdir -p ' + dst )
	os.system( 'ln -svf ' + src + ' ' + dst )


	## Other directories
	##--------------------------------------
	## link 'run' and 'site_paramfils' directories
	if example:
		os.system( 'unlink run')
		os.system( 'unlink site_paramfils')
		os.system( 'ln -sv run_EXAMPLE run')
		os.system( 'ln -sv site_paramfils_EXAMPLE site_paramfils')
	else:
		os.system( 'unlink run')
		os.system( 'unlink site_paramfils')
		os.system( 'ln -sv ../input_' + simsuite + '_sofun/run run')
		os.system( 'ln -sv ../input_' + simsuite + '_sofun/site_paramfils site_paramfils')

else:
	##--------------------------------------
	## SITE-SCALE SIMULATIONS
	##--------------------------------------

	## use same site and simulation parameter files for cnmodel and cmodel simulations
	if simsuite == 'fluxnet_fixalloc':
	  simsuite = 'fluxnet_cnmodel'

	if simsuite == 'olson_cmodel':
	  simsuite = 'olson'

	if simsuite == 'campi_cmodel':
	  simsuite = 'campi'

	os.system( 'unlink run')
	os.system( 'unlink site_paramfils')
	os.system( 'unlink input/sitedata')
	os.system( 'ln -sv ../input_' + simsuite + '_sofun/run run')
	os.system( 'ln -sv ../input_' + simsuite + '_sofun/site_paramfils site_paramfils')
	os.system( 'ln -sv ../../input_' + simsuite + '_sofun/sitedata input/sitedata')




