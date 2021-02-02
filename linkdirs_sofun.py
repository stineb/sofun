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
## - "fluxnet2015_cmodel"
##--------------------------------------------------------------------
## For global simulations, set name to 'global'.
## - "global"
##--------------------------------------------------------------------
## This links NetCDF input files from directories mirrored locally from
## /work/bstocker/labprentice/data on Imperial's HPC CX1 server into the 
## input directory structure required for SOFUN.
##--------------------------------------------------------------------
name = 'global'
# name = 'fluxnet2015'


##--------------------------------------------------------------------
## Manually edit the root directory for the local mirror of 
## the data directory (e.g., /cluster/home/bestocke/data on Euler.
##--------------------------------------------------------------------
dataroot = '/Users/benjaminstocker/data/'     # to run on Beni's iMac
# dataroot = '/cluster/work/climate/bestocke/data/'     # to run on Euler
# dataroot = '/Users/bestocke/data/'     # to run on Beni's Laptop

##--------------------------------------------------------------------
## Manually edit the root directory for SOFUN inputS 
##--------------------------------------------------------------------
inputroot = '~/sofun_inputs'
outputroot = '~/sofun_outputs'

##--------------------------------------------------------------------
## Link directories
##--------------------------------------------------------------------
## link output direcory
os.system( 'unlink output_nc' )
dir = outputroot + '/output_nc_' + name

if not(os.path.isdir(dir)):
	print('WARNING: creating output directory that has not existed before: ' + dir)
	os.system('mkdir ' + dir)

os.system( 'ln -svf ' + dir + ' output_nc'  )

## link 'run' directory (containing simulation parameter files)
os.system( 'unlink run')
if name == 'example':
	os.system( 'ln -sv run_EXAMPLE run')
else:
	os.system( 'ln -svf ' + inputroot + '/input_' + name + '_sofun/run run')


## link 'site_paramfils' directory (containing site parameter files) 
os.system( 'unlink site_paramfils')
if name == 'example':
	os.system( 'ln -sv site_paramfils_EXAMPLE site_paramfils')
else:
	os.system( 'ln -svf ' + inputroot + '/input_' + name + '_sofun/site_paramfils site_paramfils')


##--------------------------------------------------------------------
## Copy parameter files
##--------------------------------------------------------------------
os.system( 'mkdir params' )
os.system( 'cp params_std/* params' )


## link NetCDF input files for global simulations
if name == 'global':
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
	os.system( 'mkdir -p ' + dirn )
	call(['ln', '-svf', dataroot + 'co2/cCO2_rcp85_const850-1765.dat', dirn ])

	## fapar (fapar3g)
	##--------------------------------------
	dirn = 'input/global/fapar'
	os.system( 'mkdir -p ' + dirn )

	# call(['ln', '-svf', dataroot + 'fAPAR/fAPAR3g_v2/fAPAR3g_v2_1982_2016_FILLED.nc', dirn ])
	# os.system('echo \'fapar_forcing_source                    fapar3g\' >./input/dfapar_source.txt')

	# call(['ln', '-svf', dataroot + 'modis_ndvi_evi_zmaw/halfdeg/modis_vegetation__LPDAAC__v5__0.5deg_FILLED.nc', dirn ])
 #        os.system('echo \'fapar_forcing_source                    evi_modis\' >./input/dfapar_source.txt')

	# call(['ln', '-svf', dataroot + 'modis_lai_fpar_zmaw/MODIS-C006_MOD15A2__LAI_FPAR__LPDAAC__GLOBAL_0.5degree__UHAM-ICDC__2000_2018__MON__fv0.02.nc', dirn ])
	# os.system('echo \'fapar_forcing_source                    fpar_modis\' >./input/dfapar_source.txt')


	## soil
	##--------------------------------------
	dirn = 'input/global/soil'
	os.system( 'mkdir -p ' + dirn )
	call(['ln', '-svf', dataroot + 'soil/soilgrids/whc_soilgrids_halfdeg_FILLED.nc', dirn ])
	call(['ln', '-svf', dataroot + 'soil/hwsd/soil_type_hwsd_halfdeg.cdf', dirn ])

	## land cover
	##--------------------------------------
	dirn = 'input/global/landcover'
	os.system( 'mkdir -p ' + dirn )
	# call(['ln', '-svf', dataroot + 'landcover/modis_landcover_halfdeg_2010_FILLED.nc', dirn ])
	call(['ln', '-svf', dataroot + 'c4_still/final/c4_percentage.nc', dirn ])

	## nimpl predictors
	##--------------------------------------
	src = dataroot + 'nimpl_sofun_inputs/map/Final_ncfile/*'
	dst = 'input/global/nimpl'
	if not os.path.isdir( dst ):
		os.system( 'mkdir -p ' + dst )
	os.system( 'ln -svf ' + src + ' ' + dst )

	## c4 percentage
	##--------------------------------------
	src = dataroot + 'c4_still/final/c4_percentage_revlat.nc'
	dst = 'input/global/landcover'
	if not os.path.isdir( dst ):
		os.system( 'mkdir -p ' + dst )
	os.system( 'ln -svf ' + src + ' ' + dst )	

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
	dst = 'input/global/climate/ccov'
	if not os.path.isdir( dst ):
		os.system( 'mkdir -p ' + dst )

	## cloud cover
	src = dataroot + 'cru/ts_4.01/cru_ts4.01.1901.2016.cld.dat.nc'
	os.system( 'ln -svf ' + src + ' ' + dst )

	## daily minimum temperature
	src = dataroot + 'cru/ts_4.01/cru_ts4.01.1901.2016.tmn.dat.nc'
	os.system( 'ln -svf ' + src + ' ' + dst )

	## daily maximum temperature
	src = dataroot + 'cru/ts_4.01/cru_ts4.01.1901.2016.tmx.dat.nc'
	os.system( 'ln -svf ' + src + ' ' + dst )

elif name == 'example':

	# os.system( 'ln -sv ./input_EXAMPLE_sofun/sitedata input')
	here = os.getcwd()
	os.system( 'ln -svf ' + here + '/input_EXAMPLE_sofun/sitedata/ input')
	os.system( 'cp input_EXAMPLE_sofun/dfapar_source.txt input/' )

else:

	os.system( 'unlink input/global')
	os.system( 'ln -svf ' + inputroot + '/input_' + name + '_sofun/sitedata input/sitedata')


