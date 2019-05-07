import pandas
import os
import os.path
from subprocess import call

##--------------------------------------------------------------------
## Simulation suite
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
##--------------------------------------------------------------------
## For global simulations, set simsuite to 'global'.
## This links NetCDF input files from directories mirrored locally from
## /work/bstocker/labprentice/data on Imperial's HPC CX1 server into the 
## input directory structure required for SOFUN.
##--------------------------------------------------------------------
## For an example simulation (simulation name 'EXAMPLE_global'), set 
## simsuite to 'example'. This should work after cloning this repo 
## from github.
##--------------------------------------------------------------------
simsuite = 'fluxnet2015'

##--------------------------------------------------------------------
## set options
##--------------------------------------------------------------------
## standard outputs written to file (otherwise to screen)
out_to_file = True

## overwrite results
overwrite = True

## Use SWBM water balance model instead of SPLASH
swbm = False
if swbm:
    print 'WARNING: submitting jobs with SWBM option!'
else:
    print 'WARNING: submitting jobs with SPLASH option!'

## define default output variable to check for determining submission, given 'overwrite'
defaultvar = 'gpp' 

##--------------------------------------------------------------------
## set model setup, given simsuite (shorter code below)
##--------------------------------------------------------------------
do_cnmodel     = False
do_gpmodel     = False
do_pmodel      = False
do_pmodel_swbm = False
do_cmodel      = False

## C-N model setup
if simsuite == 'gcme' or simsuite == 'swissface' or simsuite == 'fluxnet_cnmodel' or simsuite == 'campi' or simsuite == 'olson':
    do_cnmodel = True

## C-model setup (fixed allocation)
if simsuite == 'fluxnet_cmodel' or simsuite == 'cmodel_test' or simsuite == 'atkin' or simsuite == 'olson_cmodel' or simsuite == 'campi_cmodel':
    do_cmodel = True

## P-model setup
if simsuite == 'fluxnet' or simsuite == 'pmodel_test' or simsuite == 'atkinfull' or simsuite == 'fluxnet2015' or simsuite == 'ameriwue':
    if swbm:
        do_pmodel_swbm = True
    else:
        do_pmodel = True

## Example global P-model simulation
if simsuite == 'example':
    do_gpmodel = True

##--------------------------------------------------------------------
## in some cases use same experiment info file
##--------------------------------------------------------------------
if simsuite == 'fluxnet_fixalloc':
    simsuite = 'fluxnet_cnmodel'

if simsuite == 'olson_cmodel':
  simsuite = 'olson'

if simsuite == 'campi_cmodel':
  simsuite = 'campi'

##--------------------------------------------------------------------
## Compile
##--------------------------------------------------------------------
if do_pmodel_swbm:
    exe = 'runpmodel_swbm'
    compiling_opt = 'pmodel_swbm'
elif do_pmodel:
    exe = 'runpmodel'
    compiling_opt = 'pmodel'
elif do_cmodel:
    exe = 'runcmodel'
    compiling_opt = 'cmodel'
elif do_cnmodel:
    exe = 'runcnmodel'
    compiling_opt = 'cnmodel'
elif do_gpmodel:
    exe = 'rungpmodel'
    compiling_opt = 'gpmodel'
else:
    print 'simsuite not valid'

if not os.path.exists( exe ):
    call(['make', compiling_opt])

##--------------------------------------------------------------------
## Get site/experiment names
##--------------------------------------------------------------------
filnam_siteinfo_csv = '../input_' + simsuite + '_sofun/experiments_' + simsuite + '_sofun.csv'

if os.path.exists( filnam_siteinfo_csv ):
    print 'reading site information file ' + filnam_siteinfo_csv + '...'
    siteinfo = pandas.read_csv( filnam_siteinfo_csv )
elif simsuite == 'example':
    print 'Executing single example simulation...'
else:
    print 'site info file does not exist: ' + filnam_siteinfo_csv
    
##--------------------------------------------------------------------
## Loop over site/experiment names and submit job for each site
##--------------------------------------------------------------------
if simsuite == 'example':

    ## Single example simulation, global, P-model only
    os.system( 'echo EXAMPLE_global | ./' + exe )

else:
    for index, row in siteinfo.iterrows():

        ## set whether simulation should be submitted (and previous output overwritten)
        outfil = 'output_nc/' + row['expname'] + '.d.' + defaultvar + '.nc'
        if overwrite:
            do_submit = True
        else: 
            if os.path.exists( outfil ):
                do_submit = False
            else:
                do_submit = True

        ## select only grasslands
        if simsuite == 'fluxnet_cmodel':
            if row['classid'] == 'GRA':
                print 'submitting job for site ' + row['mysitename'] + '...'
                os.system( 'echo ' + row['mysitename'] + '| ./' + exe + '>' + row['mysitename'] + '.out' + '2>' + row['mysitename'] + '.out' )

        ## normal case
        else:
            if do_submit:
                print 'submitting job for experiment ' + row['expname'] + '...'
                if out_to_file:
                    outfile = row['expname'] + '.out'
                    cmd = 'echo ' + row['expname'] + '| ./' + exe+ ' >' + outfile + ' 2>' + outfile
                    os.system( cmd )        
                else:
                    cmd = 'echo ' + row['expname'] + '| ./' + exe
                    os.system( cmd )   

            else:
                print 'NOT submitting for experiment ' + row['expname'] + '...'

        # raw_input('Press any key to continue')



