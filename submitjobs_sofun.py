import pandas
import os
import os.path
from subprocess import call

##--------------------------------------------------------------------
## Simulation suite
## - "swissface"
## - "fluxnet"
## - "fluxnet_cnmodel"
## - "gcme"
## - "campi"
## - "campi_cmodel"
## - "fluxnet_fixalloc"
## - "atkin"
## - "olson"
## - "olson_cmodel"
## - "gradstep"
##--------------------------------------------------------------------
simsuite = 'swissface'

##--------------------------------------------------------------------
## set options
##--------------------------------------------------------------------
## standard outputs written to file (otherwise to screen)
out_to_file = True

## overwrite results
overwrite = True

## define default output variable to check for determining submission, given 'overwrite'
defaultvar = 'npp' 

##--------------------------------------------------------------------
## set model setup, given simsuite (shorter code below)
##--------------------------------------------------------------------
do_cnmodel = False
do_pmodel  = False
do_cmodel  = False

## C-N model setup
if simsuite == 'gcme' or simsuite == 'swissface' or simsuite == 'fluxnet_cnmodel' or simsuite == 'campi' or simsuite == 'olson' or simsuite == 'gradstep':
    do_cnmodel = True

## C-model setup (fixed allocation)
if simsuite == 'fluxnet_cmodel' or simsuite == 'cmodel_test' or simsuite == 'atkin' or simsuite == 'olson_cmodel' or simsuite == 'campi_cmodel':
    do_cmodel = True

## P-model setup
if simsuite == 'fluxnet' or simsuite == 'pmodel_test':
    do_pmodel = True


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
if do_pmodel:
    exe = 'runpmodel'
    compiling_opt = 'pmodel'
elif do_cmodel:
    exe = 'runcmodel'
    compiling_opt = 'cmodel'
elif do_cnmodel:
    exe = 'runcnmodel'
    compiling_opt = 'cnmodel'
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
else:
    print 'site info file does not exist: ' + filnam_siteinfo_csv
    
##--------------------------------------------------------------------
## Loop over site/experiment names and submit job for each site
##--------------------------------------------------------------------
for index, row in siteinfo.iterrows():

    ## set whether simulation should be submitted (and previous output overwritten)
    outfil = 'output/' + row['expname'] + '.a.' + defaultvar + '.out'
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


