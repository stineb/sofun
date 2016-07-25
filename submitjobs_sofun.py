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
## - "fluxnet_fixalloc"
## - "atkin"
##--------------------------------------------------------------------
simsuite = 'atkin'

##--------------------------------------------------------------------
## set options
##--------------------------------------------------------------------
## standard outputs written to file (otherwise to screen)
out_to_file = True

##--------------------------------------------------------------------
## set model setup, given simsuite (shorter code below)
##--------------------------------------------------------------------
do_cnmodel = False
do_pmodel  = False
do_cmodel  = False
if simsuite == 'gcme' or simsuite == 'swissface' or simsuite == 'fluxnet_cnmodel' or simsuite == 'campi':
    do_cnmodel = True

if simsuite == 'fluxnet' or simsuite == 'pmodel_test':
    do_pmodel = True

if simsuite == 'fluxnet_cmodel' or simsuite == 'cmodel_test' or simsuite == 'atkin':
    do_cmodel = True


##--------------------------------------------------------------------
## in some cases use same experiment info file
##--------------------------------------------------------------------
if simsuite == 'fluxnet_fixalloc':
    simsuite = 'fluxnet_cnmodel'

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

    ## select only grasslands
    if simsuite == 'fluxnet_cmodel':
        if row['classid'] == 'GRA':
            print 'submitting job for site ' + row['mysitename'] + '...'
            os.system( 'echo ' + row['mysitename'] + '| ./' + exe + '>' + row['mysitename'] + '.out' + '2>' + row['mysitename'] + '.out' )

    ## normal case
    else:
        print 'submitting job for experiment ' + row['expname'] + '...'
        if out_to_file:
            os.system( 'echo ' + row['expname'] + '| ./' + exe+ '>' + row['expname'] + '.out' + '2>' + row['expname'] + '.out' )        
        else:
            os.system( 'echo ' + row['expname'] + '| ./' + exe )        
