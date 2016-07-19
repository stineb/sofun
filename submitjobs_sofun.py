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
##--------------------------------------------------------------------
simsuite = 'fluxnet_fixalloc'

##--------------------------------------------------------------------
## in some cases use same experiment info file
##--------------------------------------------------------------------
if simsuite == 'fluxnet_fixalloc':
    simsuite = 'fluxnet_cnmodel'

##--------------------------------------------------------------------
## Compile
##--------------------------------------------------------------------
if simsuite == 'fluxnet' or simsuite == 'pmodel_test':
    exe = 'runpmodel'
    compiling_opt = 'pmodel'
elif simsuite == 'fluxnet_cmodel' or simsuite == 'cmodel_test':
    exe = 'runcmodel'
    compiling_opt = 'cmodel'
elif simsuite == 'gcme' or simsuite == 'swissface' or simsuite == 'fluxnet_cnmodel' or simsuite == 'campi':
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
    if simsuite == 'fluxnet_cmodel':
        if row['classid'] == 'GRA':
            print 'submitting job for site ' + row['mysitename'] + '...'
            os.system( 'echo ' + row['mysitename'] + '| ./' + exe + '>' + row['mysitename'] + '.out' + '2>' + row['mysitename'] + '.out' )
    elif simsuite == 'gcme' or 'swissface' or simsuite == 'fluxnet_cnmodel' or simsuite == 'campi' or simsuite == 'fluxnet_fixalloc':
        print 'submitting job for experiment ' + row['expname'] + '...'
        os.system( 'echo ' + row['expname'] + '| ./' + exe+ '>' + row['expname'] + '.out' + '2>' + row['expname'] + '.out' )        
    else:
        print 'submitting job for site ' + row['mysitename'] + '...'
        os.system( 'echo ' + row['mysitename'] + '| ./' + exe + '>' + row['mysitename'] + '.out' + '2>' + row['mysitename'] + '.out')

