import pandas
import os
import os.path
from subprocess import call

##--------------------------------------------------------------------
## Simulation suite
##--------------------------------------------------------------------
simsuite = 'fluxnet_cmodel'
# simsuite = 'pmodel_test'

##--------------------------------------------------------------------
## Compile
##--------------------------------------------------------------------
if simsuite == 'fluxnet' or simsuite == 'pmodel_test':
    exe = 'runpmodel'
    compiling_opt = 'pmodel'
elif simsuite == 'fluxnet_cmodel' or simsuite == 'cmodel_test':
    exe = 'runcmodel'
    compiling_opt = 'cmodel'
else:
    print 'simsuite not valid'

if not os.path.exists( exe ):
    call(['make', compiling_opt])

##--------------------------------------------------------------------
## Get site names
##--------------------------------------------------------------------
filnam_siteinfo_csv = '../input_' + simsuite + '_sofun/siteinfo_' + simsuite + '_sofun.csv'

if os.path.exists( filnam_siteinfo_csv ):
    print 'reading site information file ...'
    siteinfo = pandas.read_csv( filnam_siteinfo_csv )
else:
    print 'site info file does not exist: ' + filnam_siteinfo_csv
    
##--------------------------------------------------------------------
## Loop over site names and submit job for each site
##--------------------------------------------------------------------
for index, row in siteinfo.iterrows():
    # print 'submitting job for site ' + row['mysitename'] + '...'
    # os.system( 'echo ' + row['mysitename'] + '| ./' + exe )
    if row['classid'] == 'GRA':
        print 'submitting job for site ' + row['mysitename'] + '...'
        os.system( 'echo ' + row['mysitename'] + '| ./' + exe )
