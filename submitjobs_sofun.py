import pandas
import os
import os.path

##--------------------------------------------------------------------
## Simulation suite
##--------------------------------------------------------------------
simsuite = 'fluxnet'
# simsuite = 'pmodel_test'

##--------------------------------------------------------------------
## Compile
##--------------------------------------------------------------------
if simsuite == 'fluxnet' or simsuite == 'pmodel_test':
    if not os.path.exists( 'runpmodel' ):
        call(['make', 'pmodel'])

##--------------------------------------------------------------------
## Get site names
##--------------------------------------------------------------------
filnam_siteinfo_csv = '../input_' + simsuite + '_sofun/siteinfo_' + simsuite + '_sofun.csv'

if os.path.exists( filnam_siteinfo_csv ):
    print 'reading site information file ...'
    siteinfo_sub_dens = pandas.read_csv( '../input_' + simsuite + '_sofun/siteinfo_' + simsuite + '_sofun.csv' )
else:
    print 'site info file does not exist: ' + filnam_siteinfo_csv
    
##--------------------------------------------------------------------
## Loop over site names and submit job for each site
##--------------------------------------------------------------------
siteinfo = pandas.read_csv( filnam_siteinfo_csv )
for index, row in siteinfo.iterrows():
    if row['classid'] == 'GRA':
        print 'submitting job for site ' + row['mysitename'] + '...'
        os.system( 'echo ' + row['mysitename'] + '| ./runpmodel' )



# sitenames = siteinfo['mysitename']
# for idx in sitenames:
  # print 'submitting job for site ' + idx + '...'
  # os.system( 'echo ' + idx + '| ./runpmodel' )

# ## test: CH-Oe1 only
# print 'submitting job for site ' + 'CH-Oe1' + '...'
# os.system( 'echo ' + 'CH-Oe1' + '| ./runpmodel' )
