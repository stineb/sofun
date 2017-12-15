import pandas
import os
import os.path

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
## Get site/experiment names
##--------------------------------------------------------------------
filnam_siteinfo_csv = '../input_' + simsuite + '_sofun/experiments_' + simsuite + '_sofun.csv'

if os.path.exists( filnam_siteinfo_csv ):
    print 'reading site information file ' + filnam_siteinfo_csv + ' ...'
    siteinfo = pandas.read_csv( filnam_siteinfo_csv )
elif simsuite == 'example':
    print 'Executing single example simulation...'
else:
    print 'site info file does not exist: ' + filnam_siteinfo_csv

print 'writing file sitelist.txt'
fil = open('sitelist.txt', 'w')
for index, row in siteinfo.iterrows():
	fil.write( row['expname'] + '\n' )
