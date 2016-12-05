from subprocess import call
import os
import os.path

##--------------------------------------------------------------------
## Simulation suite. Chose any of
## - "swissface"
## - "gradstep"
## - "fluxnet"
## - "fluxnet_cnmodel"
## - "gcme"
## - "campi"
## - "campi_cmodel"
## - "fluxnet_fixalloc"
## - "atkin"
## - "olson"
## - "olson_cmodel"
##--------------------------------------------------------------------
simsuite = 'gradstep'

##--------------------------------------------------------------------
## Link directories
##--------------------------------------------------------------------
## link output direcories
os.system( 'rm output' )
os.system( 'mkdir ../output_' + simsuite + '_sofun' )
call(['ln', '-svf', '../output_' + simsuite + '_sofun', 'output'])

## use same site and simulation parameter files for cnmodel and cmodel simulations
if simsuite == 'fluxnet_fixalloc':
  simsuite = 'fluxnet_cnmodel'

if simsuite == 'olson_cmodel':
  simsuite = 'olson'

if simsuite == 'campi_cmodel':
  simsuite = 'campi'

call(['ln', '-svf', '../../input_' + simsuite + '_sofun/sitedata', 'input'])
call(['ln', '-svf', '../input_' + simsuite + '_sofun/site_paramfils', '.'])
call(['ln', '-svf', '../input_' + simsuite + '_sofun/run', '.'])


