from subprocess import call
import os
import os.path

##--------------------------------------------------------------------
## Simulation suite. Chose any of
## - "swissface"
## - "fluxnet"
## - "gcme"
## - "campi"
##--------------------------------------------------------------------
simsuite = 'swissface'

##--------------------------------------------------------------------
## Link directories
##--------------------------------------------------------------------
call(['ln', '-svf', '../../input_' + simsuite + '_sofun/sitedata', 'input'])
call(['ln', '-svf', '../input_' + simsuite + '_sofun/site_paramfils', '.'])
call(['ln', '-svf', '../input_' + simsuite + '_sofun/run', '.'])
os.system( 'rm output' )
call(['ln', '-svf', '../output_' + simsuite + '_sofun', 'output'])
