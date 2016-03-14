from subprocess import call
import os
import os.path

##--------------------------------------------------------------------
## Simulation suite
##--------------------------------------------------------------------
# simsuite = 'fluxnet'
simsuite = 'pmodel_test'

##--------------------------------------------------------------------
## Link directories
##--------------------------------------------------------------------
call(['ln', '-svf', '../../input_' + simsuite + '_sofun/sitedata', 'input'])
call(['ln', '-svf', '../input_' + simsuite + '_sofun/site_paramfils', '.'])
call(['ln', '-svf', '../input_' + simsuite + '_sofun/run', '.'])
