from subprocess import call

##--------------------------------------------------------------------
## Copy standard files (*_std). Standard files are part of the 
## repository, while their copies are usually modified frequently
## without making it necessary to update respective standard files
## each time
##--------------------------------------------------------------------
call(['cp', 'Makefile_std', 'Makefile'])
