#!/bin/bash

proc_sitescale_site(){
	## Combines output files, written for each year, along time axis. For one specific site.

	## 1. argument runname/site name

	##-------------------------------------
	## Daily
	##-------------------------------------
	## GPP
	cdo -O mergetime output_nc/$1.*.d.wcont.nc   output_nc/$1.d.wcont.nc
	rm output_nc/$1.*.d.wcont.nc

	## AET
	cdo -O mergetime output_nc/$1.*.d.pet.nc   output_nc/$1.d.pet.nc
	rm output_nc/$1.*.d.pet.nc

	## PET
	cdo -O mergetime output_nc/$1.*.d.aet.nc   output_nc/$1.d.aet.nc
	rm output_nc/$1.*.d.aet.nc

	## WCONT
	cdo -O mergetime output_nc/$1.*.d.gpp.nc   output_nc/$1.d.gpp.nc
	rm output_nc/$1.*.d.gpp.nc

	## fAPAR
	cdo -O mergetime output_nc/$1.*.d.fapar.nc   output_nc/$1.d.fapar.nc
	rm output_nc/$1.*.d.fapar.nc

	## PPFD
	cdo -O mergetime output_nc/$1.*.d.ppfd.nc   output_nc/$1.d.ppfd.nc
	rm output_nc/$1.*.d.ppfd.nc

	## TEMPERATURE
	cdo -O mergetime output_nc/$1.*.d.temp.nc   output_nc/$1.d.temp.nc
	rm output_nc/$1.*.d.temp.nc


	##-------------------------------------
	## Annual
	##-------------------------------------
	## GPP
	cdo -O mergetime output_nc/$1.*.a.wcont.nc   output_nc/$1.a.wcont.nc
	rm output_nc/$1.*.a.wcont.nc

	## AET
	cdo -O mergetime output_nc/$1.*.a.pet.nc   output_nc/$1.a.pet.nc
	rm output_nc/$1.*.a.pet.nc

	## PET
	cdo -O mergetime output_nc/$1.*.a.aet.nc   output_nc/$1.a.aet.nc
	rm output_nc/$1.*.a.aet.nc

	## WCONT
	cdo -O mergetime output_nc/$1.*.a.gpp.nc   output_nc/$1.a.gpp.nc
	rm output_nc/$1.*.a.gpp.nc

	## fAPAR
	cdo -O mergetime output_nc/$1.*.a.fapar.nc   output_nc/$1.a.fapar.nc
	rm output_nc/$1.*.a.fapar.nc

	## PPFD
	cdo -O mergetime output_nc/$1.*.a.ppfd.nc   output_nc/$1.a.ppfd.nc
	rm output_nc/$1.*.a.ppfd.nc

	## TEMPERATURE
	cdo -O mergetime output_nc/$1.*.a.temp.nc   output_nc/$1.a.temp.nc
	rm output_nc/$1.*.a.temp.nc

	## ALPHA (AET/PET)
	cdo -O mergetime output_nc/$1.*.a.alpha.nc   output_nc/$1.a.alpha.nc
	rm output_nc/$1.*.a.alpha.nc	

	return 0
}

proc_global(){
	## Combines output files, written for each year, along time axis. For one specific site.

	## 1. argument runname/site name

	##-------------------------------------
	## Daily
	##-------------------------------------
	## GPP
	cdo -O mergetime output_nc/$1.*.d.gpp.nc   output_nc/$1.d.gpp.nc
	rm output_nc/$1.*.d.gpp.nc

	##-------------------------------------
	## Annual
	##-------------------------------------
	## AET
	cdo -O mergetime output_nc/$1.*.a.pet.nc   output_nc/$1.a.pet.nc
	rm output_nc/$1.*.a.pet.nc

	## PET
	cdo -O mergetime output_nc/$1.*.a.aet.nc   output_nc/$1.a.aet.nc
	rm output_nc/$1.*.a.aet.nc

	## ALPHA (AET/PET)
	cdo -O mergetime output_nc/$1.*.a.alpha.nc   output_nc/$1.a.alpha.nc
	rm output_nc/$1.*.a.alpha.nc	

	return 0
}

proc_sitescale_simsuite(){
	## Combines output files, written for each year, along time axis. For an entire simulation suite.
	## Uses ./sitelist.txt for a list of sites. Create this file using get_sitelist_simsuite.py

	python get_sitelist_simsuite.py

	sitelist=`cat sitelist.txt` 

	for idx in $sitelist
	do
		if [[ -e output_nc/$idx.d.gpp.nc ]]
		then
			echo "NetCDF output already processed for site $idx"
		else
			echo "PROCESSING SITE $idx"
			proc_sitescale_site $idx
		fi
	done

	return 0
}


