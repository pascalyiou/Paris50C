#!/bin/env bash -l
## Extracting $varname over Ile de France Region in CMIP6 simulations
## Assumes that this code is stored in:
## ${HOME}/programmes/SCRIPTS/extract_var_CMIP6_IdF_EN.sh
## Launch in batch (on IPSL spiritx server):
## sbatch --partition zen4 --time 60:00:00 ${HOME}/programmes/SCRIPTS/extract_var_CMIP6_IdF_EN.sh varname
## Adapt to your own server architecture
## This programme requires a reasonably recent version of cdo:
## https://code.mpimet.mpg.de/projects/cdo

## Pascal Yiou, LSCE, Feb. 2023, March 2023, August 2023
## This code is distributed "as is" under a CeCILL license:
## http://www.cecill.info/
## It can be downloaded and used for free for academic purposes.
## For commercial uses, please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr)

## BATCH configuration options  (on spiritx @ IPSL)
#PBS -N CMIP6_download
#PBS -l nodes=1:ppn=1
#PBS -q week
#PBS -j oe
#PBS -m abe

## Use with sbatch (on spiritx @ IPSL)
# Partition
#SBATCH --partition zen4
# Only one cpu 
#SBATCH --ntasks 1 
# asking time could be in minute 120 or 2:00:00  or 1-00:00:00(24H)
#SBATCH --time 5:00:00 
# 
# to debug script 
set -x
# purging all module to be sure to not having interferaence with current environnement
module purge
# loading only needed module for sample
module load cdo/2.0.6

## Your user name
owner=`whoami`

## A grid file encompassing the region of interest should be created.
## The pathname should be adapted for your own use.
grid=${HOME}/programmes/RStat/GREC50/GRID_WEUR
## GRID_WEUR is an ascii file that contains the list of gridpoints
##  ERA5 (remove #)
# gridtype  = lonlat
# gridsize  = 8000
# xsize     = 100
# ysize     = 80
# xname     = longitude
# xlongname = "Longitude values"
# xunits    = "degrees_east"
# yname     = latitude
# ylongname = "Latitude values"
# yunits    = "degrees_north"
# xfirst    = -4.875
# xinc      = 0.25
# yfirst    = 40.125
# yinc      = 0.25
# scanningMode = 64

#1.75,3.25,48.25,49.25
## Longitude and latitude coordinates for Ile de France. Should be adapted
## For other regions.
lon1=1.75
lon2=3.25
lat1=48.25
lat2=49.25
## Region Suffix for output files (IdF=Ile de France)
suffreg=IdF

## Output directory. Depends on user's server architecture
opath=/scratchx/${owner}/CMIP6/IDF

r0=("CMIP" "ScenarioMIP" "ScenarioMIP" "ScenarioMIP" "ScenarioMIP")
r1=("historical" "ssp126" "ssp245" "ssp370" "ssp585")
##r0=("ScenarioMIP")
##r1=("ssp534-over")

#Variable name in CMIP6
## For a test
var="tasmax"
var=$1

##fname=${HOME}/${var}_CMIP6_list.txt

## This is where CMIP6 is stored on the IPSL server. Adapt to the architecture
## of your server.
cd /bdd/CMIP6

## Search over historical and four SSP scenarios, in lists defined by r0 and r1
for i in "${!r0[@]}"; do
    echo "Starting process of ${r1[i]} files for ${var}"
## Loop on list of institutions    
    for institute in $(ls ${r0[$i]}); do
## Loop on list of models run by institutions	
        for model in $(ls ${r0[$i]}/$institute); do
	    
            #Check if model has desired scenario
            if [ ! -d ${r0[$i]}/$institute/$model/${r1[$i]} ]; then
                continue
            fi

            for run in $(ls ${r0[$i]}/$institute/$model/${r1[$i]}); do
               	cd /bdd/CMIP6 
                r=${run%i*}
                r=${r:1}

                path="${r0[i]}"/$institute/$model/"${r1[i]}"/$run/day/$var/*/latest
		#Check if path exists
                if [ ! -d $path ]; then
                    continue
                fi

                cd $path
		dumfile=${var}_${institute}_${model}_${r1[i]}_${run}.nc
		outfile=$(echo ${dumfile} | sed "s/.nc//")
		if [ ! -f ${opath}/${outfile}_${suffreg}}.nc ]; then
## Concatenation of all temporal files: sometimes, there are several files
## of a few decades in directories. The formatting is not standard in CMIP6.
		cdo cat *.nc ${opath}/${dumfile}
## Bilinear Interpolation on ERA5 grid (0.25 x 0.25) de l'IdF
## from ${grid} file defined before
		cdo remapbil,${grid} ${opath}/${dumfile} ${opath}/${var}.nc
## Daily time series of maxima over IdF region (lon & lat)
## Conversion of K to °C
		cdo subc,273.15 -fldmax -sellonlatbox,${lon1},${lon2},${lat1},${lat2} ${opath}/${var}.nc  ${opath}/${outfile}_${suffreg}.nc
		

## Remove big intermediate files
		\rm ${opath}/${var}.nc ${opath}/${dumfile}
		fi
                cd /bdd/CMIP6
            done
	    cd /bdd/CMIP6
        done
    done    
done
