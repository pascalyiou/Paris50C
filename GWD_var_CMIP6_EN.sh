#!/bin/env bash -l
## Comutation of global annual mean surface temperatures in CMIP6
## Use (on spiritx @ IPSL):
## sbatch --partition zen4 --time 60:00:00 ${HOME}/programmes/RStat/GREC50/scripts_EN/GWD_var_CMIP6_EN.sh varname
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

## Output directory (to be created before hand)
opath=/scratchx/${owner}/CMIP6/GWD

r0=("CMIP" "ScenarioMIP" "ScenarioMIP" "ScenarioMIP" "ScenarioMIP")
r1=("historical" "ssp126" "ssp245" "ssp370" "ssp585")
##r0=("ScenarioMIP")
##r1=("ssp534-over")

#Variable name in CMIP6 
var="tas"
var=$1

## This is where CMIP6 is stored on the IPSL server. Adapt to your
## own server.
cd /bdd/CMIP6

for i in "${!r0[@]}"; do
    echo "Starting process of ${r1[i]} files for ${var}"
    
    for institute in $(ls ${r0[$i]}); do
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
		if [ ! -f ${opath}/${outfile}_IdF.nc ]; then
		## Concatenation
		    cdo cat *.nc ${opath}/${dumfile}
		    cdo subc,273.15 -yearmean -fldmean ${opath}/${dumfile} ${opath}/${outfile}_GWD.nc

## Remove temporary big files
		    \rm ${opath}/${dumfile}
		fi
                cd /bdd/CMIP6
            done
	    cd /bdd/CMIP6
        done
    done    
done
