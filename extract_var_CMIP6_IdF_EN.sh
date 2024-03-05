#!/bin/env bash
## Extracting $varname over Ile de France Region in CMIP6 simulations
## Assumes that this code is stored in:
## ${HOME}/programmes/RStat/GREC50/script_EN/V2/extract_var_CMIP6_IdF_EN.sh
## Batch launch (on spiritx) with:
## sbatch --partition zen4 --ntasks=1 --cpus-per-task=10 --time 60:00:00 ${HOME}/programmes/RStat/GREC50/script_EN/V2/extract_var_CMIP6_EN_v2.sh varname operator subsca

## Pascal Yiou, LSCE, Feb. 2023, March 2023, January 2024, February 2024
## This code is distributed "as is" under a CeCILL license:
## http://www.cecill.info/
## It can be downloaded and used for free for academic purposes.
## For commercial uses, please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr)

## Environment configuration
#PBS -N CMIP6_download
#PBS -l nodes=1:ppn=1
#PBS -q week
#PBS -j oe
#PBS -m abe

## Use sbatch (on spiritx @ IPSL)
# Partition
#SBATCH --partition zen4
# Only one cpu 
#SBATCH --ntasks 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
# asking time could be in minute 120 or 2:00:00  or 1-00:00:00(24H)
#SBATCH --time 60:00:00 
# 
# to debug script 
## set -x
# purging all module to be sure to not having interferaence with current environnement
module purge
# loading only needed module for sample
module load cdo/2.0.6

## Your user name
owner=`whoami`

# slurm in environment variable SLURM_CPUS_PER_TASK
export NUM_CPUS=$SLURM_CPUS_PER_TASK

## Computation parameters
grid=${HOME}/programmes/RStat/GREC50/GRID_IDF
## GRID_IDF contains:
# gridtype  = lonlat
# gridsize  = 120
# xsize     = 12
# ysize     = 10
# xname     = longitude
# xlongname = "Longitude values"
# xunits    = "degrees_east"
# yname     = latitude
# ylongname = "Latitude values"
# yunits    = "degrees_north"
# xfirst    = 1.0
# xinc      = 0.25
# yfirst    = 48.0
# yinc      = 0.25
# scanningMode = 64
## GRID_WEUR contient:
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
opath=/scratchx/${owner}/CMIP6/${suffreg}

r0=("CMIP" "ScenarioMIP" "ScenarioMIP" "ScenarioMIP" "ScenarioMIP")
r1=("historical" "ssp126" "ssp245" "ssp370" "ssp585")
##r0=("ScenarioMIP")
##r1=("ssp534-over")

#Variable name in CMIP6 
##var="tasmax"
var=$1
oper=$2 ## mean, max or min, to use fld${oper} in cdo
subsca=$3

## Define the function to be executed in parallel
## Not used
run_parallel_task() {
    # local grid=$1
    # local opath=$2
    # local oper=$3
    local dumfile=$1
    local domfile=$2
    local outfile=$3
    cdo cat *.nc ${opath}/${dumfile}
## Interpolation sur la grille ERA5 (0.25 x 0.25) de l'IdF
    cdo remapbil,${grid} ${opath}/${dumfile} ${opath}/${domfile}
## Creation d'une serie temporelle quotidienne
    cdo subc,273.15 -fld${oper} -sellonlatbox,${lon1},${lon2},${lat1},${lat2} ${opath}/${domfile}  ${opath}/${outfile}_${oper}_IdF.nc
		
## Enlever les fichiers tres gros
    \rm ${opath}/${domfile} ${opath}/${dumfile}
}

## Alternative function for // execution
## Ajout Janvier 2024
## Not used
run_parallel_v2() {
    local datpath=$1
    local model=$2
    local r1=$3
    local run=$4
    local outfile=$5
    cd ${datpath}
    for fi in `ls *.nc`; do
	dum=$(echo ${fi} | sed "s/.nc//")_dum
	dom=$(echo ${fi} | sed "s/.nc//")_dom
	cdo remapbil,${grid} ${fi} ${opath}/${dum}.nc
	cdo -L subc,${subsca} -fld${oper} -sellonlatbox,${lon1},${lon2},${lat1},${lat2} ${opath}/${dum}.nc  ${opath}/${dom}.nc
    done
    cdo cat ${opath}/${var}_day_${model}_${r1[i]}_${run}_*_*_dom.nc ${opath}/${outfile}_${oper}_IdF.nc
## Enlever les fichiers tres gros
    \rm ${opath}/${var}_day_${model}_${r1[i]}_${run}_*_*_dom.nc ${opath}/${var}_day_${model}_${r1[i]}_${run}_*_*_dum.nc
}

## Alternative function for // execution
## // computing is done on cdo remapbil, which is computer greedy
## Ajout 18 Janvier 2024
run_sequence_v3() {
    local dumfile=$1
    local domfile=$2
    local outfile=$3
    cdo cat *.nc ${opath}/${dumfile}
## Interpolation sur la grille ERA5 (0.25 x 0.25) de l'IdF
    cdo -P ${NUM_CPUS} remapbil,${grid} ${opath}/${dumfile} ${opath}/${domfile}
## Creation d'une serie temporelle quotidienne
    cdo subc,${subsca} -fld${oper} -sellonlatbox,${lon1},${lon2},${lat1},${lat2}  ${opath}/${domfile}  ${opath}/${outfile}_${oper}_IdF.nc
		
## Remove big intermediate files
    \rm ${opath}/${domfile} ${opath}/${dumfile}
}

## Starting job loops
## Path for CMIP6 data on spiritx @ IPSL
cd /bdd/CMIP6

## Search over historical and four SSP scenarios, in lists defined by r0 and r1
ii=1
for i in "${!r0[@]}"; do
    echo "Starting process of ${r1[i]} files of ${oper} on ${var}"
    
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
		domfile=${var}_${institute}_${model}_${r1[i]}_${run}_${ii}.nc
		outfile=${var}_${institute}_${model}_${r1[i]}_${run}
##		outfile=$(echo ${dumfile} | sed "s/.nc//")
		if [ ! -f ${opath}/${outfile}_${oper}_IdF.nc ]; then
## Distribute tasks among available CPUs
		    task_id=$(( (ii - 1) % NUM_CPUS + 1 ))
		    run_sequence_v3 ${dumfile} ${domfile} ${outfile}
##		    run_parallel_task ${dumfile} ${domfile} ${outfile} &
##		    run_parallel_v2 /bdd/CMIP6/${path} ${model} ${r1[i]} ${run} ${outfile}
# 		    if (( task_id == NUM_CPUS )); then
# # Wait for tasks to complete before starting new ones
# 			wait
# 		    fi
		    (( ii++ ))
		fi
                cd /bdd/CMIP6
            done
	    cd /bdd/CMIP6
        done
    done    
done

wait
## End job
