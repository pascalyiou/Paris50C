#!/bin/env bash
## Compute ENSO and AMO indices in the CMIP6 archive
## Pascal Yiou (LSCE), Aug. 2024. Based on a script by Robin Noyelle (LSCE)

## Use (on spiritx @ IPSL):
## sbatch --partition zen4 --time 60:00:00 ${HOME}/programmes/RStat/GREC50/ENSO_AMO_CMIP6_EN.sh 
## Adapt to your own server architecture
## This programme requires a reasonably recent version of cdo:
## https://code.mpimet.mpg.de/projects/cdo

## Pascal Yiou, LSCE,  August 2024
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
## module load cdo/2.3.0

## Your user name
owner=`whoami`

## Output directory (to be created before hand)
opath=/scratchx/${owner}/CMIP6/ENSO-AMO

r0=("CMIP" "ScenarioMIP" "ScenarioMIP" "ScenarioMIP" "ScenarioMIP")
r1=("historical" "ssp126" "ssp245" "ssp370" "ssp585")
##r0=("ScenarioMIP")
##r1=("ssp534-over")

#Variable name in CMIP6 
var="tos"
##var=$1

#--- AMO lat lon boxes--
lat1=0.0
lat2=60.0
lon1=280.0
lon2=360.0

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

                path="${r0[i]}"/$institute/$model/"${r1[i]}"/$run/Omon/$var/*/latest
		#Check if path exists
                if [ ! -d $path ]; then
                    continue
                fi

                cd $path
		dumfile=${var}_${institute}_${model}_${r1[i]}_${run}.nc
		outfile=$(echo ${dumfile} | sed "s/.nc//")
		if [ ! -f ${opath}/${outfile}_ENSO-AMO.nc ]; then
## File Concatenation (creates a big file)
		    cdo -b F32 -O cat *.nc ${opath}/${dumfile}
		    
## Create monthly ENSO index (Nino 3.4)
		    ## A revoir:
		    ##https://code.mpimet.mpg.de/boards/1/topics/14037
		    # ## Climatology
		    cdo -ymonmean -sellonlatbox,190.,240.,-5.,5. ${opath}/${dumfile} ${opath}/${outfile}_Nino3.4_clim.nc
		    cdo -sub  -sellonlatbox,190.,240.,-5.,5. ${opath}/${dumfile} ${opath}/${outfile}_Nino3.4_clim.nc ${opath}/${outfile}_Nino3.4_anom.nc
		    cdo div ${opath}/${outfile}_Nino3.4_anom.nc -timstd ${opath}/${outfile}_Nino3.4_anom.nc ${opath}/${outfile}_Nino3.4.nc
		    cdo -b F32 -O fldmean ${opath}/${outfile}_Nino3.4.nc ${opath}/${outfile}_Nino3.4-index.nc
## Remove temporary big files
		    \rm ${opath}/${outfile}_Nino3.4_clim.nc ${opath}/${outfile}_Nino3.4_anom.nc
		    ## Just mean-average SST over Nino3.4 region
		    cdo -b F32 -O fldmean -sellonlatbox,190,240,-5,5 -selname,tos ${opath}/${dumfile} ${opath}/${outfile}_ENSO.nc

		    ## Create monthly AMO Index
		    ## Average over the North Atlantic
		    cdo -b F32 -O fldmean -sellonlatbox,285,353,25,60 -selname,tos ${opath}/${dumfile} ${opath}/${outfile}_NA_sst.nc
		    cdo timmean -selname,tos ${opath}/${outfile}_NA_sst.nc ${opath}/${outfile}_NA_sst_mean.nc
		    cdo sub -selname,tos ${opath}/${outfile}_NA_sst.nc ${opath}/${outfile}_NA_sst_mean.nc ${opath}/${outfile}_NA_sstano.nc
		    ## Global average SST
		    cdo -b F32 -O fldmean -sellonlatbox,0,360,-60,60 -selname,tos ${opath}/${dumfile} ${opath}/${outfile}_glo_sst.nc
		    cdo timmean -selname,tos ${opath}/${outfile}_glo_sst.nc ${opath}/${outfile}_glo_sst_mean.nc
		    cdo sub -selname,tos ${opath}/${outfile}_glo_sst.nc ${opath}/${outfile}_glo_sst_mean.nc ${opath}/${outfile}_glo_sstnano.nc
		    ## AMO=SST_NA - SST_global
		    cdo sub -selname,tos ${opath}/${outfile}_NA_sstano.nc ${opath}/${outfile}_glo_sstnano.nc ${opath}/${outfile}_AMO.nc


## Remove temporary big files
		    \rm ${opath}/${dumfile} ${opath}/${outfile}_glo_sst.nc ${opath}/${outfile}_glo_sst_mean.nc ${opath}/${outfile}_glo_sstnano.nc ${opath}/${outfile}_NA_sst.nc ${opath}/${outfile}_NA_sst_mean.nc ${opath}/${outfile}_NA_sstano.nc
		fi
                cd /bdd/CMIP6
            done
	    cd /bdd/CMIP6
        done
    done    
done
