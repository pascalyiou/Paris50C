#!/bin/sh -l
## Preparation de fichiers netcdf pour les cartes de SLP pendant
## les episodes a 50°C autour de Paris
## Pascal Yiou (LSCE), Mars 2024
## Se lance par:
## ${HOME}/programmes/RStat/GREC50/Paris50C/CMIP6_SLP_fileprep.sh &

module load cdo/2.0.6

## Your user name
owner=`whoami`

## Assumes that CMIP6 data have been downloaded there (OK for IPSL server)
CMIPdir=/bdd/CMIP6/ScenarioMIP
## Needs to be adjusted
OUTdir=/scratchx/${owner}/CMIP6/IDF

start_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\n\nStarting script at: ${start_date}\n"

pvar=psl
zvar=zg
zlev=50000
FRbox=-10,10,40,52
NAbox=-50,40,30,72

## This file contains a list of models/runs to be processed
infile=/mnt/homedafs-2.10/${owner}/CMIP6_IDF/tasmax_T48_dT20_horiz2100_CMIP6_3.txt
infile=${HOME}/programmes/RStat/GREC50/Paris50C/data/tasmax_T48_dT20_horiz2100_CMIP6_3.txt

# Verify that file exists
if [ ! -f ${infile} ]; then
    echo "Le fichier ${infile}  n'existe pas."
    exit 1
fi

# Read ${infile} line by line
while IFS= read -r line; do
    # Extraire le nom du fichier et la date
    mod1=$(echo $line | cut -d ' ' -f 2)
    mod2=$(echo $line | cut -d ' ' -f 3)
    scen=$(echo $line | cut -d ' ' -f 4)
    run=$(echo $line | cut -d ' ' -f 5)
    date=$(echo $line | cut -d ' ' -f 6)

    dirmod=${CMIPdir}/${mod1}/${mod2}/${scen}/${run}/day/${pvar}/*/latest

    # Vérifier si le repertoire existe
    if [ ! -d ${dirmod} ]; then
        echo "Le repertoire $dirmod n'existe pas."
        continue
    fi

    cd ${dirmod}

    for file in $(ls *.nc); do
	cdo  divc,100 -sellonlatbox,${NAbox} -seldate,${date} ${file} ${OUTdir}/${pvar}_${mod1}_${mod2}_${scen}_${run}_NA_${date}.nc
    done


    # Afficher un message indiquant le traitement du fichier
    echo "Extraction du champ de pression pour $date à partir de $file - OK"
done < ${infile}



end_date=`date +"%m/%d/%Y (%H:%M)"`
echo -e "\nScript finished at: ${end_date}\n"
