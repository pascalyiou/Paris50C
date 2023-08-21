## Kolmogorov-Smirnov (K-S) test of CMIP6 data with reanalyses and observations
## to determine which simulations have a PDF that is deemed reasonable
## Par Pascal Yiou, LSCE, Mars 2023, August 2023
## This code is distributed "as is" under a CeCILL license:
## http://www.cecill.info/
## It can be downloaded and used for free for academic purposes.
## For commercial uses, please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr)

##Example use:
## R CMD BATCH "--args tasmax JJA" ${HOME}/programmes/RStat/GREC50/scripts_EN/T_kstest_CMIP6.R
## This R script requires that times series of CMIP6 models have been created
## by:
## ${HOME}/programmes/RStat/GREC50/scripts_EN/extract_var_CMIP6_IdF_EN.sh varname

## This script creates two files that contain the list of CMIP6 files
## that pass a K_S test with ERA5 and EOBS

## General path parameters for input and output
SI=Sys.info()
user=SI[["user"]]
if(SI[[1]] == "Linux"){## Needs to be adapted. This is valid for spiritx @ IPSL
    Tdir=paste("/scratchx/",user,"/CMIP6/IDF/",sep="") ## Input
    OUTdir = paste("/scratchx/",user,"/CMIP6/IDF/",sep="") ## Output
}

## Requires ncdf4 and ncdf4.helpers libraries
## use:
## install.packages(c("ncdf4","ncdf4.helpers"),dependencies=TRUE)
library(ncdf4)
library(ncdf4.helpers)

T0=48

args=(commandArgs(TRUE)) ## Reads input options
print(args)
i=1
if(length(args)>0){
    varname = args[i];i=i+1
    seas =args[i];i=i+1
}else{## Default options
    varname = "tasmax"
    seas = "JJA"
}

l.seas=list(JJA=c(6,7,8),SON=c(9,10,11),DJF=c(12,1,2),MAM=c(3,4,5))

setwd(Tdir)
## List of files to test
ls.fi = system(paste("ls ",varname,"_*_historical_*_IdF.nc",sep=""),intern=TRUE)
outfit1=paste(varname,"_",seas,"_CMIP6_EOBS.txt",sep="")
outfit2=paste(varname,"_",seas,"_CMIP6_ERA5.txt",sep="")
##outfiR=paste(varname,"_",seas,"_CMIP6.Rdat",sep="")

##  ERA5 data that are provided with on the github package
## Data were obtained from the Climate Exporer
nfi=paste("/scratchx/",user,"/ERA5/iera5_tmax_daily_eu_1.75-3.25E_48.25-49.25N_n_max.nc",sep="")
nc=nc_open(nfi)
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.ERA=as.numeric(format(time.dum, "%Y%m%d"))
mm.ERA=floor(time.ERA/100) %% 100
yy.ERA=floor(time.ERA/10000)
T.ERA=ncvar_get(nc,"tmax")-273
I.ERA.seas = which(mm.ERA %in% l.seas[[seas]] & yy.ERA >= 1995 & yy.ERA <= 2014)
nc_close(nc)

##  EOBS data that are provided with on the github package
## Data were obtained from the Climate Exporer
nc=nc_open(paste("/scratchx/",user,"/EOBS/iensembles_025_tx_1.75-3.25E_48.25-49.25N_n_max.nc",sep=""))
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.EOBS=as.numeric(format(time.dum, "%Y%m%d"))
mm.EOBS=floor(time.EOBS/100) %% 100
yy.EOBS=floor(time.EOBS/10000)
T.EOBS=ncvar_get(nc,"tx")
I.EOBS.seas = which(mm.EOBS %in% l.seas[[seas]] & yy.EOBS >= 1995 & yy.EOBS <= 2014)
nc_close(nc)

## Output files
cat(file=outfit1,paste("## Models passing K-S test with EOBS\n"))
cat(file=outfit2,paste("## Models passing K-S test with ERA5\n"))
T.exc=list()
for(fi in ls.fi){
    nam=strsplit(fi,".nc")
    nc=nc_open(fi)
    time.dum = ncdf4.helpers::nc.get.time.series(nc)
    time.ref=as.numeric(format(time.dum, "%Y%m%d"))
    mm=floor(time.ref/100) %% 100
    yy=floor(time.ref/10000)
    T=ncvar_get(nc,varname)
    nc_close(nc)
## The years can be changed. This is a climatology of 20 years    
    I.seas = which(mm %in% l.seas[[seas]] & yy >= 1995 & yy <= 2014)
    if(length(I.seas)>0){
        ks.T = ks.test(T[I.seas],T.EOBS[I.EOBS.seas])
## The criterion of the ks test can be changed        
        if(ks.T$statistic <= 0.1){ 
            print(paste(fi,"EOBS"))
            cat(file=outfit1,paste(nam,"\n"),append=TRUE)
        }
        ks.T = ks.test(T[I.seas],T.ERA[I.ERA.seas])
## The criterion of the ks test can be changed        
        if(ks.T$statistic <= 0.1){
            print(paste(fi,"ERA5"))
            cat(file=outfit2,paste(nam,"\n"),append=TRUE)
        }
    }   
}
##save(file=outfiR,T.exc)

q("no")
