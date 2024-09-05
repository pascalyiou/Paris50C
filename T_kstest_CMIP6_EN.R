## Kolmogorov-Smirnov (K-S) test of CMIP6 data with reanalyses and observations
## to determine which simulations have a PDF that is deemed reasonable
## Par Pascal Yiou, LSCE, Mars 2023, August 2023, February 2024
## This code is distributed "as is" under a CeCILL license:
## http://www.cecill.info/
## It can be downloaded and used for free for academic purposes.
## For commercial uses, please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr)

##Example use:
## R CMD BATCH "--args tasmax JJA" ${HOME}/programmes/RStat/GREC50/scripts_EN/V2/T_kstest_CMIP6.R
## This R script requires that times series of CMIP6 models have been created
## by:
## ${HOME}/programmes/RStat/GREC50/scripts_EN/V2/extract_var_CMIP6_EN_v2.sh varname operator subsca

## This script creates two files that contain the list of CMIP6 files
## that pass a K_S test with ERA5 and EOBS

## General path parameters for input and output

SI=Sys.info()
user=SI[["user"]]
if(SI[[1]] == "Linux"){
    Rsource=paste("/home/",user,"/programmes/RStat/",sep="")
    Tdir=paste("/homedata/",user,"/CMIP6_IDF/",sep="")
    OUTdir = paste("/homedata/",user,"/CMIP6_IDF/",sep="") ## Needs to be adapted
}

library(ncdf4)
library(ncdf4.helpers)

varname="tasmax"
T0=48

args=(commandArgs(TRUE))
print(args)
i=1
if(length(args)>0){
    varname = args[i];i=i+1
    seas =args[i];i=i+1
}else{
    varname = "tas"
    seas = "JJA"
}

l.seas=list(JJA=c(6,7,8),SON=c(9,10,11),DJF=c(12,1,2),MAM=c(3,4,5))

setwd(Tdir)
## Liste des fichiers a tester
ls.fi = system(paste("ls ",varname,"_*_historical_*mean_IdF.nc",sep=""),
               intern=TRUE)

## Fichiers de sortie
outfit1=paste(varname,"_",seas,"_CMIP6_EOBS_2.txt",sep="")
outfit2=paste(varname,"_",seas,"_CMIP6_ERA5_2.txt",sep="")
outfit3=paste(varname,"_",seas,"_CMIP6_EOBS_m.txt",sep="")
outfit4=paste(varname,"_",seas,"_CMIP6_ERA5_m.txt",sep="")
##outfiR=paste(varname,"_",seas,"_CMIP6.Rdat",sep="")

## Donnees TX ERA5
nfi=paste("/homedata/",users,"/ERA5/iera5_tmax_daily_eu_1.75-3.25E_48.25-49.25N_n_su.nc",sep="")
nc=nc_open(nfi)
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.ERA=as.numeric(format(time.dum, "%Y%m%d"))
mm.ERA=floor(time.ERA/100) %% 100
yy.ERA=floor(time.ERA/10000)
TX.ERA=ncvar_get(nc,"tmax")
IX.ERA.seas = which(mm.ERA %in% l.seas[[seas]] & yy.ERA >= 1995 & yy.ERA <= 2014)
nc_close(nc)

## Donnees TG ERA5
nfi=paste("/homedata/",user,"/ERA5/iera5_t2m_daily_eu_1.75-3.25E_48.25-49.25N_n_su.nc",sep="")
nc=nc_open(nfi)
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.ERA=as.numeric(format(time.dum, "%Y%m%d"))
mm.ERA=floor(time.ERA/100) %% 100
yy.ERA=floor(time.ERA/10000)
TG.ERA=ncvar_get(nc,"t2m")
IG.ERA.seas = which(mm.ERA %in% l.seas[[seas]] & yy.ERA >= 1995 & yy.ERA <= 2014)
nc_close(nc)


## Donnees TX EOBS
nfi=paste("/homedata/",user,"/EOBS/iensembles_025_tx_1.75-3.25E_48.25-49.25N_n_su.nc",sep="")
nc=nc_open(nfi)
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.EOBS=as.numeric(format(time.dum, "%Y%m%d"))
mm.EOBS=floor(time.EOBS/100) %% 100
yy.EOBS=floor(time.EOBS/10000)
TX.EOBS=ncvar_get(nc,"tx")
IX.EOBS.seas = which(mm.EOBS %in% l.seas[[seas]] & yy.EOBS >= 1995 & yy.EOBS <= 2014)
nc_close(nc)

## Donnees TG EOBS
nfi=paste("/homedata/",user,"/EOBS/iensembles_025_tg_1.75-3.25E_48.25-49.25N_n.nc",sep="")
nc=nc_open(nfi)
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.EOBS=as.numeric(format(time.dum, "%Y%m%d"))
mm.EOBS=floor(time.EOBS/100) %% 100
yy.EOBS=floor(time.EOBS/10000)
TG.EOBS=ncvar_get(nc,"tg")
IG.EOBS.seas = which(mm.EOBS %in% l.seas[[seas]] & yy.EOBS >= 1995 & yy.EOBS <= 2014)
nc_close(nc)


## Donnees SAFRAN
nfi=paste("/homedata/",user,"/SAFRAN/SAFRAN_tmax_IdF_fldmean_Celsius.nc",sep="")
nc=nc_open(nfi)
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.SAFRAN=as.numeric(format(time.dum, "%Y%m%d"))
mm.SAFRAN=floor(time.SAFRAN/100) %% 100
yy.SAFRAN=floor(time.SAFRAN/10000)
T.SAFRAN=ncvar_get(nc,"tasmax")
I.SAFRAN.seas = which(mm.SAFRAN %in% l.seas[[seas]] & yy.SAFRAN >= 1995 & yy.SAFRAN <= 2014)
nc_close(nc)

T1d.ERA=tapply(TX.ERA,yy.ERA,max)
T1d.EOBS=tapply(TX.EOBS,yy.EOBS,max)
T1d.SAFRAN=tapply(T.SAFRAN,yy.SAFRAN,max)

if(varname=="tasmax"){
    T.EOBS=TX.EOBS
    I.EOBS.seas=IX.EOBS.seas
    T.ERA=TX.ERA
    I.ERA.seas=IX.ERA.seas
    oper=match.fun("max")
}
if(varname=="tas"){
    T.EOBS=TG.EOBS
    I.EOBS.seas=IG.EOBS.seas
    T.ERA=TG.ERA
    I.ERA.seas=IG.ERA.seas
    oper=match.fun("mean")
}
cat(file=outfit1,paste("## Models passing K-S test with EOBS\n"))
cat(file=outfit2,paste("## Models passing K-S test with ERA5\n"))
cat(file=outfit3,paste("## Models passing K-S test on max with EOBS\n"))
cat(file=outfit4,paste("## Models passing K-S test on max with ERA5\n"))
T.exc=list()
for(fi in ls.fi){
    nam=strsplit(fi,".nc")
    nc=nc_open(fi)
    namod=strsplit(nam[[1]],"_")[[1]][3]
    if(!(namod %in% c("UKESM1-0-LL"))){## Modeles connus pour avoir des erreurs
        time.dum = ncdf4.helpers::nc.get.time.series(nc)
        time.ref=as.numeric(format(time.dum, "%Y%m%d"))
        mm=floor(time.ref/100) %% 100
        yy=floor(time.ref/10000)
        T=ncvar_get(nc,varname)
        nc_close(nc)
        I.seas = which(mm %in% l.seas[[seas]] & yy >= 1995 & yy <= 2014)
        if(length(I.seas)>0){
 ## Test for all days of JJA       
            ks.T = ks.test(T[I.seas],TG.EOBS[I.EOBS.seas])
            if(ks.T$statistic <= 0.1){
                print(paste(fi,"EOBS"))
                cat(file=outfit1,paste(nam,"\n"),append=TRUE)
            }#end if ks.T
            ks.T = ks.test(T[I.seas],TG.ERA[I.ERA.seas])
            if(ks.T$statistic <= 0.1){
                print(paste(fi,"ERA5"))
                cat(file=outfit2,paste(nam,"\n"),append=TRUE)
            }#end if ks.T
## Test for the mean/max of JJA        
            ks.T = ks.test(tapply(T[I.seas],yy[I.seas],oper,na.rm=TRUE),
                           tapply(TG.EOBS[I.EOBS.seas],yy.EOBS[I.EOBS.seas],
                                  oper,na.rm=TRUE))
            if(ks.T$statistic <= 0.1){
                print(paste(fi,"EOBS"))
                cat(file=outfit3,paste(nam,"\n"),append=TRUE)
            }#end if ks.T
            ks.T = ks.test(tapply(T[I.seas],yy[I.seas],oper,na.rm=TRUE),
                           tapply(TG.ERA[I.ERA.seas],yy.ERA[I.ERA.seas],
                                  oper,na.rm=TRUE))
            if(ks.T$statistic <= 0.1){
                print(paste(fi,"ERA5"))
                cat(file=outfit4,paste(nam,"\n"),append=TRUE)
            }#end if ks.T
        }#end if length
    }
}
##save(file=outfiR,T.exc)
## Attention: enlever "a la main" les modeles:
## HadGEM3-GC31-HH, HadGEM3-GC31-HM, HadGEM3-GC31-LL, HadGEM3-GC31-LM, HadGEM3-GC31-MH, HadGEM3-GC31-MM, UKESM1.0-LL
q("no")
