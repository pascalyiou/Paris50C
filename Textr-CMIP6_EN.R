## Detection of days when tmax exceeds 48°C in the CMIP6 ensemble
## Par Pascal Yiou, LSCE, Mars 2023, Juillet 2023, August 2023
## This code is distributed "as is" under a CeCILL license:
## http://www.cecill.info/
## It can be downloaded and used for free for academic purposes.
## For commercial uses, please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr)

## Example use:
## R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/scripts_EN/Textr-CMIP6.R

## Requires:
## R CMD BATCH "--args tasmax JJA" ${HOME}/programmes/RStat/GREC50/scripts_EN/T_kstest_CMIP6.R
## to read: tasmax_JJA_CMIP6_ERA5.txt et tasmax_JJA_CMIP6_EOBS.txt
## which list CMIP6 models with OK simulations
## Also requires
## ${HOME}/programmes/RStat/GREC50/scripts_EN/GWD_var_CMIP6_EN.sh tas
## To compute global mean surface temperature time series for each model

## General path parameters for input and output
SI=Sys.info()
user=SI[["user"]]
if(SI[[1]] == "Linux"){
    Tdir=    paste("/scratchx/",user,"/CMIP6/IDF/",sep="") ## Input
    OUTdir = paste("/scratchx/",user,"/CMIP6/IDF/",sep="") ## Output
    GWDdir=  paste("/scratchx/",user,"/CMIP6/GWD/",sep="") ## Input GWD
}

## Requires ncdf4 and ncdf4.helpers libraries
## use:
## install.packages(c("ncdf4","ncdf4.helpers"),dependencies=TRUE)
library(ncdf4)
library(ncdf4.helpers)

args=(commandArgs(TRUE)) ## Reads input options
print(args)
i=1
if(length(args)>0){
    varname = args[i];i=i+1
## Temperature threshold to be exceeded    
    T0 =as.integer(args[i]);i=i+1
## We require a fluctuation threshold of at most "Tresh" °C
## before or after reaching 48°C
    Tthresh=as.integer(args[i]);i=i+1
## Time horizon    
    horizon=as.integer(args[i]);i=i+1
}else{
    varname = "tasmax"
## Temperature threshold to be exceeded (T0=48°C by default) 
    T0 = 48
## We require a fluctuation threshold of at most "Tresh=20" °C
## before or after reaching 48°C
    Tthresh=20
## Time horizon  (=2100 by default)  
    horizon = 2100
}

setwd(Tdir)
## List of models that pass the k-s test with ERA5
dum=readLines("tasmax_JJA_CMIP6_ERA5.txt",skip=1)
fi.ERA5=dum[2:length(dum)]
mod.ERA5=unique(unlist(strsplit(fi.ERA5,"_historical_"))[seq(1,2*length(fi.ERA5),by=2)])
## List of models that pass the k-s test with EOBS
dum=readLines("tasmax_JJA_CMIP6_EOBS.txt",skip=1)
fi.EOBS=dum[2:length(dum)]
mod.EOBS=unique(unlist(strsplit(fi.EOBS,"_historical_"))[seq(1,2*length(fi.EOBS),by=2)])
## List of models with at least one historical simulation passing at
## least one of the two tests
mod.union=union(mod.ERA5,mod.EOBS)

## List of all CMIP6 simulations
ls.fi = system(paste("ls ",varname,"_*_IdF.nc",sep=""),intern=TRUE)
outfit=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6.txt",sep="")
outfiR=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6.Rdat",sep="")

## Calculation of the # of OK simulations (passing the 2 k-s tests)
b=c()
a=strsplit(ls.fi,"_")
for(i in 1:length(a)){
    b=c(b,paste(a[[i]][1],a[[i]][2],a[[i]][3],sep="_"))
}
length(which(b %in% mod.union))

## Stat initialization for each model
## List of models
d1=strsplit(ls.fi,"_")
d2=unlist(d1)

l1=d2[seq(2,length(d2),by=6)]
l2=d2[seq(3,length(d2),by=6)]
l.mo=paste(l1,l2,sep="_")

## Correction 7/08/2023
l.mod=unique(l.mo)

CMIP.stat=list()
for(mod in l.mod){
    CMIP.stat[[mod]]=list(ntot=0,nT0=0,nensT0=0)
}

## Computation of stats of exceedances for each model
cat(file=outfit,paste("## Exceedances of",varname,"=",T0,"\n"))
T.exc=list()
for(fi in ls.fi){
    nam=strsplit(fi,".nc")[[1]]
## Extracting prefix information from simulation name
    r1=strsplit(nam,varname)[[1]][2]
    r2=strsplit(r1,"IdF")[[1]][1]
    r0=strsplit(nam,"_")
    varmod=paste(r0[[1]][1],r0[[1]][2],r0[[1]][3],sep="_")
## Does the simulation belong to the list of "chosen models"?    
    if(varmod %in% mod.union){
## Correction 7/08/2023
        modna=paste(r0[[1]][2],r0[[1]][3],sep="_")
## Count members for each model        
        CMIP.stat[[modna]]$ntot=CMIP.stat[[modna]]$ntot+1
## Read time series of model member        
        nc=nc_open(fi)
        time.dum = ncdf4.helpers::nc.get.time.series(nc)
        time.ref=as.numeric(format(time.dum, "%Y%m%d"))
        T.dum=ncvar_get(nc,varname)
        nc_close(nc)
## When does T.dum exceed T0 before the horizon        
        I.T0 = which(T.dum >= T0 & floor(time.ref/10000) <= horizon)
## Filter those days when the daily fluctuation is below Thresh        
        I0 = I.T0[T.dum[I.T0-1] >= T.dum[I.T0]-Tthresh]
        if(length(I0)>0){
## For each model: number of members at least one exceedance over T0            
            CMIP.stat[[modna]]$nensT0=CMIP.stat[[modna]]$nensT0+1
## For each model: number of exceedances of T0
            CMIP.stat[[modna]]$nT0=CMIP.stat[[modna]]$nT0+length(I0)
## Read file with global annual mean temperature (GST) for each simulation
            fiGWD=paste(GWDdir,"tas",r2,"GWD.nc",sep="")
            if(file.exists(fiGWD)){
                nc=nc_open(fiGWD)
                time.dum = ncdf4.helpers::nc.get.time.series(nc)
                yy.GWD=floor(as.numeric(format(time.dum, "%Y%m%d"))/10000)
                T.GWD=ncvar_get(nc,"tas")
                nc_close(nc)
## GST (=GWD) for events T > T0            
                GWD=c()
                for(t in I0) GWD=c(GWD,T.GWD[yy.GWD %in%
                                             floor(time.ref[t]/10000)])
                GWDdiff=GWD-T.GWD[1] ## GST increase since the first year
## Formatting output for Rstat and ascii files                
                T.exc[[nam]]=list(time=time.ref[I0], T=T.dum[I0],
                                  GWD=GWD,
                                  GWDdiff=GWDdiff,
                                  varmod=varmod,scen=r0[[1]][4],run=r0[[1]][5])
                cat(file=outfit,paste(fi,"\n"),append=TRUE)
                cat(file=outfit,time.ref[I0],append=TRUE)
                cat(file=outfit,"\n",append=TRUE)
                cat(file=outfit,GWD,append=TRUE)
                cat(file=outfit,"\n",append=TRUE)
                cat(file=outfit,GWDdiff,append=TRUE)
                cat(file=outfit,"\n",append=TRUE)
            }#end if file.exists
        }#end if length>0
    }
}#end for fi
save(file=outfiR,T.exc,CMIP.stat)

q("no")
