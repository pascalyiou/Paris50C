## Detection of days when tmax exceeds 48°C in the CMIP6 ensemble
## Par Pascal Yiou, LSCE, Mars 2023, Juillet 2023, Janvier 2024
## This code is distributed "as is" under a CeCILL license:
## http://www.cecill.info/
## It can be downloaded and used for free for academic purposes.
## For commercial uses, please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr)

## Example use:
## R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/scripts_EN/V2/Textr-CMIP6_v2_EN.R

## Requires:
## R CMD BATCH "--args tasmax JJA" ${HOME}/programmes/RStat/GREC50/scripts_EN/V2/T_kstest_CMIP6_EN.R
## to read: tasmax_JJA_CMIP6_ERA5.txt et tasmax_JJA_CMIP6_EOBS.txt
## which list CMIP6 models with OK simulations
## Also requires
## ${HOME}/programmes/RStat/GREC50/scripts_EN/V2/GMST_var_CMIP6_v2_EN.sh tas
## To compute global mean surface temperature time series for each model

SI=Sys.info()
user=SI[["user"]]
if(SI[[1]] == "Linux"){
    Tdir=    paste("/homedata/",user,"/CMIP6_IDF/",sep="") ## Input
    Rdir=    paste("/homedata/",user,"/CMIP6_IDF/",sep="") ## Input
    OUTdir = paste("/scratchx/",user,"/CMIP6/IDF/",sep="") ## Output
    GWDdir=  paste("/scratchx/",user,"/CMIP6/GWD/",sep="") ## Input GWD
}

library(ncdf4)
library(ncdf4.helpers)

args=(commandArgs(TRUE))
print(args)
i=1
if(length(args)>0){
    varname = args[i];i=i+1
## Temperature threshold to be exceeded    
    T0 =as.integer(args[i]);i=i+1
## We require a fluctuation threshold of at most "Tresh" °C
## before or after reaching 48°C
    Tthresh=as.integer(args[i]);i=i+1
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

setwd(Rdir)
## List of models that pass the k-s test with ERA5
dum=readLines("tas_JJA_CMIP6_ERA5_2.txt",skip=1)
fi.ERA5=dum[2:length(dum)]
mod.ERA5=unique(unlist(strsplit(fi.ERA5,"_historical_"))[seq(1,2*length(fi.ERA5),by=2)])
## List of models that pass the k-s test with EOBS
dum=readLines("tas_JJA_CMIP6_EOBS_2.txt",skip=1)
fi.EOBS=dum[2:length(dum)]
mod.EOBS=unique(unlist(strsplit(fi.EOBS,"_historical_"))[seq(1,2*length(fi.EOBS),by=2)])
## List of models with at least one historical simulation passing at
## least one of the two tests
mod.union=union(mod.ERA5,mod.EOBS)
a=strsplit(mod.union,"_")
u.mod.union=c()
for(i in 1:length(a)){
    u.mod.union = c(u.mod.union,paste(a[[i]][2],a[[i]][3],sep="_"))
}

## List of all available CMIP6 simulations
setwd(Tdir)
ls.fi = system(paste("ls ",varname,"_*_*_max_IdF.nc",sep=""),intern=TRUE)
outfit=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6_2.txt",sep="")
outfiR=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6_2.Rdat",sep="")

## Calculation of the # of OK simulations (passing the 2 k-s tests)
b=c()
a=strsplit(ls.fi,"_")
for(i in 1:length(a)){
    b=c(b,paste(a[[i]][2],a[[i]][3],sep="_"))
}
length(which(b %in% u.mod.union))

## Stat initialization for each model
## List of models
d1=strsplit(ls.fi,"_")
d2=unlist(d1)

l1=d2[seq(2,length(d2),by=7)]
l2=d2[seq(3,length(d2),by=7)]
l.mo=paste(l1,l2,sep="_")

## Correction 7/08/2023
l.mod=unique(l.mo)

l.mod2=unlist(strsplit(l.mod,"_"))[seq(2,length(l.mod),by=2)]

## Statistique generale sur les simulations disponibles
CMIP.glob.stat=list()
for(imod in l.mod){
    mod.na=strsplit(imod,"_")[[1]][2]
    dum.mod=list()
    for(scen in c("historical","ssp126","ssp245","ssp370","ssp585")){
        n.mod.scen=length(list.files(pattern=
                                         paste(varname,"_",imod,"_",scen,
                                               sep="")))
        dum.mod[[scen]]=n.mod.scen
    }
    CMIP.glob.stat[[mod.na]]=dum.mod
}

nsim=0
for(mod.na in sort(names(CMIP.glob.stat))){
    I.all=1
    idum=0
    for(i in 1:length(CMIP.glob.stat[[mod.na]])){
        I.all=I.all*CMIP.glob.stat[[mod.na]][[i]]
        idum=idum+CMIP.glob.stat[[mod.na]][[i]]
        }
    if(I.all > 0){
        for(i in 1:length(CMIP.glob.stat[[mod.na]]))
            cat(mod.na,"\t",names(CMIP.glob.stat[[mod.na]])[i],"\t",
                CMIP.glob.stat[[mod.na]][[i]],"\n")
        nsim=nsim+idum
        
    }
}

CMIP.stat=list()
for(mod in l.mod){
    CMIP.stat[[mod]]=list(ntot=0,nT0=0,nensT0=0)
}

##Beginning of "big loop"
## outfit=Text output
cat(file=outfit,paste("## Exceedances of",varname,"=",T0,"\n"))
GMT2000=c()
T.exc=list()
T.JJA=list()
for(fi in ls.fi){
    nam=strsplit(fi,".nc")[[1]]
## Fichier correspondant pour la moyenne annuelle globale de tas    
    r1=strsplit(nam,varname)[[1]][2]
    r2=strsplit(r1,"IdF")[[1]][1]
    r3=strsplit(r2,"max")[[1]][1]
    r0=strsplit(nam,"_")
    scen=r0[[1]][4]
    nmod=paste(r0[[1]][2],r0[[1]][3],sep="_")
    varmod=paste(r0[[1]][1],r0[[1]][2],r0[[1]][3],sep="_")
## Does the simulation belong to the list of "chosen models"?    
    if(nmod %in% u.mod.union){
## Correction 7/08/2023
        modna=paste(r0[[1]][2],r0[[1]][3],sep="_")
## Count members for each model        
        CMIP.stat[[modna]]$ntot=CMIP.stat[[modna]]$ntot+1 
        fiGWD.hist=paste(GWDdir,"tas_",r0[[1]][2],"_",r0[[1]][3],"_historical_",
                         r0[[1]][5],"_GWD.nc",sep="")
        fiGWD.scen=paste(GWDdir,"tas",r3,"GWD.nc",sep="")
## Read time series of model member        
        file.inf=file.info(fi)
        if(file.inf$size>100000){
            nc=nc_open(fi)
            time.dum = ncdf4.helpers::nc.get.time.series(nc)
            time.ref=as.numeric(format(time.dum, "%Y%m%d"))
            T.dum=ncvar_get(nc,varname)
            nc_close(nc)
## When does T.dum exceed T0 before the horizon        
            I.T0 = which(T.dum >= T0 & floor(time.ref/10000) <= horizon)
## Medianes JJA pour chaque annee
            mm=floor(time.ref/100) %% 100
            yy=floor(time.ref/10000)
            I.JJA=which(mm %in% c(6:8))
            T.q50 = tapply(T.dum[I.JJA],yy[I.JJA],median,na.rm=TRUE)
            T.JJA[[fi]] = list(filename=fi,sname=nam,modname=nmod,
                               scen=r0[[1]][4],
                               T.q50=T.q50,yy=unique(yy))
## Filter those days when the daily fluctuation is below Thresh        
            I0 = I.T0[T.dum[I.T0-1] >= T.dum[I.T0]-Tthresh]
            if(length(I0)>0){
## For each model: number of members at least one exceedance over T0            
                CMIP.stat[[modna]]$nensT0=CMIP.stat[[modna]]$nensT0+1
## For each model: number of exceedances of T0
                CMIP.stat[[modna]]$nT0=CMIP.stat[[modna]]$nT0+length(I0)
## GMST pendant les evenements  T > T0              
                if(file.exists(fiGWD.scen) & file.exists(fiGWD.hist) ){
                    nc=nc_open(fiGWD.scen)
                    time.dum = ncdf4.helpers::nc.get.time.series(nc)
                    yy.GWD=floor(as.numeric(format(time.dum, "%Y%m%d"))/10000)
                    T.GWD=ncvar_get(nc,"tas")
                    nc_close(nc)
                    nc=nc_open(fiGWD.hist)
                    time.dum = ncdf4.helpers::nc.get.time.series(nc)
                    yy.h.GWD=floor(as.numeric(format(time.dum, "%Y%m%d"))/10000)
                    T.h.GWD=ncvar_get(nc,"tas")
                    nc_close(nc)
                    GMT = c(T.h.GWD,T.GWD)
                    yy.GMT = c(yy.h.GWD,yy.GWD)
                    GMT2000=c(GMT2000,mean(GMT[yy.GMT %in% c(1950:2000)],
                                           na.rm=TRUE))
                    l1=length(T.GWD)
                    l2=length(T.h.GWD)
                    r12=range(GMT)
                    if(l1 >= 80 & l2 >=20 & r12[1] > 10 & r12[2] < 30){
## Spline smoothing (equivalent to a 30 year moving average)               
                        GMT.spl=smooth.spline(GMT,spar=0.8)
## Temperature globale pendant les evenements T > T0            
                        GWD=c()
                        for(t in I0) GWD=c(GWD,
                                           GMT.spl$y[yy.GMT %in% floor(time.ref[t]/10000)])
                        GWDdiff1=GWD-mean(GMT.spl$y[yy.GMT %in% c(1850:1900)]) ## GMST Increase since 1850-1900
                        GWDdiff2=GWD-mean(GMT.spl$y[yy.GMT %in% c(1950:2000)]) ## GMST increase since 1950-2000
## Sauvegarde du resultat                
                        T.exc[[nam]]=list(time=time.ref[I0], T=T.dum[I0],
                                          GWD=GWD,
                                          GWDdiff=GWDdiff1,GWDdiff21=GWDdiff2,
                                          varmod=varmod,scen=r0[[1]][4],
                                          run=r0[[1]][5])
                        cat(file=outfit,paste(fi,"\n"),append=TRUE)
                        cat(file=outfit,time.ref[I0],append=TRUE)
                        cat(file=outfit,"\n",append=TRUE)
                        cat(file=outfit,GWD,append=TRUE)
                        cat(file=outfit,"\n",append=TRUE)
                        cat(file=outfit,GWDdiff1,append=TRUE)
                        cat(file=outfit,"\n",append=TRUE)
                        cat(file=outfit,GWDdiff2,append=TRUE)
                        cat(file=outfit,"\n",append=TRUE)
                    }#end if
                }#end if file.exists
            }#end if length>0
        }#end if file size > 1000
    }#end if nmod in u.mod.union
}#end for fi
## outfiR: Rdata output
save(file=outfiR,T.exc,CMIP.stat,CMIP.glob.stat,GMT2000, T.JJA,ls.fi)


q("no")
