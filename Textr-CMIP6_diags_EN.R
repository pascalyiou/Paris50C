## Simple analyses of days when tmax exceeds 48°C in the CMIP6 ensemble
## Par Pascal Yiou, LSCE, Mars 2023, Juillet 2023, August 2023
## This code is distributed "as is" under a CeCILL license:
## http://www.cecill.info/
## It can be downloaded and used for free for academic purposes.
## For commercial uses, please contact Pascal Yiou (pascal.yiou@lsce.ipsl.fr)

## Example use:
## R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/scripts_EN/Textr-CMIP6_diags_EN.R
## Requires:
## R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/script_EN/Textr-CMIP6_EN.R
## and read the file:
## outfiR=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6.Rdat",sep=""

## General path parameters for input and output
SI=Sys.info()
user=SI[["user"]]
if(SI[[1]] == "Linux"){
    ERA5dir= paste("/scratchx/",user,"/ERA5/",sep="") ## Input
    Tdir=    paste("/scratchx/",user,"/CMIP6/IDF/",sep="") ## Input
    OUTdir = paste("/scratchx/",user,"/CMIP6/IDF/",sep="") ## Output
    GWDdir=  paste("/scratchx/",user,"/CMIP6/GWD/",sep="") ## Input GWD
}

## Requires ncdf4 and ncdf4.helpers libraries
## use:
## install.packages(c("ncdf4","ncdf4.helpers"),dependencies=TRUE)
library(ncdf4)
library(ncdf4.helpers)

l.seas=list(JJA=c(6,7,8),SON=c(9,10,11),DJF=c(12,1,2),MAM=c(3,4,5))

nT=100 ## Maximum number of events for each simulation
## Reads input options. Should be the same as in the call of
##  R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/script_EN/Textr-CMIP6_EN.R
args=(commandArgs(TRUE))
print(args)
i=1
if(length(args)>0){
    varname = args[i];i=i+1
    T0 =as.integer(args[i]);i=i+1
    Tthresh=as.integer(args[i]);i=i+1
    horizon=as.integer(args[i]);i=i+1
}else{
    varname = "tasmax"
    T0 = 48
    Tthresh=20
    horizon = 2100
}

## Colors for figures
l.col=list(ssp126="green",ssp245="blue",ssp370="orange",ssp585="red")

## Read GST (global surface tempetature) in ERA5
## Data is obtained from the Climate Explorer
fname=paste(ERA5dir,"era5_t2m_year_GWD.nc",sep="")
nc=nc_open(fname)
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.ERA=as.numeric(format(time.dum, "%Y%m%d"))
yy.ERA=floor(time.ERA/10000)
T.ERA=ncvar_get(nc,"t2m")
nc_close(nc)
T.ERA.late = mean(T.ERA[yy.ERA >= 1980])

## Read tmax over IdF from EOBS
## Data is obtained from the Climate Explorer
seas="JJA"
nc=nc_open(paste("/scratchx/",user,"/EOBS/iensembles_025_tx_1.75-3.25E_48.25-49.25N_n_max.nc",sep=""))
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.EOBS=as.numeric(format(time.dum, "%Y%m%d"))
mm.EOBS=floor(time.EOBS/100) %% 100
yy.EOBS=floor(time.EOBS/10000)
T.EOBS=ncvar_get(nc,"tx")
I.EOBS.seas = which(mm.EOBS %in% l.seas[[seas]] & yy.EOBS >= 1995 & yy.EOBS <= 2014)
nc_close(nc)

setwd(Tdir)
outfiR=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6.Rdat",
             sep="")
## from:
## R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/scripts_EN/Textr-CMIP6.R          
## T.exc[[nam]]=list(time=time.ref[I0], T=T.dum[I0],
##                   GWD=GWD, GWDdiff=GWDdiff,
##                   varmod=varmod,scen=r0[[1]][4],run=r0[[1]][5])
##save(file=outfiR,T.exc)
load(outfiR)

## List of simulations
lf=names(T.exc)

## List of models
lf.dum=strsplit(lf,"_")
dum=unlist(lf.dum)
mdum=t(matrix(unlist(lf.dum),nrow=6))
umod=unique(mdum[,3])
## Symbol for each model
l.symb=list()
for(i in 1:length(umod)){
    l.symb[[umod[i]]]=i
}

## File summarize results (Table 1 of paper)
fout=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6_table.txt",
           sep="")
cat(file=fout,"Data for Table 1\n")
for(ff in names(T.exc)){
    dum=strsplit(ff,"_")[[1]]
    ssp=dum[4]
    modna=dum[3]
    cat(file=fout,paste(modna,T.exc[[ff]]$scen,T.exc[[ff]]$run,
                        length(T.exc[[ff]]$time),
                        floor(min(T.exc[[ff]]$time)/10000),
                        format(min(T.exc[[ff]]$GWD),digits=3),
                        format(min(T.exc[[ff]]$GWDdiff),digits=2),
                        "\n",sep="\t"),append=TRUE)

}#end for ff


## First occurrence of T>T0 for each model
l.first=list()
for(ff in names(T.exc)){
    dum=strsplit(ff,"_")[[1]]
    ssp=dum[4]
    modna=dum[3]
    l.first[[ff]]=list(tfirst=T.exc[[ff]]$time[1],
                          GWDfirst=T.exc[[ff]]$GWD[1])
}

## First occurrence of T>T0 in each scenario
l.first.scen=list()
for(scen in c("ssp126","ssp245","ssp370","ssp585")){
    fdiag=c()
    for(ff in names(T.exc)){
        dum=strsplit(ff,"_")[[1]]
        ssp=dum[4]
        modna=dum[3]
        if(ssp == scen){
            lf=T.exc[[ff]]$time[1]
## For ssp2, 3 and 5, we eclude UKESM1-0-LL which reaches T0 too often
            if(scen=="ssp126" | (scen != "ssp126" & modna != "UKESM1-0-LL")){
                fdiag=rbind(fdiag,c(ff,modna,lf))
            }
        }
    }
    l.first.scen[[scen]]=fdiag   
}

## Result: for each scenario, the first model that reaches T0
l.f.scen=list()
for(scen in c("ssp126","ssp245","ssp370","ssp585")){
    i.f=which.min(as.numeric(l.first.scen[[scen]][,3]))
    l.f.scen[[scen]]=l.first.scen[[scen]][i.f,]
}

## Read correcponding time series, and extract summer
l.JJA.T0=list()
for(scen in c("ssp126","ssp245","ssp370","ssp585")){
    yref=floor(as.numeric(l.f.scen[[scen]][3])/10000)
    fname=paste(l.f.scen[[scen]][1],".nc",sep="")
    nc=nc_open(fname)
    time.dum = ncdf4.helpers::nc.get.time.series(nc)
    time.ref=as.numeric(format(time.dum, "%Y%m%d"))
    T.dum=ncvar_get(nc,varname)
    nc_close(nc)
    T.JJA=T.dum[time.ref >= (yref*100+6)*100+1 &
                time.ref <= (yref*100+8)*100+31]
    l.JJA.T0[[scen]]=T.JJA

}

## Read time series of tasmin, and extract summer
l.JJA.T0.tasmin=list()
for(scen in c("ssp126","ssp245","ssp370","ssp585")){
    yref=floor(as.numeric(l.f.scen[[scen]][3])/10000)
    dum=strsplit(l.f.scen[[scen]][1],paste(varname,"_",sep=""))[[1]][2]
    fname=paste("tasmin_",dum,".nc",sep="")
    nc=nc_open(fname)
    time.dum = ncdf4.helpers::nc.get.time.series(nc)
    time.ref=as.numeric(format(time.dum, "%Y%m%d"))
    T.dum=ncvar_get(nc,"tasmin")
    nc_close(nc)
    T.JJA=T.dum[time.ref >= (yref*100+6)*100+1 &
                time.ref <= (yref*100+8)*100+31]
    l.JJA.T0.tasmin[[scen]]=T.JJA

}

## Plot GST (global mean temperature) increase (GWD) when tmax exceeds T0
## for each model and each scenario
pdf(paste("T-GWD_T0",T0,"_horiz",horizon,"_CMIP6.pdf",sep=""),width=8)
layout(matrix(1:2,1,2),width=c(2,1))
par(mar=c(4,4,1,1))
plot(c(2015,2100),c(14,22),type="n",xlab="Years",ylab="T GWD [°C]")

for(ff in names(T.exc)){
    dum=strsplit(ff,"_")[[1]]
    ssp=dum[4]
    yy=unique(floor(T.exc[[ff]]$time/10000))
    yy=floor(T.exc[[ff]]$time/10000)
    if(yy[1]<=2100 & length(yy[yy<=2100])<=nT){
        GWD=T.exc[[ff]]$GWD
        GWDdiff=T.exc[[ff]]$GWDdiff
        points(yy[yy<=2100],GWD[yy<=2100],col=l.col[[ssp]],
               pch=unlist(l.symb[dum[3]]))
    }
}
abline(h=T.ERA.late,col="grey",lty=3,lwd=4)
par(mar=c(1,1,1,1))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
legend("top",legend=umod,pch=unlist(l.symb),ncol=1,cex=1,bty="n")
legend("bottom",legend=names(l.col),col=unlist(l.col),lwd=3,bty="n")
dev.off()

l.tcol=c("blue","green","grey","orange","red","brown")
"tcol"=function(T,l.col=l.tcol){
    iT=floor(T)-14
    mmiT = min(max(iT,1),length(l.col))
    return(l.col[mmiT])
}

"ssp.off"=function(ssp){
    if(ssp=="historical") return(-0.2)
    if(ssp=="ssp126") return(-0.1)
    if(ssp=="ssp245") return(0.0)
    if(ssp=="ssp370") return(0.1)
    if(ssp=="ssp585") return(0.2)
}

"ssp.sym"=function(ssp){
    if(ssp=="historical") return(21)
    if(ssp=="ssp126") return(22)
    if(ssp=="ssp245") return(23)
    if(ssp=="ssp370") return(24)
    if(ssp=="ssp585") return(25)
    }

## Plot year of first occurrence of T>T0 for each scenario and
## each model. Colors correspond to GST increase.
fout=paste("TX",T0,"_horiz",horizon,"_1st-CMIP6.pdf",sep="")
pdf(file=fout,height=10)
layout(matrix(1:4,2,2))
par(mar=c(10,4,1,1))
i=1
for(issp in c("ssp126","ssp245","ssp370","ssp585")){
    plot(c(1,length(umod)),c(2015,2100),type="n",xlab="",
         ylab=paste("1st occurrence in",issp),axes=FALSE)
    for(ff in names(T.exc)){
        dum=strsplit(ff,"_")[[1]]
        ssp=dum[4]
        modna=dum[3]
        imodna=match(modna,umod)
        yfirst=floor(l.first[[ff]]$tfirst/10000)
        colo=tcol(l.first[[ff]]$GWD)
        sspoff=ssp.off(ssp)
        if(ssp == issp){
            points(imodna,yfirst,col=colo,pch=19)
        }
    }
    axis(side=2)
    axis(side=1,at=c(1:length(umod)),labels=umod,las=3)
    box()
    if(issp == "ssp126"){
        
        legend("bottomleft",legend=c("14 < T < 15","15 < T < 16","16 < T < 17",
                         "17 < T < 18",
                         "18 < T < 19","20 < T"),col=l.tcol,lwd=3,bty="n")
    }
    legend("topleft",paste("(",letters[i],")",sep=""),bty="n")
    i=i+1
}## End for issp
dev.off()

## Plot of daily time serise for selected simulations, for four
## SSP scenarios
fout="TX_CMIP6_ssp_ts.pdf"
pdf(fout)
layout(matrix(1:4,2,2))
i=1
for(scen in c("ssp126","ssp245","ssp370","ssp585")){
    yref=floor(as.numeric(l.f.scen[[scen]][3])/10000)
    member=strsplit(l.f.scen[[scen]],"_")[[1]][5]
    par(mar=c(4,4,1,1))
    plot(l.JJA.T0[[scen]],type="n",
         xlab="Days since June 1st",ylab="TX [°C]",
         ylim=c(10,50))
    lines(T.EOBS[time.EOBS >= 20030601 & time.EOBS <= 20030831],col="red",lwd=2)
    lines(T.EOBS[time.EOBS >= 20220601 & time.EOBS <= 20220831],col="orange",lwd=2)
    lines(l.JJA.T0[[scen]],lwd=3)
    lines(l.JJA.T0.tasmin[[scen]],lwd=2,col="blue")
    abline(h=20,lty=3,col="blue")
    legend("bottomright",legend=paste("(",letters[i],")",sep=""),bty="n")
    legend("bottom",scen,bty="n")
    legend("topleft",
           legend=c(paste("TX",l.f.scen[[scen]][2],member,"(",yref,")"),
                    paste("TM",l.f.scen[[scen]][2],member,"(",yref,")"),
                    "TX EOBS (2003)",
                    "TX EOBS (2022)"),
           lwd=c(3,2,2,2),col=c("black","blue","red","orange"),bty="n")
    i=i+1
}
dev.off()

q("no")
