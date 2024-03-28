## Detection des jours ou tmax depasse 48°C dans CMIP6
## Figure drawing for the various diagnostics
## Par Pascal Yiou, LSCE, Mars 2023, Jan 2024
## Peut se lancer (par exemple) par
## R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/script_EN/V2/Textr-CMIP6_diags_v2_EN.R
## Requires:
## R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/script_EN/V2/Textr-CMIP6_v2_EN.R
SI=Sys.info()
user=SI[["user"]]
if(SI[[1]] == "Linux"){
    Tdir=    paste("/homedata/",user,"/CMIP6_IDF/",sep="") ## Input
    Rdir=    paste("/homedata/",user,"/CMIP6_IDF/",sep="") ## Input
    OUTdir = paste("/scratchx/",user,"/CMIP6/IDF/",sep="") ## Output
    GWDdir=  paste("/scratchx/",user,"/CMIP6/GWD/",sep="") ## Input GWD
    ERA5dir= paste("/scratchx/",user,"/ERA5/",sep="") ## Input ERA5
}

## Requires ncdf4 and ncdf4.helpers libraries
## use:
## install.packages(c("ncdf4","ncdf4.helpers"),dependencies=TRUE)
library(ncdf4)
library(ncdf4.helpers)
l.seas=list(JJA=c(6,7,8),SON=c(9,10,11),DJF=c(12,1,2),MAM=c(3,4,5))

nT=100 ## Nombre d'annees max ou on accepte qu'on depasse T0

## Reads input options. Should be the same as in the call of
##  R CMD BATCH "--args tasmax 48 20 2100" ${HOME}/programmes/RStat/GREC50/script_EN/V2/Textr-CMIP6_v2_EN.R
args=(commandArgs(TRUE))
print(args)
i=1
if(length(args)>0){
    varname = args[i];i=i+1
    T0 =as.integer(args[i]);i=i+1
    ## Securite: on veut une variation de moins de Tthresh (20°C) pour
    ## atteindre 48°C
    Tthresh=as.integer(args[i]);i=i+1
    horizon=as.integer(args[i]);i=i+1
}else{
    varname = "tasmax"
    T0 = 48
## Securite: on veut une variation de moins de Tthresh (10°C) pour atteindre 48°C
    Tthresh=20
    horizon = 2100
}

## Colors for figures
l.col=list(ssp126="green",ssp245="blue",ssp370="orange",ssp585="red")

## READ INPUT OBSERVATIONS
## Read GMST (global surface tempetature) in ERA5
## Data is obtained from the Climate Explorer
## The data is also in the data/ directory
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
## The data is also in the data/ directory
seas="JJA"
nc=nc_open(paste("/scratchx/",user,"/EOBS/iensembles_025_tx_1.75-3.25E_48.25-49.25N_n_max.nc",sep=""))
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.EOBS=as.numeric(format(time.dum, "%Y%m%d"))
dd.EOBS=time.EOBS %% 100
mm.EOBS=floor(time.EOBS/100) %% 100
yy.EOBS=floor(time.EOBS/10000)
T.EOBS=ncvar_get(nc,"tx")
I.EOBS.seas = which(mm.EOBS %in% l.seas[[seas]] & yy.EOBS >= 1995 & yy.EOBS <= 2014)
nc_close(nc)

nc=nc_open(paste("/scratchx/",user,"/EOBS/iensembles_025_tn_1.75-3.25E_48.25-49.25N_n_max.nc",sep=""))
time.dum = ncdf4.helpers::nc.get.time.series(nc)
time.EOBS=as.numeric(format(time.dum, "%Y%m%d"))
dd.TN.EOBS=time.EOBS %% 100
mm.TN.EOBS=floor(time.EOBS/100) %% 100
yy.TN.EOBS=floor(time.EOBS/10000)
TN.EOBS=ncvar_get(nc,"tn")
I.TN.EOBS.seas = which(mm.TN.EOBS %in% l.seas[[seas]] & yy.TN.EOBS >= 1995 & yy.TN.EOBS <= 2014)
nc_close(nc)

## Calcul du cycle saisonnier de TX et TN pour 1990-2020
I.JJA.ref=which(mm.EOBS %in% c(6:8) & yy.EOBS %in% c(1990:2020))
mmdd=mm.EOBS*100+dd.EOBS
TX.seas = tapply(T.EOBS,mmdd,mean,na.rm=TRUE)
TX.seas.spl = smooth.spline(TX.seas,spar=0.8)

I.JJA.ref=which(mm.TN.EOBS %in% c(6:8) & yy.TN.EOBS %in% c(1990:2020))
mmdd=mm.TN.EOBS*100+dd.TN.EOBS
TN.seas = tapply(TN.EOBS,mmdd,mean,na.rm=TRUE)
TN.seas.spl = smooth.spline(TN.seas,spar=0.8)

## READ RESULTS FILE
setwd(Tdir)
outfiR=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6_2.Rdat",
             sep="")
## de:  R CMD BATCH "--args tasmax 48" ${HOME}/programmes/RStat/GREC50/Textr-CMIP6.R          
##  T.exc[[nam]]=list(time=time.ref[I0], T=T.dum[I0],
##save(file=outfiR,T.exc)
load(outfiR)

## List of simulations
lf=names(T.exc)

## List of models
lf.dum=strsplit(lf,"_")
dum=unlist(lf.dum)
mdum=t(matrix(unlist(lf.dum),nrow=7))
umod=unique(mdum[,3])
## Symbole pour chaque modele
l.symb=list()
for(i in 1:length(umod)){
    l.symb[[umod[i]]]=i
}

## File for Table 2 of the paper
fout=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6_table_2.txt",
           sep="")
cat(file=fout,"Data for Table 1\n")
GWI=c() ## Global Warming Increase since 1900
GWI2=c() ## Global Warming Increase since 2000
for(ff in names(T.exc)){
    dum=strsplit(ff,"_")[[1]]
    ssp=dum[4]
    modna=dum[3]
    cat(file=fout,paste(modna,T.exc[[ff]]$scen,T.exc[[ff]]$run,
                        length(T.exc[[ff]]$time),
                        floor(min(T.exc[[ff]]$time)/10000),
                        format(min(T.exc[[ff]]$GWD),digits=3),
                        format(min(T.exc[[ff]]$GWDdiff),digits=2),
                        format(min(T.exc[[ff]]$GWDdiff21),digits=2),
                        "\n",sep="\t"),append=TRUE)
    GWI=c(GWI, T.exc[[ff]]$GWDdiff)
    GWI2=c(GWI2, T.exc[[ff]]$GWDdiff21)

}#end for ff

## Probability density of the years of TX>48°C
GWI2.dens=density(GWI2)
N.exc=length(T.exc)
N.ssp=585

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
## Pour les ssp > 1, on exclut UKESM1-0-LL qui atteint trop facilement T0
## Normalement, c'est inutile car on a enlevé UKESM1-0-LL
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
for(scen in c("ssp245","ssp370","ssp585")){
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
for(scen in c("ssp245","ssp370","ssp585")){
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

## Read pr (precip) files and extraction of summer
l.JJA.pr=list()
for(scen in c("ssp245","ssp370","ssp585")){
    yref=floor(as.numeric(l.f.scen[[scen]][3])/10000)
    dum=strsplit(l.f.scen[[scen]][1],paste(varname,"_",sep=""))[[1]][2]
    fname=paste("pr_",dum,".nc",sep="")
    nc=nc_open(fname)
    time.dum = ncdf4.helpers::nc.get.time.series(nc)
    time.pr=as.numeric(format(time.dum, "%Y%m%d"))
    pr.dum=(ncvar_get(nc,"pr")+273.15)*86400
    nc_close(nc)
    pr.JJA=pr.dum[time.pr >= (yref*100+6)*100+1 &
                time.pr <= (yref*100+8)*100+31]
    l.JJA.pr[[scen]]=pr.JJA

}

## Plot GMST (global mean temperature) increase (GWI) when tmax exceeds T0
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
## each model. Colors correspond to GMST increase.
fout=paste("TX",T0,"_horiz",horizon,"_1st-CMIP6.pdf",sep="")
pdf(file=fout,width=10)
layout(matrix(1:3,1,3))
par(mar=c(10,4,1,1))
i=1
for(issp in c("ssp245","ssp370","ssp585")){
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
            points(imodna,yfirst,col=colo,pch=19,cex=4)
##            points(imodna,yfirst,col=colo,pch=ssp.sym(ssp))
        }
    }
    axis(side=2)
    axis(side=1,at=c(1:length(umod)),labels=umod,las=3)
    box()
    if(issp == "ssp245"){
        
        legend("bottomleft",legend=c("14 < GMT < 15","15 < GMT < 16",
                                     "16 < GMT < 17",
                         "17 < GMT < 18",
                         "18 < GMT < 19","20 < GMT"),col=l.tcol,lwd=3,bty="n")
    }
    legend("topleft",paste("(",letters[i],")",sep=""),bty="n")
    i=i+1
}## End for issp
dev.off()


## Plot of daily time serise for selected simulations, for three
## SSP scenarios
fout="TX_CMIP6_ssp_ts.pdf"
pdf(fout)
layout(matrix(1:6,3,2))
i=1
for(scen in c("ssp245","ssp370","ssp585")){
    yref=floor(as.numeric(l.f.scen[[scen]][3])/10000)
    member=strsplit(l.f.scen[[scen]],"_")[[1]][5]
    par(mar=c(4,4,2,4))
    plot(l.JJA.T0[[scen]],type="n",
         xlab="Days since June 1st",ylab="TX [C]",
         ylim=c(15,50))
    axis(side=3,at=c(31,62),labels=c("Jul","Aug"))
    axis(side=4,at=seq(15,40,by=5),labels=seq(0,25,by=5))
    mtext("PR [mm/day]",side=4,line=2)
    lines(T.EOBS[time.EOBS >= 20030601 & time.EOBS <= 20030831],col="red",lwd=2)
    lines(T.EOBS[time.EOBS >= 20220601 & time.EOBS <= 20220831],col="orange",lwd=2)
    lines(l.JJA.T0[[scen]],lwd=3)
    lines(TX.seas.spl$y[as.numeric(names(TX.seas)) >= 601 &
                        as.numeric(names(TX.seas)) <= 831],col="green",lwd=3)
    lines(l.JJA.pr[[scen]]+15)
    abline(h=48,lty=3,col="red")
    abline(h=29,lty=4,col="brown")
    I=which(l.JJA.T0[[scen]]>=48)
    abline(v=I,lty="dashed",col="grey")
    legend("topright",legend=paste("(",letters[i],")",sep=""),bty="n")
    legend("bottom",scen,bty="n")
    legend("topleft",
           legend=c(paste("TX",l.f.scen[[scen]][2],member,"(",yref,")"),
                    "TX EOBS (2003)",
                    "TX EOBS (2022)",
                    "Seasonal cycle 1990-2020"),
           lwd=c(3,2,2,2,2),col=c("black","red","orange","green"),bty="n")
    i=i+1
}

for(scen in c("ssp245","ssp370","ssp585")){
    yref=floor(as.numeric(l.f.scen[[scen]][3])/10000)
    member=strsplit(l.f.scen[[scen]],"_")[[1]][5]
    par(mar=c(4,4,2,4))
    plot(l.JJA.T0.tasmin[[scen]],type="n",
         xlab="Days since June 1st",ylab="TN [C]",
         ylim=c(8,30))
    axis(side=3,at=c(31,62),labels=c("Jul","Aug"))
    lines(TN.EOBS[time.EOBS >= 20030601 & time.EOBS <= 20030831],col="red",lwd=2)
    lines(TN.EOBS[time.EOBS >= 20220601 & time.EOBS <= 20220831],col="orange",lwd=2)
    lines(l.JJA.T0.tasmin[[scen]],lwd=3)
    lines(TN.seas.spl$y[as.numeric(names(TN.seas)) >= 601 &
                        as.numeric(names(TN.seas)) <= 831],col="green",lwd=3)
    abline(h=20,lty=3,col="blue")
    legend("topright",legend=paste("(",letters[i],")",sep=""),bty="n")
    legend("bottom",scen,bty="n")
    legend("topleft",
           legend=c(paste("TN",l.f.scen[[scen]][2],member,"(",yref,")"),
                    "TN EOBS (2003)",
                    "TN EOBS (2022)",
                    "Seasonal cycle 1990-2020"),
           lwd=c(3,2,2,2,2),col=c("black","red","orange","green"),bty="n")
    i=i+1
}

dev.off()

## Figure des "codes barres" de GMST les annees où on dépasse 48°C
## Barcode figure of Delta GMST for years when T>48°C
scen.col=list(ssp245="blue",ssp370="green",ssp585="red")

fout="barcode_GWL_CMIP6.pdf"
pdf(fout)
par(mar=c(4,4,1,1))
plot(c(1,8),c(0,0.09),type="n",xlab="Delta GMST (1950-2000) [°C]",
     ylab="Proba",axes=FALSE)
axis(side=1)
axis(side=2)
box()
for(i in 1:length(T.exc)){
    abline(v=T.exc[[i]]$GWDdiff2,col=scen.col[[T.exc[[i]]$scen]])
}
lines(GWI2.dens$x,GWI2.dens$y*N.exc/N.ssp,lwd=3)
##legend("topleft",bty="n",legend="(b)")
legend("bottomleft",legend=names(scen.col),lwd=2,col=unlist(scen.col),
       bty="n")
dev.off()


q("no")

## Similar barcode figure, but with Delta GMST from 1850-1900
fout="barcode_GWL_CMIP6.pdf"
pdf(fout)
layout(matrix(1:2,2,1))
par(mar=c(4,1,1,1))
plot(c(1,8),c(0,1),type="n",xlab="Delta GMST (1850-1900) [C]",
     ylab="",axes=FALSE)
axis(side=1)
box()
for(i in 1:length(T.exc)){
    abline(v=T.exc[[i]]$GWDdiff,col=scen.col[[T.exc[[i]]$scen]])
}
legend("topleft",bty="n",legend="(a)")
legend("bottomleft",legend=names(scen.col),lwd=2,col=unlist(scen.col),
       bty="n")
plot(c(1,8),c(0,1),type="n",xlab="Delta GMST (1950-2000) [C]",
     ylab="",axes=FALSE)
axis(side=1)
box()
for(i in 1:length(T.exc)){
    abline(v=T.exc[[i]]$GWDdiff2,col=scen.col[[T.exc[[i]]$scen]])
}
legend("topleft",bty="n",legend="(b)")
legend("bottomleft",legend=names(scen.col),lwd=2,col=unlist(scen.col),
       bty="n")
dev.off()
