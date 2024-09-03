## ENSO and AMO states during tmax events in CMIP6
## Requires ${HOME}/programmes/RStat/GREC50/ENSO_AMO_CMIP6_EN.sh to compute
## those indices
## Par Pascal Yiou, LSCE, Sept. 2024
SI=Sys.info()
user=SI[["user"]]
if(SI[[1]] == "Linux"){
    Tdir=paste("/homedata/",user,"/CMIP6_IDF/",sep="") ## Needs to be adapted
    ENSOdir=paste("/scratchx/",user,"/CMIP6/ENSO-AMO/",sep="") ## Needs to be adapted
   OUTdir = paste("/scratchx/",user,"/CMIP6/ENSO-AMO/",sep="") ## Needs to be adapted
}

library(ncdf4)
library(ncdf4.helpers)
l.seas=list(JJA=c(6,7,8),SON=c(9,10,11),DJF=c(12,1,2),MAM=c(3,4,5))

varname="tasmax"
T0=48
nT=100 ## Nombre d'annees max ou on accepte qu'on depasse T0

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
## Securite: on veut une variation de moins de Tresh (10°C) pour atteindre 48°C
    Tthresh=20
    horizon = 2100
}

l.col=list(ssp126="green",ssp245="blue",ssp370="orange",ssp585="red")

setwd(Tdir)
outfiR=paste(varname,"_T",T0,"_dT",Tthresh,"_horiz",horizon,"_CMIP6_2.Rdat",
             sep="")
## de:  R CMD BATCH "--args tasmax 48" ${HOME}/programmes/RStat/GREC50/Textr-CMIP6.R          
##  T.exc[[nam]]=list(time=time.ref[I0], T=T.dum[I0],
##save(file=outfiR,T.exc)
load(outfiR)

## Liste des simulations
lf=names(T.exc)

## Liste des modeles avec depassement de tX>50C
lf.dum=strsplit(lf,"_")
dum=unlist(lf.dum)
mdum=t(matrix(unlist(lf.dum),nrow=7))
umod=unique(mdum[,3])
## Symbole pour chaque modele
l.symb=list()
for(i in 1:length(umod)){
    l.symb[[umod[i]]]=i
}

## Liste des fichiers ENSO et AMO pour lesquels il y a un dépassement
## de seuil avant 2100
setwd(ENSOdir)
lf.AMO=list()
lf.ENSO=list()
for(mod in umod){
    f.AMO = system(paste("ls tos_*_",mod,"_*_AMO.nc",sep=""),intern=TRUE)
    f.ENSO = system(paste("ls tos_*_",mod,"_*_Nino3.4-index.nc",sep=""),
                    intern=TRUE)
    lf.AMO[[mod]]=f.AMO
    lf.ENSO[[mod]]=f.ENSO
}

## Read AMO and ENSO data from CMIP6 simulations
l.AMO=list()
l.ENSO=list()
for(mod in umod){
    ## Lecture de l'AMO
    for(isim in 1:length(lf.AMO[[mod]])){
        nasim=lf.AMO[[mod]][isim]
        pref.nasim=strsplit(nasim,".nc")[[1]]
        nc=nc_open(nasim)
        time.dum = ncdf4.helpers::nc.get.time.series(nc)
        time.AMO=as.numeric(format(time.dum, "%Y%m%d"))
        mm=floor(time.AMO/100) %% 100
        yy=floor(time.AMO/10000)
        AMO=ncvar_get(nc,"tos")
        nc_close(nc)
        AMOJJA=tapply(AMO[mm %in% c(6,7,8)],yy[mm %in% c(6,7,8)],mean)
        l.AMO[[pref.nasim]]=cbind(unique(yy),AMOJJA)
    }
    ## Lecture de l'I NIno3.4
     for(isim in 1:length(lf.ENSO[[mod]])){
        nasim=lf.ENSO[[mod]][isim]
        pref.nasim=strsplit(nasim,".nc")[[1]]
        nc=nc_open(nasim)
        time.dum = ncdf4.helpers::nc.get.time.series(nc)
        time.ENSO=as.numeric(format(time.dum, "%Y%m%d"))
        mm=floor(time.ENSO/100) %% 100
        yy=floor(time.ENSO/10000)
        ENSO=ncvar_get(nc,"tos")
        nc_close(nc)
        ENSOJJA=tapply(ENSO[mm %in% c(6,7,8)],yy[mm %in% c(6,7,8)],mean)
        l.ENSO[[pref.nasim]]=cbind(unique(yy),ENSOJJA)
    } 
}

## Scatter plot of Nino3.4 vs AMO indices
fout="ENSO-AMO_CMIP6.pdf"
pdf(fout)
par(mar=c(4,4,1,1))
plot(c(0,5),c(-3,3),type="n",xlab="AMO index (JJA)",ylab="Nino3.4 index (JJA)")
for(i in 1:length(l.AMO)){
    nn=strsplit(names(l.AMO)[i],"tos_")[[1]][2]
    nn=strsplit(nn,"_AMO")[[1]]
    points(l.AMO[[i]][,2],l.ENSO[[i]][,2],cex=0.1)
}
for(i in 1:length(l.AMO)){
    nn=strsplit(names(l.AMO)[i],"tos_")[[1]][2]
    nn=strsplit(nn,"_AMO")[[1]]
    nn.conv=paste("tasmax_",nn,"_max_IdF",sep="")
    if(nn.conv %in% names(T.exc)){
        TT=T.exc[[nn.conv]]$time
        ii=which(l.AMO[[i]][,1] %in% floor(TT/10000))
        points(l.AMO[[i]][ii,2],l.ENSO[[i]][ii,2],pch=19,col="red",cex=2)
    }
}
dev.off()

q("no")
