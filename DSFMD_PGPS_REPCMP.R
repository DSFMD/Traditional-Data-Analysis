# DSFMD - Dynamic Systems Framework for Movement Data
# Example implementation for calculation for comparing repeated flights
# Code by Ingo Schiffner 2018

#set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(sp) #spatial transformation library
library(data.table) # better than data frame
library(useful) # cart2pol
library(ARTool) # aligned rank transform anova
library(car) # for leveneTest
library(emmeans) # lsmeans post hoc comparison
library(plotly) # interactive plots
library(compiler) # byte compiler

#sub/helper functions
source('../SHARED/GetPGPSData.R') # read GPS data from csv file
source("../SHARED/CalSimil.R") # comapre and estimate similarity indices
source("../SHARED/CalEff.R") # calculate efficiency
source("../SHARED/CalNN.R") # calculate nearest neighbour distance
source("../SHARED/CalAngPErr.R") # calculate angular prediction error
source("../SHARED/RankHist.R") # histogram based ranking
source("../SHARED/CAlNMI.R") # calculate normalied mutual information 

#compile sub/helper functions
GetPGPSData <- cmpfun(GetPGPSData)
CalSimil <- cmpfun(CalSimil)
CalEff <- cmpfun(CalEff)
CalNN <- cmpfun(CalNN)
CalAngPErr <- cmpfun(CalAngPErr)
RankHist <- cmpfun(RankHist) 
CalNMI <- cmpfun(CalNMI)

#projection string for conversion to wgs84 latitude and longitude format
wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#projection string for Universal Transverse Mercator (UTM) format
utm<-"+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#select data folder
cdir<-choose.dir(default = "", caption = "Select folder")
SubDirNames <- list.files(path = cdir, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)

#process timing
ptm <- proc.time()

#load metadata
mdir<-dirname(cdir)
loc_dt <- read.table(paste(mdir,"/meta_data/sites.txt",sep=''), header=T, sep=",", quote='"', dec=",", na.strings="", colClasses="character", fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)
loc_dt$lat <- as.numeric(loc_dt$lat)
loc_dt$lon <- as.numeric(loc_dt$lon)
coordinates(loc_dt)<-~lon+lat
proj4string(loc_dt)<-CRS(wgs84)
loc_dt<-spTransform(loc_dt, utm)

#get data from all subfolders
gps_dt <-NULL
for (csdir in SubDirNames)
{
  print(paste("reading data from", csdir))
  csdir_full <- paste(cdir,'\\',csdir,sep="")
  FileNames <- list.files(path = csdir_full, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  
  #site index
  sdx <- which(loc_dt$Site == csdir)
  
  #process all files in subfolder
  for (i in FileNames)
  {
    t_dat <- GetPGPSData(paste(csdir_full,'\\',i ,sep=""))
    
    #limit to flight data
    fdx <- which(t_dat$SPEED>5)[1]
    if (!is.na(fdx))
    {
      dx <- t_dat$X-loc_dt@coords[1,1]
      dy <- t_dat$Y-loc_dt@coords[1,2]
      de <- (dx^2 + dy^2)^0.5
      hdx<-which (de<50)[1]
      if (is.na(hdx))
      {
        hdx<-nrow(t_dat)
      }
      t_dat<-t_dat[fdx:hdx,]
      
      spl<-unlist(strsplit(i,'_'))
      t_dat$DATE <- t_dat$TIME[1]
      tmp <-  unlist(strsplit(spl[2],'.csv'))
      if (nchar(tmp)==1)
      {
        tmp<-unlist(strsplit(spl[3],'.csv'))
      }
      t_dat$ID <- tmp
      t_dat$SITE <- loc_dt$Site[sdx]
      gps_dt <- rbind(gps_dt,t_dat)
    }
  }
}
gps_dt<-as.data.table(gps_dt)
rm(t_dat,csdir, csdir_full, FileNames, loc_dt, sdx, SubDirNames, spl,utm,wgs84)
rm(fdx,hdx,dx,dy,de)

#determine main release site
tmp <- gps_dt[,6:8]
tmp <- as.data.table(lapply(tmp, as.factor))
tmp$CNT<-1
tmp<-aggregate(CNT~SITE+DATE,tmp,max)
scnt<-aggregate(CNT~SITE,tmp,sum)
mdx <- which(scnt$CNT==max(scnt$CNT))
sv <- unique(gps_dt$SITE)

mres<-NULL
ores<-NULL
tmp <- subset(gps_dt,SITE==scnt$SITE[mdx])
bv <-unique(tmp$ID)
for (i in bv)
{
  stmp <- subset(tmp,ID==i)
  dv <- unique(stmp$DATE)
  #make sure data is ordered
  dv<-dv[order(dv)]
  for (j in 1:length(dv))
  {
    print(paste("processing data from", i, "track No.", j))
    #get data for current track
    cxy <- subset(stmp,DATE==dv[j])
    cxy <- cxy[order(cxy$TIME),]
    cxy <- cxy[,c(2,3)]
    
    if (j!=length(dv)){
      
      #get data for random track 
      rbx <- sample(which(bv!=i))[1]
      rxy <- subset(tmp,ID==bv[rbx])
      rdx <- sample(unique(rxy$DATE))[1]
      rxy <- subset(rxy,DATE==rdx)
      rxy <- rxy[order(rxy$TIME),]
      rxy <- rxy[,c(2,3)]
      
      #get data for next track
      nxy <- subset(stmp,DATE==dv[j+1])
      nxy <- nxy[order(nxy$TIME),]
      nxy <- nxy[,c(2,3)]
      
      nres<-CalSimil(cxy,nxy)
      rres<-CalSimil(cxy,rxy)
      mres <- rbind(mres,c(i,j,nres,rres[2:6]))
      
      
    }else{
      #compare to offsite releases
      ov <- which(sv!=sv[mdx])
      for (h in ov)
      {
        otmp <- subset(gps_dt,SITE==sv[h])
        
        #get data for offsite track
        oxy <- subset(otmp, ID==i)
        oxy <- oxy[order(oxy$TIME),]
        oxy <- oxy[,c(2,3)]
        if (nrow(oxy)!=0)
        {
          #get data for random offsite track
          oid <- unique(otmp$ID)
          rxy <- subset(otmp, ID==oid[sample(which(oid!=i))[1]])
          rxy <- rxy[order(rxy$TIME),]
          rxy <- rxy[,c(2,3)]
          
          nres <-CalSimil(cxy,oxy)
          rres <-CalSimil(cxy,rxy)
          ores <- rbind(ores,c(i,sv[h],nres,rres[2:6]))
        }
      }
    }
  }
}
#clear workspace
rm(cxy,nxy,oxy,rxy,nres,rres,scnt,stmp,otmp,tmp,gps_dt)
rm(i,j,h,oid,mdx,rbx,rdx,bv,dv,sv,ov)

mres <- as.data.table(mres)
colnames(mres)<-c('ID','FNO','CE','NE','dE','NND','APE','NMI','RE','RdE','RNND','RAPE','RNMI')
mres[,3:13] = lapply(mres[,3:13], as.numeric)

#stats
tstat <- function (x,y)
{
  print(mean(x,na.rm=T))
  print(mean(y,na.rm=T))
  t.test(x,y,paired=T)
}
tstat(mres$dE, mres$RdE)
tstat(mres$NND, mres$RNND)
tstat(mres$NMI, mres$RNMI)
tstat(mres$APE, mres$RAPE)

ores <- as.data.table(ores)
colnames(ores)<-c('ID','FNO','CE','NE','dE','NND','APE','NMI','RE','RdE','RNND','RAPE','RNMI')
ores[,3:13] = lapply(ores[,3:13], as.numeric)

tstat(ores$dE, ores$RdE)
tstat(ores$NND, ores$RNND)
tstat(ores$NMI, ores$RNMI)
tstat(ores$APE, ores$RAPE)

#make summary graphs  
sem <- function(x) {
  sd(x)/sqrt(length(x))
}

ldx<-unique(ores$FNO)[1]
lv<-which(ores$FNO==ldx)
mres_add <-ores[lv,]
mres_add$FNO <- as.character(max(as.numeric(mres$FNO))+1)
mres_add$FNO <- as.numeric(mres_add$FNO)
#efficiency
EFF_avg<-aggregate(data=rbind(mres,mres_add), CE~FNO, FUN=mean)
EFF_sem<-aggregate(data=rbind(mres,mres_add), CE~FNO, FUN=sem)
plot_ly(EFF_avg, x=~FNO, y=~CE, type='scatter', mode='lines',line=list(color='rgba(0,100,80,1)'),name='Average') %>%
  add_trace(x=~FNO, y=~CE-EFF_sem$CE,type='scatter',mode='lines',line=list(color='rgba(0,100,80,1)'),name='Lower') %>%
  add_trace(x=~FNO, y=~CE+EFF_sem$CE,type='scatter',mode='lines',fill = 'tonexty', fillcolor='rgba(0,100,80,0.15)', line=list(color='rgba(0,100,80,1)'),name='Upper')%>%
  layout(xaxis = list(title="Flight No."),
         yaxis = list(title="Efficiency",range=c(0,1)))

m <- art(CE ~ as.factor(FNO) + (1|ID), data=rbind(mres,mres_add))
leveneTest(CE ~ as.factor(FNO), data=rbind(mres,mres_add))
anova(m)
marginal<-emmeans(artlm(m, "as.factor(FNO)"),~ as.factor(FNO))
pairs(marginal,adjust="tukey")
Sum = cld(marginal,alpha = 0.05,adjust = "tukey")         ###  Tukey-adjusted comparisons
Sum
cor.test(as.numeric(EFF_avg$FNO),EFF_avg$CE)
rm(ldx,lv,mres_add,EFF_avg,EFF_sem)

NN_avg<-aggregate(data=mres, NND~FNO, FUN=mean)
NN_sem<-aggregate(data=mres, NND~FNO, FUN=sem)
plot_ly(NN_avg, x=~FNO, y=~NND, type='scatter', mode='lines',line=list(color='rgba(0,100,80,1)'),name='Average') %>%
  add_trace(x=~FNO, y=~NND-NN_sem$NND,type='scatter',mode='lines',line=list(color='rgba(0,100,80,1)'),name='Lower') %>%
  add_trace(x=~FNO, y=~NND+NN_sem$NND,type='scatter',mode='lines',fill = 'tonexty', fillcolor='rgba(0,100,80,0.15)', line=list(color='rgba(0,100,80,1)'),name='Upper')%>%
  layout(xaxis = list(title="Flight No."),
         yaxis = list(title="Nearest Neighbour (m)",range=c(0,600)))

m <- art(NND ~ as.factor(FNO) + (1|ID), data=mres)
leveneTest(NND ~ as.factor(FNO), data=mres)
anova(m)
marginal<-emmeans(artlm(m, "as.factor(FNO)"),~ as.factor(FNO))
pairs(marginal,adjust="tukey")
Sum = cld(marginal,alpha = 0.05,adjust = "tukey")         ###  Tukey-adjusted comparisons
Sum
cor.test(as.numeric(NN_avg$FNO),NN_avg$NND)
rm(NN_avg,NN_sem)

NMI_avg<-aggregate(data=mres, NMI~FNO, FUN=mean)
NMI_sem<-aggregate(data=mres, NMI~FNO, FUN=sem)
plot_ly(NMI_avg, x=~FNO, y=~NMI, type='scatter', mode='lines',line=list(color='rgba(0,100,80,1)'),name='Average') %>%
  add_trace(x=~FNO, y=~NMI-NMI_sem$NMI,type='scatter',mode='lines',line=list(color='rgba(0,100,80,1)'),name='Lower') %>%
  add_trace(x=~FNO, y=~NMI+NMI_sem$NMI,type='scatter',mode='lines',fill = 'tonexty', fillcolor='rgba(0,100,80,0.15)', line=list(color='rgba(0,100,80,1)'),name='Upper')%>%
  layout(xaxis = list(title="Flight No."),
         yaxis = list(title="Normalized Mutual Information",range=c(0,1)))

m <- art(NMI ~ as.factor(FNO) + (1|ID), data=na.omit(mres))
leveneTest(NMI ~ as.factor(FNO), data=mres)
anova(m)
marginal<-emmeans(artlm(m, "as.factor(FNO)"),~ as.factor(FNO))
pairs(marginal,adjust="tukey")
Sum = cld(marginal,alpha = 0.05,adjust = "tukey")         ###  Tukey-adjusted comparisons
Sum
cor.test(as.numeric(NMI_avg$FNO),NMI_avg$NMI)
rm(NMI_avg,NMI_sem)

APE_avg<-aggregate(data=mres, APE~FNO, FUN=mean)
APE_sem<-aggregate(data=mres, APE~FNO, FUN=sem)
plot_ly(APE_avg, x=~FNO, y=~APE, type='scatter', mode='lines',line=list(color='rgba(0,100,80,1)'),name='Average') %>%
  add_trace(x=~FNO, y=~APE-APE_sem$APE,type='scatter',mode='lines',line=list(color='rgba(0,100,80,1)'),name='Lower') %>%
  add_trace(x=~FNO, y=~APE+APE_sem$APE,type='scatter',mode='lines',fill = 'tonexty', fillcolor='rgba(0,100,80,0.15)', line=list(color='rgba(0,100,80,1)'),name='Upper')%>%
  layout(xaxis = list(title="Flight No."),
         yaxis = list(title="Angular Prediction Error",range=c(0,90)))

m <- art(APE ~ as.factor(FNO) + (1|ID), data=mres)
leveneTest(APE ~ as.factor(FNO), data=mres)
anova(m)
marginal<-emmeans(artlm(m, "as.factor(FNO)"),~ as.factor(FNO))
pairs(marginal,adjust="tukey")
Sum = cld(marginal,alpha = 0.05,adjust = "tukey")         ###  Tukey-adjusted comparisons
Sum
cor.test(as.numeric(APE_avg$FNO),APE_avg$APE)
rm(APE_avg,APE_sem)
rm(m,marginal,Sum)

#write rep results
NND_avg<-aggregate(data=mres, NND~ID, FUN=mean)
NMI_avg<-aggregate(data=mres, NMI~ID, FUN=mean)
APE_avg<-aggregate(data=mres, APE~ID, FUN=mean)
REP_M_RES <- merge(NND_avg,NMI_avg)
REP_M_RES <- merge(REP_M_RES,APE_avg)
fn <- paste(dirname(cdir),"/res/REP_M_RES.csv",sep='')
write.table(REP_M_RES, file = fn, append = T, quote = F, sep = ";",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"))
print(ptm-proc.time())
rm(list=ls())
