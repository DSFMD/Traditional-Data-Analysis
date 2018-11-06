# DSFMD - Dynamic Systems Framework for Movement Data
# Example implementation for calculation of traditional measurements for GPS data
# Code by Ingo Schiffner 2018

#set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(sp) #spatial transformation library
library(data.table) # better than data frame
library(compiler) # byte compiler

#sub/helper functions
source('../SHARED/GetPGPSData.R') # read GPS data from csv file
source('../SHARED/GetVVB.R') # calculate vanishing bearings at specified intervals
source('../SHARED/CalAng.R') # convert cartesian to polar coordinates
source('../SHARED/CalEff.R') # calculate efficiency
source('../SHARED/ConTime.R') # convert time (s) to HH:MM:SS format

#compile sub/helper functions for better performance
GetPGPSData <- cmpfun(GetPGPSData)
GetVVB <- cmpfun(GetVVB)
CalAng <- cmpfun(CalAng)
CalEff <- cmpfun(CalEff)
ConTime <- cmpfun(ConTime)

vv <- c(500,1000,1500,2000,2500) # take vanishing bearings at following intervals

#projection string for conversion to wgs84 latitude and longitude format
wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#projection string for Universal Transverse Mercator (UTM) format
utm<-"+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#select data folder
cdir<-choose.dir(default = "", caption = "Select folder")
SubDirNames <- list.files(path = cdir, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)

#process timing
ptm <- proc.time()

#load meta data
mdir<-dirname(cdir)
loc_dt <- read.table(paste(mdir,"/meta_data/sites.txt",sep=''), header=T, sep=",", quote='"', dec=",", na.strings="", colClasses="character", fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)
loc_dt$lat <- as.numeric(loc_dt$lat)
loc_dt$lon <- as.numeric(loc_dt$lon)
coordinates(loc_dt)<-~lon+lat
proj4string(loc_dt)<-CRS(wgs84)
loc_dt<-spTransform(loc_dt, utm)

#process subfolders
dir.create(paste(mdir,"/res/",sep=''))
for (csdir in SubDirNames)
{
  print(paste("processing data from", csdir))
  csdir_full <- paste(cdir,'\\',csdir,sep="")
  isdir <- !file_test("-f", csdir_full) # returns TRUE
  FileNames <- list.files(path = csdir_full, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  
  sdx <- which(loc_dt$Site == csdir)
  
  #process all files in subfolder
  for (cfile in FileNames)
  {
    
    #get info from file name
    IDV <-unlist(strsplit(cfile,'.csv'))
    IDV <-unlist(strsplit(IDV,'_'))
    if (length(IDV)<3) {IDV <- c(IDV,'C')}
    
    #read data
    fullfile <- paste(csdir_full,'\\',cfile ,sep="")
    t_dat <- GetPGPSData(fullfile)
    
    #limit to flight data
    rp <- which(t_dat$SPEED>5)[1]
    if (!is.na(rp))
    {
      
      #calculation of distance and direction to home
      dx <- t_dat$X[2:nrow(t_dat)] - loc_dt@coords[1,1]
      dy <- t_dat$Y[2:nrow(t_dat)] - loc_dt@coords[1,2]
      dst_h <- (dx^2+dy^2)^0.5 
      
      #first point at loft
      hp <- min(which(dst_h <= 50))
      if (is.infinite(hp)){
        hp <- nrow(t_dat)
      }
      t_dat<-t_dat[rp:hp,]
      
      #calculate virtual vanishing bearings
      vvb_res<-GetVVB(t_dat,sdx,loc_dt,vv)
      
      #calculate efficiency, duration, N pauses, pause duration.
      ce <- CalEff(t_dat$X,t_dat$Y)[[3]]
      
      #duration
      dur<-ConTime(t_dat$TIME[nrow(t_dat)]-t_dat$TIME[1])

      #mean speed
      sv <- which(t_dat$SPEED>15)
      spd <- round(mean(t_dat$SPEED[sv]),digits=1) 
      
      #flight duration
      dt <- diff(t_dat$TIME)
      fdur<-ConTime(sum(dt[sv],na.rm=T))
      
      #pause duration 
      pv <- which(t_dat$SPEED<=15)
      pdur<-ConTime(sum(dt[pv],na.rm=T))
      
      #number of pauses
      rp <- rep(0,nrow(t_dat))
      rp[pv]<-1
      rlp <- rle(rp)
      np <- length(which(rlp$values==1 & rlp$length>30))
      
      #write traditional results
      vdx <- which(vv==2500)
      res_strng <- paste(c(paste(IDV,collapse=';'),csdir,vvb_res[vdx],ce,spd,dur,fdur,pdur,np),collapse=";")
      fn <- paste(mdir,"/res/TRD_RES.csv",sep='')
      if (!file.exists(fn))
      {
        write.table(t(c('DATE','ID','EXP','SITE','VB','CE','SPD','DUR','FDUR','PDUR','NP')), file = fn, append = F, quote = F, sep = ";",
                    eol = "\n", na = "NA", dec = ".", row.names = F,
                    col.names = F, qmethod = c("escape", "double"))
      }
      write.table(res_strng, file = fn, append = T, quote = F, sep = ";",
                  eol = "\n", na = "NA", dec = ".", row.names = F,
                  col.names = F, qmethod = c("escape", "double"))
      
      #write vvb results
      fn <- paste(mdir,"/res/VVB_RES.csv",sep='')
      if (!file.exists(fn))
      {
        write.table(t(c('DATE','ID','EXP','SITE',vv)), file = fn, append = F, quote = F, sep = ";",
                    eol = "\n", na = "NA", dec = ".", row.names = F,
                    col.names = F, qmethod = c("escape", "double"))
      }
      vvb_strng <- paste(c(paste(IDV,collapse=';'),csdir,t(vvb_res)),collapse=";")
      write.table(vvb_strng, file = fn, append = T, quote = F, sep = ";",
                  eol = "\n", na = "NA", dec = ".", row.names = F,
                  col.names = F, qmethod = c("escape", "double"))
    }
  }
}
print(ptm-proc.time())
rm(list=ls())

