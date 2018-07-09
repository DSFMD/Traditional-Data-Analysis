# DSFMD - Dynamic Systems Framework for Movement Data
# example implementation for calculation of traditional measurements for gps data
# Code by Ingo Schiffner 2018

#Traditional Analysis
source('GetVVB.R')
vdx <- c(500,1000,1500,2000,2500) # take vanishing bearings at following intervals

#sub/helper functions
source('CalAng.R')

library(sp) #spatial transformation library
library(data.table) # better than data frame

#projection string for conversion to wgs84 latitude and longitude format
wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#projection string for Universal Transverse Mercator (UTM) format
utm<-"+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#select data folder
cdir<-choose.dir(default = "", caption = "Select folder")
SubDirNames <- list.files(path = cdir, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)

#load meta data
loc_dt <- read.table("meta_data/sites.txt", header=T, sep=",", quote='"', dec=",", na.strings="", colClasses="character", fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)
loc_dt$lat <- as.numeric(loc_dt$lat)
loc_dt$lon <- as.numeric(loc_dt$lon)
coordinates(loc_dt)<-~lon+lat
proj4string(loc_dt)<-CRS(wgs84)
loc_dt<-spTransform(loc_dt, utm)

#bird meta data
bid_dt <- as.data.table(read.table("meta_data/birds.txt", header=T, sep=",", quote='"', dec=",", na.strings="", colClasses="character", fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE))

#process subfolders
gps_data_dt <- NULL
ref_data_dt <- NULL

dir.create('res')
write.table(t(c('FID',vdx)), file = "res//vvb.csv", append = F, quote = F, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"))

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
    fullfile <- paste(csdir_full,'\\',cfile ,sep="")
    t_dat <- as.data.table(read.table(fullfile, header=T, sep=",", quote='"', dec=".", na.strings="", colClasses="character", fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)) 
    t_dat<-t_dat[,1:6]
    colnames(t_dat)<-c('DATE','TIME','LAT','LON','ALT','SPEED')
    t_dat$LAT<-as.numeric(t_dat$LAT)
    t_dat$LON<-as.numeric(t_dat$LON)
    t_dat$ALT<-as.numeric(t_dat$ALT)
    t_dat$SPEED<-round(as.numeric(t_dat$SPEED)/1000, digits=1) 
    t_dat <- subset(t_dat,t_dat$SPEED>0)
    
    #merge date and time
    t_dat$TIME <- as.numeric(as.POSIXct(paste(t_dat$DATE,t_dat$TIME),"%Y/%m/%d %H:%M:%S", tz='UTC'))
    t_dat$DATE<-NULL
    
    #convert to utm
    sdf <- data.frame(t_dat$LAT,t_dat$LON)
    colnames(sdf) <- c("lat","lon")
    coordinates(sdf)<-~lon+lat
    proj4string(sdf)<-CRS(wgs84)
    sdf<-spTransform(sdf, utm)
    t_dat$LAT<-sdf@coords[,1]
    t_dat$LON<-sdf@coords[,2]
    colnames(t_dat)[2:3]<-c('X','Y')
    
    #calculate virtual vanishing bearings
    vvb_res<-GetVVB(t_dat,sdx,loc_dt,vdx)
    
    #write vvb results
    write.table(t(c(cfile,vvb_res)), file = "res//vvb.csv", append = T, quote = F, sep = ",",
                eol = "\n", na = "NA", dec = ".", row.names = F,
                col.names = F, qmethod = c("escape", "double"))
  }
}

