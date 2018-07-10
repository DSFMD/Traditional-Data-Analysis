# DSFMD - Dynamic Systems Framework for Movement Data
# example for using leaflet for interactive display of gps tracks
# Code by Ingo Schiffner 2018

library(leaflet) # required for interactive map display
library(data.table) # better than data frame
library(htmlwidgets) # for saving html widgets

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

#bird meta data
bid_dt <- as.data.table(read.table("meta_data/birds.txt", header=T, sep=",", quote='"', dec=",", na.strings="", colClasses="character", fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE))

#process subfolders
gps_data_dt <- NULL
ref_data_dt <- NULL

dir.create('plot')
for (csdir in SubDirNames)
{
  print(paste("processing data from", csdir))
  csdir_full <- paste(cdir,'\\',csdir,sep="")
  isdir <- !file_test("-f", csdir_full) # returns TRUE
  FileNames <- list.files(path = csdir_full, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)
  
  sdx <- which(loc_dt$Site == csdir)
  
  ll<-leaflet() %>% addTiles() %>%
    #load options for background maps
    addProviderTiles("CartoDB.Positron", group="Positron") %>%
    addProviderTiles("Stamen.Toner", group = "Toner") %>%
    addProviderTiles('Esri.WorldImagery', group = "Satellite") %>%
    addProviderTiles('Esri.WorldShadedRelief', group = "Relief")
  
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
    
    #convert to spatial data frmae
    sdf <- data.frame(t_dat$LAT,t_dat$LON)
    colnames(sdf) <- c("lat","lon")
    coordinates(sdf)<-~lon+lat
    proj4string(sdf)<-CRS(wgs84)
    
    ll<-addPolylines(ll,data=sdf@coords, group = cfile ,color='#FF0000', weight = 2, opacity = 1)
    
  }
  #layers control elements
  ll<-addLayersControl(ll,baseGroups = c("Positron", "Toner", "Satellite", "Relief"), overlayGroups = FileNames,
                       options = layersControlOptions(collapsed = FALSE))
  saveWidget(ll, file=paste('.\\plot\\',csdir,".html",sep=''))
}
