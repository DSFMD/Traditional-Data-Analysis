# DSFMD - Dynamic Systems Framework for Movement Data
# Example for using leaflet for interactive display of gps tracks
# Code by Ingo Schiffner 2018

#set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(leaflet) # required for interactive map display
library(sp) # spatial data frame
library(data.table) # better than data frame
library(htmlwidgets) # for saving html widgets

#sub/helper functions
source('../SHARED/GetPGPSData.R') # read GPS data from csv file

#compile sub/helper functions
GetPGPSData <- cmpfun(GetPGPSData)

#projection string for conversion to wgs84 latitude and longitude format
wgs84<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

#projection string for Universal Transverse Mercator (UTM) format
utm<-"+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#select data folder
cdir<-choose.dir(default = "", caption = "Select folder")
SubDirNames <- list.files(path = cdir, pattern = NULL, all.files = FALSE, full.names = FALSE, recursive = FALSE, ignore.case = FALSE)

vdx <- c(2500)
IDVC<-NULL

#load meta data
mdir<-dirname(cdir)
loc_dt <- read.table(paste(mdir,"/meta_data/sites.txt",sep=''), header=T, sep=",", quote='"', dec=",", na.strings="", colClasses="character", fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)
loc_dt$lat <- as.numeric(loc_dt$lat)
loc_dt$lon <- as.numeric(loc_dt$lon)
coordinates(loc_dt)<-~lon+lat
proj4string(loc_dt)<-CRS(wgs84)

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
    addProviderTiles('Esri.WorldImagery', group = "Satellite")
  
  #colors for different tracks
  factpal <- colorFactor(topo.colors(length(FileNames)), FileNames, na.color = "#808080")
  
  #process all files in subfolder
  for (cfile in FileNames)
  {
    #read data
    fullfile <- paste(csdir_full,'\\',cfile ,sep="")
    t_dat <- GetPGPSData(fullfile)
    
    #convert to spatial data frmae
    sdf <- data.frame(t_dat$X,t_dat$Y)
    colnames(sdf) <- c("lon","lat")
    coordinates(sdf)<-~lon+lat
    proj4string(sdf)<-CRS(utm)
    sdf<-spTransform(sdf, wgs84)
    fdx <- which(FileNames==cfile)
    ll<-addPolylines(ll,data=sdf@coords, group=cfile ,color=factpal(cfile), weight = 2, opacity = 1)
    
  }
  #layers control elements
  ll<-addCircleMarkers(ll,data=loc_dt[sdx,],
                     radius = 10,
                     weight = 3,
                     color = "#FFFFFF",
                     stroke = T,
                     fillColor = "#00FF00",
                     fillOpacity = 1,
                     group ="Release Site") %>%
    
    addCircleMarkers(data=loc_dt[1,],
                     radius = 10,
                     weight = 3,
                     color = "#FFFFFF",
                     stroke = T,
                     fillColor = "#FF0000",
                     fillOpacity = 1,
                     group ="Home Loft") %>%
    addCircles(data=loc_dt[sdx,],
                     radius = vdx,
                     weight = 3,
                     color = "#000000",
                     stroke = T,
                     fillColor = "#0000FF",
                     fillOpacity = 0,
                     group ="Vanishing Bearing") %>%
    addLayersControl(baseGroups = c("Positron", "Toner", "Satellite", "Relief"), overlayGroups = c("Release Site","Home Loft","Vanishing Bearing", FileNames),
                       options = layersControlOptions(collapsed = FALSE))
  
  #save plot
  dir.create(paste(mdir,"/plot/",sep=''))
  f<-paste(mdir,'/plot/',csdir,".html",sep='')
  saveWidget(ll,file.path(normalizePath(dirname(f)),basename(f)))
}
rm(list=ls())
