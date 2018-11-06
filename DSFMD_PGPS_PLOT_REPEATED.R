# DSFMD - Dynamic Systems Framework for Movement Data
# Example for using leaflet to plot repeated flights (including offsite releases)
# Code by Ingo Schiffner 2018

#set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#libraries
library(sp) #spatial transformation library
library(data.table) # better than data frame
library(leaflet) # required for interactive map display
library(htmlwidgets) # for saving html widgets
library(compiler) # byte compiler for faster execution of functiuons
library(mapview) # convert maps to png

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

#load meta data
mdir<-dirname(cdir)
loc_dt <- read.table(paste(mdir,"/meta_data/sites.txt",sep=''), header=T, sep=",", quote='"', dec=",", na.strings="", colClasses="character", fill=TRUE, strip.white=TRUE, blank.lines.skip=TRUE)
loc_dt$lat <- as.numeric(loc_dt$lat)
loc_dt$lon <- as.numeric(loc_dt$lon)
coordinates(loc_dt)<-~lon+lat
proj4string(loc_dt)<-CRS(wgs84)
loc_dt<-spTransform(loc_dt, utm)

gps_dt <-NULL
#get data from all subfolders
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
rm(t_dat,csdir, csdir_full, FileNames, sdx, SubDirNames, spl, dx, dy, de, fdx, hdx)

#determine main release site
tmp <- gps_dt[,6:8]
tmp <- as.data.table(lapply(tmp, as.factor))
tmp$CNT<-1
tmp<-aggregate(CNT~SITE+DATE,tmp,max)
scnt<-aggregate(CNT~SITE,tmp,sum)
mdx <- which(scnt$CNT==max(scnt$CNT))
rm(scnt)
sv <- unique(gps_dt$SITE)
bv <-unique(gps_dt$ID)

GV <- c('Start','Mid','End','Offsite')
CV <- c('#00FF00','#0000FF','#FF0000','#FFFFFF')
dir.create(paste(mdir,"/plot/",sep=''))
loc_dt<-spTransform(loc_dt, wgs84)

for (i in bv)
{
  
  ll<-leaflet(options = leafletOptions(zoomControl = FALSE)) %>% addTiles() %>%
    #load options for background maps
    addProviderTiles('Esri.WorldImagery', group = "Satellite")
  
  #get data
  tmp <- subset(gps_dt,ID==i)
  dv <- unique(tmp$DATE)
  
  #make sure data is ordered
  dv<-dv[order(dv)]
  mc <- length(dv)-length(sv)+1-3
  
  #reorder data
  v1 <- 1:3
  v2 <- 4:mc
  v3 <- (mc+1):(mc+3)
  vo <- (mc+4):length(dv)
  dxv <- c(v2,vo,v1,v3)
  
  for (j in dxv)
  {
    #get data for current track
    t_dat <- subset(tmp,DATE==dv[j])
    
    #convert to spatial data frmae
    sdf <- data.frame(t_dat$X,t_dat$Y)
    colnames(sdf) <- c("lon","lat")
    coordinates(sdf)<-~lon+lat
    proj4string(sdf)<-CRS(utm)
    sdf<-spTransform(sdf, wgs84)
    
    #choose group
    if (j<=3){
      gdx <- 1
      pv <-300
    } else if (j>mc) {
      gdx <- 3
      pv <-300
    } else {
      gdx <- 2
      pv <-200
    }
    if (t_dat$SITE!=sv[mdx]) { gdx <- 4 }
      
    ll<-addPolylines(ll,data=sdf@coords, group=GV[gdx] ,color=CV[gdx], weight=2, opacity=1, options = list(zIndex = pv))
  }
  #layers control elements
  sdx <- which(loc_dt$Site==sv[mdx])
  ll<-addCircleMarkers(ll,data=loc_dt[sdx,],
                       radius = 10,
                       weight = 3,
                       color = "#FFFFFF",
                       opacity = 1,
                       stroke = T,
                       fillColor = "#00FF00",
                       fillOpacity = 1,
                       group ="Release Site") %>%
    
    addCircleMarkers(data=loc_dt[1,],
                     radius = 10,
                     weight = 3,
                     color = "#FFFFFF",
                     opacity = 1,
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
               group ="Vanishing Bearing")

    # addLayersControl(overlayGroups = c("Release Site","Home Loft","Vanishing Bearing", GV),
    #                  options = layersControlOptions(collapsed = FALSE))
  
  #save plot
  #f<-paste(mdir,'/plot/',i,"_rep.html",sep='')
  #saveWidget(ll,file.path(normalizePath(dirname(f)),basename(f)))
  f<-paste(mdir,'/plot/',i,"_rep.png",sep='')
  mapshot(ll, url = NULL, file = f, remove_url = TRUE)
}
rm(list=ls())
