# function to calculate virtual vanishing bearings
# with...
# t_dat data.frame with X and Y coordinates (UTM)
# sdx location index
# loc_dt data table containing the location of all release sites
# vdx indices at which vanishing bearings are to be calculated
# code by Dr. Ingo Schiffner 2018
GetVVB <- function (t_dat,sdx,loc_dt,vdx)
{
  #calculation of distance and direction from release site
  dx <- t_dat$X[2:nrow(t_dat)] - loc_dt@coords[sdx,1]
  dy <- t_dat$Y[2:nrow(t_dat)] - loc_dt@coords[sdx,2]
  dst_s <- (dx^2+dy^2)^0.5 
  dir_s <- round(mapply(CalAng, dx, dy, dst_s),digits=1)
  vvb_res <- NULL
  for (i in vdx)
  {
    cdx <- which(dst_s >= i)
    fdx <- min(cdx)
    vvb_res <- rbind(vvb_res,dir_s[fdx])
  }
  return (vvb_res)
}