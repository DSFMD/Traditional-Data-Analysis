# DSFMD ###########################################################################################################
# function to calculate angles between pairs of points in space
# INPUT:
# dx vector of delta x (~latitude;~cosine) 
# dy vector of delta y (~longitude;~sine) 
# edst vector of euclidean distance
# OUPUT:
# angle vector of angles
# Original code by Ingo Schiffner 2005-2011 coverted to R by Ingo Schiffner 2018
CalAng <- function (dx,dy,edst)
{
  #special cases
  if (dx==0 | dy==0)
  {
    if (dx == 0 & dy > 0){
      return(90)
    } else if (dx == 0 & dy < 0) {
      return(270)
    } else if (dy == 0 & dx >= 0) {
      return(360)
    } else if (dy == 0 & dx < 0) {
      return(180)
    }
  }
  
  #convert to degrees
  angle <- atan((dy/edst)/(dx/edst)) * (180/pi)
  
  if (dx > 0 & dy > 0) { # Quadrant 1, upper right
    angle <- angle
    #convert 0 degrees to 360 degrees
    if (angle == 0) 
    { 
      angle = 360 
    }
  } else if (dx < 0) {  # Quadrant 2 and 3, lower right and lower left
    angle <- 180 + angle
  } else if (dx > 0 & dy < 0) {  # Quadrant 4, upper left
    angle <- 360 + angle 
  }
  return(angle)
}