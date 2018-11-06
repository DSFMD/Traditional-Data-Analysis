# Traditional-Data-Analysis
Traditional Methods for Analysis of Animal Navigation</br>
<b><i>requires the files in the "shared" repository</i></b></br>
Gps data should be sorted by release site with subdirectories for each site. Code assumes that the file format is .csv, as produced by 
IgotU GPS recorders. Format should be: Date, Time, Latitude, Longitude, Altitude, Speed, Course, Type, Distance, Essential. File naming convention is YYMMDD_XX-YYYYY_EXP (i.e. DATE (in reverse order) BIRD ID and type of EXPERIMENT). You will also need to set up a folder meta_data containing a file called sites. txt, file needs to contain name and the wgs84 coordinates(lat/lon) for the release sites e.g.:</br></br>
Site,lat,lon</br>
HOME,53.216808,-4.173131</br>
EXAMPLE,53.062065,-4.253445</br>

<i>NOTE: first entry must be HOME</i>
  </br>
  <ul>
  <li><b>Plotting:</b>
<li>DSFMD_PGPS_PLOT_GROUP.R -  Example for using leaflet for interactive display of GPS tracks by Group
<li>DSFMD_PGPS_PLOT_INDIVIDUAL.R -  Example for using leaflet for interactive display of GPS tracks by Individual
<li>DSFMD_PGPS_PLOT_REPEATED.R -  Example for using leaflet to plot repeated flights (including offsite releases)
  </ul>
  </br>
  <ul>
  <li><b>Basic Measurements:</b>
<li>DSFMD_PGPS_TRAD.R - Example implementation for calculation of traditional measurements for GPS data (e.g. efficiency, flight speed vanishing bearings)
<li>DSFMD_PGPS_REPCMP.R - Example implementation for calculation for comparing repeated flights (inlcuding nearest neighbour distance, normalized mutual information, angular prediction error)
</ul>
