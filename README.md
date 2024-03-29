# Building Land Use Regression model - an example

## Overview of the project

- https://adithi-spatial-learning.netlify.app/, made using [Quarto](https://quarto.org/). 
- A portion of LUR project data analysis is shown here as an example.

## Walk through the code

- https://walk-through-lur.netlify.app/, made using [Quarto](https://quarto.org/). 
- The resultant dataframe's 5 observations are shown. 

## Folder overview 
- `R` folder  
  - `functions.R` all custom made functions are defined here
- `index.qmd` is an example document on how to run the Land Use Regression model
- `img` folder has images used in here
- `data` folder has data
  - `training_data\yearly_avg_2018_2021.csv` and `Delhi_code_spatial_info\Delhi_site_lat_long_info.csv` these are the measurement sites in Delhi and the file used for the training, validation, and prediction, this file needs to have a column called `CODE` to uniquely identify each point for the analysis; this is usually a `multipoint` object
  - `airport\airport_point.csv` latitude, longitude data of airport runway in Delhi
  - `AOD\interpolated\Delhi_AOD_map_2019.tif` a raster file of aerosol  optical depth of the year 2019 in Delhi, this is resampled using bilinear interpolation to 50 m
  - `DEM\Delhi_SRTM_Elevation.tif` a raster file of the digital elevation model of Delhi
  - `direction_of_effect\parameters_pm.csv` is the file for the desired direction of effect, the columns in it are `param` - the parameter / variable name, and the direction of effect for them either positive (> 0) or negative (< 0) or no effect with respect to the response variable. 
  - `industries\Delhi_Industry_Details.csv` locations of the various industries in Delhi 
  - `landuse\lulc_vars.csv` area covered by various land use type buffers around the measurement sites in Delhi derived using [Earth Engine](https://earthengine.google.com/)
  - `population\pop_vars.csv` population data in the buffers around the measurement sites in Delhi derived using Earth Engine
  - `railways\Delhi_railways.shp` railway shapefile of Delhi


## Methodology 

#### Eeftens, M., Beelen, R., De Hoogh, K., Bellander, T., Cesaroni, G., Cirach, M., . & Hoek, G. (2012). Development of land use regression models for PM2.5, PM2.5 absorbance, PM10 and PMcoarse in 20 European study areas; results of the ESCAPE project. Environmental Science & Technology, 46(20), 11195-11205. https://doi.org/10.1021/es301948k

## Variables / Parameters and their buffers 

![\label{fig:parameters}](img/parameters_2.jpg)
    
## [Land use parameters](https://developers.google.com/earth-engine/datasets/catalog/ESA_WorldCover_v100) 

class - renamed_as 

10 - tree_cover

20 - shrubland

30 - shrubland

40 - cropland

50 - builtup

60 - bare_land

70 - snow_ice

80 - per_water_bodies

90 - shrubland

95 - mangroves

100 - moss_lichen

## Column names 

```{r}
"CODE"                    ,   "PM2.5"                      
"lat"                     ,   "long"                       
"tree_cover_buffer_100"   ,   "tree_cover_buffer_1000"     
"tree_cover_buffer_300"   ,   "tree_cover_buffer_500"      
"tree_cover_buffer_5000"  ,   "shrubland_buffer_100"       
"shrubland_buffer_1000"   ,   "shrubland_buffer_300"       
"shrubland_buffer_500"    ,   "shrubland_buffer_5000"      
"cropland_buffer_100"     ,   "cropland_buffer_1000"       
"cropland_buffer_300"     ,   "cropland_buffer_500"        
"cropland_buffer_5000"    ,   "builtup_buffer_100"         
"builtup_buffer_1000"     ,   "builtup_buffer_300"         
"builtup_buffer_500"      ,   "builtup_buffer_5000"        
"bare_land_buffer_100"    ,   "bare_land_buffer_1000"      
"bare_land_buffer_300"    ,   "bare_land_buffer_500"       
"bare_land_buffer_5000"   ,   "per_water_bod_buffer_100"   
"per_water_bod_buffer_1000",   "per_water_bod_buffer_300"   
"per_water_bod_buffer_500",   "per_water_bod_buffer_5000"  
"rail_buffer_1000"        ,   "rail_buffer_500"            
"rail_buffer_5000"        ,   "pop_buffer_300m"            
"pop_buffer_500m"         ,   "pop_buffer_1000m"           
"pop_buffer_5000m"        ,   "inverse_distance_industries"
"inverse_distance_airport",   "elevation"                  
"aod"  

```

## Other similar models 

#### Machine learning model:
  - [`Random forest`](https://link.springer.com/article/10.1023/a:1010933404324?utm_source=getftr&utm_medium=getftr&utm_campaign=getftr_pilot) using [`ranger`](https://cran.r-project.org/web/packages/ranger/ranger.pdf)

#### Geostatistical model:
  - [`Geographically weighted regression`](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1538-4632.1996.tb00936.x) using [`GWmodel`](https://cran.r-project.org/web/packages/GWmodel/GWmodel.pdf)

## Session information and package versions 

![\label{fig:session_info}](img/session_info.jpg)

## Earth Engine parameter download links

- [LULC predictors extraction](https://code.earthengine.google.com/17ed646fc4afed19d58b6d3c8ebd9f5d) 
- [NDVI predictors extraction](https://code.earthengine.google.co.in/bab3d8ccce61b72e4719f5b9b855031d)
- https://code.earthengine.google.co.in/c0c232549638701710e359d057af9fa4
- https://code.earthengine.google.co.in/3bb33602b139cbd8559c2b3d31fa503c
- https://code.earthengine.google.co.in/4ea4a284dabdcad8d1a69892d2d4b9d1
- https://code.earthengine.google.co.in/b01c5360bdb674104a0c377681c55bf4
- https://code.earthengine.google.co.in/fc9544c89948a6416f2e049cfd00bd95
- https://code.earthengine.google.co.in/581db9ee4a89363de6847cf12a2faaf3
- https://code.earthengine.google.co.in/4d1f67bd2f1715287b84aa65653cec3f
- https://code.earthengine.google.co.in/bb1889f92a3490729bb58ab0e32903a4
- https://code.earthengine.google.co.in/192b69b5716f853733cbe431ebec8b3b
- https://code.earthengine.google.co.in/81286c7bfe5460abe1511723d5503fcc
- https://code.earthengine.google.co.in/994671b6b43489fa6909bf8ba71fd094
- https://code.earthengine.google.co.in/0180b79190550b390882cf933a53908a

## Resources 

- [Spatial Data Science in R](https://r-spatial.org/book/)

## LICENSE

[Click here.](https://github.com/adithirgis/code_examples/blob/main/LICENSE)

## Suggestions / comments 

- Please create a pull request in this repository. 
