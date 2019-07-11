#############################
# Script to collect and manipulate input data for protein mapping project
# Script is divided into different sections for each data source
# URL, paper, etc for each data source is noted
#############################

#######
# Read me
# This short intro section defines the map projections
# resolution the data is sampled to. 
# I have no preference on the projection we use, and 
# am sampling to wgs84 in this script because it 
# is one of the more commonly used projections
# For resolution, we`ll need to make a decision based
# on the resolution of the different data sets we`re
# using.
# 
# Projection and resolution can be changed by changing
# the lines of code in this section. 
# 
# Projections can be defined in two ways.
# First, using a proj4 string.
# Second, using the script crs("+init epsg:projection number")
#######


#####
### Installing libraries
library(raster)
library(rgdal)
library(gdalUtils)
library(maptools)
library(dplyr)
library(countrycode)
library(ncdf4)
library(reshape2)
library(tidyr)

#####
### Defining projection

###
# First method, using proj4 string
crs.out = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

###
# Second method, using epsg number
# Change the numerical code to use a different epsg code
# Numbers associated with proj4 strings can be found
# by searching the website "spatialreference.org"
# For a specific epsg numerical code
crs.out <- crs("+init=EPSG:4326")


#####
### Defining resolution
# Doing this by creating an empty raster with the
# projection, resolution, and geographic extent we're using
# Note that resolution will be automatically calculated based on 
# 1) the dimensions of the raster
# 2) The extent of raster

# For now, creating a raster that is at 0.5 by 0.5 degree resolution
# We'll want to use something more resolute than this when running analyses
# But setting this resolution for now because resampling and reprojecting
# to more resolute rasters
# 1) Takes a while, and
# 2) Can result in fairly large rasters (few hundred megabytes) that aren't easy to email

###
# First step is to create an empty matrix with the number of cells we want
# Current creating a matrix with 360 rows and 720 columns, filled entirely with NAs 

template.matrix <- matrix(nrow = 1800/.83333, ncol = 3600/.83333, NA)

# Now converting into a raster with the pre-defined CRS
template.raster <- raster(template.matrix, 
                       crs = crs.out)

# And updating extent
# Extent is updated in the order (xmin, xmax, ymin, ymax)
extent(template.raster) <- c(-180, 180, # updating x values
                             -90, 90) # updating y values

#####
# Data set for travel time to urban areas
# Data is from Weiss et al 2018, Nature
# Paper URL is: https://www.nature.com/articles/nature25181
# Data URL is: https://map.ox.ac.uk/research-project/accessibility_to_cities/
# Data file is the "Accessability to cities" zip file
# Apologies for the file size...it's a bit large...
#####

###
# Importing data
raw.travel <- raster("/Users/maclark/Desktop/Mapping Food Consumption/Raw data/Accessibility to cities/accessibility_to_cities_2015_v1.0/accessibility_to_cities_2015_v1.0.tif")
# Removing negative values
raw.travel[raw.travel < 0 & !is.na(raw.travel)] <- NA

###
# Reprojecting and updating resolution
# Good news! A single function does all of this for us
raw.travel <- projectRaster(from = raw.travel,
                            to = template.raster,
                            method = 'bilinear',
                            na.rm = TRUE)

# Note that data cleaning will be done later
# Updating data in rasters takes a very very long time
# But updating data when it's in a data set (a vector) takes much less time
# For this data set, -9999 indicates NA values

#####
# Data for total GDP and per capita GDP
# Data is from Kummu et al 2018, Scientific Data
# Paper URL is: https://www.nature.com/articles/sdata20184
# Data URL is: https://datadryad.org/resource/doi:10.5061/dryad.dk1j0
# Data set is "GDP_PPP_30arcsec" file
#####

# Slightly different format for importing here
# The raw ".nc" file has multiple layers
# Layer 1 = 1990 data, layer 2 = 2000 data, layer 3 = 2015 data
# Need to specify which layer to import, by using band = "layer number"

###
# 2015 data

raw.gdp2015 <- raster("/Users/maclark/Desktop/Mapping Food Consumption/Raw data/GDP Data/GDP_PPP_30arcsec_v2.nc",
                  band = 3) # Specifying layer number, this imports the 2015 data
# Reprojecting and resampling
raw.gdp2015 <- projectRaster(from = raw.gdp2015,
                         to = template.raster,
                         method = 'bilinear',
                         na.rm = TRUE)

# Reprojecting and resampling interpolates values from cells, so this no longer sums to total GDP
# Value in raster is now GDP produced per square km, rather than total GDP in the cell
# As such, need to multiply values in each cell by the cells area

# Creating a raster that contains the area of each cell
area.raster <- area(raw.gdp2015)

# Now have two options to adjust
# 1) Multiply the rasters
# This works fine on small rasters, but I would not recommend doing this with large rasters
# (raster math in R is not particulary efficient)
# 2) Extract values, multiply the values, and then reinput them into the raster
# Raw data contains GDP by cell
# Reprojecting interpolates values based on neighboring cells
# Which means reprojected raster needs to be adjusted by cell area

# Going to do this by extracting balues and then multiplying the vectors
# Because this is a fair bit quicker than multiplying rasters

# Extracting values
gdp.data <- getValues(raw.gdp2015)
area.values <- getValues(area.raster)

# Multiplying vector
gdp.total.data <- gdp.data * area.values

# Recreating raster
gdp.2015.raster <- raster(matrix(gdp.total.data,
                                 nrow = template.raster@nrows,
                                 ncol = template.raster@ncols,
                                 byrow = TRUE),
                          crs = crs.out)
# And adjusting extent
extent(gdp.2015.raster) <- extent(template.raster)


###
# 2000 data
raw.gdp2000 <- raster("/Users/maclark/Desktop/Mapping Food Consumption/Raw data/GDP Data/GDP_PPP_30arcsec_v2.nc",
                      band = 2) # Specifying layer number, this imports the 2000 data
# Reprojecting and resampling
raw.gdp2000 <- projectRaster(from = raw.gdp2000,
                             to = template.raster)

# And adjusting by cell size, as above
# Extracting values
gdp.data <- getValues(raw.gdp2000)

# Multiplying vector
gdp.total.data <- gdp.data * area.values

# Recreating raster
gdp.2000.raster <- raster(matrix(gdp.total.data,
                                 nrow = template.raster@nrows,
                                 ncol = template.raster@ncols,
                                 byrow = TRUE),
                          crs = crs.out)
# And adjusting extent
extent(gdp.2000.raster) <- extent(template.raster)

###
# 1990 data
raw.gdp1990 <- raster("/Users/maclark/Desktop/Mapping Food Consumption/Raw data/GDP Data/GDP_PPP_30arcsec_v2.nc",
                      band = 1) # Specifying layer number, this imports the 1990 data
# Reprojecting and resampling
raw.gdp1990 <- projectRaster(from = raw.gdp1990,
                             to = template.raster)

# And adjusting by cell size, as above
# Extracting values
gdp.data <- getValues(raw.gdp1990)

# Multiplying vector
gdp.total.data <- gdp.data * area.values

# Recreating raster
gdp.1990.raster <- raster(matrix(gdp.total.data,
                                 nrow = template.raster@nrows,
                                 ncol = template.raster@ncols,
                                 byrow = TRUE),
                          crs = crs.out)
# And adjusting extent
extent(gdp.1990.raster) <- extent(template.raster)


#####
# Now getting population maps
# Note that this is population per grid cell, not population density
# Using the same maps used in Kummu et al 2018 (where GDP data comes from)
# This way we are consistent across data sets
# Dataset is from GHS Population Grid, Gridded Human Settlement layerStats
# For more info, website is: http://drdsi.jrc.ec.europa.eu/dataset/ghs-population-grid-derived-from-gpw4-multitemporal-1975-1990-2000-2015
# Data for 2015 is from: http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GPW4_GLOBE_R2015A/GHS_POP_GPW42015_GLOBE_R2015A_54009_1k/V1-0/
# Data for 2000 is from: http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GPW4_GLOBE_R2015A/GHS_POP_GPW42000_GLOBE_R2015A_54009_1k/V1-0/
# Data for 1990 is from: http://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GPW4_GLOBE_R2015A/GHS_POP_GPW41990_GLOBE_R2015A_54009_1k/V1-0/
# And manipulating this data takes...a while
# These are the 1km x 1km resolution maps, there are also 250m x 250m maps
#####

###
# 2015 data
# Importing
raw.popdensity.2015 <- raster("/Users/maclark/Desktop/Mapping Food Consumption/Raw data/Pop Density Data/GHS_POP_GPW42015_GLOBE_R2015A_54009_1k_v1_0/GHS_POP_GPW42015_GLOBE_R2015A_54009_1k_v1_0.tif")

# Sanity check on total population
# This should sum to ~7.3 billion, or global population in 2015
# tmp.vector <- getValues(raw.popdensity.2015)
# sum(tmp.vector, na.rm = TRUE)/1e9


# Reprojecting and resampling
raw.popdensity.2015 <- projectRaster(raw.popdensity.2015,
                                     template.raster,
                                     method = 'bilinear')

# Reprojecting and resampling interpolates values from cells, so this no longer sums to global population
# Value in raster is now population density per km, rather than population in the cell
# As such, need to multiply values in each cell by the cells area

# Creating a raster that contains the area of each cell
area.raster <- area(raw.popdensity.2015)

# Now have two options to adjust
# 1) Multiply the rasters
# This works fine on small rasters, but I would not recommend doing this with large rasters
# (raster math in R is not particulary efficient)
# 2) Extract values, multiply the values, and then reinput them into the raster

# We're going with option 2 because it's a fair bit faster
# Extracting values
pop.values <- getValues(raw.popdensity.2015)
area.values <- getValues(area.raster)

# Multiplying values
pop.total.values <- pop.values * area.values

# And recreating the raster
pop.2015.raster <- raster(matrix(pop.total.values,
                                 nrow = template.raster@nrows,
                                 ncol = template.raster@ncols,
                                 byrow = TRUE),
                          crs = crs.out)
# And adjusting extent
extent(pop.2015.raster) <- extent(template.raster)



###
# 2000 data
# Importing
raw.popdensity.2000 <- raster("/Users/maclark/Desktop/Mapping Food Consumption/Raw data/Pop Density Data/GHS_POP_GPW42000_GLOBE_R2015A_54009_1k_v1_0/GHS_POP_GPW42000_GLOBE_R2015A_54009_1k_v1_0.tif")

# Sanity check, should sum to ~6.1 billion
# tmp.vector <- getValues(raw.popdensity.2000)
# sum(tmp.vector, na.rm = TRUE)
# Reprojecting and resampling
raw.popdensity.2000 <- projectRaster(raw.popdensity.2000,
                                     template.raster)

# Same as above, this raster now contains population density, not population
# Need to adjust

# Extracting pop values
pop.values <- getValues(raw.popdensity.2000)
# Multiplying values
pop.total.values <- pop.values * area.values

# And recreating the raster
pop.2000.raster <- raster(matrix(pop.total.values,
                                 nrow = template.raster@nrows,
                                 ncol = template.raster@ncols,
                                 byrow = TRUE),
                          crs = crs.out)
# And adjusting extent
extent(pop.2000.raster) <- extent(template.raster)


###
# 1990 data
# Importing
raw.popdensity.1990 <- raster("/Users/maclark/Desktop/Mapping Food Consumption/Raw data/Pop Density Data/GHS_POP_GPW41990_GLOBE_R2015A_54009_1k_v1_0/GHS_POP_GPW41990_GLOBE_R2015A_54009_1k_v1_0.tif")

# Sanity check, should sum to ~5.3 billion
# tmp.vector <- getValues(raw.popdensity.1990)
# sum(tmp.vector, na.rm = TRUE)
# Reprojecting and resampling
raw.popdensity.1990 <- projectRaster(raw.popdensity.1990,
                                     template.raster)

# And as with above, need to adjust for population density vs population size per cell
# Extracting pop values
pop.values <- getValues(raw.popdensity.1990)
# Multiplying values
pop.total.values <- pop.values * area.values

# And recreating the raster
pop.1990.raster <- raster(matrix(pop.total.values,
                                 nrow = template.raster@nrows,
                                 ncol = template.raster@ncols,
                                 byrow = TRUE),
                          crs = crs.out)
# And adjusting extent
extent(pop.1990.raster) <- extent(template.raster)

sum(pop.total.values, na.rm = TRUE)/1e9

###
# Converting 0s to NAs
pop.2015.raster[pop.2015.raster == 0] <- NA
pop.2000.raster[pop.2000.raster == 0] <- NA
pop.1990.raster[pop.1990.raster == 0] <- NA

###
# Getting maps of per capita GDP
pc.gdp.2015 <- gdp.2015.raster / pop.2015.raster
pc.gdp.2000 <- gdp.2000.raster / pop.2000.raster
pc.gdp.1990 <- gdp.1990.raster / pop.1990.raster




#####
# Map of primary religion
# Data is from: http://worldmap.harvard.edu/data/geonode:wrd_province_religion_qg0
# Data download is upper right corner of the web page
# Original project is: https://worldreligiondatabase.org/
# The meta data for this one is a mess
# as in they don't seem to have a file that IDs different numeric values with different primary religion values
# So making one ourselves
# To do this, need to download shapefile
# Extract province IDs and numeric values from the shapefile
# Download the .csv file, which indicates province id, primary relgion, and % primarily religion
# Clean the .csv file manually, because it isn't formatted properly
#####

###
# Loading shapefile
religion.shape <- readOGR("/Users/maclark/Desktop/Mapping Food Consumption/Raw data/Religion/wrd_province_religion_qg0/",
                            layer = 'wrd_province_religion_qg0')

# Extracting province IDs
religion.tif.province_id <- rasterize(religion.shape,
                                      template.raster,
                                      field = religion.shape$ProvinceID)



###
# Manually creating meta data
# Will be used to later convert data into a usable format for regressions
# Extracting data manually
religion.data <- data.frame(Province_ID = religion.shape$ProvinceID,
                            Primary_Religion = religion.shape$HighRelig,
                            Percent_Primary = religion.shape$ReligCatP,
                            Percent_Christian = religion.shape$ChrCatP) # I really have no idea why there's a data set for christians but not any of the other major religions

# Getting unique matches with primary religion, percent_primary, and percent_christian
religion.unique <- unique(dplyr::select(religion.data,
                                        Primary_Religion,
                                        Percent_Primary,
                                        Percent_Christian))

# Updating format of numeric columns
# Specifically converting from percents to numbers
# First for percent primary
religion.unique$Percent_Primary_lower <- as.numeric(gsub("-.*",
                                                         "",
                                                         religion.unique$Percent_Primary))
religion.unique$Percent_Primary_upper <- gsub(".*-",
                                              "",
                                              religion.unique$Percent_Primary)
religion.unique$Percent_Primary_upper <- as.numeric(gsub("%",
                                                         "",
                                                         religion.unique$Percent_Primary_upper))
# Now taking average of the range
religion.unique$Percent_Primary_Numeric <-
  (religion.unique$Percent_Primary_upper + 
  religion.unique$Percent_Primary_lower) / 2

# Now repeating for percent christian
religion.unique$Percent_Christian_lower <- as.numeric(gsub("-.*",
                                                         "",
                                                         religion.unique$Percent_Christian))
religion.unique$Percent_Christian_upper <- gsub(".*-",
                                              "",
                                              religion.unique$Percent_Christian)
religion.unique$Percent_Christian_upper <- as.numeric(gsub("%",
                                                         "",
                                                         religion.unique$Percent_Christian_upper))
# Now taking average of the range
religion.unique$Percent_Christian_Numeric <-
  (religion.unique$Percent_Christian_upper + 
     religion.unique$Percent_Christian_lower) / 2

# Merging this smaller data frame back into the data frame that contains a single row for each province
religion.data <- left_join(religion.data,
                           dplyr::select(religion.unique,
                                         Primary_Religion,
                                         Percent_Primary,
                                         Percent_Christian,
                                         Percent_Primary_Numeric,
                                         Percent_Christian_Numeric))

### 
# We now have a data frame that contains every combination of primary religion, 
# percent of population that practices primary religion,
# and percent of population that is christian
# 
# We also have a data frame that contains the same data for each province
# 
# Now need to use these data sets to create maps that contain numeric values for 
# show primary religion, and percent of population that practices primary religion 
# and christianity
###

### 
# Creating unique numeric index of the religions
# Need to do this because R doesn't allow for categorical variables in rasters
religion.index <- data.frame(Religion = unique(religion.data$Primary_Religion),
                             Religion_Number = 1:length(unique(religion.data$Primary_Religion)))


###
# Now creating new maps manually
# This takes a few steps
# 1) Convert the province id data into a data frame
# 2) Merge data the religion data to the province id data
# 2.etc) e.g. add the numeric religion identifiers into the data frame
# 3) Convert the data frames back into maps 
# 3 notes) This is optional because we'll eventually convert all the rasters into a
# data frame for the regressions. Do this if you want to see the maps, don't do this
# if you want to save time. This section is commented out on this script
# Doing step 3 is pretty quick if the rasters are coarse, but could take a few mins
# At finer resolutions


###
# 1) Convert the province id data into a data frame
# This one is quick (at least with coarse rasters)
religion.data.frame <- data.frame(Province_ID = getValues(religion.tif.province_id))

###
# 2) Merging in the religion data
# This is also quick with coarse rasters
# First adding the primary, percent primary, and christian population data
religion.data.frame <- left_join(religion.data.frame,
                                 religion.data)

# Second adding in the numeric religion id values
religion.data.frame <- left_join(religion.data.frame,
                                 dplyr::rename(religion.index,
                                               Primary_Religion = Religion))


###
# 3) Recreating maps

# First for the religion numeric map
# Creating map
religion.numeric.map <- raster(matrix(religion.data.frame$Religion_Number,
                                      nrow = template.raster@nrows,
                                      ncol = template.raster@ncols,
                                      byrow = TRUE),
                               crs = crs.out)
# Updating extent
extent(religion.numeric.map) <- extent(template.raster)

# Second for the percent of population that practices the primary religion
percent_primary.map <- raster(matrix(religion.data.frame$Percent_Primary_Numeric,
                                     nrow = template.raster@nrows,
                                     ncol = template.raster@ncols,
                                     byrow = TRUE),
                              crs = crs.out)
# Updating extent
extent(percent_primary.map) <- extent(template.raster)

# Third for percent of population that is christian
percent_christian.map <- raster(matrix(religion.data.frame$Percent_Christian_Numeric,
                                       nrow = template.raster@nrows,
                                       ncol = template.raster@ncols,
                                       byrow = TRUE),
                                crs = crs.out)
# Updating extent
extent(percent_christian.map) <- extent(template.raster)


#####
# Creating a map that contains ISO3 Numeric values by country
# Data is from: http://thematicmapping.org/downloads/world_borders.php
# Note that this map currently doesn't separate Sudan and South Sudan
#####

# importing shape file with world borders
country.map <- 
  raster("/Users/maclark/Desktop/WorldMap_With_Lesotho/Country.Raster.Lesotho.Mollweide.tif")

# And converting from a shape file into a raster
country.map <-
  projectRaster(country.map,
                template.raster,
                method = 'ngb')

#####
# Getting coordinates of cells
#####
temp.raster.coords <- template.raster
temp.raster.coords[is.na(template.raster)] <- 0
cell.coords <- rasterToPoints(temp.raster.coords)
cell.coords <- as.data.frame(cell.coords)

###
# Soil suitability overall
soil.suit <- raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/Soil.Suitability.Mollweide.1.5km.tif")
crs(soil.suit) <- crs(raster("/Users/maclark/Desktop/Mike's Files/Global Mollweide 1.5km/MollweideCountryID_1.5km.tif"))
soil.suit <- projectRaster(soil.suit,
                           template.raster,
                           method = 'bilinear',
                           na.rm = TRUE)

###
# Soil suitability for each crop type
suit.files <- list.files(path = "/Users/maclark/Desktop/Soil Suitability Rasters",
                         full.names = TRUE)
rast.names <- gsub("/Users/maclark/Desktop/Soil Suitability Rasters/","",suit.files)

for(i in 1:length(suit.files)) {
  tmp.suit.rast <- raster(suit.files[i])
  tmp.suit.raster <- projectRaster(tmp.suit.rast,
                                   template.raster,
                                   method = 'bilinear',
                                   na.rm = TRUE)
  tmp.suit.raster <- squish(tmp.suit.raster,c(0,9))
  assign(rast.names[i], tmp.suit.raster)
}

###
# Adding other data
# Can merge these in with the data set below
# Single value by country
# % employed in ag - from world bank
# % GDP from ag - from world bank
# Continent/world region (can do this using countrycode on a data frame)
# % income spent on food - this is missing a huge amount of countries
# avg food demand by country 
# GINI coefficient - from world bank
# Quality of transit network
# Social stability



#####
# Have all the data we need
# Now making a big data frame that contains all of this data
# Data frame will be used to run the regressions
#####

# Making data frame
regression.dat <- data.frame(ISO3N = getValues(country.map),
                             Primary_Religion = religion.data.frame$Primary_Religion,
                             Percent_Primary_Religion = religion.data.frame$Percent_Primary_Numeric,
                             Percent_Christian = religion.data.frame$Percent_Christian_Numeric,
                             GDP_PPP_2015 = getValues(gdp.2015.raster),
                             GDP_PPP_2000 = getValues(gdp.2000.raster),
                             GDP_PPP_1990 = getValues(gdp.1990.raster),
                             Population_2015 = getValues(pop.2015.raster),
                             Population_2000 = getValues(pop.2000.raster),
                             Population_1990 = getValues(pop.1990.raster),
                             Travel_time_to_cities = getValues(raw.travel),
                             Soil_suit_overall = getValues(soil.suit),
                             Soil_suit_cereals = getValues(Cereals_Suitability.tif),
                             Soil_suit_fruits = getValues(Fruits_Suitability.tif),
                             Soil_suit_veges = getValues(Veges_Suitability.tif),
                             Soil_suit_nuts_seeds = getValues(`Nuts and Seeds_Suitability.tif`),
                             Soil_suit_pulses = getValues(Pulses_Suitability.tif),
                             Soil_suit_roots = getValues(Roots_Suitability.tif),
                             Soil_suit_chocolate = getValues(Chocolate_Suitability.tif),
                             Soil_suit_coffee = getValues(Coffee_Suitability.tif),
                             Cell_Coord_x = cell.coords$x,
                             Cell_Coord_y = cell.coords$y)

# Adding GDP per capita
regression.dat$PerCapita_GDP_2015 <- regression.dat$GDP_PPP_2015 / regression.dat$Population_2015
regression.dat$PerCapita_GDP_2000 <- regression.dat$GDP_PPP_2000 / regression.dat$Population_2000
regression.dat$PerCapita_GDP_1990 <- regression.dat$GDP_PPP_1990 / regression.dat$Population_1990

# Adjusting for cells where pop == 0
regression.dat$PerCapita_GDP_2015[regression.dat$Population_2015 == 0] <- 0
regression.dat$PerCapita_GDP_2000[regression.dat$Population_2000 == 0] <- 0
regression.dat$PerCapita_GDP_1990[regression.dat$Population_1990 == 0] <- 0


# Adding country ISO3C values
regression.dat$ISO3C <- countrycode(origin = 'iso3n',
                                    destination = 'iso3c',
                                    regression.dat$ISO3N)

regression.dat$ISO3C[regression.dat$ISO3N == 736] <- 'SDN'

# Adding continents and world regions to the data frame
# These will give a warning about codes not being matched
# The codes are 530 and 736
# These correspond with Netherlands Antilles and Sudan, respectively
regression.dat$Continent <- countrycode(origin = 'iso3c',
                                        destination = 'continent',
                                        regression.dat$ISO3C)

regression.dat$World_Region <- countrycode(origin = 'iso3c',
                                           destination = 'region',
                                           regression.dat$ISO3C)

write.csv(regression.dat,
          "/Users/maclark/Desktop/Mapping Food Consumption/Managed Regression Data/Dat 24June2019.csv",
          row.names = FALSE)

#####
# Break in code
# Importing regression dat that has DHS data incorporated
#####

###
# Importing dat
regression.dat1 <-
  read_csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/DHS data aggregated food groups with cells with survey year 26June2019.csv")
regression.dat2 <-
  read_csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed Regression Data/DHS data with cells 24June2019.csv")

###
# Getting median survey year by country
survey_year_country <-
  regression.dat1[!is.na(regression.dat1$Survey_Year),]
survey_year_country <- survey_year_country %>% mutate(Survey_Year = as.numeric(Survey_Year)) %>% group_by(ISO3C) %>% summarise(Survey_Year = median(Survey_Year))



###
# Adding in political stability data
gini.dat <- read.csv('/Users/maclark/Desktop/Mapping Food Consumption/GINI by Country/API_SI.POV.GINI_DS2_en_csv_v2_10576671.csv',
                     skip = 4)
gini.dat$gini.min <- NA
gini.dat$gini.max <- NA
gini.dat$year.min = NA
gini.dat$year.max = NA
gini.dat$Survey_Year = NA

gini.dat <- gini.dat[gini.dat$Country.Code %in% survey_year_country$ISO3C,]

# Getting last year of data for each country
for(i in 1:nrow(gini.dat)) {
  gini.dat$Survey_Year[i] = survey_year_country$Survey_Year[survey_year_country$ISO3C == gini.dat$Country.Code[i]][1]
  if(!is.na(gini.dat[i,grep(survey_year_country$Survey_Year[survey_year_country$ISO3C == gini.dat$Country.Code[i]][1],names(gini.dat))])) {
    k = survey_year_country$Survey_Year[survey_year_country$ISO3C == gini.dat$Country.Code[i]][1]
    gini.dat[i,c('gini.min','gini.max')] <- gini.dat[i,grep(k, names(gini.dat))]
    gini.dat[i,c('year.min','year.max')] <- k
  }
  k = survey_year_country$Survey_Year[survey_year_country$ISO3C == gini.dat$Country.Code[i]][1]
  while(is.na(gini.dat$gini.min[i]) & k > 1961) {
    
    gini.dat$gini.min[i] <- gini.dat[i,grep(k,names(gini.dat))]
    gini.dat$year.min[i] <- k
    k = k-1
  }
  k = survey_year_country$Survey_Year[survey_year_country$ISO3C == gini.dat$Country.Code[i]][1]
  while(is.na(gini.dat$gini.max[i]) & k < 2019) {
    
    gini.dat$gini.max[i] <- gini.dat[i,grep(k,names(gini.dat))]
    gini.dat$year.max[i] <- k
    k = k+1
  }
}

# Taking average
gini.dat$avg_gini <-
  gini.dat$gini.min * (gini.dat$year.max - gini.dat$Survey_Year)/(gini.dat$year.max - gini.dat$year.min) +
  gini.dat$gini.max * (gini.dat$Survey_Year - gini.dat$year.min)/(gini.dat$year.max - gini.dat$year.min)

gini.dat$avg_gini[is.na(gini.dat$avg_gini)] <- gini.dat$gini.min[is.na(gini.dat$avg_gini)]

###
# Now repeating for the other countries
gini.dat.others <- read.csv('/Users/maclark/Desktop/Mapping Food Consumption/GINI by Country/API_SI.POV.GINI_DS2_en_csv_v2_10576671.csv',
                     skip = 4)
gini.dat.others <- gini.dat.others[!(gini.dat.others$Country.Code %in% gini.dat$Country.Code),]
gini.dat.others$gini.min <- NA
gini.dat.others$gini.max <- NA
gini.dat.others$year.min = NA
gini.dat.others$year.max = NA
gini.dat.others$Survey_Year = NA
gini.dat.others$avg_gini = NA
for(i in 1:nrow(gini.dat.others)) {
  k = 2018
  while(is.na(gini.dat.others$avg_gini[i]) & k >= 1961) {
    gini.dat.others$avg_gini[i] <-
      gini.dat.others[i,grep(k,names(gini.dat.others))]
    gini.dat.others$Survey_Year[i] <-
      k
    k = k-1
  }
}

gini.dat <- rbind(gini.dat,gini.dat.others)
gini.dat <- unique(gini.dat)

###
# Repeating with GDP dat
gdp.ag <- read.csv('/Users/maclark/Desktop/Mapping Food Consumption/Percent GDP from Ag/API_NV.AGR.TOTL.ZS_DS2_en_csv_v2_10580827.csv',
                     skip = 4)


gdp.ag$gdp.min <- NA
gdp.ag$gdp.max <- NA
gdp.ag$year.min = NA
gdp.ag$year.max = NA
gdp.ag$Survey_Year = NA

gdp.ag <- gdp.ag[gdp.ag$Country.Code %in% survey_year_country$ISO3C,]

# Getting last year of data for each country
for(i in 1:nrow(gdp.ag)) {
  gdp.ag$Survey_Year[i] = survey_year_country$Survey_Year[survey_year_country$ISO3C == gdp.ag$Country.Code[i]][1]
  if(!is.na(gdp.ag[i,grep(survey_year_country$Survey_Year[survey_year_country$ISO3C == gdp.ag$Country.Code[i]][1],names(gdp.ag))])) {
    k = survey_year_country$Survey_Year[survey_year_country$ISO3C == gdp.ag$Country.Code[i]][1]
    gdp.ag[i,c('gdp.min','gdp.max')] <- gdp.ag[i,grep(k, names(gdp.ag))]
    gdp.ag[i,c('year.min','year.max')] <- k
  }
  k = survey_year_country$Survey_Year[survey_year_country$ISO3C == gdp.ag$Country.Code[i]][1]
  while(is.na(gdp.ag$gdp.min[i]) & k > 1961) {
    
    gdp.ag$gdp.min[i] <- gdp.ag[i,grep(k,names(gdp.ag))]
    gdp.ag$year.min[i] <- k
    k = k-1
  }
  k = survey_year_country$Survey_Year[survey_year_country$ISO3C == gdp.ag$Country.Code[i]][1]
  while(is.na(gdp.ag$gdp.max[i]) & k < 2019) {
    
    gdp.ag$gdp.max[i] <- gdp.ag[i,grep(k,names(gdp.ag))]
    gdp.ag$year.max[i] <- k
    k = k+1
  }
}

# Taking average
gdp.ag$avg_gdp <-
  gdp.ag$gdp.min * (gdp.ag$year.max - gdp.ag$Survey_Year)/(gdp.ag$year.max - gdp.ag$year.min) +
  gdp.ag$gdp.max * (gdp.ag$Survey_Year - gdp.ag$year.min)/(gdp.ag$year.max - gdp.ag$year.min)

gdp.ag$avg_gdp[is.na(gdp.ag$avg_gdp)] <- gdp.ag$gdp.min[is.na(gdp.ag$avg_gdp)]

# Now for other countries
gdp.ag.others <- read.csv('/Users/maclark/Desktop/Mapping Food Consumption/Percent GDP from Ag/API_NV.AGR.TOTL.ZS_DS2_en_csv_v2_10580827.csv',
                   skip = 4)


gdp.ag.others$gdp.min <- NA
gdp.ag.others$gdp.max <- NA
gdp.ag.others$year.min = NA
gdp.ag.others$year.max = NA
gdp.ag.others$Survey_Year = NA
gdp.ag.others$avg_gdp = NA

gdp.ag.others <- gdp.ag.others[!(gdp.ag.others$Country.Code %in% survey_year_country$ISO3C),]
for(i in 1:nrow(gdp.ag.others)) {
  k = 2018
  while(is.na(gdp.ag.others$avg_gdp[i]) & k >= 1961) {
    gdp.ag.others$avg_gdp[i] <-
      gdp.ag.others[i,grep(k,names(gdp.ag.others))]
    gdp.ag.others$Survey_Year[i] <-
      k
    k = k-1
  }
}

gdp.dat <- rbind(gdp.ag,gdp.ag.others)
gdp.dat <- unique(gdp.dat)


###
# Repeating for labor in ag
labor.ag <- read.csv('/Users/maclark/Desktop/Mapping Food Consumption/Percent Employed by Ag/API_SL.AGR.EMPL.ZS_DS2_en_csv_v2_10576798.csv',
                     skip = 4)


labor.ag$labor.min <- NA
labor.ag$labor.max <- NA
labor.ag$year.min = NA
labor.ag$year.max = NA
labor.ag$Survey_Year = NA

labor.ag <- labor.ag[labor.ag$Country.Code %in% survey_year_country$ISO3C,]

# Getting last year of data for each country
for(i in 1:nrow(labor.ag)) {
  labor.ag$Survey_Year[i] = survey_year_country$Survey_Year[survey_year_country$ISO3C == labor.ag$Country.Code[i]][1]
  if(!is.na(labor.ag[i,grep(survey_year_country$Survey_Year[survey_year_country$ISO3C == labor.ag$Country.Code[i]][1],names(labor.ag))])) {
    k = survey_year_country$Survey_Year[survey_year_country$ISO3C == labor.ag$Country.Code[i]][1]
    labor.ag[i,c('labor.min','labor.max')] <- labor.ag[i,grep(k, names(labor.ag))]
    labor.ag[i,c('year.min','year.max')] <- k
  }
  k = survey_year_country$Survey_Year[survey_year_country$ISO3C == labor.ag$Country.Code[i]][1]
  while(is.na(labor.ag$labor.min[i]) & k > 1961) {
    
    labor.ag$labor.min[i] <- labor.ag[i,grep(k,names(labor.ag))]
    labor.ag$year.min[i] <- k
    k = k-1
  }
  k = survey_year_country$Survey_Year[survey_year_country$ISO3C == labor.ag$Country.Code[i]][1]
  while(is.na(labor.ag$labor.max[i]) & k < 2019) {
    
    labor.ag$labor.max[i] <- labor.ag[i,grep(k,names(labor.ag))]
    labor.ag$year.max[i] <- k
    k = k+1
  }
}

# Taking average
labor.ag$avg_labor <-
  labor.ag$labor.min * (labor.ag$year.max - labor.ag$Survey_Year)/(labor.ag$year.max - labor.ag$year.min) +
  labor.ag$labor.max * (labor.ag$Survey_Year - labor.ag$year.min)/(labor.ag$year.max - labor.ag$year.min)

labor.ag$avg_labor[is.na(labor.ag$avg_labor)] <- labor.ag$labor.min[is.na(labor.ag$avg_labor)]

# Other countries
labor.ag.others <- read.csv('/Users/maclark/Desktop/Mapping Food Consumption/Percent Employed by Ag/API_SL.AGR.EMPL.ZS_DS2_en_csv_v2_10576798.csv',
                     skip = 4)


labor.ag.others$labor.min <- NA
labor.ag.others$labor.max <- NA
labor.ag.others$year.min = NA
labor.ag.others$year.max = NA
labor.ag.others$Survey_Year = NA
labor.ag.others$avg_labor = NA

labor.ag.others <- labor.ag.others[!(labor.ag.others$Country.Code %in% survey_year_country$ISO3C),]
labor.ag.others <- labor.ag.others[!(labor.ag.others$Country.Code %in% survey_year_country$ISO3C),]
for(i in 1:nrow(labor.ag.others)) {
  k = 2018
  while(is.na(labor.ag.others$avg_labor[i]) & k >= 1961) {
    labor.ag.others$avg_labor[i] <-
      labor.ag.others[i,grep(k,names(labor.ag.others))]
    labor.ag.others$Survey_Year[i] <-
      k
    k = k-1
  }
}

labor.ag <- rbind(labor.ag,labor.ag.others)
labor.ag <- unique(labor.ag)

###
# Merging in
country.dat <-
  full_join(dplyr::select(gini.dat,
                          Country.Name,
                          Country.Code,
                          gini = avg_gini),
            dplyr::select(gdp.ag,
                          Country.Name,
                          Country.Code,
                          prop.gdp.ag = avg_gdp))

country.dat <-
  full_join(country.dat,
            dplyr::select(labor.ag,
                          Country.Name,
                          Country.Code,
                          prop.labor.ag = avg_labor))


###
# Getting FAO data for consumption of major food groups
fao.dat <- read.csv('/Users/maclark/Desktop/Mapping Food Consumption/FAOSTAT_data_6-24-2019.csv',
                    stringsAsFactors = FALSE)
fao.dat$Item[fao.dat$Item %in% c('Treenuts','Oilcrops')] <- 'Nuts_and_Seeds'
fao.dat.sum <- fao.dat %>% filter(!is.na(Value)) %>% group_by(Country.Code,Item,Year) %>% summarise(cons = sum(Value))

final.fao.dat <- c()
# Looping through countries to match survey year
for(i in 1:nrow(survey_year_country[!is.na(survey_year_country$ISO3C),])) {
  if(survey_year_country$Survey_Year[i] >= 2011) {
    tmp.fao.dat <- fao.dat.sum %>% filter(Year >= 2009) %>% filter(Country.Code == survey_year_country$ISO3C[i]) %>% group_by(Country.Code,Item) %>% summarise(cons_kg.yr = mean(cons))
  }
  if(survey_year_country$Survey_Year[i] < 2011) {
    tmp.fao.dat <- fao.dat.sum %>% filter(Year >= (survey_year_country$Survey_Year[i]-2)) %>% filter(Year <= (survey_year_country$Survey_Year[i]+2)) %>% filter(Country.Code == survey_year_country$ISO3C[i]) %>% group_by(Country.Code,Item) %>% summarise(cons_kg.yr = mean(cons))
  }
  final.fao.dat <- rbind(final.fao.dat,
                         tmp.fao.dat)
}


###
# Converting to wide format
# Reshaping long to wide
reshape.wide <- function(dat,
                         column.split,
                         column.data) {
  # getting list of unique items in the column to split
  split.list <- sort(unique(dat[,column.split]))
  
  # Looping through these
  for(i in 1:length(split.list)) {
    # Getting temporary data set with only the indicator
    dat.tmp <- dat[dat[,column.split] == split.list[[i]],]
    # Changing column name
    names(dat.tmp)[names(dat.tmp) == column.data] <- split.list[[i]]
    # Dropping column
    col.indicator <- which(names(dat.tmp) == column.split)
    dat.tmp <- dat.tmp[,-col.indicator]
    
    # And making master file
    if(i == 1) {
      master.tmp <- dat.tmp
    } else{
      master.tmp <- full_join(master.tmp,
                              dat.tmp)
    }
  }
  return(master.tmp)
}

fao.dat.wide <-
  reshape.wide(as.data.frame(final.fao.dat),
               column.split = 'Item',
               column.data = 'cons_kg.yr')

###
# Getting data for the out of sample countries
# Looping through countries to match survey year
fao.dat.sum <- fao.dat.sum[!(fao.dat.sum$Country.Code %in% survey_year_country$ISO3C),]
fao.dat.others <- fao.dat.sum %>% filter(Year >= 2009) %>% group_by(Country.Code, Item) %>% summarise(cons_kg.yr = mean(cons))
fao.dat.others.wide <- 
  reshape.wide(as.data.frame(fao.dat.others),
               column.split = 'Item',
               column.data = 'cons_kg.yr')

fao.dat.wide <- rbind(fao.dat.wide,fao.dat.others.wide)


###
# Merging in
regression.dat2 <- as.data.frame(dplyr::select(regression.dat1,Cell_Order,ISO3C))
regression.dat2$gini = NA
regression.dat2$prop.gdp.ag = NA
regression.dat2$prop.labor.ag = NA
for(i in 1:nrow(country.dat)) {
  regression.dat2$gini[regression.dat2$ISO3C == country.dat$Country.Code[i] & !is.na(regression.dat2$ISO3C)] <- country.dat$gini[i]
  regression.dat2$prop.gdp.ag[regression.dat2$ISO3C == country.dat$Country.Code[i] & !is.na(regression.dat2$ISO3C)] <- country.dat$prop.gdp.ag[i]
  regression.dat2$prop.labor.ag[regression.dat2$ISO3C == country.dat$Country.Code[i] & !is.na(regression.dat2$ISO3C)] <- country.dat$prop.labor.ag[i]
}

###
# Merging in fao dat
regression.dat2 <-
  left_join(regression.dat2,
            dplyr::rename(fao.dat.wide,
                          ISO3C = Country.Code))

regression.dat <-
  cbind(regression.dat1,
        regression.dat2[,3:ncol(regression.dat2)])

names(regression.dat)[34:44] <- paste("Cons_kg.yr",
                                      names(regression.dat)[34:44])

#####
# Writing the data frame
#####

write.csv(regression.dat,
          "/Users/maclark/Desktop/Mapping Food Consumption/Managed Regression Data/Global Regression Dat 27June2019 No DHS Data.csv",
          row.names = FALSE)


