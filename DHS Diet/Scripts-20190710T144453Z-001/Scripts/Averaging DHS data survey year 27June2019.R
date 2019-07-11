#####
# Merging cluster data with gps data
#####

###
# Installing libraries

library(plyr)
library(dplyr)
library(raster)
library(readr)

rm(list = ls())

###
# Cluster dat
cluster.dat <- 
  read.csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/Survey_Year.csv",
           stringsAsFactors = FALSE)

###
# GPS dat
gps.dat <- 
  read_csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/Stacked_GPS_Data.csv")
# Getting country
gps.dat$Country <-
  substr(gps.dat$Case_ID,
         start = 1,
         stop = 2)

###
# Merging
gps.cluster <-
  left_join(gps.dat,
            dplyr::select(cluster.dat,
                          Year,
                          Country,
                          Cluster_Num = Cluster.number,
                          Survey_Year = Year))

#####
# Mapping for the hell of it
# Even though there's a fair bit of data that isn't matching
#####

###
# Making a template raster
template.raster <-
  raster(matrix(0, 
                ncol = 4320,
                nrow = 2160),
         crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
extent(template.raster) <-
  c(-180, 180,
    -90, 90)


###
# Making a data frame with lat and long coordinates
tmp.df <-
  rasterToPoints(template.raster)

tmp.df <-
  data.frame(tmp.df)

###
# Rounding so it merges with lat and long coordinates from the surveys
tmp.df$x <-
  round(tmp.df$x / .0833333) * .0833333

tmp.df$y <-
  round(tmp.df$y / .0833333) * .0833333

# And rounding again to two digits
tmp.df$x <-
  round(tmp.df$x, digits = 3)
tmp.df$y <-
  round(tmp.df$y, digits = 3)

# Setting order of cells
tmp.df$Cell_Order <-
  1:nrow(tmp.df)


tmp.df2 <- tmp.df

###
# And doing the same with lat and long coordinates in the surveys
gps.cluster$Lat <-
  round(gps.cluster$Lat / .0833333) * .0833333
gps.cluster$Lon <-
  round(gps.cluster$Lon / .0833333) * .0833333

# And rounding again to two digits
gps.cluster$Lat <-
  round(gps.cluster$Lat, digits = 3)
gps.cluster$Lon <-
  round(gps.cluster$Lon, digits = 3)

###
# Merging DHS data with gps data
map.gps.cluster <-
  left_join(dplyr::select(tmp.df,
                          Lon = x,
                          Lat = y,
                          Cell_Order),
            gps.cluster)

###
# Writing csv
# write.csv(gps.cluster,
#           "/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/Survey_Year_AveragedbyCell.csv",
#           row.names = FALSE)

###
# And taking average by cell order
###

###
# Saving temporary file
cell.orders <-
  data.frame(Cell_Order = 1:9331200)

###
# Getting average year by cell
tmp.df <- map.gps.cluster[!is.na(map.gps.cluster$Survey_Year),]

year.cell <-
  tmp.df %>% group_by(Lon,Lat,Cell_Order) %>% summarise(Survey_Year = median(Survey_Year))

year.cell.sum <- year.cell %>% group_by(Survey_Year) %>% summarise(count = n())



###
# Adding regression dat
# This is a huge file
regression.dat <-
  read_csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed Regression Data/Global Regression Dat 27June2019 No DHS Data.csv")

###
# Need to add lat and lon estimates
regression.dat <- left_join(regression.dat,
                        dplyr::select(tmp.df2,
                                      Lat = y,
                                      Lon = x,
                                      Cell_Order))

###
# Then merge in DHS data
# First by aggregate food groups
# Then by individual food groups
dhs.dat.ag <- read.csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/DHS_Regression_Dat_AggregatedFoodGroups_Averaged_by_Cell 23June2019.csv",
                       stringsAsFactors = FALSE)
dhs.dat.ag$ISO3C <- countrycode(dhs.dat.ag$Country,origin = 'iso2c',destination = 'iso3c')
dhs.dat.ag$ISO3C[dhs.dat.ag$ISO3C == 'IA'] <- 'IND'
dhs.dat.ag$ISO3C[dhs.dat.ag$ISO3C == 'BU'] <- 'BDI'
dhs.dat.ag$ISO3C[dhs.dat.ag$ISO3C == 'DR'] <- 'DOM'
dhs.dat.ag$ISO3C[dhs.dat.ag$ISO3C == 'MB'] <- 'MDA'
dhs.dat.ag$ISO3C[dhs.dat.ag$ISO3C == 'NM'] <- 'NAM'

###
# First trying to merge entire data set
regression.dat1 <- 
  left_join(regression.dat,
            dhs.dat.ag[,c(4,5,8:25)])

###
# Few logic checks to make sure these are merging appropriately
# Looks like we're missing a few
sum(!is.na(dhs.dat.ag$vegetables))
sum(!is.na(dhs.dat.ag$fruits))
sum(!is.na(dhs.dat.ag$meat_fish))

sum(!is.na(regression.dat1$vegetables))
sum(!is.na(regression.dat1$fruits))
sum(!is.na(regression.dat1$meat_fish))

###
# Limiting to only cells that are within DHS survey countries
regression.dat1 <- regression.dat %>% filter(!is.na(ISO3N)) %>% filter(ISO3C %in% unique(dhs.dat.ag$ISO3C))

regression.dat1 <- 
  left_join(regression.dat,
            dhs.dat.ag[,c(4,5,8:25)])
  

cell.orders1 <-
  cbind(cell.orders,
        regression.dat)

# cell.orders1 <- left_join(regression.dat,
#                           dplyr::select(year.cell,
#                                         Cell_Order,
#                                         Survey_Year))

cell.orders1 <- 
  dplyr::select(cell.orders1,
                -Lat)
cell.orders1 <- 
  dplyr::select(cell.orders1,
                -Lon)

###
# And writing data for all cells globally
# write.csv(cell.orders1,
#           "/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/DHS data aggregated food groups with cells with survey year 26June2019.csv",
#           row.names = FALSE)

gps.cluster1 <-
  left_join(gps.cluster,
            dplyr::select(tmp.df,
                          Lat,
                          Lon,
                          Cell_Order))

###
# And dropping rows that don't contain any data
# Doing this to reduce file size for the regressions
cell.orders2 <-
  cell.orders1[cell.orders1$Cell_Order %in% gps.cluster1$Cell_Order,]

###
# And saving
write.csv(cell.orders2,
          "/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/DHS_Data_By_Cell_with_Regression_Dat_AggregatedFoodGroups_SurveyYear_26June2019.csv",
          row.names = FALSE)

###
# Number of fruit, because why not
fruit.raster <-
  raster(matrix(as.numeric(cell.orders1$Survey_Year),
                nrow = template.raster@nrows,
                ncol = template.raster@ncols,
                byrow = TRUE),
         crs = crs(template.raster))
extent(fruit.raster) <-
  extent(template.raster)

fruit.raster <-
  raster(matrix(cell.orders$Gave.child.tea.or.coffee,
                nrow = template.raster@nrows,
                ncol = template.raster@ncols,
                byrow = TRUE),
         crs = crs(template.raster))
extent(fruit.raster) <-
  extent(template.raster)
