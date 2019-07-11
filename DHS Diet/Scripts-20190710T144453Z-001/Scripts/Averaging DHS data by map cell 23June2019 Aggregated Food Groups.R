#####
# Merging cluster data with gps data
#####

###
# Installing libraries

library(plyr)
library(dplyr)
library(raster)

###
# Cluster dat
cluster.dat <- 
  read.csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed Regression Data/DHS_averaged_by_cluster_num_AggregateFoodGroups_23June2019.csv",
           stringsAsFactors = FALSE)

###
# GPS dat
gps.dat <- 
  read.csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/Stacked_GPS_Data.csv",
           stringsAsFactors = FALSE)
# Getting country
gps.dat$Country <-
  substr(gps.dat$Case_ID,
         start = 1,
         stop = 2)

###
# Merging
gps.cluster <-
  left_join(gps.dat,
            dplyr::rename(cluster.dat,
                          Cluster_Num = Cluster.Num))

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
write.csv(gps.cluster,
          "/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/DHS_Regression_Dat_AggregatedFoodGroups_Averaged_by_Cell 23June2019.csv",
          row.names = FALSE)

###
# And taking average by cell order
###

###
# Saving temporary file
cell.orders <-
  data.frame(Cell_Order = 1:9331200)

###
# Looping through columns to get avearges by cluster
for(i in 1:22) {
  ###
  # Getting select columns
  tmp.dat <-
    map.gps.cluster[,c(3,(9 + (i - 1) * 2), (10 + (i - 1) * 2))]
  # Removing NAs
  tmp.dat <-
    tmp.dat[!is.na(tmp.dat[,2]),]
  # Adding count variable
  tmp.dat$Count <-
    1
  
  # Multiplying by number of households in the cluster
  tmp.dat[,2] <-
    tmp.dat[,2] *
    tmp.dat[,3]
  
  # Summing by cluster number
  tmp.dat <-
    rowsum(tmp.dat,
           group = tmp.dat$Cell_Order,
           na.rm = TRUE)
  
  # Getting average by cluster
  tmp.dat[,2] <-
    tmp.dat[,2] /
    tmp.dat[,3]
  tmp.dat[,1] <-
    tmp.dat[,1] /
    tmp.dat$Count
  # And resetting cell order by dividing by count
  # Renaming column for merging
  names(tmp.dat)[2] <-
    "TMP"
  names(tmp.dat)[3] <-
    "TMP.2"
  
  ###
  # Merging based on our case ids
  cell.orders <-
    left_join(cell.orders,
              dplyr::select(tmp.dat,
                            Cell_Order,
                            TMP,
                            TMP.2))
  
  # Updating column name
  names(cell.orders)[ncol(cell.orders) - 1] <-
    names(map.gps.cluster)[(9 + (i - 1) * 2)]
  
  names(cell.orders)[ncol(cell.orders)] <-
    names(map.gps.cluster)[(10 + (i - 1) * 2)]
}

###
# Adding regression dat
# This is a huge file
regression.dat <-
  read_csv("/Users/maclark/Desktop/Mapping Food Consumption/Managed Regression Data/Dat 24June2019.csv")

cell.orders1 <-
  cbind(cell.orders,
        regression.dat)

###
# And writing data for all cells globally
write.csv(cell.orders1,
          "/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/DHS data aggregated food groups with cells 24June2019.csv",
          row.names = FALSE)

gps.cluster1 <-
  left_join(gps.cluster,
            dplyr::rename(tmp.df,
                          Lon = x,
                          Lat = y))

###
# And dropping rows that don't contain any data
# Doing this to reduce file size for the regressions
cell.orders2 <-
  cell.orders1[cell.orders1$Cell_Order %in% gps.cluster1$Cell_Order,]

###
# And saving
write.csv(cell.orders2,
          "/Users/maclark/Desktop/Mapping Food Consumption/Managed_DHS_Data/DHS_Data_By_Cell_with_Regression_Dat_AggregatedFoodGroups_24June2019.csv",
          row.names = FALSE)

###
# Number of fruit, because why not
fruit.raster <-
  raster(matrix(as.numeric(cell.orders1$Soil_suit_cereals),
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
