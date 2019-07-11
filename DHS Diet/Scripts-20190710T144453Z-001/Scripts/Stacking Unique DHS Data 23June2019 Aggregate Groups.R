#####
# Looking at similarity of household IDs across the data sets
#####
rm(list = ls())


###
# Importing files
br.dat <-
  read.csv("/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/BR_Managed.csv",
           stringsAsFactors = FALSE)

cr.dat <- 
  read.csv("/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/CR_Managed.csv",
           stringsAsFactors = FALSE)

ir.dat <- 
  read.csv("/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/IR_Managed.csv",
           stringsAsFactors = FALSE)
kr.dat <- 
  read.csv("/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/KR_Managed.csv",
           stringsAsFactors = FALSE)
mr.dat <- 
  read.csv("/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/MR_Managed.csv",
           stringsAsFactors = FALSE)
  
###
# Unique case ids
br.ids <-
  unique(br.dat$Case.Identification)

tmp.ids <-
  cr.dat$Case.Identification[cr.dat$Case.Identification %in% br.ids]

tmp.ids <-
  unique(tmp.ids)

cr.dat.tmp <-
  cr.dat[cr.dat$Case.Identification %in% tmp.ids[1:10],]

br.dat.tmp <-
  br.dat[br.dat$Case.Identification %in% tmp.ids[1:10],]

###
# Limiting to household IDs that are not in other files
stacked.dat <-
  rbind(br.dat,
        cr.dat,
        ir.dat,
        kr.dat,
        mr.dat)

stacked.dat <-
  unique(stacked.dat)

nrow(stacked.dat)

# 
# 
# unique.ids <-
#   unique(stacked.dat$Case.Identification)
# 
# # ir.dat <-
# #   dplyr::rename(ir.dat,
# #                 Case.Identification = Case.Identification)
# 
# stacked.dat <-
#   rbind(stacked.dat,
#         ir.dat[!(ir.dat$Case.Identification %in% unique.ids),])
# 
# unique.ids <-
#   unique(stacked.dat$Case.Identification)
# # 
# # kr.dat <-
# #   dplyr::rename(kr.dat,
# #                 Case.Identification = Case.Identification)
# 
# stacked.dat <-
#   rbind(stacked.dat,
#         kr.dat[!(kr.dat$Case.Identification %in% unique.ids),])
# 
# unique.ids <-
#   unique(stacked.dat$Case.Identification)
# 
# # mr.dat <-
# #   dplyr::rename(mr.dat,
# #                 Case.Identification = Case.Identification)
# 
# stacked.dat <-
#   rbind(stacked.dat,
#         mr.dat[!(mr.dat$Case.Identification %in% unique.ids),])
# 
# stacked.dat <-
#   unique(stacked.dat)
# 
# ###
# # Writing file
# write.csv(stacked.dat,
#           "/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/Stacked_Data.csv",
#           row.names = FALSE)

#####
# Converting character strings to numbers
#####
# sort(unique(stacked.dat$Gave.child.plain.water))
# sort(unique(stacked.dat$NA...Gave.child.sugar.water))
# sort(unique(stacked.dat$Gave.child.juice))
# sort(unique(stacked.dat$Gave.child.baby.formula))

cols.character <-
  as.numeric(which(sapply(stacked.dat,
             is.numeric) == FALSE))
library(profvis)

for(i in 3:length(cols.character)) {
  print(names(stacked.dat)[cols.character[i]])
  print(sort(unique(stacked.dat[,cols.character[i]])))
  
  pause(0.5)
}

stacked.dat.tmp <-
  stacked.dat

for(i in 3:length(cols.character)) {
  # Saving vector
  tmp.vector <-
    stacked.dat.tmp[,cols.character[i]]
  
  # Updating Yes and Nos to 1s and 0s, respectively
  tmp.vector[tmp.vector == 'Yes'] <- 1
  tmp.vector[tmp.vector == 'No'] <- 0
  
  stacked.dat.tmp[,cols.character[i]] <-
    tmp.vector
}

# And converting columns to numeric format
for(i in 6:75) {
  stacked.dat.tmp[,i] <-
    as.numeric(stacked.dat.tmp[,i])
}

gen.dat <-
  stacked.dat.tmp[,c(2:5,76)]

# and dropping columns that contain no information
for(i in 6:75) {
  tmp.vector <- stacked.dat.tmp[,i]
  if(length(tmp.vector[is.na(tmp.vector)]) != length(tmp.vector)) {
    gen.dat$tmp <-
      tmp.vector
    names(gen.dat)[ncol(gen.dat)] <-
      names(stacked.dat.tmp)[i]
  }
}

#####
# And now getting averages by cluster number
#####

library(plyr)
library(dplyr)

###
# Saving data frame
gen.dat1 <-
  gen.dat

###
# Getting unique case id and merging in
case.ids <-
  data.frame(Case.Identification = gen.dat1$Case.Identification,
             Country = gen.dat1$Country,
             Cluster.Num = gen.dat1$Cluster.number)

case.ids <-
  unique(case.ids)

case.ids$Our_ID <-
  1:nrow(case.ids)

###
# Merging our case ids back in
gen.dat1 <-
  left_join(gen.dat1,
            dplyr::select(case.ids,
                          Case.Identification,
                          Our_ID))


###
# Looping through food groups
food.groups <- c('fruits','vegetables','dairy','eggs','meat_fish','nuts_pulses','sweets','grains_tubers','animal_source_foods','fruits_vegetables')
fruit.list <- c('aprictos..peaches','any.other.fruits')
vegetable.list <- c('gave.child.any.dark.green.leafy.vegetables')
dairy.list <- c('tinned..powdered.or.fresh.milk','cheese..yogurt..other.milk.products','gave.child.yogurt')
egg.list <- c('gave.child.eggs')
meat_fish.list <- c('gave.child.liver','gave.child.fish')
nuts_pulses.list <- c('food.made.from.beans..peas')
sweets.list <- c('chocolates..sweets..candies..pastries')
grains_tubers.list <- c('gave.child.bread','gave.child.potatoes')
animal_source_foods.list <- c(dairy.list,egg.list,meat_fish.list)
fruits_vegetables.list <- c(fruit.list,vegetable.list)
food.groups.list <- list(fruit.list,vegetable.list,dairy.list,egg.list,meat_fish.list,nuts_pulses.list,sweets.list,grains_tubers.list,animal_source_foods.list,fruits_vegetables.list)
names(food.groups.list) <- food.groups

gen.dat2 <-
  dplyr::select(gen.dat1,
                Case.Identification,
                Country.code.and.phase,
                Cluster.number,
                Household.number,
                Country,
                Our_ID)

gen.dat3 <- gen.dat1
gen.dat1 <- gen.dat3

#
for(i in 1:10) {
  col.index <- c()
  for(j in 1:length(food.groups.list[[i]])) {
    col.index <- c(col.index,
                  grep(food.groups.list[[i]][j],names(gen.dat1),ignore.case=TRUE))
  }
  if(length(col.index) > 1) {
    # Getting sum across the food groups
    tmpsum <- rowSums(gen.dat1[,c(col.index)],na.rm=TRUE)
    # tmpsum[tmpsum>1]<-1
    # Reinserting NAs
    for(k in 1:length(col.index)) {
      gen.dat1[is.na(gen.dat1[,col.index[k]]),col.index[k]] <- 11
    }
    tmpsumna <- rowSums(gen.dat1[,c(col.index)],na.rm=TRUE)
    tmpsum[tmpsumna == c(length(col.index)*11)] <- NA
    for(k in 1:length(col.index)) {
      gen.dat1[gen.dat1[,col.index[k]] == 11,col.index[k]] <- NA
    }
    
    gen.dat2$tmp <- tmpsum
    names(gen.dat2)[ncol(gen.dat2)] <- names(food.groups.list)[i]
  }
  if(length(col.index) == 1) {
    gen.dat2$tmp <- gen.dat1[,col.index]
    names(gen.dat2)[ncol(gen.dat2)] <- names(food.groups.list)[i]
  }
}

###
# And looping through columns to get summary by household
for(i in 7:15) {
  ###
  # Getting select columns
  tmp.dat <-
    gen.dat2[,c(6,i)]
  # Removing NAs
  tmp.dat <-
    tmp.dat[!is.na(tmp.dat[,2]),]
  # Adding a count variable
  tmp.dat$Count = 1
  # Summing by household
  tmp.dat <-
    rowsum(tmp.dat,
           group = tmp.dat$Our_ID,
           na.rm = TRUE)
  # Getting average by household
  tmp.dat[,2] <-
    tmp.dat[,2] /
    tmp.dat$Count
  tmp.dat[,1] <-
    tmp.dat[,1] /
    tmp.dat$Count
  # Renaming column for merging
  names(tmp.dat)[2] <-
    "TMP"
  
  ###
  # Merging based on our case ids
  case.ids <-
    left_join(case.ids,
              dplyr::select(tmp.dat,
                            Our_ID,
                            TMP))
  
  # Updating column name
  names(case.ids)[ncol(case.ids)] <-
    names(gen.dat2)[i]
}


###
# Writing csv
write.csv(case.ids,
          "/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/DHS_Averaged_by_Household_AggregateFoodGroups_23June2019.csv",
          row.names = FALSE)


###
# And repeating the same for cluster numbers
cluster.nums <-
  data.frame(Cluster.Num = case.ids$Cluster.Num,
             Country = case.ids$Country)

cluster.nums <- 
  unique(cluster.nums)

cluster.nums$Our.Cluster.Num <-
  1:nrow(cluster.nums)

###
# Merging back in
case.ids <-
  left_join(case.ids,
            cluster.nums)

###
# Saving temporary file
gen.dat1 <-
  case.ids

###
# Looping through columns to get avearges by cluster
for(i in 5:13) {
  ###
  # Getting select columns
  tmp.dat <-
    gen.dat1[,c(14,i)]
  # Removing NAs
  tmp.dat <-
    tmp.dat[!is.na(tmp.dat[,2]),]
  # Adding a count variable
  tmp.dat$Count = 1
  # Summing by cluster number
  tmp.dat <-
    rowsum(tmp.dat,
           group = tmp.dat$Our.Cluster.Num,
           na.rm = TRUE)
  # Getting average by household
  tmp.dat[,2] <-
    tmp.dat[,2] /
    tmp.dat$Count
  tmp.dat[,1] <-
    tmp.dat[,1] /
    tmp.dat$Count
  # Renaming column for merging
  names(tmp.dat)[2] <-
    "TMP"
  
  ###
  # Merging based on our case ids
  cluster.nums <-
    left_join(cluster.nums,
              dplyr::select(tmp.dat,
                            Our.Cluster.Num,
                            TMP,
                            Count))
  
  # Updating column name
  names(cluster.nums)[ncol(cluster.nums) - 1] <-
    names(gen.dat1)[i]
  
  names(cluster.nums)[ncol(cluster.nums)] <-
    paste("Num_Households_Responded_",
          names(gen.dat1)[i],
          sep = "")
}

###
# Writing csv
write.csv(cluster.nums,
          "/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/DHS_averaged_by_cluster_num_AggregateFoodGroups_23June2019.csv",
          row.names = FALSE)

