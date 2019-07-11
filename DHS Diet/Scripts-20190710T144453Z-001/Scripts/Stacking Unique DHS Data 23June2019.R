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

# Checking to see unique values in each column
cols.character <-
  as.numeric(which(sapply(stacked.dat,
             is.numeric) == FALSE))
library(profvis)

for(i in 3:length(cols.character)) {
  print(names(stacked.dat)[cols.character[i]])
  print(sort(unique(stacked.dat[,cols.character[i]])))
  
  pause(0.5)
}

# Saving data frame just in case...
stacked.dat.tmp <-
  stacked.dat

# Converting responses reported as Yes/No to 1/0s
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

# And converting columns to numeric format so that we can do math on them
for(i in 6:75) {
  stacked.dat.tmp[,i] <-
    as.numeric(stacked.dat.tmp[,i])
}

# Getting subset of data, will be used later
# this contains data on how to identify the different surveys
# 1) case identification
# 2) country code and phase
# 3) cluster number
# 4) household number
# 5) country
gen.dat <-
  stacked.dat.tmp[,c(2:5,76)]

# and dropping columns that contain no information (e.g. all responses are blanks)
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
# Creating our own identification
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
# And looping through columns to get summary by household
for(i in 6:27) {
  ###
  # Getting select columns
  tmp.dat <-
    gen.dat1[,c(28,i)]
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
    names(gen.dat1)[i]
}


###
# Writing csv
write.csv(case.ids,
          "/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/DHS_Averaged_by_Household_23June2019.csv",
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
for(i in 5:26) {
  ###
  # Getting select columns
  tmp.dat <-
    gen.dat1[,c(27,i)]
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
          "/Volumes/Citadel/Oxford/Research Projects/Protein Mapping/Managed_DHS_Data/DHS_averaged_by_cluster_num_23June2019.csv",
          row.names = FALSE)

