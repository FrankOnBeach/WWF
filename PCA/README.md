# Principal Component Analysis(PCA) on DHS Food Balance Sheet Data with Socioeconomic Data.

For this project, we use PCA to reduce the dimensionality of our dataset, and we are interested to see if there is any natural grouping among our dataset that aligns with their attributes. In this case: continent level, region level, and income group level.  On top of PCA, we run Kmeans to see if Kmeans clustering agrees with the clustering and see if machine learning can pick up anything else that is shared between countries. 

We retreieved region and income level information from UN's country classification dataset. 
https://datahelpdesk.worldbank.org/knowledgebase/articles/906519

https://nbviewer.jupyter.org/github/FrankOnBeach/WWF/blob/master/PCA/PCA%20Analysis%20.ipynb for plotly embedded Jupyter Notebook. 

## Analysis Part 1

### Analysis 1

In our second analysis, we are interested to see if countries cluster together based solely on protein food balance sheet data. This incluses Milk, Bovine meat, Sheep and goat meat, Pigmeat, Poultry, Eggs, Fish, Maize, Rice, Wheat, Rye, Barley, Pulses, Starchy roots on continent region and income group level. 


### Analysis 2

In our second analysis, we are interested to see if countries cluster together based on protein food balance sheet data and some selected socioeconomic variables. We group protein foods into groups and introduce three new variables which are Average GDP, Urban and rural ratio(urban population/rural population), and permenent crops ratio(permenent crops area/total country area). This result into 'Bovine_Sheep_Pigmeat_Poultry','Cereals','Milk_Egg','Fish',' Av 2011+2012+2013 extracted from WDI (constant 2010 US$)','Urban_Rural_ratio','Permenent_Crops_Ratio'as variables.

### Analysis 3

In our third analysis, we are interested to see if countries cluster together based on protein food balance sheet data and child health data. We group protein foods into groups and introduce two new variables which are children under 5 years old stunting rate and overweight rate. This data is extracted from http://apps.who.int/gho/data/node.main.CHILDSTUNTED?lang=en.  This result into 'Bovine_Sheep_Pigmeat_Poultry','Cereals','Milk_Egg','Fish','stunting(percentage)', 'overweight(percentage)' as variables. 

### Analysis 4

Based on previous analyses, countries' natural clustering aligns best with the in come group level. The plots show that low-middle income are often clustered together by k means. We are interested to see if a lower number of clusters(3) will yield a better result. 

## Analysis Part 2 

In the second part of this project, we've decided to explore some other variables. 

### Analysis 1 

For this analysis, we are using the production level data from 2011 to 2013. 

Production Level on the following products:

Barley and products Bovine Meat Maize and products Mutton & Goat Meat Pigmeat Poultry Meat Pulses, Other and products Rice (Milled Equivalent) Rye and products Wheat and products Eggs + (Total) Fish, Seafood + (Total) Milk - Excluding Butter + (Total) Starchy Roots + (Total)

We have also tried to group these products to see if produce a better result. This result into these categories.Fish and Seafood,Bovine_Sheep_Pigmeat_Poultry,Milk_Egg, Cereal.

Unit is 1000 tonns.

## Datasets Inventory

- **master_final_merged_df_protein_gdp.csv**:
  - This is the master data file containing food balance sheet data and socioecnomoic data. 
- **CLASS.xls**:
  - Source file. This file contains country classification data from UN. https://datahelpdesk.worldbank.org/knowledgebase  
- **NUTRITION_HA_2.csv**:
  - Source file. This file contains stunting percentage of children under 5 years old. This dataset is from WHO.  
- **NUTRITION_WH2.csv.csv**:
  - Source file. This file contains overweight percentage of children under 5 years old. This dataset is from WHO.  
  
### Result Folder
- **pca_protein.csv**:
  - Result file. This file includes Principal componenet 1 and 2 computed using only protein foods data, the data used in computing pca and the clusters K means assigned the datapoint to. 
- **pca_socio.csv**:
  - Result file. This file includes Principal componenet 1 and 2 computed using grouped protein foods data along with socioeconomic data, the data used in computing pca,  and the clusters K means assigned the datapoint to.
- **pca_health.csv**:
  - Result file. This file includes Principal componenet 1 and 2 computed using grouped protein foods data along with health data, the data used in computing pca,  and the clusters K means assigned the datapoint to.
/articles/906519 Using 2012 data
- **production_avg_2011_2013.csv**:
  - Intermediate file. This is the pivoted and averaged dataset for production level. From year 2011 to 2013. Added Country Code and Continent for easier mergeing. Source from FAO food balance sheets. 

### part 2 analysis 2
- **agri_employment.csv**:
  - This is the percentage of agricultural employment file.  https://data.worldbank.org/indicator/SL.AGR.EMPL.ZS (modeled estimate)
- **agri_forestry_fish_gdp.csv**:
  - This is the percentage of agricultural contribution to GDP.  https://data.worldbank.org/indicator/NV.AGR.TOTL.ZS 
- **unemployment.csv**:
  - This is the unemployment rate.  https://data.worldbank.org/indicator/SL.UEM.TOTL.ZS (modeled ILO)
- **gdp_constant_2010.csv**:
  - This is the GDP file.  https://data.worldbank.org/indicator/NY.GDP.MKTP.CD (Constant 2010 US dollars)
  
  
  
