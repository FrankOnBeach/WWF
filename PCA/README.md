# Principal Component Analysis(PCA) on DHS Food Balance Sheet Data with Socioeconomic Data.

For this project, we use PCA to reduce the dimensionality of our dataset, and we are interested to see if there is any natural grouping among our dataset that aligns with their attributes. In this case: continent level, region level, and income group level.  On top of PCA, we run Kmeans to see if Kmeans clustering agrees with the clustering and see if machine learning can pick up anything else that is shared between countries. 

We retreieved region and income level information from UN's country classification dataset. 
https://datahelpdesk.worldbank.org/knowledgebase/articles/906519

https://nbviewer.jupyter.org/github/FrankOnBeach/WWF/blob/master/PCA/PCA%20Analysis%20.ipynb for plotly embedded Jupyter Notebook. 

## Analysis 1

In our second analysis, we are interested to see if countries cluster together based solely on protein food balance sheet data. This incluses Milk, Bovine meat, Sheep and goat meat, Pigmeat, Poultry, Eggs, Fish, Maize, Rice, Wheat, Rye, Barley, Pulses, Starchy roots on continent region and income group level. 


## Analysis 2

In our second analysis, we are interested to see if countries cluster together based on protein food balance sheet data and some selected socioeconomic variables. We group protein foods into groups and introduce three new variables which are Average GDP, Urban and rural ratio(urban population/rural population), and permenent crops ratio(permenent crops area/total country area). This result into 'Bovine_Sheep_Pigmeat_Poultry','Cereals','Milk_Egg','Fish',' Av 2011+2012+2013 extracted from WDI (constant 2010 US$)','Urban_Rural_ratio','Permenent_Crops_Ratio'as variables.

## Analysis 3

In our third analysis, we are interested to see if countries cluster together based on protein food balance sheet data and child health data. We group protein foods into groups and introduce two new variables which are children under 5 years old stunting rate and overweight rate. This data is extracted from http://apps.who.int/gho/data/node.main.CHILDSTUNTED?lang=en.  This result into 'Bovine_Sheep_Pigmeat_Poultry','Cereals','Milk_Egg','Fish','stunting(percentage)', 'overweight(percentage)' as variables. 

## Analysis 4

Based on previous analyses, countries' natural clustering aligns best with the in come group level. The plots show that low-middle income are often clustered together by k means. We are interested to see if a lower number of clusters(3) will yield a better result. 


## Datasets Inventory
- **master_final_merged_df_protein_gdp.csv**:
  - This is the master data file containing food balance sheet data and socioecnomoic data. 
- **pca_protein.csv**:
  - Result file. This file includes Principal componenet 1 and 2 computed using only protein foods data, the data used in computing pca and the clusters K means assigned the datapoint to. 
- **pca_socio.csv**:
  - Result file. This file includes Principal componenet 1 and 2 computed using grouped protein foods data along with socioeconomic data, the data used in computing pca,  and the clusters K means assigned the datapoint to.
- **pca_health.csv**:
  - Result file. This file includes Principal componenet 1 and 2 computed using grouped protein foods data along with health data, the data used in computing pca,  and the clusters K means assigned the datapoint to.
- **CLASS.xls**:
  - Source file. This file contains country classification data from UN. https://datahelpdesk.worldbank.org/knowledgebase/articles/906519 Using 2012 data
- **NUTRITION_HA_2.csv**:
  - Source file. This file contains stunting percentage of children under 5 years old. This dataset is from WHO.  
- **NUTRITION_WH2.csv.csv**:
  - Source file. This file contains overweight percentage of children under 5 years old. This dataset is from WHO.  
  
