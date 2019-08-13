# Machine Learning On DHS Children's Nutrition Data

## Project Overview
In this project, we grouped the DHS Children's Nutrition data by cell. Each cell represents households in the area.  We were interested to see if we are able to predict the precentage of each item is given to children of the households in each cell.

The predictor variables we used in this analysis include: 
- dhs_predictor_socio_cols = ['ari_employment value (AVG 2011-2013)',
                            'agri_gdp_percentage(AVG from 2011 to 2013)',
                           'AVG GDP 2011-2013 Constant 2010']

- dhs_predictor_consump_cols = ['Cons_kg.yr Cereals - Excluding Beer', 'Cons_kg.yr Eggs'
                        ,'Cons_kg.yr Fish, Seafood'
                       ,'Cons_kg.yr Meat', 'Cons_kg.yr Milk - Excluding Butter'
                       ,'Cons_kg.yr Nuts_and_Seeds', 'Cons_kg.yr Pulses'
                       ,'Cons_kg.yr Starchy Roots']

- dhs_predictor_env_cols = ['Agricultural Emissions (CO2eq)'
                           ,'Protected_Terrestrial'
                          , 'Forest_change' 
                          , 'Redlist'
                         ,'Freshwater']

- dhs_predictor_production_cols = ['Fish, Seafood', 'Bovine_Sheep_Pigmeat_Poultry', 'Milk_Egg', 'Cereals']
                    
- dhs_predictor_other_cols = ['Population_2015','Primary_Religion', 'Travel_time_to_cities','Soil_suit_overall', 'Soil_suit_cereals']

And our outcome variables are: eggs,dairy,meat_fish, nuts_pulses, grains_tubers, animal_source_food.

We tried different binning strategies for our outcome variables to see which was more appropriate while with better accuracy. We started by cutting our outcome variables into bins with width of 0.25. This resulted in 4 extremely unequal bins. We also tried customized bins as 0,0.1,0.2,0.5. While this resulted in a more evenly distributed classes, the accuracy dropped sigificantly for all models. 

Different models were used in this project including K Nearest Neighbors, Decision Tree, Random Forest. To further improve the accuracy we also tried enemble methods such as bagging and boosting. The best performaing model in all variables is XGBoost. 

We achieved an accuracy score around 50% when testing the model on the test set, also a similar core while using cross validation. We further attempted to improve the model by fine tuning the hyperparameters. This helped us slightly increase the accuracy. 

## Future Work

So far, we've been only working on the grouped dataset. It would be interesting to work on the individual data. This would potentially give us a better result since the outcome variable will be binary as in whether the household feed the child with this item. Secondly, we can use other information from the survey as new features so that we are not limited to regional and global datasets. 

The caviet here is that the individual data is too large to be processed on local computers. We should use clusters on a cloud computing service such as aws with the help of spark. 




## Datasets Inventory

