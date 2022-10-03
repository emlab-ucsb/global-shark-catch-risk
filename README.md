# Description

This repository contains code used in the paper: Burns, E.S., Bradley, D., Thomas, L.R. (in prep). Global hotspots of shark interactions with industrial longline fisheries. 

For any questions, comments, or concerns, please contact Echelle Burns [echelle_burns@uscb.edu](echelle_burns@ucsb.edu).

This repository will be citeable at Zenodo.

<insert doi>

Final datasets will be uploaded to Dryad.

<insert doi>

# Instructions

The order of running scripts should be as follows: 

- The first scripts run should be from the `src/data-wrangling` folder. The scripts within this folder are labeled starting at `01`. Scripts with the same numbers can be run in any order. For example, it makes no difference whether you clean the ICCAT or IOTC data first, but all RFMO cleaning (scripts starting with `01` should be completed before moving on to `02_combine_all_rfmo_data.Rmd`). 
- The second scripts run should be from the `src/models` folder. The scripts within this folder are labeled starting at `01`.
- Additional scripts for figure and table creation in `src/figures` and `src/tables`, respectively, should be run after data wrangling and models, but order of scripts within these folders does not matter.

# Folder Schema

Please ensure that your folder schema is the same as described below to ensure that all codes run appropriately. Several datasets were pulled from independent sources for this modeling effort and can be found in `data-updated/[relevant-dataset]/inputs` folder. All data within `intermediate` or `output` folders are generated by the codes provided. 

## Overview
+ `data-updated`
  - `chlorophyll-a`
    - `intermediates`
    - `outputs`
  - `ex-vessel-prices`
  - `global-fishing-watch`
    - `inputs`
    - `outputs`
  - `iucn-sdm-data`
    - `intermediates`
    - `outputs`
  - `mapping-templates`
  - `model-data`
    - `inputs`
    - `outputs`
  - `rfmo-boundaries`
    - `RFB_IATTC`
    - `RFB_ICCAT`
    - `RFB_IOTC`
    - `RFB_WCPFC`
  - `rfmo-data`
    - `inputs`
    - `intermediates`
    - `outputs`
  - `sea-surface-height`
    - `inputs`
    - `intermediates`
    - `outputs`
  - `sea-surface-temperature`
    - `intermediates`
    - `outputs`
  - `species-information`
+ `figures`
  - `final`
  - `intermediates`
  - `supplemental`
+ `src`
  - `data-wrangling`
  - `figures`
  - `functions`
  - `models`
  - `tables`
+ `tables`
  - `supplemental`
  
## Detailed list
+ `data-updated`: should include all datasets used and created for this project
  - `chlorophyll-a`: chlorophyll-a collated and cleaned for model use
    - `intermediates`: chlorophyll-a data collected for each year from ERDDAP (id: pmlEsaCCI42OceanColorMonthly)
      - `chla_<YY>.csv`: where <YY> represents a 4-digit year from 2012-2020
    - `outputs`: gridded chlorophyll-a data at a 1x1 and 5x5 degree spatial scale, generated from `src/data-wrangling/04_generate_grids.Rmd`
      - `binned_global_chla_1x1.csv`: gridded chlorophyll-a data at a 1x1 degree spatial scale
      - `binned_global_chla_5x5.csv`: : gridded chlorophyll-a data at a 5x5 degree spatial scale
  - `ex-vessel-prices`: ex-vessel prices for fishery caught species from 1976-2019
    - `exvessel_price_database_1976_2019.csv`: ex-vessel price data gathered from [Melnychuk et al. 2016](https://doi.org/10.1093/icesjms/fsw169) and updated to 2019 using methods described in the public-facing [github repo](https://github.com/SFG-UCSB/price-db-sfg) associated with the Melnychuk et al. 2016 paper
  - `global-fishing-watch`: remotely-sensed fishing effort collected through a collaboration with [Global Fishing Watch](https://globalfishingwatch.org/)
    - `inputs`: data downloaded directly via Global Fishing Watch queries
      - `gfw_effort_ll_ps_1x1_2012_2021.csv`: Global Fishing Watch effort (kwh) data collected for purse seines and longlines at a 1x1 degree scale from 2012-2021
      - `gfw_tonnage-groups-v20220803.csv`: specifies the gross tonnage groupings for individual vessel classes
    - `outputs`: gridded fishing effort (kwh) data at a 1x1 and 5x5 degree spatial scale, generated from `src/data-wrangling/04_generate_grids.Rmd`
      - `binned_global_gfw_1x1.csv`: gridded fishing effort (kwh) data at a 1x1 degree spatial scale
      - `binned_global_gfw_5x5.csv`: gridded fishing effort (kwh) data at a 5x5 degree spatial scale
    - `iucn-sdm-data`: species distribution data for all shark species downloaded from [IUCN](https://www.iucnredlist.org/resources/spatial-data-download) spatial data and mapping resources
      - `intermediates`: species distribution data for each species with presence re-categorized from a 1-5 scale to a 0-1 scale
        - `iucn_sid_<id>_wgs.csv`: where <id> represents a unique species identification number, as assigned by the IUCN
      - `outputs`: species distribution data for all shark species at a 1x1 and 5x5 degree spatial scale, generated from `src/data-wrangling/04_generate_grids.Rmd`
        - `global_combined_iucn_sdm_1x1.csv`: gridded species distribution data at a 1x1 degree spatial scale
        - `global_combined_iucn_sdm_5c5.csv`: gridded species distribution data at a 5x5 degree spatial scale
    - `mapping-templates`: geoTif files for plotting
      - `land_low_res_moll.tif`: used as a basemap for land designation for all figures
    - `model-data`: data generated and produced by the random forest machine learning models to predict shark catch risk
      - `inputs`: data generated for use in the machine learning model
        - `all-rfmo-models`: data generated for each RFMO random forest model that includes chlorophyll-a, sea surface temperature, sea surface height, ex-vessel prices, Global Fishing Watch, and RFMO collected datasets, generated by `src/data-wrangling/05_all_rfmo_model_data_cleaning_collation.Rmd`
          - `<rfmo>_ll_data_<res>_count_kwh`: where <rfmo> is each RFMO for which we run the models and <res> refers to a 1x1 or 5x5 degree spatial scale; a dataset for testing the models using Global Fishing Watch effort using just reported shark counts (no data reported as metric tonnes)
          - `<rfmo>_ll_data_<res>_count_hooks`: where <rfmo> is each RFMO for which we run the models and <res> refers to a 1x1 or 5x5 degree spatial scale; a dataset for testing the models using RFMO reported effort using just reported shark counts (no data reported as metric tonnes)
          - `<rfmo>_ll_data_<res>_mt_to_count_kwh`: where <rfmo> is each RFMO for which we run the models and <res> refers to a 1x1 or 5x5 degree spatial scale; a dataset for testing the models using Global Fishing Watch effort using reported shark counts and data reported as metric tonnes and converted to counts
          - `<rfmo>_ll_data_<res>_mt_to_count_hooks`: where <rfmo> is each RFMO for which we run the models and <res> refers to a 1x1 or 5x5 degree spatial scale; a dataset for testing the models using RFMO reported effort using reported shark counts and data reported as metric tonnes and converted to counts
          - `<rfmo>_ll_data_<res>_tuna_hooks`: where <rfmo> is each RFMO for which we run the models and <res> refers to a 1x1 or 5x5 degree spatial scale; a dataset with tuna catch and effort reported by RFMO
      - `outputs`: model results for each RFMO
        - `all-rfmo-models`: model results for each RFMO, generated by `src/models/01_all_rfmo_models.Rmd`
          - `<res>_count_all_rfmos_ll_effort_results.csv`: where <res> refers to a 1x1 or 5x5 degree spatial scale; results from the first round of machine learning models to determine which effort (Global Fishing Watch or RFMO reported effort) metric performed best at predicting shark catch (count only) for that resolution 
          - `<res>_mt_to_count_all_rfmos_ll_effort_results.csv`: where <res> refers to a 1x1 or 5x5 degree spatial scale; results from the first round of machine learning models to determine which effort (Global Fishing Watch or RFMO reported effort) metric performed best at predicting shark catch (count and metric tonnes converted to count) for that resolution 
          - `<rfmo>_ll_models_other_results.csv`: where <rfmo> refers to each RFMO for which we run the models; results from the second round of machine learning models to determine which additional parameters (e.g., sea surface height, ex-vessel prices) result in the best performing model
          - `<rfmo>_ll_untuned_final_predict.csv`: where <rfmo> refers to each RFMO for which we run the models; the final, global prediction from the trained machine learning model for that RFMO
          - `<rfmo>_ll_untuned_model.rds`: where <rfmo> refers to each RFMO for which we run the models; the final fitted classification and regression model, the model prediction on the test dataset, the metrics for the prediction on the test dataset, the final global prediction on the novel dataset, the metrics for the final global prediction on the novel dataset, and the feature importances for the model
    - `rfmo-boundaries`: boundaries of the RFMOs for plotting from the [FAO](https://data.apps.fao.org/map/catalog/srv/eng/catalog.search#/home)
      - `RFB_IATTC`: shapefile for the IATTC boundary
      - `RFB_ICCAT`: shapefile for the ICCAT boundary
      - `RFB_IOTC`: shapefile for the IOTC boundary
      - `RFB_WCPFC`: shapefile for the WCPFC boundary
    - `rfmo-data`: data collected from each RFMO for use in the model
      - `inputs`: data downloaded directly from each RFMO
        - `iattc`: data collected from [IATTC's public domain data](https://www.iattc.org/en-US/Data/Public-domain)
          - `PublicLLSharkMt.csv`: longline data for sharks reported in metric tonnes
          - `PublicLLSharkNum.csv`: longline data for sharks reported in counts
          - `PublicLLTunaBillfishMt.csv`: longline data for tunas and tuna-like species reported in metric tonnes
          - `PublicLLTunaBillfishNum.csv`: longline data for tunas and tuna-like species reported in number
          - `PublicPSBillfishFlag.csv`: purse seine data for billfish by flag (not used in current models, but used in data collation)
          - `PublicPSBillfishSetType.csv`: purse seine data for billfish by set type (not used in current models, but used in data collation)
          - `PublicPSSharkFlag.csv`: purse seine data for sharks by flag (not used in current models, but used in data collation)
          - `PublicPSSharkSetType.csv`: purse seine data for sharks by set type (not used in current models, but used in data collation)
          - `PublicPSTunaFlag.csv`: purse seine data for tunas by flag (not used in current models, but used in data collation)
          - `PublicPSTunaSetType.csv`: purse seine data for tunas by set type (not used in current models, but used in data collation)
          - `PublicSizePSBillfish.csv`: purse seine size data for billfish (not used in current models, but used in data collation)
          - `flag_codes.csv`: a reference to convert flag codes to country names accessible [here](https://www.iban.com/country-codes)
          - `ASFIS_sp_2019.txt`: a reference to convert species codes into species names, available from the [FAO](https://www.fao.org/fishery/en/collection/asfis/en)
        - `iccat`: data collected from [ICCAT's public domain data](https://www.iccat.int/en/accesingdb.html) - Task 2 catch/effort
          - `t2ce_20220131web.mdb`: web database that was downloaded from the public domain data
          - `t2ce_20220131web.csv`: web database that was converted to a csv for usability
          - `CODES_EffortTypes.xls`: a reference to convert effort type codes to effort types
          - `CODES_Flags-Fleets.xlsx`: a reference to convert flag codes to countries
          - `CODES_Gears.xls`: a reference to convert gear codes into gears
          - `CODES_Other.xls`: a reference to convert other codes into text
          - `CODES_SamplingAreas.xls`: a reference to convert areas into spatial coordinates
          - `CODES_Species.xlsx`: a reference to convert species codes to species name
          - `CODES_SquareTypes.xls`: a reference to convert spatial resolution codes to spatial resolutions
          - `CODES_TimePeriods.xls`: a reference to convert time period codes into time periods
          - `Species.xlsx`: a secondary reference to convert species to species name
          - `ASFIS_sp_2020.xlsx`: a secondary reference to convert species codes into species names, available from the [FAO](https://www.fao.org/fishery/en/collection/asfis/en)
        - `iotc`: data collected from [IOTC's public domain data](https://iotc.org/data/datasets/latest/CEAll)
          - `IOTC-2020-WPEB16-DATA12_CE.xlsx`: data collected from observers, primarily for purse seine (not used in current models, but used in data collation)
          - `IOTC-DATASETS-2022-02-23-CELongline_1950-2020.csv`: longline catch data from 1950-2020
          - `IOTC-DATASETS-2022-02-23-CEOther_1950-2020.csv`: other gear catch data from 1950-2020
          - `IOTC-DATASETS-2022-02-23-CESurface_1950-2020.csv`: surface gear catch data from 1950-2020
          - `ASFIS_sp_2020.xlsx`: a reference to convert species codes into species names, available from the [FAO](https://www.fao.org/fishery/en/collection/asfis/en)
        - `wcpfc`: data collected from WCPFC's public domain data
          - `WCPFC_L_PUBLIC_BY_FLAG_YR.csv`: longline catch data reported by year and flag for target species' [public domain data](https://www.wcpfc.int/public-domain)
          - `WCPFC_S_PUBLIC_BY_FLAG_YEAR.csv`: purse seine catch data reported by year and flag for target species' [public domain data](https://www.wcpfc.int/public-domain) (not used in current models, but used in data collation)
          - `BDEP Tables (MASTER - 27 July 2021).xlsx`: longline and purse seine [bycatch public domain data](https://www.wcpfc.int/public-domain-bycatch)
          - `ASFIS_sp_2020.xlsx`: a reference to convert species codes into species names, available from the [FAO](https://www.fao.org/fishery/en/collection/asfis/en)
      - `intermediates`: cleaned files originating from the input files
        - `iattc`: cleaned IATTC files that have the same data as the original input files but have been structured with consistent column names, generated by `src/data-wrangling/01_rfmo_cleaning_iattc.Rmd`
          - `PublicPSTunaSetType.csv`: originates from `../../inputs/iattc/PublicPSTunaSetType.csv`
          - `PublicPSTunaFlag.csv`: originates from `../../inputs/iattc/PublicPSTunaFlag.csv`
          - `PublicPSSharkSetType_n.csv`: originates from `../../inputs/iattc/PublicPSSharkSetType_n.csv`
          - `PublicPSSharkSetType_mt.csv`: originates from `../../inputs/iattc/PublicPSSharkSetType_mt.csv`
          - `PublicPSSharkFlag_n.csv`: originates from `../../inputs/iattc/PublicPSSharkFlag_n.csv`
          - `PublicPSSharkFlag_mt.csv`: originates from `../../inputs/iattc/PublicPSSharkFlag_mt.csv`
          - `PublicPSBillfishSetType.csv`: originates from `../../inputs/iattc/PublicPSBillfishSetType.csv`
          - `PublicPSBillfishFlag.csv`: originates from `../../inputs/iattc/PublicPSBillfishFlag.csv`
          - `PublicLLTunaBillfishNum.csv`: originates from `../../inputs/iattc/PublicLLTunaBillfishNum.csv`
          - `PublicLLTunaBillfishMt.csv`: originates from `../../inputs/iattc/PublicLLTunaBillfishMt.csv`
          - `PublicLLSharkNum.csv`: originates from `../../inputs/iattc/PublicLLSharkNum.csv`
          - `PublicLLSharkMt.csv`: originates from `../../inputs/iattc/PublicLLSharkMt.csv`
        - `iccat`: cleaned ICCAT files that have the same data as the original input files, but have been structured with consistent column names, generated by `src/data-wrangling/01_rfmo_cleaning_iccat.Rmd`
          - `t2ce_web.csv`: originates from `../../inputs/iccat/t2ce_20220131web.csv`
        - `iotc`: cleaned IOTC files that have the same data as the original input files, but have been structured with consistent column names, generated by `src/data-wrangling/01_rfmo_cleaning_iotc.Rmd`
          - `interactions.csv`: originates from `../../inputs/iotc/IOTC-2020-WPEB16-DATA12_CE.xlsx`
          - `CESurface.csv`: originates from `../../inputs/iotc/IOTC-DATASETS-2022-02-23-CESurface_1950-2020.csv`
          - `CEOther.csv`: originates from `../../inputs/iotc/IOTC-DATASETS-2022-02-23-CEOther_1950-2020.csv`
          - `CELongline.csv`: originates from `../../inputs/iotc/IOTC-DATASETS-2022-02-23-CELongline_1950-2020.csv`
        - `wcpfc`: cleaned WCPFC files that have the same data as the original input files, but have been structured with consistent column names, generated by `src/data-wrangling/01_rfmo_cleaning_wcpfc.Rmd`
          - `purseseine_flagyear.csv`: originates from `../../inputs/wcpfc/WCPFC_S_PUBLIC_BY_FLAG_YEAR.csv`
          - `longline_flagyear.csv`: originates from `../../inputs/wcpfc/WCPFC_L_PUBLIC_BY_FLAG_YR.csv`
          - `bycatch_seine.csv`: originates from `../../inputs/wcpfc/BDEP Tables (MASTER - 27 July 2021).xlsx`
          - `bycatch_longline.csv`: originates from `../../inputs/wcpfc/BDEP Tables (MASTER - 27 July 2021).xlsx`
      - `outputs`: collated outputs by RFMO and all RFMOs combined
        - `wcpfc-all.csv`: all individual WCPFC files combined, generated by `src/data-wrangling/01_rfmo_cleaning_wcpfc.Rmd`
        - `iccat-all.csv`: all individual ICCAT files combined, generated by `src/data-wrangling/01_rfmo_cleaning_iccat.Rmd`
        - `iotc-all.csv`: all individual IOTC files combined, generated by `src/data-wrangling/01_rfmo_cleaning_iotc.Rmd`
        - `iattc-all.csv`: all individual IATTC files combined, generated by `src/data-wrangling/01_rfmo_cleaning_iattc.Rmd`
        - `all_data.csv`: data combined from all RFMOs to be used in the machine learning models, generated by `src/data-wrangling/02_combine_all_rfmo_data.Rmd`
    - `sea-surface-height`: sea surface height data collated and cleaned for model use
      - `inputs`: sea surface height data gathered from [MEaSUREs Gridded Sea Surface Height Anomalies Version 1812](ttps://doi.org/10.5067/SLREF-CDRV2)
        - `ssh_grids_v1812_<dateid>.nc`: where <dateid> is a unique combination of date and time; netcdf files gathered directly from the MEaSUREs database
      - `intermediates`: sea surface height data converted to .csv files and collated by mean, generated from `src/data-wrangling/04_generate_grids.Rmd`
        - `ssh_grids_v1812_<dateid>.csv`: where <dateid> is a unique combination of date and time; csv files converted from netcdf files
        - `mean_ssh_<YY>.csv`: where <YY> is a 4 digit year; the mean sea surface height for a particular year
      - `outputs`: gridded sea surface height data at a 1x1 and 5x5 degree spatial scale, generated from `src/data-wrangling/04_generate_grids.Rmd`
        - `binned_global_ssh_1x1.csv`: gridded sea surface height data at a 1x1 degree spatial scale
        - `binned_global_ssh_5x5.csv`: gridded sea surface height data at a 5x5 degree spatial scale
    - `sea-surface-temperature`: sea surface temperature collated and cleaned for model use
      - `intermediates`: sea surface temperature data collected for each year from ERDDAP (id: erdHadISST)
        - `sst_<YY>.csv`: where <YY> represents a 4-digit year from 2012-2020
      - `outputs`: gridded sea surface temperature data at a 1x1 and 5x5 degree spatial scale, generated from `src/data-wrangling/04_generate_grids.Rmd`
        - `binned_global_sst_1x1.csv`: gridded sea surface temperature data at a 1x1 degree spatial scale
        - `binned_global_sst_5x5.csv`: : gridded sea surface temperature data at a 5x5 degree spatial scale
    - `species-information`: various files that provide information on species groups, species weight (mt) to count conversions, etc. 
      - `species_groups.xlsx`: a list of species scientific names, common names, and species group for categorization
      - `species_list_for_sdms.csv`: a list of species for which we should grab species distribution data for modeling
      - `species_list.csv`: the full list of species, including all possible species within "nei" categories
      - `spp_list_fishing_thread.csv`: a list of species that are threatened by fishing
      - `spp_weight_count_conversion.csv`: the conversion metrics we use to calculate count from metric tonnes using methods described in [Worm et al., 2013](https://doi.org/10.1016/j.marpol.2012.12.034)
+ `figures`: output figures
  - `final`: figures used in the submitted manuscript
  - `intermediates`: figures generated during intermediate steps of the modeling process
  - `supplemental`: figures used in the supplemental data in the submitted manuscript
+ `src`: all scripts used to clean, model, and analyze data
  - `data-wrangling`: scripts for data cleaning and collation
    - `01_rfmo_cleaning_<rfmo>.Rmd`: where <rfmo> refers to each RFMO for which we run the models; script to pull raw data, clean and collate
    - `02_combine_all_rfmo_data.Rmd`: script that combines and harmonizes all RFMO data
    - `03_species_list_for_sdms.Rmd`: script that gathers relevant species for pulling IUCN species distribution data
    - `04_generate_grids.Rmd`: converts non-RFMO spatial data (chlorophyll-a, sea surface temperature, sea surface height, Global Fishing Watch effort, IUCN species distribution models) into a gridded format to match with the RFMO data
    - `05_all_rfmo_model_data_cleaning_collation.Rmd`: generates datasets for each RFMO to be used in model training
  - `figures`: all scripts used to generate figures for the submitted manuscript
    - `figure_<x>.R`: where <x> corresponds to a figure number in the submitted manuscript; script for creating individual figures
    - `plot_dfaults.R`: script that provides datasets and themes used across all figures
    - `supplemental_scaled_predicted_catch.R`: script that generates a figure for the predicted catch where each RFMO is independently scaled
    - `supplemental_workflow.R`: script that generates a figure for the result of assumptions we made during model creation and data manipulation processes
  - `functions`: functions used throughout
    - `func-all_rfmo_effort_models.R`: function `all_rfmo_effort_models()` used in `src/models/01_all_rfmo_models.Rmd` to test the predictive power of different effort sources
    - `func-all_rfmo_other_models.R`: function `all_rfmo_other_models()` used in `src/models/01_all_rfmo_models.Rmd` to test the predictive power of different predictor variables, once the effort source had been chosen
    - `func-all_rfmo_untuned_models.R`: function `all_rfmo_untuned_models()` used in `src/models/01_all_rfmo_models.Rmd` to generate the final models and predictions used for manuscript results
    - `func-get_names.R`: function `get_names()` used in `src/data-wrangling/03_species_list_for_sdms.Rmd` to gather names of species within particular species groups
    - `func-train_test_split.R`: function `train_test_split()` used in `all_rfmo_effort_models()`, `all_rfmo_other_models()`, and `all_rfmo_untuned_models()` to split the data into a training and testing dataset spatially and temporally
  - `tables`: all scripts used to generate tables
    - `table_1.R`: script used to generate table 1 in the submitted manuscript
    - `supplemental_raw_catch_rfmo_species.R`: script used to calculate the raw proportions of shark catch by each RFMO
    - `supplemental_percent_zero_catch_predicted_nonzero.R`: script used to determine which cells were reported as non-zero values but predicted as 0 values for gut-checking 
+ `tables`: output tables
  - `supplemental`: tables used in the supplemental data in the submitted manuscript

# R Version
All code was run using RStudio: 2022.07.1+554 for macOS and R version 4.2.1

# Required Libraries
+ Data Ingestion, Cleaning, Harmonization, and Organization
  - `tidyverse` (version 1.3.2)
  - `readxl` (version 1.4.1)
  - `reshape` (version 0.8.9)
  - `here` (version 1.0.1)
  - `rfishbase` (version 4.0.0)
  - `googledrive` (version 2.0.0)
  - `readr` (version 2.1.2)
  - `vroom` (version 1.5.7)
+ Geospatial Data Ingestion and Manipulation
  - `sf` (version 1.0-8)
  - `rdgal` (version 1.5-32)
  - `raster` (version 3.5-29)
  - `fasterize` (version 1.0.3)
  - `geosphere` (version 1.5-14)
+ Data Visualization
  - `scales` (version 1.2.1)
  - `knitr` (version 1.40)
  - `RColorBrewer` (version 1.1-3)
  - `cowplot` (version 1.1.1)
  - `paletteer` (version 1.4.1)
  - `tmap` (version 3.3-3)
+ Data Analysis
  - `pscl` (version 1.5.1)
  - `tidymodels` (version 1.0.0)
  - `broom.mixed` (version 0.2.9.4)
  - `MultivariateRandomForest` (version 1.1.5)
