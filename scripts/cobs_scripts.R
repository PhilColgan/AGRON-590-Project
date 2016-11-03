## Loading data and setting working directory

install.packages(c("dplyr" , "ggplot2"))
getwd()
cobs_data<- read.csv("data/KBase_MGRast_Metadata_9May2013_EMB.csv")

##deleting blank columns


## changing heading names

colnames(cobs_data)<-c("sample_Id" , "sample_month" , "sample_year" , "crop" , "sample_block" , "agg_frac" , "MGRAST_Id" , "agrochem_addition" , "crop_rot" , "land_use" , "veg_class" , "veg_class_meth" , "drain_class" , "extreme_event" , "FAO_class" , "fire_hist" , "soil_hor" , "soil_hor_meth" , "link_soil_method" , "soil_tax" , "soil_tax_meth" , "MGRAST_Id" , "micro_bm" , "micro_bm_meth" , "misc_param" , "pH" , "pH_meth" , "dna_mix" , "land_use_pre" , "land_use_pre_meth" , "sample_position" , "salinity_meth" , "sample_wt_dna" , "siev_size" , "slope_aspect" , "slope_grad" , "soil_type" , "soil_type_meth" , "store_cond" , "texture" , "texture_meth" , "till" , "total_N" , "total_N_meth" , "total_OC_meth" , "total_OC" , "soil_water" , "soil_water_meth" , "total_C" , "misc_param_1" , "MBN_dry" , "MBN_applied" , "Ext_C_dry" , "Ext_C_applied" , "Ext_C_N_dry" , "Ext_N_applied" , "Bulk_dense" , "Ext_P_dry" , "AMF_col" , "AP_act" , "BG_act" , "BX_act" , "CB_act" , "NAG_act" , "Sum_C_act" , "MBC_dry" , "MBC_applied" , "MBC_MBN_meth" , "Ext_C_Ext_N_meth" , "AMF_col_meth" , "root_bm" , "root_dep" , "AMF_col_bm" , "MBC:MBN" , "MWD" , "agg_frac_prop" , "N2O_2011" , "CH4_2011" , "N2O_2012" , "CO2_2011" , "CO2_2012")
cobs_data
