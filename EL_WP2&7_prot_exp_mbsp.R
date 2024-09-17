###Multivariate Bayesian Model with Shrinkage Priors
### Zhiyang Wang
###May 3rd, 2024

library(fastDummies)
library(MBSP)
library(foreign)

set.seed(050320241)

im_dataset_FT12<-read.dta("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/imputated_FT12_dataset.dta")
#load imputated data from the first stage of ExWAS
#i save the data in dta form, therefore i used read.dta, you could change to read_csv, etc. 
## change FT12 to WALNUT in coding for better smoothing usage. 

#------------ MBSP--------------------------------------------------------------------------------------------------------#
row.names(im_dataset_FT12)<-im_dataset_FT12$person_nb ##person_nb is the participant ID in FT12 

im_dataset_FT12_bmi<-na.omit(im_dataset_FT12)
##i set a separated protein, exposure amd confounder, because a part of participants missing information in BMI 

protein <- as.matrix(im_dataset_FT12[,c("P05155", 	"P01008", 	"Q9UK55", 	"P04180", 	"P02787", 	"O95445", 	"P08571", 	
                                        "P06681", 	"P00734", 	"P00751", 	"P00742", 	"P08185", 	"P10909", 	"P11279", 	
                                        "Q14126", 	"P09172", 	"P23142", 	"Q9Y6R7", 	"P07911", 	"Q86TH1", 	"P0C0L5", 	
                                        "P40189", 	"Q13449", 	"P28799", 	"Q12805", 	"Q9BY67")])
protein_bmi <- as.matrix(im_dataset_FT12_bmi[,c("P05155", 	"P01008", 	"Q9UK55", 	"P04180", 	"P02787", 	"O95445", 	"P08571", 	
                                        "P06681", 	"P00734", 	"P00751", 	"P00742", 	"P08185", 	"P10909", 	"P11279", 	
                                        "Q14126", 	"P09172", 	"P23142", 	"Q9Y6R7", 	"P07911", 	"Q86TH1", 	"P0C0L5", 	
                                        "P40189", 	"Q13449", 	"P28799", 	"Q12805", 	"Q9BY67")])
                                     
#for WALNUTs, please change the list to c("P05155", 	"P01008", 	"P25786", 	"Q9UK55", 	"P04180", 	"P02787", 	"O95445", 	"P08571", 	"P06681", 	
##                                    "P51693", 	"P00734", 	"P00751", 	"P00742", 	"P08185", 	"P13611", 	"O95274", 	"P10909", 	"P11279", 	
##                                    "Q14126", 	"P09172", 	"P23142", 	"Q9Y6R7", 	"Q99983", 	"P07911", 	"Q86TH1", 	"P0C0L5", 	"P40189", 	
##                                    "Q13449", 	"P28799", 	"Q12805", 	"Q9BY67", 	"P07858")
                                     
# i have inquire the package developer, there is no need for scaling 

im_dataset_FT12<-dummy_columns(im_dataset_FT12, select_colums<-c('occup_m', 'occup_f'), remove_first_dummy=TRUE)
im_dataset_FT12<-im_dataset_FT12[, -which(names(im_dataset_FT12) %in% c( 'occup_m', 'occup_f'))] 
# Here is convert non-ordered categorical exposures to dummy variables
row.names(im_dataset_FT12)<-im_dataset_FT12$person_nb
exposure<-im_dataset_FT12[,c("built_ghs_300",	"built_ghs_500",	"dist_green_clc",	"size_green_clc",	"dist_blue_clc",	"size_blue_clc",	"pop_wp_100",	
                        "pop_wp_300",	"pop_wp_500",	"ndvi_5yrs_all_100",	"ndvi_5yrs_greenest_100",	"ndvi_5yrs_all_300",	"ndvi_5yrs_greenest_300",	
                        "ndvi_5yrs_all_500",	"ndvi_5yrs_greenest_500",	"msavi_5yrs_all_100",	"msavi_5yrs_greenest_100",	"msavi_5yrs_all_300",	
                        "msavi_5yrs_greenest_300",	"msavi_5yrs_all_500",	"msavi_5yrs_greenest_500",	"treecover_100",	"treecover_300",	
                        "treecover_500",	"ne_dem_100",	"ne_dem_300",	"ne_dem_500",	"ne_slo_100",	"ne_slo_300",	"ne_slo_500",	"plyVideogame_ch_cr",	
                        "alc_ch_cr",	"cob_f_fr",	"edu_f_fr",	"cob_m_mr",	"edu_m_mr",	"preg_smk_m_mr", "hh_adultsnr_pr",	
                        "famsize_child_pr",	"smk_exp_ch_fr",	"age0",	"age10",	"age20",	"age30",	"age40",	"age50",	"age60",	"age70",	"age80",	
                        "age90",	"over66",	"under18",	"background",	"edu_high", "occup_m_2", "occup_m_3", "occup_f_2", "occup_f_3","BC","NO2","O3","PM25")] 
# Here is to choose exposures, you will get a matrix with 62 exposure (58 non-dummies and 4 dummies)

## keep binary variable (0,1) and numeric style
exposure$cob_m_mr<-(as.numeric(exposure$cob_m_mr)-1)
exposure$cob_f_fr<-(as.numeric(exposure$cob_f_fr)-1)
exposure$alc_ch_cr<-(as.numeric(exposure$alc_ch_cr)-1)
exposure$preg_smk_m_mr<-(as.numeric(exposure$preg_smk_m_mr)-1)
exposure$smk_exp_ch_fr<-(as.numeric(exposure$smk_exp_ch_fr)-1)


exposure$plyVideogame_ch_cr<-as.numeric(exposure$plyVideogame_ch_cr)
exposure$edu_m_mr<-as.numeric(exposure$edu_m_mr)
exposure$edu_f_fr<-as.numeric(exposure$edu_f_fr)

exposure<-as.matrix(sapply(exposure,as.numeric))  


im_dataset_FT12_bmi<-dummy_columns(im_dataset_FT12_bmi, select_colums<-c('occup_m', 'occup_f'), remove_first_dummy=TRUE)
im_dataset_FT12_bmi<-im_dataset_FT12_bmi[, -which(names(im_dataset_FT12_bmi) %in% c( 'occup_m', 'occup_f'))] 
# Here is convert non-ordered categorical exposures to dummy variables
row.names(im_dataset_FT12_bmi)<-im_dataset_FT12_bmi$person_nb
exposure_bmi<-im_dataset_FT12_bmi[,c("built_ghs_300",	"built_ghs_500",	"dist_green_clc",	"size_green_clc",	"dist_blue_clc",	"size_blue_clc",	"pop_wp_100",	
                             "pop_wp_300",	"pop_wp_500",	"ndvi_5yrs_all_100",	"ndvi_5yrs_greenest_100",	"ndvi_5yrs_all_300",	"ndvi_5yrs_greenest_300",	
                             "ndvi_5yrs_all_500",	"ndvi_5yrs_greenest_500",	"msavi_5yrs_all_100",	"msavi_5yrs_greenest_100",	"msavi_5yrs_all_300",	
                             "msavi_5yrs_greenest_300",	"msavi_5yrs_all_500",	"msavi_5yrs_greenest_500",	"treecover_100",	"treecover_300",	
                             "treecover_500",	"ne_dem_100",	"ne_dem_300",	"ne_dem_500",	"ne_slo_100",	"ne_slo_300",	"ne_slo_500",	"plyVideogame_ch_cr",	
                             "alc_ch_cr",	"cob_f_fr",	"edu_f_fr",	"cob_m_mr",	"edu_m_mr",	"preg_smk_m_mr", "hh_adultsnr_pr",	
                             "famsize_child_pr",	"smk_exp_ch_fr",	"age0",	"age10",	"age20",	"age30",	"age40",	"age50",	"age60",	"age70",	"age80",	
                             "age90",	"over66",	"under18",	"background",	"edu_high", "occup_m_2", "occup_m_3", "occup_f_2", "occup_f_3","BC","NO2","O3","PM25")] 
# Here is to choose exposures, you will get a matrix with 62 exposure (54 non-dummies and 4 dummies)

## keep binary variable (0,1) and numeric style
exposure_bmi$cob_m_mr<-(as.numeric(exposure_bmi$cob_m_mr)-1)
exposure_bmi$cob_f_fr<-(as.numeric(exposure_bmi$cob_f_fr)-1)
exposure_bmi$alc_ch_cr<-(as.numeric(exposure_bmi$alc_ch_cr)-1)
exposure_bmi$preg_smk_m_mr<-(as.numeric(exposure_bmi$preg_smk_m_mr)-1)
exposure_bmi$smk_exp_ch_fr<-(as.numeric(exposure_bmi$smk_exp_ch_fr)-1)


exposure_bmi$plyVideogame_ch_cr<-as.numeric(exposure_bmi$plyVideogame_ch_cr)
exposure_bmi$edu_m_mr<-as.numeric(exposure_bmi$edu_m_mr)
exposure_bmi$edu_f_fr<-as.numeric(exposure_bmi$edu_f_fr)

exposure_bmi<-as.matrix(sapply(exposure_bmi,as.numeric))  


confounder_withoutbmi <- as.matrix(im_dataset_FT12[, c("age", "sex")])
confounder_withbmi <- as.matrix(im_dataset_FT12_bmi[, c("age", "sex", "bmi_ya")])
#Here is to get two matrix for confounder (one is with BMI, one is without BMI)


#-------------------------------------------------------------------------------------------------------------------------------------------#
#Now, we will try MBSP models with 9 different sets of parameter and adjusted confounders including BMI
#-------------------------------------------------------------------------------------------------------------------------------------------#
#-------- HORSESHOES PRIOR u=0.5, a=0.5
#-- Try tau=1e-8 (very small: 10^{-8})
mbsp_model_bmi_FT12_1 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=0.5, a=0.5, tau=1e-8, model_criteria=TRUE)
# variables selected by the MBSP model with horseshoes prior and when tau is very small
mbsp_model_bmi_FT12_1$active_predictors

#-- Try tau=1e4 (very large: 10^4)
mbsp_model_bmi_FT12_2 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=0.5, a=0.5, tau=1e4, model_criteria=TRUE)
# variables selected by the MBSP model with horseshoes prior and when tau is very large
mbsp_model_bmi_FT12_2$active_predictors

#-------- NEG priors u=1, a=0.1
#-- Try tau=1e-8 (very small: 10^{-8})
mbsp_model_bmi_FT12_3 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=1, a=0.1, tau=1e-8, model_criteria=TRUE)
# variables selected by the MBSP model with NEG prior and when tau is very small
mbsp_model_bmi_FT12_3$active_predictors

#-- Try tau=1e4 (very large: 10^4)
mbsp_model_bmi_FT12_4 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=1, a=0.1, tau=1e4, model_criteria=TRUE)
# variables selected by the MBSP model with NEG prior and  when tau is very large
mbsp_model_bmi_FT12_4$active_predictors

#-------- Strawderman-berger priors (u=1, a=0.5)
#-- Try tau=1e-8 (very small: 10^{-8})
mbsp_model_bmi_FT12_5 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=1, a=0.5, tau=1, model_criteria=TRUE)
# variables selected by the MBSP model with trawderman-berger priors and when tau is very small
mbsp_model_bmi_FT12_5$active_predictors

#-- Try tau=1e4 (very large: 10^4)
mbsp_model_bmi_FT12_6 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=1, a=0.5, tau=1e4, model_criteria=TRUE)
# variables selected by the MBSP model with trawderman-berger priors and when tau is very large
mbsp_model_bmi_FT12_6$active_predictors

#-------- Use 'default' tau of tau=1/(p*sqrt(n*log(n)) but change u and a
# Corresponds to NEG prior tau=1/(p*sqrt(n*log(n)) but change u and a
mbsp_model_bmi_FT12_7 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=1, a=0.1, model_criteria=TRUE)
# variables selected by the MBSP model when NEG prior is used. 
# Gives a less sparse model, due to NEG having lighter tails. 
mbsp_model_bmi_FT12_7$active_predictors

#-------- Corresponds to horseshoe prior with default tau=1/(p*sqrt(n*log(n))
mbsp_model_bmi_FT12_8 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=0.5, a=0.5, model_criteria=TRUE)
# variables selected by the MBSP model
mbsp_model_bmi_FT12_8$active_predictors

#-------- Corresponds to Strawderman-berger prior (u=1, a=0.5)  with default tau=1/(p*sqrt(n*log(n))
mbsp_model_bmi_FT12_9 <- MBSP(Y=protein_bmi,X=exposure_bmi,confounders=confounder_withbmi, u=1, a=0.5, model_criteria=TRUE)
# variables selected by the MBSP model
mbsp_model_bmi_FT12_9$active_predictors

#--------Now, lets sreen the DIC and WAIC of each model, there is no cutoff for these two fit indeces
mbsp_model_bmi_FT12_1$DIC # HORSESHOES PRIOR u=0.5, a=0.5  (very small: 10^{-8})
mbsp_model_bmi_FT12_2$DIC # HORSESHOES PRIOR u=0.5, a=0.5 (very large: 10^4)
mbsp_model_bmi_FT12_3$DIC # NEG prior u=1, a=0.1  (very small: 10^{-8})
mbsp_model_bmi_FT12_4$DIC # NEG prior u=1, a=0.1 (very large: 10^4)
mbsp_model_bmi_FT12_5$DIC # Strawderman-berger prior (u=1, a=0.5)  (very small: 10^{-8})
mbsp_model_bmi_FT12_6$DIC # Strawderman-berger prior (u=1, a=0.5) (very large: 10^4)
mbsp_model_bmi_FT12_7$DIC # NEG prior tau=1/(p*sqrt(n*log(n)) 
mbsp_model_bmi_FT12_8$DIC # horseshoe prior with default tau=1/(p*sqrt(n*log(n))
mbsp_model_bmi_FT12_9$DIC # Strawderman-berger prior (u=1, a=0.5)  with default tau=1/(p*sqrt(n*log(n))

#--------based on the model which have capture any exposures (model$active_predictors is not "integer(0)"), please select the model with lowest DIC
mbsp.model_bmi_FT12 <- #none model based on FT12 select any predictor, so i left blank here, please put your selected model here
productCI_bmi_FT12 <- mbsp.model_bmi_FT12$B_CI_lower * mbsp.model_bmi_FT12$B_CI_upper

# Display the active predictors (correspond to row indices in B)
colnames(exposure)[mbsp.model_bmi_FT12$active_predictors]
list_active_pred_bmi_FT12 <- which(apply((!productCI_bmi_FT12<0),1,any))
for (i in list_active_pred_bmi_FT12)
{
  cols_active_pred_bmi_FT12 <- which(!productCI_bmi_FT12[i,] < 0)
  mbsp.model_bmi_FT12$B_est[i,cols_active_pred_bmi_FT12] <- 0
}

# Set to zero the coefficients of non selected predictors according to CIs
mbsp.model_bmi_FT12$B_est[-list_active_pred_bmi_FT12,] <- matrix(0,dim(mbsp.model_bmi_FT12$B_est[-list_active_pred_bmi_FT12,])[1],dim(mbsp.model_bmi_FT12$B_est[-list_active_pred_bmi_FT12,])[2])
mbsp.model_bmi_FT12$B_est
rownames(mbsp.model_bmi_FT12$B_est) <- colnames(exposure)
colnames(mbsp.model_bmi_FT12$B_est) <- colnames(protein_scaled)

#-- generate table for ploting --#
# Function for identifying the breaks for generating N equal size groups of beta values
quantile_breaks <- function(xs, n = 8) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks_neg_bmi_FT12 <- quantile_breaks(mbsp.model_bmi_FT12$B_est[which(mbsp.model_bmi_FT12$B_est<0)], n = 3)
mat_breaks_pos_bmi_FT12 <- quantile_breaks(mbsp.model_bmi_FT12$B_est[which(mbsp.model_bmi_FT12$B_est>0)], n = 3)
mat_breaks_bmi_FT12 <- c(mat_breaks_neg_bmi_FT12,-0.00000000000000000000001,0.00000000000000000000001,mat_breaks_pos_bmi_FT12)
table(cut(as.numeric(mbsp.model_bmi_FT12$B_est), breaks=mat_breaks_bmi_FT12))

# Keep only selected predictors
sel_exp_bmi_FT12 <- rownames(mbsp.model_bmi_FT12$B_est)[list_active_pred_bmi_FT12]
mbsp.model_filter_bmi_FT12 <- as.data.frame(mbsp.model$B_est[sel_exp_bmi_FT12,])


###Walnut could save this dataframe to me and send to me
write.csv(mbsp.model_filter_bmi_FT12 , "//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/mbsp.model_filter_bmi_FT12.csv")
#-------------------------------------------------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------------------------------------------------------------------#
#Now, we will try MBSP models with 9 different sets of parameter and adjusted confounders without BMI
#-------------------------------------------------------------------------------------------------------------------------------------------#
#-------- HORSESHOES PRIOR u=0.5, a=0.5
#-- Try tau=1e-8 (very small: 10^{-8})
mbsp_model_nobmi_FT12_1 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=0.5, a=0.5, tau=1e-8, model_criteria=TRUE)
# variables selected by the MBSP model with horseshoes prior and when tau is very small
mbsp_model_nobmi_FT12_1$active_predictors

#-- Try tau=1e4 (very large: 10^4)
mbsp_model_nobmi_FT12_2 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=0.5, a=0.5, tau=1e4, model_criteria=TRUE)
# variables selected by the MBSP model with horseshoes prior and when tau is very large
mbsp_model_nobmi_FT12_2$active_predictors

#-------- NEG priors u=1, a=0.1
#-- Try tau=1e-8 (very small: 10^{-8})
mbsp_model_nobmi_FT12_3 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=1, a=0.1, tau=1e-8, model_criteria=TRUE)
# variables selected by the MBSP model with NEG prior and when tau is very small
mbsp_model_nobmi_FT12_3$active_predictors

#-- Try tau=1e4 (very large: 10^4)
mbsp_model_nobmi_FT12_4 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=1, a=0.1, tau=1e4, model_criteria=TRUE)
# variables selected by the MBSP model with NEG prior and  when tau is very large
mbsp_model_nobmi_FT12_4$active_predictors

#-------- Strawderman-berger priors (u=1, a=0.5)
#-- Try tau=1e-8 (very small: 10^{-8})
mbsp_model_nobmi_FT12_5 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=1, a=0.5, tau=1, model_criteria=TRUE)
# variables selected by the MBSP model with trawderman-berger priors and when tau is very small
mbsp_model_nobmi_FT12_5$active_predictors

#-- Try tau=1e4 (very large: 10^4)
mbsp_model_nobmi_FT12_6 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=1, a=0.5, tau=1e4, model_criteria=TRUE)
# variables selected by the MBSP model with trawderman-berger priors and when tau is very large
mbsp_model_nobmi_FT12_6$active_predictors

#-------- Use 'default' tau of tau=1/(p*sqrt(n*log(n)) but change u and a
# Corresponds to NEG prior tau=1/(p*sqrt(n*log(n)) but change u and a
mbsp_model_nobmi_FT12_7 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=1, a=0.1, model_criteria=TRUE)
# variables selected by the MBSP model when NEG prior is used. 
# Gives a less sparse model, due to NEG having lighter tails. 
mbsp_model_nobmi_FT12_7$active_predictors

#-------- Corresponds to horseshoe prior with default tau=1/(p*sqrt(n*log(n))
mbsp_model_nobmi_FT12_8 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=0.5, a=0.5, model_criteria=TRUE)
# variables selected by the MBSP model
mbsp_model_nobmi_FT12_8$active_predictors

#-------- Corresponds to Strawderman-berger prior (u=1, a=0.5)  with default tau=1/(p*sqrt(n*log(n))
mbsp_model_nobmi_FT12_9 <- MBSP(Y=protein,X=exposure,confounders=confounder_withoutbmi, u=1, a=0.5, model_criteria=TRUE)
# variables selected by the MBSP model
mbsp_model_nobmi_FT12_9$active_predictors

#--------Now, lets sreen the DIC and WAIC of each model, there is no cutoff for these two fit indeces
mbsp_model_nobmi_FT12_1$DIC # HORSESHOES PRIOR u=0.5, a=0.5  (very small: 10^{-8})
mbsp_model_nobmi_FT12_2$DIC # HORSESHOES PRIOR u=0.5, a=0.5 (very large: 10^4)
mbsp_model_nobmi_FT12_3$DIC # NEG prior u=1, a=0.1  (very small: 10^{-8})
mbsp_model_nobmi_FT12_4$DIC # NEG prior u=1, a=0.1 (very large: 10^4)
mbsp_model_nobmi_FT12_5$DIC # Strawderman-berger prior (u=1, a=0.5)  (very small: 10^{-8})
mbsp_model_nobmi_FT12_6$DIC # Strawderman-berger prior (u=1, a=0.5) (very large: 10^4)
mbsp_model_nobmi_FT12_7$DIC # NEG prior tau=1/(p*sqrt(n*log(n)) 
mbsp_model_nobmi_FT12_8$DIC # horseshoe prior with default tau=1/(p*sqrt(n*log(n))
mbsp_model_nobmi_FT12_9$DIC # Strawderman-berger prior (u=1, a=0.5)  with default tau=1/(p*sqrt(n*log(n))

#--------based on the model which have capture any exposures (model$active_predictors is not "integer(0)"), please select the model with lowest DIC
mbsp.model_nobmi_FT12 <- #none model based on FT12 select any predictor, so i left blank here, please put your selected model here
  productCI_nobmi_FT12 <- mbsp.model_nobmi_FT12$B_CI_lower * mbsp.model_nobmi_FT12$B_CI_upper

# Display the active predictors (correspond to row indices in B)
colnames(exposure)[mbsp.model_nobmi_FT12$active_predictors]
list_active_pred_nobmi_FT12 <- which(apply((!productCI_nobmi_FT12<0),1,any))
for (i in list_active_pred_nobmi_FT12)
{
  cols_active_pred_nobmi_FT12 <- which(!productCI_nobmi_FT12[i,] < 0)
  mbsp.model_nobmi_FT12$B_est[i,cols_active_pred_nobmi_FT12] <- 0
}

# Set to zero the coefficients of non selected predictors according to CIs
mbsp.model_nobmi_FT12$B_est[-list_active_pred_nobmi_FT12,] <- matrix(0,dim(mbsp.model_nobmi_FT12$B_est[-list_active_pred_nobmi_FT12,])[1],dim(mbsp.model_nobmi_FT12$B_est[-list_active_pred_nobmi_FT12,])[2])
mbsp.model_nobmi_FT12$B_est
rownames(mbsp.model_nobmi_FT12$B_est) <- colnames(exposure)
colnames(mbsp.model_nobmi_FT12$B_est) <- colnames(protein_scaled)

#-- Plotting outputs: Heatmap of beta coefficients --#
# Function for identifying the breaks for generating N equal size groups of beta values
quantile_breaks <- function(xs, n = 8) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks_neg_nobmi_FT12 <- quantile_breaks(mbsp.model_nobmi_FT12$B_est[which(mbsp.model_nobmi_FT12$B_est<0)], n = 3)
mat_breaks_pos_nobmi_FT12 <- quantile_breaks(mbsp.model_nobmi_FT12$B_est[which(mbsp.model_nobmi_FT12$B_est>0)], n = 3)
mat_breaks_nobmi_FT12 <- c(mat_breaks_neg_nobmi_FT12,-0.00000000000000000000001,0.00000000000000000000001,mat_breaks_pos_nobmi_FT12)
table(cut(as.numeric(mbsp.model_nobmi_FT12$B_est), breaks=mat_breaks_nobmi_FT12))

# Keep only selected predictors
sel_exp_nobmi_FT12 <- rownames(mbsp.model_nobmi_FT12$B_est)[list_active_pred_nobmi_FT12]
mbsp.model_filter_nobmi_FT12 <- as.data.frame(mbsp.model$B_est[sel_exp_nobmi_FT12,])


###Walnut could save this dataframe to me and send to me
write.csv(mbsp.model_filter_nobmi_FT12 , "//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/mbsp.model_filter_nobmi_FT12.csv")
#-------------------------------------------------------------------------------------------------------------------------#

save.image("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/FT12_mbsp.Rdata") 
##save your results, and dont upload this, it contain your sensitive information. It will save your time not to perform again


