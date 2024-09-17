##title: ExWAS for WP2&7 protein-exposure-pfactor study in Equal-Life
##Coder: Zhiyang Wang
##Date: 04.30.2024

library(rexposome)
library(readr)
library(foreign)
library(haven)

set.seed(04292023) # change whatever data

##Here is the instruction of the rexposome package: https://www.bioconductor.org/packages/release/bioc/vignettes/rexposome/inst/doc/exposome_data_analysis.html
##set for exposome dataset
##The merged dataset for protein, exposures and covariate is already processed and cleaned in Stata
## change FT12 to WALNUT in coding for better smoothing usage. 

#---------data preparation----------------------------------------------------------------------------------------------------------------#
description_FT12 <- read.csv("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/description.csv", header = TRUE)
exposure_FT12 <- read.csv("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/exposure.csv", header = TRUE)
phenotype_FT12 <- read.csv("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/phenotype.csv", header = TRUE)
phenotype_FT12_2 <- read.csv("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/phenotype.csv", header = TRUE)

rownames(description_FT12) <- description_FT12[, 1]
description_FT12 <- description_FT12[, -1]
rownames(phenotype_FT12) <- phenotype_FT12[, 1]
phenotype_FT12 <- phenotype_FT12[, -1]
rownames(exposure_FT12) <- exposure_FT12[, 1]
exposure_FT12 <- exposure_FT12[, -1]

exp_FT12 <- loadExposome(exposures = exposure_FT12, description = description_FT12, phenotype = phenotype_FT12,
                         description.famCol = "Family",exposures.asFactor = 5)
exp_FT12@featureData@data[[".type"]][24]<-"numeric"
##manually set variable "hh_adultspr_pr_11" as numeric, it is the 26th variable in the alphabet order (A-Z).
## If in WALNUTs, hh_adultspr_pr_11 is more than 5 level, you could just delete this line of code

##imputation
exp_FT12 <- imputation(exp_FT12, messages = TRUE)

##generate imputated exposure set
imputation_FT12 <- expos(exp_FT12)
imputation_FT12$person_nb<-rownames(imputation_FT12)
im_dataset_FT12 <- merge(phenotype_FT12_2, imputation_FT12 ,by = 'person_nb')

write.dta(im_dataset_FT12, "//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/imputated_FT12_dataset.dta") 
##generate for an imputated dataset, which will be ready for another analysis 
##(write.dta is for stata, you could choose write.excel, write.csv, etc, suit yourself)
#-------------------------------------------------------------------------------------------------------------------------#

#------------ ExWAS without BMI--------------------------------------------------------------------------------------------------------#

proteins<-c("P05155", 	"P01008", 	"Q9UK55", 	"P04180", 	"P02787", 	"O95445", 	"P08571", 	
            "P06681", 	"P00734", 	"P00751", 	"P00742", 	"P08185", 	"P10909", 	"P11279", 	
            "Q14126", 	"P09172", 	"P23142", 	"Q9Y6R7", 	"P07911", 	"Q86TH1", 	"P0C0L5", 	
            "P40189", 	"Q13449", 	"P28799", 	"Q12805", 	"Q9BY67")

#for WALNUTs, please change the list to c("P05155", 	"P01008", 	"P25786", 	"Q9UK55", 	"P04180", 	"P02787", 	"O95445", 	"P08571", 	"P06681", 	
##                                    "P51693", 	"P00734", 	"P00751", 	"P00742", 	"P08185", 	"P13611", 	"O95274", 	"P10909", 	"P11279", 	
##                                    "Q14126", 	"P09172", 	"P23142", 	"Q9Y6R7", 	"Q99983", 	"P07911", 	"Q86TH1", 	"P0C0L5", 	"P40189", 	
##                                    "Q13449", 	"P28799", 	"Q12805", 	"Q9BY67", 	"P07858")


predictor_nobmi<-c("age", "sex")

for (i in proteins){
  setwd("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/exwas_FT12_air/withoutBMI/")# create a file to store exwas results without adjusting bmi first
  formula <- formula(paste(i, "~", paste(predictor_nobmi, collapse = "+")))  
  
  assign(paste(i, "_nobmi_FT12_result", sep = ""), exwas(exp_FT12, formula = formula, family = "gaussian", robust = TRUE))
  result <- data.frame(extract(get(paste(i, "_nobmi_FT12_result", sep = ""))))
  result$exposure <- row.names(result)
  result$protein <- i
  file_name <- paste("FT12_result_nobmi_", i, ".csv",sep = "")
  write_csv(result, file_name, col_names = TRUE)}
Q9Y6R7_nobmi_FT12_result@effective # get your p-value for ExWAS without covariate of BMI
#-------------------------------------------------------------------------------------------------------------------------#




#------------ ExWAS with BMI--------------------------------------------------------------------------------------------------------#
im_dataset_FT12 <- read.dta("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/imputated_FT12_dataset.dta")
im_dataset_FT12_bmi<-na.omit(im_dataset_FT12) ##here is to create an exposome dataset after removing missing in BMI, if you dont have any missing in BMI, you could skep it
phenotype_bmi_FT12<-im_dataset_FT12_bmi[, c("person_nb", "sex", "age", "bmi_ya", 
                                            "P05155", 	"P01008", 	"Q9UK55", 	"P04180", 	"P02787", 	"O95445", 	"P08571", 	
                                            "P06681", 	"P00734", 	"P00751", 	"P00742", 	"P08185", 	"P10909", 	"P11279", 	
                                            "Q14126", 	"P09172", 	"P23142", 	"Q9Y6R7", 	"P07911", 	"Q86TH1", 	"P0C0L5", 	
                                            "P40189", 	"Q13449", 	"P28799", 	"Q12805", 	"Q9BY67")] ##one participant id, 3 covariates, and proteins were selected
exposure_bmi_FT12<-im_dataset_FT12_bmi[, -which(names(-im_dataset_FT12_bmi) %in% c("sex", "age", "bmi_ya", 
                                            "P05155", 	"P01008", 	"Q9UK55", 	"P04180", 	"P02787", 	"O95445", 	"P08571", 	
                                            "P06681", 	"P00734", 	"P00751", 	"P00742", 	"P08185", 	"P10909", 	"P11279", 	
                                            "Q14126", 	"P09172", 	"P23142", 	"Q9Y6R7", 	"P07911", 	"Q86TH1", 	"P0C0L5", 	
                                            "P40189", 	"Q13449", 	"P28799", 	"Q12805", 	"Q9BY67"))]

rownames(phenotype_bmi_FT12) <- phenotype_bmi_FT12[, 1]
phenotype_bmi_FT12 <- phenotype_bmi_FT12[, -1]
rownames(exposure_bmi_FT12) <- exposure_bmi_FT12[, 1]
exposure_bmi_FT12 <- exposure_bmi_FT12[, -1]

exp_bmi_FT12 <- loadExposome(exposures = exposure_bmi_FT12, description = description_FT12, phenotype = phenotype_bmi_FT12,
                         description.famCol = "Family",exposures.asFactor = 5)
exp_bmi_FT12@featureData@data[[".type"]][24]<-"numeric"

predictor<-c("age", "sex", "bmi_ya") 

for (i in proteins){
  setwd("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/exwas_FT12_air/withBMI/")# create a file to store exwas results adjusting bmi first
  formula <- formula(paste(i, "~", paste(predictor, collapse = "+")))  
  assign(paste(i, "_FT12_result", sep = ""), exwas(exp_bmi_FT12, formula = formula, family = "gaussian", robust = TRUE)) # robust SE for cluster RCT or twin design
  result <- data.frame(extract(get(paste(i, "_FT12_result", sep = ""))))
  result$exposure <- row.names(result)
  result$protein <- i
  file_name <- paste("FT12_result_", i, ".csv",sep = "")
  write_csv(result, file_name, col_names = TRUE)}
Q9Y6R7_FT12_result@effective # get your p-value for ExWAS including covariate of BMI
#-------------------------------------------------------------------------------------------------------------------------#


save.image("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/FT12_ExWAS.Rdata") ##save your results, and dont upload this, it contain your sensitive information
load("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/FT12_ExWAS.Rdata")



