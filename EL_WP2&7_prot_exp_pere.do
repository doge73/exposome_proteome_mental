use "\\ad.helsinki.fi\home\z\zhiywang\Desktop\SDQ_exp_protein\imputated_FT12_dataset.dta", clear
//import the imputated dataset contain all the information 
//Stata also allowed other file types


merge 1:1 person_nb using "\\ad.helsinki.fi\home\z\zhiywang\Desktop\SDQ_exp_protein\pfactor_age14_FT12.dta" 
// import and merge the dataset with p-factor information (if your p-factor dataset is not in dta file, please convert to dta at first)
drop if _merge==2
drop _merge
merge 1:1 person_nb using "\\ad.helsinki.fi\home\z\zhiywang\Desktop\SDQ_exp_protein\MPNI_interview_age.dta", keepusing (int14ag2)
drop if _merge==2

drop if p_factor ==.
//i drop participant without p-factor information. If there is no missing in WALNUTs of P-factor, just skip
//check the protocol how to clean the p-factor variable


replace int14ag2 =14.15 if int14ag2 ==. //imputed by median

///pecent of excess risk explained (pere) in FinnTwin12 by generalized linear regression, there is five exposure-protein sets 

glm p_factor O95445 i.alc_ch_cr sex int14ag2 , robust //exposure only model
//if exposure is a categorical variables, add i. as prefix
glm p_factor O95445  sex int14ag2 , robust
glm p_factor i.alc_ch_cr sex int14ag2 , robust

glm p_factor P0C0L5 hh_adultsnr_pr sex int14ag2 , robust
glm p_factor P0C0L5 sex int14ag2 , robust
glm p_factor hh_adultsnr_pr sex int14ag2 , robust

glm p_factor P05155 i.plyVideogame_ch_cr sex int14ag2 , robust
glm p_factor P05155 sex int14ag2 , robust
glm p_factor i.plyVideogame_ch_cr sex int14ag2 , robust

glm p_factor P40189 ndvi_5yrs_greenest_500 sex int14ag2 , robust
glm p_factor P40189 sex int14ag2 , robust
glm p_factor ndvi_5yrs_greenest_500 sex int14ag2 , robust



///pecent of excess risk explained (pere) in WALNUTs by generalized linear regression, there is three exposure-protein sets 
glm p_factor P28799 treecover_300 sex age, robust   // exposure and protein model
glm p_factor P28799  sex age, robust                // protein only model
glm p_factor treecover_300 sex age, robust          // exposure only model

glm p_factor P28799 treecover_500 sex age, robust   // exposure and protein model
glm p_factor P28799  sex age, robust                // protein only model
glm p_factor treecover_500 sex age, robust          // exposure only model

glm p_factor Q99983 i.Occup_f_pr sex age, robust    // exposure and protein model
glm p_factor Q99983  sex age, robust                // protein only model
glm p_factor i.Occup_f_pr sex age, robust           // exposure only model
//if exposure is a categorical variables, add i. as prefix
