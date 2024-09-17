use "\\ad.helsinki.fi\home\z\zhiywang\Desktop\SDQ_exp_protein\imputated_FT12_dataset.dta", clear


merge 1:1 person_nb using "\\ad.helsinki.fi\home\z\zhiywang\Desktop\SDQ_exp_protein\FT exposures\harm_ft12_exp.dta", keepusing (fam_nb)
drop if _merge==2
drop _merge


xtset fam_nb

xtreg O95445 sex age alc_ch_cr,fe 
xtreg P05155 sex age i.plyVideogame,fe
xtreg P40189 ndvi_5yrs_greenest_500 sex age,fe


glm O95445 sex age alc_ch_cr, robust
glm P05155 sex age i.plyVideogame, robust
glm P40189 ndvi_5yrs_greenest_500 sex age, robust
