rm(list=ls())
library(haven)
library(OpenMx)
library(psych); library(polycor)
set.seed(225242421)
# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA
# Load Data
imputated_twin_CTCT <- read_dta("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/imputated_FT12_CTCT.dta")
dim(imputated_twin_CTCT)
describe(imputated_twin_CTCT[,1:12], skew=F)

#######CTCT for both continuous variable############################################
# Select Variables for Analysis
vars <- c('ndvi_5yrs_greenest_500','P40189') # list of variables names
nv <- 2 # number of variables
ntv <- nv*2 # number of total variables
selVars <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")
# Select Data for Analysis
mzData <- subset(imputated_twin_CTCT, ZYG1==1, selVars)
dzData <- subset(imputated_twin_CTCT, ZYG1==2, selVars)
# Generate Descriptive Statistics
colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)
cov(mzData,use="complete")
cov(dzData,use="complete")
# Set Starting Values
svMe <- c(15,5) # start value for means
svVa <- c(.5,.8) # start value for variances
lbVa <- .0001 # lower bound for variances
# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL
# Create Algebra for expected Mean Matrices
meanMZ <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c('mean1', 'mean2', 'mean1', 'mean2'), name="meanMZ" )
meanDZ <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c('mean1', 'mean2', 'mean1', 'mean2'), name="meanDZ" )
# Create Algebra for expected Variance/Covariance Matrices
covMZ <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, 
                   values=c(.5, .3, .2, .2,
                            .3, .5, .2, .2,
                            .2, .2, .5, .3,
                            .2, .2, .3, .5), 
                   lbound=0,
                   labels=c('V1', 'Cov12', 'MZCov1', 'MZCross',
                            'Cov12', 'V2', 'MZCross', 'MZCov2',
                            'MZCov1', 'MZCross', 'V1', 'Cov12',
                            'MZCross', 'MZCov2', 'Cov12', 'V2'), name="covMZ", byrow = T )
covDZ <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, 
                   values=c(.5, .3, .1, .1,
                            .3, .5, .1, .1,
                            .1, .1, .5, .3,
                            .1, .1, .3, .5),  
                   lbound=0,
                   labels=c('V1', 'Cov12', 'DZCov1', 'DZCross',
                            'Cov12', 'V2', 'DZCross', 'DZCov2',
                            'DZCov1', 'DZCross', 'V1', 'Cov12',
                            'DZCross', 'DZCov2', 'Cov12', 'V2'), name="covDZ", byrow = T )
# Create Data Objects for Multiple Groups
dataMZ <- mxData( observed=mzData, type="raw" )
dataDZ <- mxData( observed=dzData, type="raw" )
# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal( covariance="covMZ", means="meanMZ", dimnames=selVars )
expDZ <- mxExpectationNormal( covariance="covDZ", means="meanDZ", dimnames=selVars )
funML <- mxFitFunctionML()
# Create Model Objects for Multiple Groups
modelMZ <- mxModel( meanMZ, covMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ <- mxModel( meanDZ, covDZ, dataDZ, expDZ, funML, name="DZ" )
multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
# Create Confidence Interval Objects
ciCov <- mxCI( c('MZ.covMZ','DZ.covDZ') )
ciMean <- mxCI( c('MZ.meanMZ','DZ.meanDZ') )
# Build Saturated Model with Confidence Intervals
modelSAT <- mxModel( "twoSATc", modelMZ, modelDZ, multi, ciCov, ciMean )
# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL
# Run Saturated Model
mxOption(NULL, "Default optimizer", "SLSQP")
fitSAT <- mxRun( modelSAT, intervals=F )
sumSAT <- summary( fitSAT )
sumSAT
# ----------------------------------------------------------------------------------------------------------------------
fitSAT$output$matrices$MZ.covMZ
fitSAT$output$matrices$DZ.covDZ

cov2cor(fitSAT$output$matrices$MZ.covMZ)
cov2cor(fitSAT$output$matrices$DZ.covDZ)
####################################################################################
