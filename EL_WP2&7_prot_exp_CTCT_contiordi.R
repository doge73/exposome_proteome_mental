rm(list=ls())
library(haven)
library(OpenMx)
library(psych); library(polycor)
source("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/miFunctions.R")
set.seed(52420242)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA
# Load Data
imputated_twin_CTCT <- read_dta("//ad.helsinki.fi/home/z/zhiywang/Desktop/SDQ_exp_protein/imputated_FT12_CTCT.dta")
#######CTCT for one continuous and one ordinal variables############################################

imputated_twin_CTCT$plyVideogame_ch_cr1<-as.factor(imputated_twin_CTCT$plyVideogame_ch_cr1)
imputated_twin_CTCT$plyVideogame_ch_cr2<-as.factor(imputated_twin_CTCT$plyVideogame_ch_cr2)

vars <- c('P05155') # list of continuous variables names
nvc <- 1 # number of continuous variables
ntvc <- nvc*2 # number of total continuous variables
conVars <- paste(vars,c(rep(1,nvc),rep(2,nvc)),sep="")
# Select Ordinal Variables
nth <- 3 # number of thresholds
vars <- c('plyVideogame_ch_cr') # list of ordinal variables names
nvo <- 1 # number of ordinal variables
ntvo <- nvo*2 # number of total ordinal variables
ordVars <- paste(vars,c(rep(1,nvo),rep(2,nvo)),sep="")
ordData <- imputated_twin_CTCT
wtquant <- quantile(ordData[,c('plyVideogame_ch_cr1','plyVideogame_ch_cr2')],(0:(nth+1))/(nth+1),na.rm=TRUE)
###for (i in c('plyVideogame_ch_cr1','plyVideogame_ch_cr2')) { ordData[,i] <- cut(ordData[,i], breaks=wtquant, labels=c(0:nth)) }
# Select Variables for Analysis
vars <- c('P05155','plyVideogame_ch_cr') # list of variables names
nv <- nvc+nvo # number of variables
ntv <- nv*2 # number of total variables
selVars <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")
# Select Data for Analysis
mzData <- subset(ordData, ZYG1==1, selVars)
dzData <- subset(ordData, ZYG1==2, selVars)
mzDataF <- cbind(mzData[,conVars],mxFactor( x=mzData[,ordVars], levels=c(1:4)) )
dzDataF <- cbind(dzData[,conVars],mxFactor( x=dzData[,ordVars], levels=c(1:4)) )

# Set Starting Values
frMV <- c(TRUE,FALSE) # free status for variables
frCvD <- valDiagLU(frMV,T,T,ntv) # free status for diagonal, lower & upper elements of covariance matrix
frCv <- matrix(as.logical(frCvD),4)
svMe <- c(15,0) # start value for means
svVa <- c(.5,1) # start value for variances
lbVa <- .0001 # lower bound for variances
svLTh <- 0 # start value for first threshold
svITh <- 1 # start value for increments
svTh <- matrix(rep(c(svLTh,(rep(svITh,nth-1)))),nrow=nth,ncol=nv) # start value for thresholds
lbTh <- matrix(rep(c(-3,(rep(0.001,nth-1))),nv),nrow=nth,ncol=nv) # lower bound for thresholds

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL
# Create Algebra for expected Mean & Threshold Matrices
meanMZ <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=frMV, values=svMe, labels=c('mean1', 'mean2', 'mean1', 'mean2'), name="meanMZ" )
meanDZ <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=frMV, values=svMe, labels=c('mean1', 'mean2', 'mean1', 'mean2'), name="meanDZ" )
thinMZ <- mxMatrix( type="Full", nrow=nth, ncol=ntvo, free=TRUE, values=svTh, lbound=lbTh, labels=c('thre1', 'thre2', 'thre3'), name="thinMZ" )
thinDZ <- mxMatrix( type="Full", nrow=nth, ncol=ntvo, free=TRUE, values=svTh, lbound=lbTh, labels=c('thre1', 'thre2', 'thre3'), name="thinDZ" )
inc <- mxMatrix( type="Lower", nrow=nth, ncol=nth, free=FALSE, values=1, name="inc" )
threMZ <- mxAlgebra( expression= inc %*% thinMZ, name="threMZ" )
threDZ <- mxAlgebra( expression= inc %*% thinDZ, name="threDZ" )
# Create Algebra for expected Covariance Matrices
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
dataMZ <- mxData( observed=mzDataF, type="raw" )
dataDZ <- mxData( observed=dzDataF, type="raw" )
# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal( covariance="covMZ", means="meanMZ", dimnames=selVars, thresholds="threMZ", threshnames=ordVars )
expDZ <- mxExpectationNormal( covariance="covDZ", means="meanDZ", dimnames=selVars, thresholds="threDZ", threshnames=ordVars )
funML <- mxFitFunctionML()
# Create Model Objects for Multiple Groups
modelMZ <- mxModel( meanMZ, covMZ, thinMZ, inc, threMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ <- mxModel( meanDZ, covDZ, thinDZ, inc, threDZ, dataDZ, expDZ, funML, name="DZ" )
multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
# Create Confidence Interval Objects
ciCor <- mxCI( c('MZ.covMZ','DZ.covDZ' ))
ciThre <- mxCI( c('MZ.threMZ','DZ.threDZ' ))
# Build Saturated Model with Confidence Intervals
modelSAT <- mxModel( "twoSATj", modelMZ, modelDZ, multi, ciCor, ciThre )
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
########################################################################################

