
# This script was written to perform Multiple Linear Regression on internode morphology data.
# Robert Anderson - February, 2014.

#Setting the working directory and installing the required packages
setwd("C:/Users/Robert/Google Drive/Mullet Lab/Projects/Stem Growth Projects/PhyB & Diameter Growth")
getwd()
library("car", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")

# Optional: Cross validation of the regression model
#install.packages("cvTools")
library("cvTools", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")


multiLR <- read.delim("C:/Users/Robert/Google Drive/Mullet Lab/Projects/Stem Growth Projects/PhyB & Diameter Growth/MLR-input_100M-regress.txt")
View(multiLR)

VolumeFIT = lm(log.INT.volume ~ log.cell.height + log.trans.area +	log.num.long +	log.num.trans, data=multiLR)
sum.VolumeFIT = summary(VolumeFIT)
sum.VolumeFIT

LengthFIT = lm(log.INT.length ~ log.cell.height + log.trans.area +  log.num.long +	log.num.trans, data=multiLR)
sum.LengthFIT = summary(LengthFIT)
sum.LengthFIT

DiameterFIT = lm(log.INT.diameter ~ log.cell.height + log.trans.area +  log.num.long +  log.num.trans, data=multiLR)
sum.DiameterFIT = summary(DiameterFIT)
sum.DiameterFIT


# -----( test whether we can remove certain variables from the model )---------

VolumeFIT.sub1 = lm(log.INT.volume ~ log.trans.area +  log.num.long +	log.num.trans, data=multiLR)
anova(VolumeFIT.sub1,VolumeFIT)
# Above shows that cell height can be dropped from this dataset.

VolumeFIT.sub2 = lm(log.INT.volume ~ log.cell.height +  log.num.long +  log.num.trans, data=multiLR)
anova(VolumeFIT.sub2,VolumeFIT)
# Above shows that cell transverse area cannot be dropped from this dataset.

LengthFIT.sub1 = lm(log.INT.length ~ log.cell.height +  log.num.long, data=multiLR)
anova(LengthFIT.sub1,LengthFIT)
# Above shows that cell trans area and cell trans number can be dropped from this dataset.

DiameterFIT.sub1 = lm(log.INT.diameter ~ log.trans.area +  log.num.long +  log.num.trans, data=multiLR)
anova(DiameterFIT.sub1,DiameterFIT)
# Above shows that cell height can be dropped from this dataset.

DiameterFIT.sub2 = lm(log.INT.diameter ~ log.cell.height + log.trans.area + log.num.trans, data=multiLR)
anova(DiameterFIT.sub2,DiameterFIT)
# Above shows that cell number in the longitudinal plane cannot be dropped from this dataset.


# ------( Regression Diagnostics )---------------------

VolumeFIT.resid = residuals(VolumeFIT)
LengthFIT.resid = residuals(LengthFIT)
DiameterFIT.resid = residuals(DiameterFIT)

Volume.resid.table = as.data.frame(VolumeFIT.resid)
Length.resid.table = as.data.frame(LengthFIT.resid)
Diameter.resid.table = as.data.frame(DiameterFIT.resid)

qqnorm(VolumeFIT.resid, ylab = "100M & SM100  INT Volume MLR Residuals")
qqline(VolumeFIT.resid)

qqnorm(LengthFIT.resid, ylab = "100M & SM100  INT Length MLR Residuals")
qqline(LengthFIT.resid)

qqnorm(DiameterFIT.resid, ylab = "100M & SM100  INT Diameter MLR Residuals")
qqline(DiameterFIT.resid)

shapiro.test(VolumeFIT.resid)
shapiro.test(LengthFIT.resid)
shapiro.test(DiameterFIT.resid)

hist(VolumeFIT.resid, xlab = "100M & SM100  INT Volume MLR residuals", main = " ")
hist(LengthFIT.resid, xlab = "100M & SM100  INT Length MLR residuals", main = " ")
hist(DiameterFIT.resid, xlab = "100M & SM100  INT Diameter MLR residuals", main = " ")

# ---( Variance-inflation factors)

Volume.vif <- vif(lm(log.INT.volume ~ log.cell.height + log.trans.area +  log.num.long +	log.num.trans, data=multiLR))
Volume.vif
Length.vif <- vif(lm(log.INT.length ~ log.cell.height + log.trans.area +  log.num.long +  log.num.trans, data=multiLR))
Length.vif
Diameter.vif <- vif(lm(log.INT.diameter ~ log.cell.height + log.trans.area +  log.num.long +  log.num.trans, data=multiLR))
Diameter.vif

#  *Note,  tolerance = 1/vif  A tolerance<0.01 (vif>100) means trait should be removed from model


# ---( RMSE cross validation )

volume.CV = cvFit(VolumeFIT, data = multiLR, y = multiLR$log.INT.volume, cost = rmspe, K = 5, R = 10)
volume.CV
length.CV = cvFit(LengthFIT, data = multiLR, y = multiLR$log.INT.length, cost = rmspe, K = 5, R = 10)
length.CV
diameter.CV = cvFit(DiameterFIT, data = multiLR, y = multiLR$log.INT.diameter, cost = rmspe, K = 5, R = 10)
diameter.CV


# Values for the full model
volume.CV.full = rmspe(multiLR$log.INT.volume, predict(VolumeFIT), includeSE=TRUE)
volume.CV.full
length.CV.full = rmspe(multiLR$log.INT.length, predict(LengthFIT), includeSE=TRUE)
length.CV.full
diameter.CV.full = rmspe(multiLR$log.INT.diameter, predict(DiameterFIT), includeSE=TRUE)
diameter.CV.full

#library("DAAG", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")
#CVlm(df=multiLR, form.lm=formula(GlucFIT2), m=5, plotit="Observed")

#install.packages("caret", dependencies = c("Depends", "Suggests"))
#library("caret", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")
#library("mlbench", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")
#dat <- data.frame(multiLR)
#set.seed(107)


