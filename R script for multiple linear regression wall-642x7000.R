
# This is the Multiple Regression for comp*conversion in the BTx642*Tx7000 mapping data
# Written Feb, 2014 by Robert Anderson for the BTx642*Tx7000 comp*conversion publication

# Set directory, load MLR library
setwd("C:/Users/Robert/Google Drive/Mullet Lab/Projects/QTL Mapping projects/BTx642 x Tx7000")
getwd()
library("car", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")

# Optional: Input standardization, as in Gelman, 2008.
library("arm", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")
# --Note: To use unstandardized values, simply remove "standardize()" from the GlucFIT and PentFIT equations.

# Optional: Cross validation of the regression model
#install.packages("cvTools")
library("cvTools", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")


multiLR <- read.delim("C:/Users/Robert/Google Drive/Mullet Lab/Projects/QTL Mapping projects/BTx642 x Tx7000/642x7000_compxconver_R-input.txt")
View(multiLR)
GlucFIT = standardize(lm(yield.glucose ~ stem.length +  int.diameter +  wall.density + percent.wall +  percent.lignin.in.wall + grams.cellulose + pca.pellet + fa.pellet, data=multiLR))
sum.GlucFIT = summary(GlucFIT)
sum.GlucFIT
PentFIT = standardize(lm(yield.pentose ~ stem.length +  int.diameter +  wall.density +  + percent.wall + percent.lignin.in.wall + grams.xylan + pca.pellet + fa.pellet, data=multiLR))
sum.PentFIT = summary(PentFIT)
sum.PentFIT

# -----( test whether we can remove certain variables from the model )---------

GlucFIT.sub1 = standardize(lm(yield.glucose ~ dtf + stem.length +  int.diameter +  wall.density +  percent.wall +  percent.lignin.in.wall + grams.cellulose + pca.pellet + fa.pellet, data=multiLR))
anova(GlucFIT.sub1,GlucFIT)
# Above shows that _.

PentFIT.sub1 = standardize(lm(yield.pentose ~ dtf + stem.length +  int.diameter +  wall.density +  percent.wall +  percent.lignin.in.wall + grams.xylan + pca.pellet + fa.pellet, data=multiLR))
anova(PentFIT.sub1,PentFIT)
# Above shows that _.


# ------( Regression Diagnostics )---------------------

GlucFIT.resid = residuals(GlucFIT)
PentFIT.resid = residuals(PentFIT)

Gluc.resid.table = as.data.frame(GlucFIT.resid)
Pent.resid.table = as.data.frame(PentFIT.resid)

qqnorm(GlucFIT.resid, ylab = "Glucose yield MLR Residuals")
qqline(GlucFIT.resid)

qqnorm(PentFIT.resid, ylab = "Pentose Yield MLR Residuals")
qqline(PentFIT.resid)

shapiro.test(GlucFIT.resid)
shapiro.test(PentFIT.resid)

hist(GlucFIT.resid, xlab = "Glucose yield MLR residuals", main = " ")
hist(PentFIT.resid, xlab = "Pentose yield MLR residuals", main = " ")

# ---( Variance-inflation factors)

GlucFIT.vif <- vif(standardize(lm(yield.glucose ~ stem.length +  int.diameter +  wall.density + percent.wall +	percent.lignin.in.wall + grams.cellulose + pca.pellet + fa.pellet, data=multiLR)))
GlucFIT.vif
PentFIT.vif <- vif(standardize(lm(yield.pentose ~ stem.length +  int.diameter +  wall.density +	percent.wall + percent.lignin.in.wall + grams.xylan + pca.pellet + fa.pellet, data=multiLR)))
PentFIT.vif

GlucFIT.dtf.vif <- vif(standardize(lm(yield.glucose ~ dtf + stem.length +  int.diameter +  wall.density + percent.wall + percent.lignin.in.wall + grams.cellulose + pca.pellet + fa.pellet, data=multiLR)))
GlucFIT.dtf.vif
PentFIT.dtf.vif <- vif(standardize(lm(yield.pentose ~ dtf + stem.length +  int.diameter +  wall.density + percent.wall + percent.lignin.in.wall + grams.xylan + pca.pellet + fa.pellet, data=multiLR)))
PentFIT.dtf.vif

#  *Note,  tolerance = 1/vif  A tolerance<0.01 (vif>100) means trait should be removed from model


# ---( RMSE cross validation )

# Removing standardization, since this affected variable length (NULL values)
GlucFIT2 = lm(yield.glucose ~ dtf + stem.length +  int.diameter +  wall.density + percent.wall + percent.lignin.in.wall + grams.cellulose + pca.pellet + fa.pellet, data=multiLR)
PentFIT2 = lm(yield.pentose ~ dtf + stem.length +  int.diameter +  wall.density + percent.wall + percent.lignin.in.wall + grams.xylan + pca.pellet + fa.pellet, data=multiLR)

Gluc.CV = cvFit(GlucFIT2, data = multiLR, y = multiLR$yield.glucose, cost = rmspe, K = 5, R = 10)
Gluc.CV
Pent.CV = cvFit(PentFIT2, data = multiLR, y = multiLR$yield.pentose, cost = rmspe, K = 5, R = 10)
Pent.CV

# Values for the full model
Gluc.CV.full = rmspe(multiLR$yield.glucose, predict(GlucFIT2), includeSE=TRUE)
Gluc.CV.full
Pent.CV.full = rmspe(multiLR$yield.pentose, predict(PentFIT2), includeSE=TRUE)
Pent.CV.full

#library("DAAG", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")
#CVlm(df=multiLR, form.lm=formula(GlucFIT2), m=5, plotit="Observed")

#install.packages("caret", dependencies = c("Depends", "Suggests"))
#library("caret", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")
#library("mlbench", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")
#dat <- data.frame(multiLR)
#set.seed(107)


