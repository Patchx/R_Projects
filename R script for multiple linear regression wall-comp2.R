
# This is the Multiple Regression for comp*conversion in the E-SAP data
# Written Feb, 2014 by Robert Anderson for the E-SAP comp*conversion publication

# Set directory, load MLR library
setwd("C:/Users/Robert/Google Drive/Mullet Lab/Projects/QTL Mapping projects/DGA - ESP panel/ESP - 2012")
getwd()
library("car", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")

# Optional: Input standardization, as in Gelman, 2008.
library("arm", lib.loc="C:/Revolution/R-Enterprise-6.1/R-2.14.2/library")
# --Note: To use unstandardized values, simply remove "standardize()" from the GlucFIT and PentFIT equations.

multiLR <- read.delim("C:/Users/Robert/Google Drive/Mullet Lab/Projects/QTL Mapping projects/DGA - ESP panel/ESP - 2012/ESP-multilinear-vegeta-for_R.txt")
View(multiLR)
GlucFIT = standardize(lm(Yield.Glucose ~ Cell.wall.density + Lignin.in.wall +  Galact.Arabinan.in.Xylan +	Conversion.weight +	grams.cellulose +	X3INT.Diameter +	X3INT.Length, data=multiLR))
sum.GlucFIT = summary(GlucFIT)
sum.GlucFIT
PentFIT = standardize(lm(Yield.Pentose ~ Cell.wall.density + Lignin.in.wall +  Galact.Arabinan.in.Xylan +  Conversion.weight +  grams.Xylan +	X3INT.Diameter +	X3INT.Length, data=multiLR))
sum.PentFIT = summary(PentFIT)
sum.PentFIT

# -----( test whether we can remove certain variables from the model )---------

GlucFIT.sub1 = lm(Yield.Glucose ~ Cell.wall.density + Lignin.in.wall +  Galact.Arabinan.in.Xylan +  Conversion.weight +	grams.cellulose +	X3INT.Diameter, data=multiLR)
anova(GlucFIT.sub1,GlucFIT)
# Above shows that int length is important to the Gluc model, apart from density effects.

GlucFIT.sub2 = lm()
anova(GlucFIT.sub2,GlucFIT)
# Above shows that int diameter is important to the Gluc model, apart from density effects.

PentFIT.sub1 = lm()
anova(PentFIT.sub1,PentFIT)
# Above shows that int length and diameter can be removed from the Pent model.

GlucFIT.sub3 = lm()
anova(GlucFIT.sub3,GlucFIT)
# Above shows that cell/xylan/lignin in the wall can be removed from the model as long as the grams of cellulose or xylan are included in the Gluc model.

PentFIT.sub2 = lm()
anova(PentFIT.sub2,PentFIT)
# Above shows that cell/xylan/lignin in the wall can be removed from the model as long as the grams of cellulose or xylan are included in the Pent model.

GlucFIT.sub4 = lm()
anova(GlucFIT.sub4,GlucFIT)
# Above shows that wall density affects conversion, apart from the relationship with INT length and diameter.
  
PentFIT.sub3 = lm()
anova(PentFIT.sub3,PentFIT)
# Above shows that wall density affects conversion, apart from the relationship with INT length and diameter.

GlucFIT.sub5 = lm()
anova(GlucFIT.sub5,GlucFIT)
# Above shows the GlucFIT partial F-test for removing all traits that have vif > 10, except Lignin-in-wall, cellulose-grams and xylan-grams.


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

GlucFIT.vif <- vif(lm(Yield.Glucose ~ Cell.wall.density + Lignin.in.wall +  Galact.Arabinan.in.Xylan +	Conversion.weight +	grams.cellulose +	X3INT.Diameter +	X3INT.Length, data=multiLR))
GlucFIT.vif
PentFIT.vif <- vif(lm(Yield.Pentose ~ Cell.wall.density + Lignin.in.wall +  Galact.Arabinan.in.Xylan +  Conversion.weight +  grams.Xylan +	X3INT.Diameter +	X3INT.Length, data=multiLR))
PentFIT.vif


#  *Note,  tolerance = 1/vif  A tolerance<0.01 (vif>100) means trait should be removed from model

shapiro.test(GlucFIT.resid)
shapiro.test(PentFIT.resid)

