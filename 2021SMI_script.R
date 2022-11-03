getwd()

library(assertr)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(boot)
library(mvnormtest)
library(ggplot2)
library(magrittr)
library(data.table)
library(smatr)

theme_set(theme_bw())

###############################################################
### Load in our main data sets ################################
###############################################################

SMI_zoos <- read.csv("2021SMI.csv", header = TRUE)
head(SMI_zoos)
view(SMI_zoos) #contains wild frogs caught only in 2021
dim(SMI_zoos) #110 rows, 9 columns
SMI_zoos$EB <- as.factor(SMI_zoos$EB) #as factor. Two levels: 0=EB 1=not (OK)
SMI_zoos$pop <- as.factor(SMI_zoos$pop) #factor with 6 levels
SMI_zoos$locale <- as.factor(SMI_zoos$locale) #factor with 2 levels (Wild, Zoo)
#yr_birth <- as.Date(SMI_zoos$birth_yr, format = "%Y") # make new column or set birth_yr as numeric?
str(SMI_zoos)
summary(SMI_zoos)

mass_zoos <- read.csv("2021mass.csv", header = TRUE)
head(mass_zoos)
view(mass_zoos)
dim(mass_zoos) #64 rows, 6 columns
str(mass_zoos) #set EB as factor
mass_zoos$EB <- as.factor(mass_zoos$EB) #as factor. Two levels: 0=EB 1=not (OK)
summary(mass_zoos)

wild_SMI <- read.csv("FULLwildSMI.csv", header = TRUE)
head(wild_SMI)
view(wild_SMI) # data from all wild frogs caught from 2011-2022
dim(wild_SMI) #502 rows, 6 columns
wild_SMI$frog_id <- as.factor(wild_SMI$frog_id) #as factor. 490 levels (i.e. individual frogs)
wild_SMI$pop <- as.factor(wild_SMI$pop) #as factor. 5 levels: "CH" "MS" "MT" "MV" "ST"
str(wild_SMI)
summary(wild_SMI)

prepo <-read.csv("prepost_SMI.csv", header = TRUE) #have now added calculated svl (using linear regressions below)
head(prepo)
view(prepo)
dim(prepo) #87 rows, 12 columns
prepo$pop <- as.factor(prepo$pop)
prepo$EB <- as.factor(prepo$EB)
prepo$birth_yr <- as.factor(prepo$birth_yr) #9 levels
str(prepo)
summary(prepo)

######################################################################
### Do some visualization of datasets ################################
######################################################################
hist(SMI_zoos$age_yrs, main = "OSF age") #looking for bell curve to indicate normalcy of data
hist(SMI_zoos$mass, main = "OSF mass")
hist(SMI_zoos$SVL, main = "OSF SVL")
#none looking super normal
hist(mass_zoos$mass_2021, main = "OSF mass 2021")
hist(mass_zoos$mass_2022, main = "OSF mass 2022") #looks fairly normal
#wild data
hist(wild_SMI$avg_mass, main = "wild mass")
hist(wild_SMI$avg_SVL, main = "wild SVL") #looks fairly normal

ggdensity(SMI_zoos$mass, main = "2021 Mass", xlab = "mass (g)")
ggdensity(wild_SMI$avg_mass, main = "mass for wild frogs", xlab = "mass (g)")
ggdensity(wild_SMI$avg_SVL, main = "SVL for wild frogs", xlab = "length (mm)")

#Shapiro-wilks test for normalcy 
shapiro.test(SMI_zoos$age_yrs)
shapiro.test(SMI_zoos$mass)
shapiro.test(SMI_zoos$SVL)

ggqqplot(SMI_zoos$age_yrs, na.rm = TRUE)
ggqqplot(SMI_zoos$mass) #looks the closest to normal distribution
ggqqplot(SMI_zoos$SVL)
#pre-post data
ggqqplot(prepo$mass_pre20, main = "qq-plot mass pre-brum 2020", na.rm = TRUE)
ggqqplot(prepo$svl_pre20, na.rm = TRUE)
ggqqplot(prepo$mass_post21, na.rm = TRUE)
ggqqplot(prepo$svl_post21, na.rm = TRUE)
ggqqplot(prepo$mass_pre_21, main = "qq-plot mass pre-brum 2021", na.rm = TRUE)
ggqqplot(prepo$svl_pre21, na.rm = TRUE) #these actually all look close to normal if not normal

plot(mass ~ SVL, data = SMI_zoos, main = "mass by SVL")
plot(svl_pre20 ~ mass_pre20, data = prepo, main = "svl by mass pre-brumation 2020")
plot(svl_post21 ~ mass_post21, data = prepo, main = "svl by mass post-brumation 2021")
plot(svl_pre21 ~ mass_pre_21, data = prepo, main = "svl by mass pre-brumation 2021")

#check mass_zoos data
ggdensity(mass_zoos$mass_2022, na.rm = TRUE)
ggdensity(mass_zoos$mass_2021, na.rm = TRUE)

shapiro.test(mass_zoos$mass_2022) #nearly but not quite normal - non-parametric distribution
shapiro.test(mass_zoos$mass_2021)

#shapiro for wild data
shapiro.test(wild_SMI$avg_mass)
shapiro.test(wild_SMI$avg_SVL)
#shapiro for pre-post data
shapiro.test(prepo$mass_pre20) #normal
shapiro.test(prepo$length_pre20)
shapiro.test(prepo$svl_pre20) #add calculated svl values to check for normalcy
shapiro.test(prepo$mass_post21)
shapiro.test(prepo$length_post21) #same as below
shapiro.test(prepo$svl_post21) #normal
shapiro.test(prepo$mass_pre_21) #normal
shapiro.test(prepo$svl_pre21) #nearly normal but not quite

###### some potential subsets of data that I might need? from confidence interval calculations ###
#filter out NAs 
# age_na <- SMI_zoos %>% 
#   filter(!is.na(age_yrs))
# post_2021 <- prepo %>% 
#   filter(!is.na(mass_post21))
# pre_2021 <- prepo %>% 
#   filter(!is.na(svl_pre21))
# po_2021 <- prepo %>% 
#   filter(!is.na(svl_post21))


##############################################################
####linear regressions to predict SVL for VanAqua SUL values #
##############################################################
#### using length measurements taken by Pourya - Oct 2022 ####
conversion <- read.csv("VA_lengths_oct2022.csv", header = TRUE)
#add columns for average SUL and SVL
conversion$avgSUL <- rowMeans(conversion[ , c(3:5)], na.rm=TRUE)
conversion$avgSVL <- rowMeans(conversion[ , c(6:8)], na.rm = TRUE)
conversion <- conversion %>% 
  relocate(avgSUL, .after = SUL3) %>% 
  relocate(avgSVL, .after = SVL3)
view(conversion)
conversion$frog <- as.factor(conversion$frog)
str(conversion)
summary(conversion)

###### Linear regression of SVL by SUL
linreg <- lm(avgSVL ~ avgSUL, data = conversion)
summary(linreg)

plot(avgSVL ~ avgSUL, data = conversion, main = "Oct 2022 SUL by SVL") # n=24
abline(linreg)
plot(linreg$residuals, pch = 16, col = "red", main = "SVL ~ SUL residuals") #plot the residuals - do they look random?
plot(cooks.distance(linreg), pch = 16, col = "blue", main = "SVL ~ SUL cooks distance") 

# looks like one outlier, maybe two
identify_outliers(conversion, avgSUL) #frog 10 and 23 are outliers according to SUL
identify_outliers(conversion, avgSVL) #frog 23 is outlier according to SVL
# none are extreme outliers and I have verified they are not errors in recording - frog 23 is just very small
################################################### will leave the outliers because of sample size ##########

######## Create equation of model fit using estimates: 
############# SVL = 1.1957693x - 10.7090655, R2 = 0.6163136 (moderate fit), n=24, p=3.3613e-06 ########
#######################################################################################################

## try again with frog 10 and 23 removed as outliers
conv_outliers <- conversion %>% 
  filter(frog != "10")%>% 
  filter(frog != "23")
view(conv_outliers)
lm.out <- lm(avgSVL ~ avgSUL, data = conv_outliers)
summary(lm.out)
plot(avgSVL ~ avgSUL, data = conv_outliers, main = "Oct 2022 SUL by SVL") # n=24
abline(lm.out)
##### This looks worse. more spread in data. R2 and p have gotten worse
##### stick with outliers included. not actual outliers just limited sample size

#subset prepo data frame by populations: GVZ, VA
prepo.GVZ <- subset(prepo, pop == "GVZ")
str(prepo.GVZ)
prepo.VA <- subset(prepo, pop == "VA")
str(prepo.VA)

str(wild_SMI) #n=502
#subset wild data by populations
wild_SMI_CH <- subset(wild_SMI, pop =="CH") #n=7
wild_SMI_MS <- subset(wild_SMI, pop =="MS") #n=155
wild_SMI_MT <- subset(wild_SMI, pop =="MT") #n=1
wild_SMI_MV <- subset(wild_SMI, pop =="MV") #n=336
wild_SMI_ST <- subset(wild_SMI, pop =="ST") #n=3

########################################################################################
########### Linear regressions to find best fit for calculating bsma #####################
########################################################################################

### Mass (dependant) vs length (independent)
############################ Pre-brumation (i.e. Nov) 2020
plot(prepo$svl_pre20, prepo$mass_pre20, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(prepo$pop)]),
     xlab = "snout-vent length (mm)", ylab = "mass (g)")
legend("bottomright", title = "Population", c("GVZoo", "VanAqua"), fill = c('blue', 'green'), cex = 0.8)
abline(lm(mass_pre20 ~ svl_pre20, data = prepo), col = "black")
lm(mass_pre20 ~ svl_pre20, data = prepo)
summary(lm(mass_pre20 ~ svl_pre20, data = prepo))
abline(lm(mass_pre20 ~ svl_pre20, data = prepo.GVZ), col = "blue")
summary(lm(mass_pre20 ~ svl_pre20, data = prepo.GVZ))
abline(lm(mass_pre20 ~ svl_pre20, data = prepo.VA), col = "green")
summary(lm(mass_pre20 ~ svl_pre20, data = prepo.VA))
legend("topleft", title = "Regression lines", c("n=87, Adjusted R-squared: 0.8325, p-value: 2.2e-16",
                                                "n=42, Adjusted R-squared: 0.6823, p-value: 9.996e-12",
                                                "n=45, Adjusted R-squared: 0.8839, p-value: 2.2e-16"),
       fill = c('black', 'blue', 'green'), cex = 0.8)
title("Pre-brumation (Nov) 2020")

#Log-transformation
prepo$log.svl_pre20 = log(prepo$svl_pre20)
prepo$log.mass_pre20 = log(prepo$mass_pre20)
prepo.GVZ$log.svl_pre20 = log(prepo.GVZ$svl_pre20)
prepo.VA$log.svl_pre20 = log(prepo.VA$svl_pre20)
prepo.GVZ$log.mass_pre20 = log(prepo.GVZ$mass_pre20)
prepo.VA$log.mass_pre20 = log(prepo.VA$mass_pre20)

### Log Mass vs Log Length for 
############################ Pre-brumation (i.e. Nov) 2020
plot(prepo$log.svl_pre20, prepo$log.mass_pre20, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(prepo$pop)]),
     xlab = "Log transformed snout-vent length", ylab = "Log transformed mass")
legend("topleft", title = "Population", c("GVZoo", "VanAqua"), fill = c('blue', 'green'), cex = 0.8)
abline(lm(log.mass_pre20 ~ log.svl_pre20, data = prepo), col = "black")
lm(log.mass_pre20 ~ log.svl_pre20, data = prepo)
summary(lm(log.mass_pre20 ~ log.svl_pre20, data = prepo))
abline(lm(log.mass_pre20 ~ log.svl_pre20, data = prepo.GVZ), col = "blue")
summary(lm(log.mass_pre20 ~ log.svl_pre20, data = prepo.GVZ))
abline(lm(log.mass_pre20 ~ log.svl_pre20, data = prepo.VA), col = "green")
summary(lm(log.mass_pre20 ~ log.svl_pre20, data = prepo.VA))
legend("bottomright", title = "Regression lines", c("n=87, Adjusted R-squared: 0.9227, p-value: 2.2e-16",
                                                    "n=42, Adjusted R-squared: 0.6915, p-value: 5.527e-12",
                                                    "n=45, Adjusted R-squared: 0.9515, p-value: 2.2e-16"),
       fill = c('black', 'blue', 'green'), cex = 0.8)
title("Pre-brumation (Nov) 2020")

### Mass (dependant) vs length (independent)
############################ Post-brumation (i.e. Mar) 2021
plot(prepo$svl_post21, prepo$mass_post21, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(prepo$pop)]),
     xlab = "snout-vent length (mm)", ylab = "mass (g)")
legend("bottomright", title = "Population", c("GVZoo", "VanAqua"), fill = c('blue', 'green'), cex = 0.8)
abline(lm(mass_post21 ~ svl_post21, data = prepo), col = "black")
lm(mass_post21 ~ svl_post21, data = prepo)
summary(lm(mass_post21 ~ svl_post21, data = prepo))
abline(lm(mass_post21 ~ svl_post21, data = prepo.GVZ), col = "blue")
summary(lm(mass_post21 ~ svl_post21, data = prepo.GVZ))
abline(lm(mass_post21 ~ svl_post21, data = prepo.VA), col = "green")
summary(lm(mass_post21 ~ svl_post21, data = prepo.VA))
legend("topleft", title = "Regression lines", c("n=55, Adjusted R-squared: 0.4238, p-value: 4.462e-08",
                                                "n=40, Adjusted R-squared: 0.2616, p-value: 0.0004406",
                                                "n=15, Adjusted R-squared: 0.7145, p-value: 4.422e-05"),
       fill = c('black', 'blue', 'green'), cex = 0.6)
title("Post-brumation (Mar) 2021")

#Log-transformation
prepo$log.svl_post21 = log(prepo$svl_post21)
prepo$log.mass_post21 = log(prepo$mass_post21)
prepo.GVZ$log.svl_post21 = log(prepo.GVZ$svl_post21)
prepo.VA$log.svl_post21 = log(prepo.VA$svl_post21)
prepo.GVZ$log.mass_post21 = log(prepo.GVZ$mass_post21)
prepo.VA$log.mass_post21 = log(prepo.VA$mass_post21)
wild_SMI$log.avg_mass = log(wild_SMI$avg_mass)
wild_SMI$log.avg_SVL = log(wild_SMI$avg_SVL)

### Log Mass vs Log Length for 
############################ Post-brumation (i.e. Mar) 2021
plot(prepo$log.svl_post21, prepo$log.mass_post21, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(prepo$pop)]),
     xlab = "Log transformed snout-vent length", ylab = "Log transformed mass")
legend("bottomright", title = "Population", c("GVZoo", "VanAqua"), fill = c('blue', 'green'), cex = 0.8)
abline(lm(log.mass_post21 ~ log.svl_post21, data = prepo), col = "black")
lm(log.mass_post21 ~ log.svl_post21, data = prepo)
summary(lm(log.mass_post21 ~ log.svl_post21, data = prepo))
abline(lm(log.mass_post21 ~ log.svl_post21, data = prepo.GVZ), col = "blue")
summary(lm(log.mass_post21 ~ log.svl_post21, data = prepo.GVZ))
abline(lm(log.mass_post21 ~ log.svl_post21, data = prepo.VA), col = "green")
summary(lm(log.mass_post21 ~ log.svl_post21, data = prepo.VA))
legend("topleft", title = "Regression lines", c("n=55, Adjusted R-squared: 0.4555, p-value: 9.62e-09",
                                                "n=40, Adjusted R-squared: 0.2462, p-value: 0.000668",
                                                "n=15, Adjusted R-squared: 0.7341, p-value: 2.761e-05"),
       fill = c('black', 'blue', 'green'), cex = 0.6)
title("Post-brumation (Mar) 2021")


# Mass (dependant) vs length (independent)
############################ Pre-brumation (i.e. Nov) 2021
plot(prepo$svl_pre21, prepo$mass_pre_21, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(prepo$pop)]),
     xlab = "snout-vent length (mm)", ylab = "mass (g)")
legend("topleft", title = "Population", c("GVZoo", "VanAqua"), fill = c('blue', 'green'), cex = 0.8)
abline(lm(mass_pre_21 ~ svl_pre21, data = prepo), col = "black")
lm(mass_pre_21 ~ svl_pre21, data = prepo)
summary(lm(mass_pre_21 ~ svl_pre21, data = prepo))
abline(lm(mass_pre_21 ~ svl_pre21, data = prepo.GVZ), col = "blue")
summary(lm(mass_pre_21 ~ svl_pre21, data = prepo.GVZ))
abline(lm(mass_pre_21 ~ svl_pre21, data = prepo.VA), col = "green")
summary(lm(mass_pre_21 ~ svl_pre21, data = prepo.VA))
legend("bottomright", title = "Regression lines", c("n=75, Adjusted R-squared: 0.3753, p-value: 3.116e-09",
                                                    "n=40, Adjusted R-squared: 0.0288, p-value: 0.1502",
                                                    "n=35, Adjusted R-squared: 0.6113, p-value: 1.78e-08"),
       fill = c('black', 'blue', 'green'), cex = 0.7)
title("Pre-brumation (Nov) 2021")

#Log-transformation
prepo$log.svl_pre21 = log(prepo$svl_pre21)
prepo$log.mass_pre_21 = log(prepo$mass_pre_21)
prepo.GVZ$log.svl_pre21 = log(prepo.GVZ$svl_pre21)
prepo.VA$log.svl_pre21 = log(prepo.VA$svl_pre21)
prepo.GVZ$log.mass_pre_21 = log(prepo.GVZ$mass_pre_21)
prepo.VA$log.mass_pre_21 = log(prepo.VA$mass_pre_21)

# Log Mass vs Log Length for 
############################ Pre-brumation (i.e. Nov) 2021
plot(prepo$log.svl_pre21, prepo$log.mass_pre_21, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(prepo$pop)]),
     xlab = "Log transformed snout-vent length", ylab = "Log transformed mass")
legend("topleft", title = "Population", c("GVZoo", "VanAqua"), fill = c('blue', 'green'), cex = 0.8)
abline(lm(log.mass_pre_21 ~ log.svl_pre21, data = prepo), col = "black")
lm(log.mass_pre_21 ~ log.svl_pre21, data = prepo)
summary(lm(log.mass_pre_21 ~ log.svl_pre21, data = prepo))
abline(lm(log.mass_pre_21 ~ log.svl_pre21, data = prepo.GVZ), col = "blue")
summary(lm(log.mass_pre_21 ~ log.svl_pre21, data = prepo.GVZ))
abline(lm(log.mass_pre_21 ~ log.svl_pre21, data = prepo.VA), col = "green")
summary(lm(log.mass_pre_21 ~ log.svl_pre21, data = prepo.VA))
legend("bottomright", title = "Regression lines", c("n=75, Adjusted R-squared: 0.4983, p-value: 9.185e-13",
                                                    "n=40, Adjusted R-squared: 0.02358, p-value: 1716",
                                                    "n=35, Adjusted R-squared: 0.7226, p-value: 6.32e-11"),
       fill = c('black', 'blue', 'green'), cex = 0.6)
title("Pre-brumation (Nov) 2021")

##################################### The pre-brum 2020 regression is the best (highest R-squared) for GVZoo
### Will use this regression to calculate GVZoo's own bsma #################################################
#################### Will use larger VA dataset to calculate bsma for that population ######################


# Mass (dependant) vs length (independent)
############################ Wild Frog Data - full dataset (from 2011-2022)
plot(avg_SVL~avg_mass, data = wild_SMI)
plot(wild_SMI$avg_SVL, wild_SMI$avg_mass, pch = 16, cex = 1.0, ####### This isn't working ???? ################
     col = (c('blue', 'green', 'purple', 'orange', 'red')[as.numeric(wild_SMI$pop)]), #something up with this line
     xlab = "snout-vent length (mm)", ylab = "mass (g)")
legend("bottomright", title = "Population", c("Chaplin", "Maria Slough", "Mountain Slough", "Morris Valley", "Semmihault"),
       fill = c('blue', 'green', 'purple', 'orange', 'red'), cex = 0.6)
abline(lm(avg_mass ~ avg_SVL, data = wild_SMI), col = "black")
lm(avg_mass ~ avg_SVL, data = wild_SMI)
summary(lm(avg_mass ~ avg_SVL, data = wild_SMI))
abline(lm(avg_mass ~ avg_SVL, data = wild_SMI_CH), col = "blue")
summary(lm(avg_mass ~ avg_SVL, data = wild_SMI_CH))
abline(lm(avg_mass ~ avg_SVL, data = wild_SMI_MS), col = "green")
summary(lm(avg_mass ~ avg_SVL, data = wild_SMI_MS))
abline(lm(avg_mass ~ avg_SVL, data = wild_SMI_MV), col = "orange")
summary(lm(avg_mass ~ avg_SVL, data = wild_SMI_MV))
abline(lm(avg_mass ~ avg_SVL, data = wild_SMI_ST), col = "red")
summary(lm(avg_mass ~ avg_SVL, data = wild_SMI_ST))
#no abline for MT pop because only 1 observation
legend("topright", title = "Regression lines", c("n=502, Adjusted R-squared: 0.5608, p-value: 2.2e-16",
                                                 "n=7, Adjusted R-squared: 0.5825, p-value: 0.0281",
                                                 "n=155, Adjusted R-squared: 0.3152, p-value: 5.7e-15",
                                                 "n=336, Adjusted R-squared: 0.7011, p-value: 2.2e-16",
                                                 "n=3, Adjusted R-squared: 0.9991, p-value: 0.0133"),
       fill = c('black', 'blue', 'green', 'orange', 'red'), cex = 0.45)
title("Full Wild M by L")

#Log-transformation
wild_SMI$log.avg_SVL = log(wild_SMI$avg_SVL)
wild_SMI$log.avg_mass = log(wild_SMI$avg_mass)
wild_SMI_CH$log.avg_SVL = log(wild_SMI_CH$avg_SVL)
wild_SMI_MS$log.avg_SVL = log(wild_SMI_MS$avg_SVL)
wild_SMI_MV$log.avg_SVL = log(wild_SMI_MV$avg_SVL)
wild_SMI_ST$log.avg_SVL = log(wild_SMI_ST$avg_SVL)
wild_SMI_CH$log.avg_mass = log(wild_SMI_CH$avg_mass)
wild_SMI_MS$log.avg_mass = log(wild_SMI_MS$avg_mass)
wild_SMI_MV$log.avg_mass = log(wild_SMI_MV$avg_mass)
wild_SMI_ST$log.avg_mass = log(wild_SMI_ST$avg_mass)

# Log Mass vs Log Length for 
############################ Wild Frog Data - full dataset (from 2012-2022*)
plot(wild_SMI$log.avg_SVL, wild_SMI$log.avg_mass, pch = 16, cex = 1.0,
     col = (c('blue', 'green', 'purple', 'orange', 'red')[as.numeric(wild_SMI$pop)]),
     xlab = "snout-vent length (mm)", ylab = "mass (g)")
legend("topleft", title = "Population", c("Chaplin", "Maria Slough", "Mountain Slough", "Morris Valley", "Semmihault"),
       fill = c('blue', 'green', 'purple', 'orange', 'red'), cex = 0.6)
abline(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI), col = "black")
lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI)
summary(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI))
abline(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI_CH), col = "blue")
summary(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI_CH))
abline(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI_MS), col = "green")
summary(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI_MS))
abline(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI_MV), col = "orange")
summary(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI_MV))
abline(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI_ST), col = "red")
summary(lm(log.avg_mass ~ log.avg_SVL, data = wild_SMI_ST))
#no abline for MT pop because only 1 observation
legend("bottomright", title = "Regression lines", c("n=502, Adjusted R-squared: 0.5812, p-value: 2.2e-16",
                                                    "n=7, Adjusted R-squared: 0.5935, p-value: 0.0261",
                                                    "n=155, Adjusted R-squared: 0.3641, p-value: 2.2e-16",
                                                    "n=336, Adjusted R-squared: 0.7080, p-value: 2.2e-16",
                                                    "n=3, Adjusted R-squared: 0.9982, p-value: 0.0191"),
       fill = c('black', 'blue', 'green', 'orange', 'red'), cex = 0.5)
title("Log of M by L Full Wild dataset")
########## Also missing data points here. need to come back to and figure out what's going wrong #####

###########################################################
## VA SMI data organized by age
###########################################################
VA_byage_wide <- read.csv("VA_SMI_byage.csv", header = TRUE)
#transform from wide to long 
VA_age <- melt(setDT(VA_byage_wide), id=1:8, measure=patterns("^mass", "^sul"),
               value.name=c("mass", "sul"), variable.name="ageof", na.rm = TRUE)

VA_age <- VA_age %>% 
  mutate(age = case_when(ageof == "1" ~ "0", ageof == "2" ~ "1", ageof == "3" ~ "2", ageof == "4" ~ "3",
                         ageof == "5" ~ "4", ageof == "6" ~ "5", ageof == "7" ~ "6", ageof == "8" ~ "7",
                         ageof == "9" ~ "8", ageof == "10" ~ "9", ageof == "11" ~ "10")) %>%
  relocate(age, .after = EB) #need to make way cleaner way of doing this. set eqn instead of each value
VA_age$ageof <- NULL #drop ageof column now that we have proper age
view(VA_age)
VA_age$frog_id <- as.factor(VA_age$frog_id) #factor 144 levels (individual frogs)
VA_age$EB <- as.factor(VA_age$EB) #as factor. Two levels: 0=EB 1=not (OK)
VA_age$pop <- as.factor(VA_age$pop) #factor with 1 level (all from VA)
VA_age$source <- as.factor(VA_age$source) #factor with 13 levels
VA_age$birth_yr <- as.factor(VA_age$birth_yr) #13 levels
# VA_age$age <- as.factor(VA_age$age) #11 levels - do we actually want age as a factor? always?
str(VA_age)
dim(VA_age) #444 rows, 11 columns
summary(VA_age)

#convert sul values into SVL using previous regression equation
#SVL = 1.1957x - 10.709, R2 = 0.6163
SULconvert <- #set function for SUL to SVL conversion
  function(x) {
    1.1957693 * x - 10.7090655
  }

VA_age$svl <- SULconvert(VA_age$sul) #new column of converted svl values

#visualize distribution
hist(VA_age$mass, main = "VA mass")
hist(VA_age$svl, main = "VA svl")
ggdensity(VA_age$mass, main = "VA mass", xlab = "mass (g)")
ggdensity(VA_age$svl, main = "VA svl", xlab = "snout-vent length (mm)")
ggqqplot(VA_age$mass, na.rm = TRUE)
ggqqplot(VA_age$svl, na.rm = TRUE)
#calculate normality by shapiro (H0 = normal)
shapiro.test(VA_age$mass) #significant = not normal
shapiro.test(VA_age$svl) #significant = not normal
######neither are normal but have large sample size (n=444 total observations)
######samples are not independent because we have multiple (repeated) measurements per individual frog

#####subset data for adults (>2 years old)
#SMI should standardize for age but want to look at normality of data as is
VA_adults <- VA_age %>% 
  filter(age != 0 & age != 1) #exclude frogs at 0 and 1 yrs
view(VA_adults)
str(VA_adults) #n=263, still have all  individual frogs though

VA_adults2 <- VA_adults %>% #exclude 0-2 years old. technically "maturity" is reached after 2yrs
  filter(age != 2)
view(VA_adults2)
str(VA_adults2) #n=168 obs

#view their distribution
hist(VA_adults$mass)
hist(VA_adults$svl)
hist(VA_adults2$mass)
hist(VA_adults2$svl) #these all look considerably more normal

### Mass (dependant) vs length (independent)
############################ VA full dataset
#subset by status (EB or OK)
VA_EB <- subset(VA_age, EB == "0")
str(VA_EB) #n=158
VA_OK <- subset(VA_age, EB == "1")
str(VA_OK) #n=276

plot(VA_age$svl, VA_age$mass, pch = 16, cex = 1.3,
     col = (c('red', 'blue')[as.numeric(VA_age$EB)]),
     xlab = "snout-vent length (mm)", ylab = "mass (g)")
legend("bottomright", title = "status", c("Egg bound", "Other"), fill = c('red', 'blue'), cex = 0.8)
abline(lm(mass ~ svl, data = VA_age), col = "black")
lm(mass ~ svl, data = VA_age)
summary(lm(mass ~ svl, data = VA_age))
abline(lm(mass ~ svl, data = VA_EB), col = "red")
summary(lm(mass ~ svl, data = VA_EB))
abline(lm(mass ~ svl, data = VA_OK), col = "blue")
summary(lm(mass ~ svl, data = VA_OK))
legend("topleft", title = "Regression lines", c("n=434, Adjusted R-squared: 0.7607, p-value: 2.2e-16", 
                                                "n=158, Adjusted R-squared: 0.7774, p-value: 2.2e-16",
                                                "n=276, Adjusted R-squared: 0.7577, p-value: 2.2e-16"),
       fill = c('black', 'red', 'blue'), cex = 0.65)
title("MbyL for full VA dataset")

#Log-transformation
VA_age$log.svl = log(VA_age$svl)
VA_age$log.mass = log(VA_age$mass)
VA_EB$log.svl = log(VA_EB$svl)
VA_OK$log.svl = log(VA_OK$svl)
VA_EB$log.mass = log(VA_EB$mass)
VA_OK$log.mass = log(VA_OK$mass)

### Log Mass vs Log Length for 
############################ Pre-brumation (i.e. Nov) 2020
plot(VA_age$log.svl, VA_age$log.mass, pch = 16, cex = 1.3,
     col = (c('red', 'blue')[as.numeric(VA_age$EB)]),
     xlab = "Log transformed snout-vent length", ylab = "Log transformed mass")
legend("topleft", title = "Status", c("Egg bound", "Other"), fill = c('red', 'blue'), cex = 0.8)
abline(lm(log.mass ~ log.svl, data = VA_age), col = "black")
lm(log.mass ~ log.svl, data = VA_age)
summary(lm(log.mass ~ log.svl, data = VA_age))
abline(lm(log.mass ~ log.svl, data = VA_EB), col = "red")
summary(lm(log.mass ~ log.svl, data = VA_EB))
abline(lm(log.mass ~ log.svl, data = VA_OK), col = "blue")
summary(lm(log.mass ~ log.svl, data = VA_OK))
legend("bottomright", title = "Regression lines", c("n=434, Adjusted R-squared: 0.906, p-value: 2.2e-16",
                                                    "n=158, Adjusted R-squared: 0.9264, p-value: 2.2e-16",
                                                    "n=276, Adjusted R-squared: 0.897, p-value: 2.2e-16"),
       fill = c('black', 'red', 'blue'), cex = 0.7)
title("Log of MbyL for VA dataset")

##############################adults vs young
#already have VA_adults2 for >2 year old frogs
str(VA_adults2) #n=167
VA_0 <- subset(VA_age, age == "0")
str(VA_0) #n=75
VA_1 <- subset(VA_age, age == "1")
str(VA_1) #n=106
VA_2 <- subset(VA_age, age == "2")
str(VA_2) #n=95

VA_age <- VA_age %>% 
  mutate(maturity = case_when(age > 2 ~ "adult",
                              age == 0  ~ "0",
                              age == 1 ~ "1",
                              age == 2 ~ "2")) %>% 
  relocate(maturity, .after = age)
VA_age$maturity <- as.factor(VA_age$maturity)
str(VA_age) #4 levels: 0, 1, 2, adult (>2 yrs old)

plot(VA_age$svl, VA_age$mass, pch = 16, cex = 1.3,
     col = (c('orange', 'green', 'purple', 'blue')[as.numeric(VA_age$maturity)]),
     xlab = "snout-vent length (mm)", ylab = "mass (g)")
legend("bottomright", title = "Age", c("0yo", "1yo", "2yo", "adult (>2yo)"), 
       fill = c('orange', 'green', 'purple', 'blue'), cex = 0.7)
abline(lm(mass ~ svl, data = VA_age), col = "black")
lm(mass ~ svl, data = VA_age)
summary(lm(mass ~ svl, data = VA_age))
abline(lm(mass ~ svl, data = VA_0), col = "orange")
summary(lm(mass ~ svl, data = VA_0))
abline(lm(mass ~ svl, data = VA_1), col = "green")
summary(lm(mass ~ svl, data = VA_1))
abline(lm(mass ~ svl, data = VA_2), col = "purple")
summary(lm(mass ~ svl, data = VA_2))
abline(lm(mass ~ svl, data = VA_adults2), col = "blue")
summary(lm(mass ~ svl, data = VA_adults2))
legend("topleft", title = "Regression lines", c("n=444, Adjusted R-squared: 0.7607, p-value: 2.2e-16", 
                                                "n=168, Adjusted R-squared: 0.7582, p-value: 2.2e-16",
                                                "n=75, Adjusted R-squared: 0.7657, p-value: 2.2e-16",
                                                "n=106, Adjusted R-squared: 0.7274, p-value: 2.2e-16",
                                                "n=95, Adjusted R-squared: 0.3358, p-value: 2.2e-16"),
       fill = c('black', 'orange', 'green', 'purple', 'blue'), cex = 0.6)
title("MbyL by age at VA")

#######Log transformation
VA_0$log.svl = log(VA_0$svl)
VA_0$log.mass = log(VA_0$mass)
VA_1$log.svl = log(VA_1$svl)
VA_1$log.mass = log(VA_1$mass)
VA_2$log.svl = log(VA_2$svl)
VA_2$log.mass = log(VA_2$mass)
VA_adults2$log.svl = log(VA_adults2$svl)
VA_adults2$log.mass = log(VA_adults2$mass)
#######Plot Log M by L
plot(VA_age$log.svl, VA_age$log.mass, pch = 16, cex = 1.3,
     col = (c('orange', 'green', 'purple', 'blue')[as.numeric(VA_age$maturity)]),
     xlab = "Log transformed snout-vent length", ylab = "log transformed mass")
legend("topleft", title = "Age", c("0yo", "1yo", "2yo", "adult (>2yo)"), 
       fill = c('orange', 'green', 'purple', 'blue'), cex = 0.7)
abline(lm(log.mass ~ log.svl, data = VA_age), col = "black")
lm(log.mass ~ log.svl, data = VA_age)
summary(lm(log.mass ~ log.svl, data = VA_age))
abline(lm(log.mass ~ log.svl, data = VA_0), col = "orange")
summary(lm(log.mass ~ log.svl, data = VA_0))
abline(lm(log.mass ~ log.svl, data = VA_1), col = "green")
summary(lm(log.mass ~ log.svl, data = VA_1))
abline(lm(log.mass ~ log.svl, data = VA_2), col = "purple")
summary(lm(log.mass ~ log.svl, data = VA_2))
abline(lm(log.mass ~ log.svl, data = VA_adults2), col = "blue")
summary(lm(log.mass ~ log.svl, data = VA_adults2))
legend("bottomright", title = "Regression lines", c("n=444, Adjusted R-squared: 0.9060, p-value: 2.2e-16", 
                                                    "n=168, Adjusted R-squared: 0.7699, p-value: 2.2e-16",
                                                    "n=75, Adjusted R-squared: 0.8214, p-value: 2.2e-16",
                                                    "n=106, Adjusted R-squared: 0.8183, p-value: 2.2e-16",
                                                    "n=95, Adjusted R-squared: 0.4102, p-value: 2.2e-16"),
       fill = c('black', 'orange', 'green', 'purple', 'blue'), cex = 0.6)
title("Log of MbyL by age at VA")

#########################################################################################
# Do regression of full 2021 snapshot of zoos and wild ##################################
# Need bsma for TZ pop still so we can calculate SMI for 2021 ###########################
#########################################################################################

str(SMI_zoos)
#subset SMI_zoos data frame by populations: CH, GVZ, MS, MV, TZ, VA
SMI_zoos.CH <- subset(SMI_zoos, pop == "CH")
str(SMI_zoos.CH) #n=7
SMI_zoos.GVZ <- subset(SMI_zoos, pop == "GVZ")
str(SMI_zoos.GVZ) #n=41
SMI_zoos.MS <- subset(SMI_zoos, pop == "MS")
str(SMI_zoos.MS) #n=25
SMI_zoos.MV <- subset(SMI_zoos, pop == "MV")
str(SMI_zoos.MV) #n=13
SMI_zoos.TZ <- subset(SMI_zoos, pop == "TZ")
str(SMI_zoos.TZ) #n=9
SMI_zoos.VA <- subset(SMI_zoos, pop == "VA")
str(SMI_zoos.VA) #n=15
#subset by locale (zoos vs wild)
SMI_zoos.z <- subset(SMI_zoos, locale == "Zoo")
str(SMI_zoos.z) #n=65
SMI_zoos.w <- subset(SMI_zoos, locale == "Wild")
str(SMI_zoos.w) #n=45

### Mass (dependant) vs length (independent)
############################ 2021 snapshot of all zoos and wild populations (wild from all years)
plot(SMI_zoos$SVL, SMI_zoos$mass, pch = 16, cex = 1.3,
     col = (c('blue', 'purple', 'green', 'orange', 'yellow', 'red')[as.numeric(SMI_zoos$pop)]),
     xlab = "snout-vent length (mm)", ylab = "mass (g)")
legend("bottomright", title = "Population", c("Chaplin", "GVZoo", "Maria Slough", "Morris Valley", "TZoo", "VanAqua"), 
       fill = c('blue', 'purple', 'green', 'orange', 'yellow', 'red'), cex = 0.6)
abline(lm(mass ~ SVL, data = SMI_zoos), col = "black")
lm(mass ~ SVL, data = SMI_zoos)
summary(lm(mass ~ SVL, data = SMI_zoos))
abline(lm(mass ~ SVL, data = SMI_zoos.CH), col = "blue")
summary(lm(mass ~ SVL, data = SMI_zoos.CH))
abline(lm(mass ~ SVL, data = SMI_zoos.GVZ), col = "purple")
summary(lm(mass ~ SVL, data = SMI_zoos.GVZ))
abline(lm(mass ~ SVL, data = SMI_zoos.MS), col = "green")
summary(lm(mass ~ SVL, data = SMI_zoos.MS))
abline(lm(mass ~ SVL, data = SMI_zoos.MV), col = "orange")
summary(lm(mass ~ SVL, data = SMI_zoos.MV))
abline(lm(mass ~ SVL, data = SMI_zoos.TZ), col = "yellow")
summary(lm(mass ~ SVL, data = SMI_zoos.TZ))
abline(lm(mass ~ SVL, data = SMI_zoos.VA), col = "red")
summary(lm(mass ~ SVL, data = SMI_zoos.VA))
legend("topleft", title = "Regression lines", c("n=7, Adjusted R-squared: 0.5825, p-value: 0.02806",
                                                "n=41, Adjusted R-squared: 0.252, p-value: 0.000489",
                                                "n=25, Adjusted R-squared: 0.6651, p-value: 4.14e-07",
                                                "n=13, Adjusted R-squared: 0.6776, p-value: 0.0003331",
                                                "n=9, Adjusted R-squared: 0.4312, p-value: 0.03256",
                                                "n=15, Adjusted R-squared: 0.7145, p-value: 4.422e-05"),
       fill = c('blue', 'purple', 'green', 'orange', 'yellow', 'red'), cex = 0.6)
title("2021 spring Mass by SVL")

#Log-transformation
SMI_zoos$log.SVL = log(SMI_zoos$SVL)
SMI_zoos$log.mass = log(SMI_zoos$mass)
SMI_zoos.CH$log.SVL = log(SMI_zoos.CH$SVL)
SMI_zoos.CH$log.mass = log(SMI_zoos.CH$mass)
SMI_zoos.GVZ$log.SVL = log(SMI_zoos.GVZ$SVL)
SMI_zoos.GVZ$log.mass = log(SMI_zoos.GVZ$mass)
SMI_zoos.MS$log.SVL = log(SMI_zoos.MS$SVL)
SMI_zoos.MS$log.mass = log(SMI_zoos.MS$mass)
SMI_zoos.MV$log.SVL = log(SMI_zoos.MV$SVL)
SMI_zoos.MV$log.mass = log(SMI_zoos.MV$mass)
SMI_zoos.TZ$log.SVL = log(SMI_zoos.TZ$SVL)
SMI_zoos.TZ$log.mass = log(SMI_zoos.TZ$mass)
SMI_zoos.VA$log.SVL = log(SMI_zoos.VA$SVL)
SMI_zoos.VA$log.mass = log(SMI_zoos.VA$mass)

# Log Mass vs Log Length
############################ 2021 snapshot of all zoos and wild populations
plot(SMI_zoos$log.SVL, SMI_zoos$log.mass, pch = 16, cex = 1.3,
     col = (c('blue', 'purple', 'green', 'orange', 'yellow', 'red')[as.numeric(SMI_zoos$pop)]),
     xlab = "Log transformed snout-vent length", ylab = "Log transformed mass (g)")
legend("bottomright", title = "Population", c("Chaplin", "GVZoo", "Maria Slough", "Morris Valley", "TZoo", "VanAqua"), 
       fill = c('blue', 'purple', 'green', 'orange', 'yellow', 'red'), cex = 0.6)
abline(lm(log.mass ~ log.SVL, data = SMI_zoos), col = "black")
lm(log.mass ~ log.SVL, data = SMI_zoos)
summary(lm(log.mass ~ log.SVL, data = SMI_zoos))
abline(lm(log.mass ~ log.SVL, data = SMI_zoos.CH), col = "blue")
summary(lm(log.mass ~ log.SVL, data = SMI_zoos.CH))
abline(lm(log.mass ~ log.SVL, data = SMI_zoos.GVZ), col = "purple")
summary(lm(log.mass ~ log.SVL, data = SMI_zoos.GVZ))
abline(lm(log.mass ~ log.SVL, data = SMI_zoos.MS), col = "green")
summary(lm(log.mass ~ log.SVL, data = SMI_zoos.MS))
abline(lm(log.mass ~ log.SVL, data = SMI_zoos.MV), col = "orange")
summary(lm(log.mass ~ log.SVL, data = SMI_zoos.MV))
abline(lm(log.mass ~ log.SVL, data = SMI_zoos.TZ), col = "yellow")
summary(lm(log.mass ~ log.SVL, data = SMI_zoos.TZ))
abline(lm(log.mass ~ log.SVL, data = SMI_zoos.VA), col = "red")
summary(lm(log.mass ~ log.SVL, data = SMI_zoos.VA))
legend("topleft", title = "Regression lines", c("n=7, Adjusted R-squared: 0.5935, p-value: 0.02613",
                                                "n=41, Adjusted R-squared: 0.2388, p-value: 0.0007024",
                                                "n=25, Adjusted R-squared: 0.6846, p-value: 2.048e-07",
                                                "n=13, Adjusted R-squared: 0.6446, p-value: 0.0005799",
                                                "n=9, Adjusted R-squared: 0.4231, p-value: 0.03438",
                                                "n=15, Adjusted R-squared: 0.7345, p-value: 2.73e-05"),
       fill = c('blue', 'purple', 'green', 'orange', 'yellow', 'red'), cex = 0.5)
title("Log of MbyL spring 2021")
## the un-transformed regression is the best R-squared (higher value is best) for TZoo (though still not great)
## will use it for calculating TZ's bsma later on


################################################################################
################################################################################
########## Calculate scaled mass index (SMI) ###################################
################################################################################

# Scaled mass index, Pieg & Green 2009

# Scaled mass index (SMI): ^Mi = Mi (Lo/Li)^bsma
# where: ^Mi is scaled mass: predicted body mass for individual i when the linear body measure is standardized to L0.
#         Mi is mass of individual i
#         Lo is an arbitraty value of length (e.g. arithmetric mean of L for the population of study)
#         Li is the linear body measurement of individual i
#         bsma is the scaling exponent: calculated indirectly by dividing the slope from an OLS regression (bOLS) 
#by the Pearson’s correlation coefficient r (LaBarbera 1989), or directly using online software (Bohonak and van der Linde 2004).
# calculate a different bsma for each population (VA, TZ, GVZ, wild) to account for different growth patterns in each

# calculate bsma manually (described in Pieg & Green 2009 p.1886)
#for Wild frogs - using wild allometric growth as standard for scaling exponent
lm(wild_SMI$log.avg_mass ~ wild_SMI$log.avg_SVL) #using log
slope <- 2.312869 #from variable above
cor.test(wild_SMI$log.avg_SVL, wild_SMI$log.avg_mass, method = c("pearson")) #using log again 
pcorr <- 0.762914 #from above (Pearson's product-moment correlation)
bsma.w <- slope / pcorr #manual calculation of bsma (scaling exponent)
bsma.w #3.03162

#calculate bsma again but this time automated using smatr() - Standardized Major Axis
library(smatr)
sma(wild_SMI$log.avg_mass ~ wild_SMI$log.avg_SVL) # the slope of the estimate is the bsma 
bsma.w <- 3.0316258

########################################################################################
#### Each population needs its own bsma - calcaulted based on their best regression ####
########################################################################################

##calculate bsma for GVZoo based on pre20 data - the best regression for GVZ
##### each population have their own scaling exponent - VA will be calculated later
sma(prepo.GVZ$log.mass_pre20 ~ prepo.GVZ$log.svl_pre20)
bsma.GVZ <- 3.183815

#calculate Lo
Lo <- mean(wild_SMI$avg_SVL, na.rm = TRUE)
Lo #70.820

Lo.GVZ <- mean(prepo.GVZ$svl_pre20, na.rm = TRUE)
Lo.GVZ #67.71429

#######create function for ScaledMassIndex using ^Mi = Mi ((Lo/Li)^bsma)
##using bsma from MbyL regression of wild frog data to make SMI function for wild populations
ScaledMassIndex.w <-
  function(x, y) {
    y * ((Lo / x) ^ bsma.w)
  }

wild_SMI$SMI <- ScaledMassIndex.w(wild_SMI$avg_SVL, wild_SMI$avg_mass)
head(wild_SMI)

#### make SMI equation using GVZ bsma - but keep same Lo because this is arbitrary chosen length
ScaledMassIndex.gvz <-
  function(x, y) {
    y * ((Lo / x) ^ bsma.GVZ)
  }

##### SMI equation for TZ population using un-transformed regression from SMI_zoos dataset
sma(SMI_zoos.TZ$mass ~ SMI_zoos.TZ$SVL)
bsma.tz <- 2.058488

ScaledMassIndex.tz <-
  function(x, y) {
    y * ((Lo / x) ^ bsma.tz)
  }

### calculate bsma for VA from VA_age dataset - best MbyL regression for VA data (most full dataset)
sma(VA_age$log.mass ~ VA_age$log.svl) #use log #the slope of the estimate is the bsma 
VA.bsma <- 2.80075
#check if different if just using VA adults
sma(VA_adults2$log.mass ~ VA_adults2$log.svl)
VA.a.bsma <- 2.50028
#wild bsma was 3.0316
#bsma from VA is closer to this and closer to slope of 3 in general (3 supposed to be best?)
########################
#calculate Lo
Lo.VA <- mean(VA_age$svl, na.rm = TRUE)
Lo.VA #56.4495
Lo.VA.a <- mean(VA_adults2$svl, na.rm = TRUE)
Lo.VA.a #69.7527
#wild Lo was 70.820 - but the Lo should be arbitrary

########Create new SMI functions for new bsma using old Lo (from wild population)
ScaledMassIndex.VA <- function(x, y) {
  y * ((Lo / x) ^ VA.bsma)
}

ScaledMassIndex.va <- function(x, y) {
  y * ((Lo / x) ^ VA.a.bsma)
}

##### Calculate and add SMI to VA dataset #####
VA_age$SMI.w <- ScaledMassIndex.w(VA_age$svl, VA_age$mass) #this is using wild bsma
VA_age$SMI <- ScaledMassIndex.VA(VA_age$svl, VA_age$mass) #using VA bsma
VA_age$SMI.va <- ScaledMassIndex.va(VA_age$svl, VA_age$mass) #using VA adults only bsma
head(VA_age)

#calculate SMI for VA adult subset
VA_adults$SMI <- ScaledMassIndex.VA(VA_adults$svl, VA_adults$mass)
VA_adults2$SMI <- ScaledMassIndex.VA(VA_adults2$svl, VA_adults2$mass)

#plot SMI in VA_age dataset
######################################################
#GRAPHING#############################################
ggplot(data = VA_age, aes(x = maturity, y = SMI, fill=factor(maturity, labels = c("0yo", "1yo", "2yo",
                                                                                  "adult (>2yo)"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Age (years)", y = "Scaled Mass Index (g)", title = "SMI at VanAqua (2012-2020)")+
  theme_classic()+
  scale_fill_manual(values = c("orange", "green", "purple", "blue"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0yo", "1" = "1yo", "2" = "2yo", "adult" = "adult (>2yo)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

ggplot(data = VA_age, aes(x = as.factor(age), y = SMI, fill=factor(age)))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Age (years)", y = "Scaled Mass Index (g)", title = "SMI at VanAqua (2012-2020)")+
  theme_classic()+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
## graphing above with SMI.w actually looks more standardized than SMI
## SMI (using own bsma) changes slightly over age. should be fairly consistent. SMI.va shows huge variation
## could reflect actual SMI changes in frogs as they age at VA though
## Stick with using full VA dataset for bsma and SMI calculations

ggboxplot(VA_age %>% filter(!is.na(maturity)), x = "EB", y = "SMI", facet.by = "maturity")+
  labs(x = "Status", y = "Scaled Mass Index (g)", title = "SMI by status at VanAqua")+
  scale_x_discrete(labels=c("0" = "egg bound", "1" = "OK"))

#############################################################################################
### Calculate SMI for 2021 data using Lo and bsma(s) defined previously - unique bsma per pop
#############################################################################################
## use subset data and combine back together at end 
####################################################################

SMI_zoos.w$SMI <- ScaledMassIndex.w(SMI_zoos.w$SVL, SMI_zoos.w$mass) 
view(SMI_zoos.w)
SMI_zoos.GVZ$SMI <- ScaledMassIndex.gvz(SMI_zoos.GVZ$SVL, SMI_zoos.GVZ$mass)
view(SMI_zoos.GVZ)
SMI_zoos.TZ$SMI <- ScaledMassIndex.tz(SMI_zoos.TZ$SVL, SMI_zoos.TZ$mass)
view(SMI_zoos.TZ)
SMI_zoos.VA$SMI <- ScaledMassIndex.VA(SMI_zoos.VA$SVL, SMI_zoos.VA$mass)
view(SMI_zoos.VA)
###### rbind all these subsets with their unique SMI values - such that each SMI is calculated with own bsma
# SMI_zoos.w needs log columns
SMI_zoos.w$log.SVL = log(SMI_zoos.w$SVL)
SMI_zoos.w$log.mass = log(SMI_zoos.w$mass)

full_21 <- rbind(SMI_zoos.w, SMI_zoos.GVZ, SMI_zoos.TZ, SMI_zoos.VA)
view(full_21)

######################################################################################################
# plot SMI ####################################################################
ggplot(data = full_21, aes(x = pop, y = SMI, fill=factor(pop, labels = c("Chaplin", "GVZoo", "Maria Slough",
                                                                         "Morris Valley", "TZoo", "VanAqua"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Population", y = "Scaled Mass Index (g)", title = "SMI in spring 2021")+
  theme_classic()+
  scale_fill_manual(values = c("blue", "purple", "green", "orange", "yellow", "red"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("CH" = "Chaplin", "GVZ" = "GVZoo", "MS" = "Maria Slough",
                            "MV" = "Morris Valley", "TZ" = "TZoo", "VA" = "VanAqua"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 11))

#wild vs zoos 
#create new column to group three wild populations together but separate from zoos
full_21 <- full_21 %>% 
  mutate(population = case_when(pop == "CH" ~ "wild",
                                pop == "MS" ~ "wild",
                                pop == "MV" ~ "wild",
                                pop == "GVZ" ~ "GVZ",
                                pop == "VA" ~ "VA",
                                pop == "TZ" ~ "TZ")) %>% 
  relocate(population, .after = 3)
full_21$population <- as.factor(full_21$population) #4 levels, GVZ TZ VA wild
view(full_21)

ggplot(data = full_21, aes(x = population, y = SMI, fill=factor(population, labels = c("GVZoo", 
                                                                                       "TZoo", "VanAqua", "Wild"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Population", y = "Scaled Mass Index (g)", title = "SMI in spring 2021")+
  theme_classic()+
  scale_fill_manual(values = c("purple", "yellow", "red", "blue"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("wild" = "Wild", "GVZ" = "GVZoo", "TZ" = "TZoo", "VA" = "VanAqua"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 13))

ggboxplot(full_21, x = "EB", y = "SMI", facet.by = "population")+
  labs(x = "Status", y = "Scaled Mass Index (g)", title = "2021 SMI in egg bound vs OK OSF")+
  scale_x_discrete(labels=c("0" = "egg bound", "1" = "other females"))+
  theme(text = element_text(family = "Arial"))

####################################################################################################
###############subset prepo to just mass and svl and transform to long / stacked form ##############
####################################################################################################
prepo.long <- prepo[, c(1:5, 7:8, 10:12)]
prepo.l <- melt(setDT(prepo.long), id=1:4, measure=patterns("^mass", "^svl"),
                value.name=c("mass", "svl"), variable.name="when", na.rm = TRUE)
prepo.l <- prepo.l %>% 
  mutate(season = case_when(when == "1" ~ "pre20", when == "2" ~ "post21", when == "3" ~ "pre21")) %>%
  relocate(season, .after = EB) 
prepo.l$when <- NULL #drop column now that we have proper season
view(prepo.l)

################# subset out GVZoo and VA into own dataframes to calculate separate SMI's with own bsma ########
prepo.gvz <- subset(prepo.l, pop == "GVZ")
prepo.va <- subset(prepo.l, pop == "VA")

prepo.gvz$SMI <- ScaledMassIndex.gvz(prepo.gvz$svl, prepo.gvz$mass)
prepo.va$SMI <- ScaledMassIndex.VA(prepo.va$svl, prepo.va$mass)
view(prepo.gvz) #n=122
view(prepo.va) #n=95
###### combine these two datasets, now with their proper respective SMI's
combine_prepo <- rbind(prepo.gvz, prepo.va)

#add in wild SMI - need to make wild dataset that includes EB and season
wild_SMI$EB <- 1
wild_SMI <- wild_SMI %>% 
  relocate(EB, .after = pop) %>% 
  rename(mass = avg_mass,
         svl = avg_SVL,
         log.mass = log.avg_mass,
         log.svl = log.avg_SVL)
view(wild_SMI)

wild.move <- wild_SMI[, c(1:3, 5:7, 10)]
wild.move$season <- "postbrum" 
wild.move <- wild.move %>% 
  relocate(season, .after = EB) %>% 
  relocate(birth_yr, .after = pop)

###want to add TZ data for spring 2021 as well. Grab TZ data from SMI_zoos.TZ
pp_TZ <- SMI_zoos.TZ[, c(1:2, 4, 7:9, 12)] 
#add column to match other datasets 
pp_TZ$season <- "post21"
pp_TZ <- pp_TZ %>% 
  rename(svl = SVL) %>%
  relocate(season, .after = EB)
view(pp_TZ)

#combine dataframes - all have 8 columns
full_pp <- rbind(combine_prepo, wild.move, pp_TZ)
view(full_pp)
#add some other columns to make more cohesive for grouping
full_pp <- full_pp %>% 
  mutate(when = case_when(grepl("pre", season) ~ "fall", grepl("post", season) ~ "spring")) %>% 
  relocate(when, .after = season)
#column for wild vs all zoos
full_pp <- full_pp %>% 
  mutate(population = case_when(pop == "GVZ" ~ "GVZ", pop == "TZ" ~ "TZ", pop == "VA" ~ "VA", pop == "CH" ~ "wild",
                                pop == "ST" ~ "wild", (grepl("M", pop) ~ "wild"))) %>% 
  relocate(population, .after = pop)
view(full_pp)
full_pp$population <- as.factor(full_pp$population)
full_pp$season <- factor(full_pp$season, levels = c("pre20", "post21", "postbrum", "pre21"), ordered = TRUE)
str(full_pp)
dim(full_pp) #728 rows, 10 columns

###### Plot all together now! ##############################################
###########################################################################
library(EnvStats) # for adding sample size
ggplot(data = full_pp, aes(x = population, y = SMI, fill=population))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~season)+
  labs(x = "Population", y = "Scaled Mass Index (g)", title = "pre and post-brumation SMI from 2020-2021")+
  theme_classic()+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))+
  stat_n_text() #sample size - figure out how to put it in the actual box

ggplot(data = full_pp, aes(x = population, y = SMI, fill=population))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~when)+
  labs(x = "Population", y = "Scaled Mass Index (g)", title = "Seasonal SMI from 2020-2021")+
  theme_classic()+
  scale_fill_brewer(palette = "Set2")+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))+
  stat_n_text()
# see difference in fall and spring SMI between wild and zoos - look for significant differences here

# facet by fall or spring to show difference in EB vs other in each season (using fill)
ggplot(data = full_pp, aes(x = population, y = SMI, fill=EB))+
  geom_boxplot(show.legend = TRUE)+
  facet_grid(~when)+
  labs(x = "Population", y = "Scaled Mass Index (g)", title = "SMI over seasons")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
  theme(legend.title = element_blank())+
  scale_fill_discrete(name = "Status", labels = c("0" = "egg bound", "1" = "other"))+
  scale_fill_manual(values = c("#DC7466", "#9C8C78"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
# this is actually a great visualization of EB frogs in each season
##### need to try only including frogs who became egg bound in that year - i.e. 2021 (exclude 2022 EB mortalities)

### use combine_prepo to compare EB vs OK in GVZ vs VA over fall '20 and spring '21
### this still includes both their 2022 EB deaths as well... need to remove or look at without ?
str(combine_prepo) #217 obs from just GVZ and VA
combine_prepo$season <- as.factor(combine_prepo$season)
combine_prepo <- combine_prepo %>% 
  mutate(when = case_when(grepl("pre", season) ~ "fall", grepl("post", season) ~ "spring")) %>% 
  relocate(when, .after = season)
view(combine_prepo)
## Plot again with EB vs OK
ggplot(data = combine_prepo, aes(x = EB, y = SMI, fill=EB))+
  geom_boxplot(alpha = 0.8, show.legend = FALSE)+
  facet_grid(when~pop)+
  labs(x = "Status", y = "Scaled Mass Index (g)")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels = c("0" = "egg bound", "1" = "other"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#####mean SMI per population
mean(wild_SMI$SMI, na.rm = TRUE) #43.5887
#subset wild populations
wildCH_SMI <- subset(wild_SMI, pop == "CH")
str(wildCH_SMI) #n=7
wildMS_SMI <- subset(wild_SMI, pop == "MS")
wildMT_SMI <- subset(wild_SMI, pop == "MT")
wildMV_SMI <- subset(wild_SMI, pop == "MV")
wildST_SMI <- subset(wild_SMI, pop == "ST")
#means for wild pops 
mean(wildCH_SMI$SMI, na.rm = TRUE) #40.19464
mean(wildMS_SMI$SMI, na.rm = TRUE) #42.23196
mean(wildMT_SMI$SMI, na.rm = TRUE) #51.71032
mean(wildMV_SMI$SMI, na.rm = TRUE) #44.23784
mean(wildST_SMI$SMI, na.rm = TRUE) #46.19538

#according to Stevens 2013 (Intermediate Statistics: A Modern Approach), normality is not required (for ANOVA) for cases of large data sets (n>30)
#Kruskal-Wallis is an alternative option, compares medians of pops instead of means- does not require normality or homoscedasticity

#do a Levene or Bartlett test first before ANOVA to check for homoscedasticity (equality of variance) leveneTest()

####assume normality for mass and length data because most samples sizes are >30 and half the data are normal, half are not
#check homogeneity of variances (as assumption of ANOVA)


################################################################################
## create new dataframe combining VA and wild full datasets for SMI over time ##
################################################################################

adults <- VA_adults2[, c('frog_id', 'pop', 'EB', 'birth_yr', 'age', 'mass', 'svl', 'log.mass', 'log.svl', 'SMI')]
view(adults)
#combine dataframes
adult <- rbind(adults, wild_SMI)
view(adult)
str(adult) #670 rows, 10 columns - check: combined VA_adults2 (n=168) + wild_SMI (n=502) = GOOD
#recalculate SMI and try SMI using VA.bsma also
adult$SMI.w <- ScaledMassIndex.w(adult$svl, adult$mass)
adult$SMI.VA <- ScaledMassIndex.VA(adult$svl, adult$mass)
adult$SMI.va <- ScaledMassIndex.va(adult$svl, adult$mass)
view(adult)

########### Plot the SMI
ggplot(data = adult, aes(x = pop, y = SMI, fill=factor(pop, labels = c("VA", "CH", "MS",
                                                                       "MT", "MV", "ST"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Population", y = "Scaled Mass Index (g)", title = "SMI in adult VA vs wild over ~10 yrs")+
  theme_classic()+
  scale_fill_manual(values = c("green", "orange", "yellow", "purple", "blue", "brown"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("VA" = "VanAqua", "CH" = "Chaplin", "MS" = "Maria Slough", "MT" = "Mountain",
                            "MV" = "Morris Valley", "ST" = "Semmihault"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#create new column to group wild populations together
adult <- adult %>% 
  mutate(population = case_when(pop == "CH" ~ "wild",
                                pop == "MS" ~ "wild",
                                pop == "MV" ~ "wild",
                                pop == "MT" ~ "wild",
                                pop == "ST" ~ "wild",
                                pop == "VA" ~ "VA",)) %>% 
  relocate(population, .after = 2)
adult$population <- as.factor(adult$population)
str(adult)

#plot again but comparing full overall wild with VA
ggplot(data = adult, aes(x = population, y = SMI, fill=factor(population, labels = c("VA", "wild"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Population", y = "Scaled Mass Index (g)", title = "SMI in adult VA vs overall wild")+
  theme_classic()+
  scale_fill_manual(values = c("green", "orange"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("VA" = "Vancovuer Aquarum", "wild" = "Combined Wild"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
# Interpretation: VA SMI is significantly** higher than combined wild SMI over past 10 years
# confirms single year (2021) snapshot

############ NOTE: VanAqua SMI is based on fall measurements while wild is from breeding season (spring) #####

##### facet with status
ggboxplot(adult, x = "EB", y = "SMI", facet.by = "population")+ #I like this graph better than the one below
  labs(x = "Status", y = "Scaled Mass Index (g)", title = "SMI by status")+
  scale_x_discrete(labels=c("0" = "egg bound", "1" = "OK"))
## again this confirms the single year (2021) snapshot of VA SMI data - showing no significant differences in EB vs other

###############################################################
####### look at (fall) mass changes with age in TZ and VA
###############################################################
mass_age.wide <- read.csv("fall_mass_byage.csv", header = TRUE)
head(mass_age.wide)
view(mass_age.wide)
#need to transform into long form
mass_age <- melt(setDT(mass_age.wide), id=1:9, measure=patterns("^mass"),
              value.name=c("mass"), variable.name="ageof", na.rm = TRUE)
mass_age <- mass_age %>% 
  mutate(age_mass = case_when(ageof == 'mass0' ~ '0',
                              ageof == 'mass1' ~ '1',
                              ageof ==  'mass2' ~ '2',
                              ageof == 'mass3' ~ '3',
                              ageof ==  'mass4' ~ '4',
                              ageof ==  'mass5' ~ '5',
                              ageof ==  'mass6' ~ '6',
                              ageof == 'mass7' ~ '7',
                              ageof ==  'mass8' ~ '8',
                              ageof ==  'mass9' ~ '9',
                              ageof ==  'mass10' ~ '10',
                              ageof == 'mass11' ~ '11',
                              ageof == 'mass12' ~ '12')) %>% 
  relocate(age_mass, .after = EB)
mass_age$ageof <- NULL  #can now drop the ageof column since we have proper age       
view(mass_age)        
dim(mass_age) #566 rows, 11 columns
 
mass_age$ID <- as.factor(mass_age$ID) #convert a bunch of chr to factor
mass_age$pop <- as.factor(mass_age$pop)
mass_age$source <- as.factor(mass_age$source)
mass_age$birth_site <- as.factor(mass_age$birth_site) #12 levels
mass_age$birth_type <- as.factor(mass_age$birth_type) #2 levels 
mass_age$EB <- as.factor(mass_age$EB) #as factor. Two levels: 0=EB 1=not (OK)
mass_age$age_mass <- as.factor(mass_age$age_mass)
mass_age$age_mass <- factor(mass_age$age_mass, 
                            levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"),
                            ordered = TRUE) #double check the level order
str(mass_age)
summary(mass_age)

#visualize normality
hist(mass_age$mass) #bimodal thing again from 0-1 year olds
hist(mass_age$age_yrs)
ggqqplot(mass_age$mass, na.rm = TRUE)
shapiro.test(mass_age$mass) #significant = not-normal
# subset out adults (>2 years old at mass)
mass_adult <- mass_age %>% 
  filter(age_mass != 0 & age_mass != 1 & age_mass != 2)
view(mass_adult)
dim(mass_adult) #273 rows, 11 columns (down from 566)
#order the levels of age_mass correctly
mass_adult$age_mass <- factor(mass_adult$age_mass, 
                              levels = c("3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

#check if normality is any different for adults
hist(mass_adult$mass) #looks good
shapiro.test(mass_adult$mass) # NORMAL! p-value = 0.06646

#plot mass by age
ggplot(data = mass_age, aes(x = age_mass, y = mass, fill=factor(pop)))+ #can swap out mass_age for mass_adult to compare
  geom_boxplot()+
  labs(x = "Age", y = "Mass (g)", title = "mass of adult OSF with age")+
  theme_classic()+
  scale_fill_discrete(name = "Population",
                      breaks = c("TZ", "VA"),
                      labels = c("TZ", "VA"))+
   scale_fill_manual(name = "Population", 
                     values = c("blue", "green"))+
  theme(legend.position = "top")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

####### facet grid to show TZ vs VA with also status? #########
#filter out NA for status
mass.clean <- mass_age %>% #filter out the NA values for status
  filter(!is.na(EB))

ggplot(data = mass.clean, aes(x = age_mass, y = mass, fill=EB))+
  geom_boxplot()+
  facet_grid(~pop)+
  labs(x = "Age (yrs)", y = "Mass (g)", title = "mass changes with age and status")+
  theme_classic()+
  scale_fill_discrete(name = "Population",
                      breaks = c("TZ", "VA"),
                      labels = c("TZ", "VA"))+
  scale_fill_manual(name = "Status", 
                    values = c("red", "purple"))+
  theme(legend.position = "top")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
### I would like to look at above graph but with EB/OK and VA/TZ plotted on some axis - diff colours
# add column to split into groups

## subset the data
mass.VA <- subset(mass.clean, pop == "VA") #  n=556
mass.VA.EB <- subset(mass.VA, EB == "0")
mass.VA.OK <- subset(mass.VA, EB == "1")
str(mass.VA.EB) # n=167
str(mass.VA.OK) # n=303
#do same with TZ
mass.TZ <- subset(mass.clean, pop == "TZ")
mass.TZ.EB <- subset(mass.TZ, EB == "0")
mass.TZ.OK <- subset(mass.TZ, EB == "1")
str(mass.TZ.EB) # n=12
str(mass.TZ.OK) # n=74
## add columns to the subsets for EB + pop
mass.VA <- mass.VA %>% 
  mutate(status = case_when(EB == "0" ~ "EB.va", EB == "1" ~ "OK.va")) %>% 
  relocate(status, .after = EB)
mass.TZ <- mass.TZ %>% 
  mutate(status = case_when(EB == "0" ~ "EB.tz", EB == "1" ~ "OK.tz")) %>% 
  relocate(status, .after = EB)
#combine two datasets with new columns
EB_comp <- rbind(mass.TZ, mass.VA)
view(EB_comp) #566 obs again
#try graph again but now group by status
#####################################################
# library(devtools)
# devtools::install_github('Mikata-Project/ggthemr')
# library(ggthemr)
# ggthemr('dust')

######### PLOT for CHS
ggplot(data = EB_comp, aes(x = age_mass, y = mass, fill=status))+
  geom_boxplot()+
  labs(x = "Age (yrs)", y = "Mass (g)", title = "Mass Changes with Age", subtitle = "Vancouver Aquarium and Toronto Zoo")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
  scale_fill_manual(name = "Pop Status", 
                    labels = c("EB.va" = "EB at VA", "OK.va" = "OK at VA", "EB.tz" = "EB at TZ", "OK.tz" = "OK at TZ"),
                    values = c("#663300", "#CC0000", "#CC9933","#FF9999"))+
  theme(legend.position = "right")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#### trying to figure out how to make this a scatter plot and add curved lines? 
## or even just leave as box plots but show all four facets
plot(mass.clean$age_mass, mass.clean$mass, pch = 16, cex = 1.3,
     col = (c('red', 'blue')[as.numeric(mass.clean$pop)]),
     xlab = "Age (yrs)", ylab = "Mass (g)")
# legend("bottomright", title = "status", c("Egg bound", "Other"), fill = c('red', 'blue'), cex = 0.8)
# abline(lm(mass ~ svl, data = VA_age), col = "black")
# lm(mass ~ svl, data = VA_age)
# summary(lm(mass ~ svl, data = VA_age))
# abline(lm(mass ~ svl, data = VA_EB), col = "red")
# summary(lm(mass ~ svl, data = VA_EB))
# abline(lm(mass ~ svl, data = VA_OK), col = "blue")
# summary(lm(mass ~ svl, data = VA_OK))
# legend("topleft", title = "Regression lines", c("n=434, Adjusted R-squared: 0.7607, p-value: 2.2e-16", 
#                                                 "n=158, Adjusted R-squared: 0.7774, p-value: 2.2e-16",
#                                                 "n=276, Adjusted R-squared: 0.7577, p-value: 2.2e-16"),
#        fill = c('black', 'red', 'blue'), cex = 0.65)
# title("MbyL for full VA dataset")
##########################################

ggboxplot(mass_age, x = "pop", y = "mass", facet.by = "EB")+
  labs(x = "Population", y = "Mass (g)", title = "Fall mass by status")+
  scale_x_discrete(labels=c("VA" = "VanAqua", "TZ" = "TZoo"))


ggboxplot(mass.clean, x = "age_mass", y = "mass", facet.by = "EB")+
  labs(x = "Age (years)", y = "Mass (g)", title = "Fall mass by status - combined TZ + VA", 
       subtitle = "0 = egg bound, 1 = other")

ggboxplot(mass.clean, x = "EB", y = "mass", facet.by = "age_mass")+
  labs(x = "Status", y = "Mass (g)", title = "TZ + VA fall mass by status")+
  scale_x_discrete(labels=c("0" = "EB", "1" = "OK"))

#Regression of M by Age
#subset by zoo
ma_VA <- subset(mass_age, pop == "VA") #n=480
ma_TZ <- subset(mass_age, pop == "TZ") #n=86

plot(as.numeric(mass_age$age_mass), mass_age$mass, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(mass_age$pop)]),
     xlab = "Age (yrs)", ylab = "Mass (g)")
legend("topleft", title = "Population", c("Toronto Zoo", "Vancouver Aquarium"), 
       fill = c('blue', 'green'), cex = 0.7)
abline(lm(mass ~ age_mass, data = mass_age), col = "black")
lm(mass ~ age_mass, data = mass_age)
summary(lm(mass ~ age_mass, data = mass_age))
abline(lm(mass ~ age_mass, data = ma_TZ), col = "blue")
summary(lm(mass ~ age_mass, data = ma_TZ))
abline(lm(mass ~ age_mass, data = ma_VA), col = "green")
summary(lm(mass ~ age_mass, data = ma_VA))
legend("bottomright", title = "Regression lines", c("n=566, Adjusted R-squared: 0.5867, p-value: 2.2e-16", 
                                                "n=86, Adjusted R-squared: 0.2408, p-value: 0.000471",
                                                "n=480, Adjusted R-squared: 0.6118, p-value: 2.2e-16"),
       fill = c('black', 'blue', 'green'), cex = 0.6)
title("Mass by age")

#try with adult only data
maa_VA <- subset(mass_adult, pop == "VA") #n=195
maa_TZ <- subset(mass_adult, pop == "TZ") #n=78

plot(as.numeric(mass_adult$age_mass), mass_adult$mass, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(mass_adult$pop)]),
     xlab = "Age (yrs)", ylab = "Mass (g)")
legend("topleft", title = "Population", c("Toronto Zoo", "Vancouver Aquarium"), 
       fill = c('blue', 'green'), cex = 0.6)
abline(lm(mass ~ age_mass, data = mass_adult), col = "black")
lm(mass ~ age_mass, data = mass_adult)
summary(lm(mass ~ age_mass, data = mass_adult))
abline(lm(mass ~ age_mass, data = maa_TZ), col = "blue")
summary(lm(mass ~ age_mass, data = maa_TZ))
abline(lm(mass ~ age_mass, data = maa_VA), col = "green")
summary(lm(mass ~ age_mass, data = maa_VA))
legend("bottomright", title = "Regression lines", c("n=273, Adjusted R-squared: 0.126, p-value: 9.905e-07", 
                                                    "n=78, Adjusted R-squared: 0.1301, p-value: 0.02194",
                                                    "n=195, Adjusted R-squared: 0.2036, p-value: 4.825e-08"),
       fill = c('black', 'blue', 'green'), cex = 0.6)
title("Mass by age - adults only")
####can we plot a more curving line to this rather than straight line of best fit? may not be linear relationship

#######Log transformation
mass_age$log.mass = log(mass_age$mass)
ma_TZ$log.mass = log(ma_TZ$mass)
ma_VA$log.mass = log(ma_VA$mass)
mass_age$log.age = log(as.numeric(mass_age$age_mass))
ma_TZ$log.age = log(as.numeric(ma_TZ$age_mass))
ma_VA$log.age = log(as.numeric(ma_VA$age_mass))

plot(mass_age$log.age, mass_age$log.mass, pch = 16, cex = 1.3,
     col = (c('blue', 'green')[as.numeric(mass_age$pop)]),
     xlab = "Log transformed age", ylab = "Log transformed mass")
legend("topleft", title = "Population", c("Toronto Zoo", "Vancouver Aquarium"), 
       fill = c('blue', 'green'), cex = 0.6)
abline(lm(log.mass ~ log.age, data = mass_age), col = "black")
lm(log.mass ~ log.age, data = mass_age)
summary(lm(log.mass ~ log.age, data = mass_age))
abline(lm(log.mass ~ log.age, data = ma_TZ), col = "blue")
summary(lm(log.mass ~ log.age, data = ma_TZ))
abline(lm(log.mass ~ log.age, data = ma_VA), col = "green")
summary(lm(log.mass ~ log.age, data = ma_VA))
legend("bottomright", title = "Regression lines", c("n=566, Adjusted R-squared: 0.659, p-value: 2.2e-16", 
                                                    "n=86, Adjusted R-squared: 0.2095, p-value: 5.587e-06",
                                                    "n=480, Adjusted R-squared: 0.6845, p-value: 2.2e-16"),
       fill = c('black', 'blue', 'green'), cex = 0.8)
title("Log of mass by age")

################################################################################################
####### Let's do some ANOVAs and/or Kruskal-Wallis tests to see if differences are significant
################################################################################################
####### Compare SMI among populations (using full_21 dataset and looking btw 'populations')

#Hypotheses
#Null: There is no difference in the means of factor A, B, C, D
#Alternative: There is a difference in the means of factor A, B, C, D

#Assumptions:
#Randomly sampled
#Independently sampled - this is true because single year snapshot so only one measurement per frog
#Normality - my variables going in (mass & svl) were not normal but need to check for calculated SMI
####ANOVA is "robust" against non-normality when large sample sizes are at play. 
####Test normality with shapiro-wilks, q-q plots and histogram
#Equal variance - test with bartlett

#Design:
#Balanced = sample sizes of all groups are the same
#Unbalanced = sample sizes of all groups are not the same
##If unbalanced, need to use alternative SS calculation (type IV)

######### Case 1: SMI ~ population #######
######## snapshot of 2021 spring SMI #####
# Continuous dependent variable: SMI
# Categorical independent variables: GVZ, TZ, VA, Wild
#Unbalanced design, unequal sample sizes

##################################################################
# Test assumptions
hist(full_21$SMI) #looking for bell shape
# Histogram looks skewed, long right tail

qqnorm(full_21$SMI) #looking for straight line
#looks quite sloped

#Shapiro Wilks test: null hypothesis = variance of data is normally distributed
#check by population
shapiro.test(subset(full_21, population == "GVZ")$SMI) # p-value: 2.918e-05
shapiro.test(subset(full_21, population == "TZ")$SMI) # p-value: 0.7690
shapiro.test(subset(full_21, population == "VA")$SMI) # p-value: 0.00131
shapiro.test(subset(full_21, population == "wild")$SMI) # p-value: 0.5217
# TZ and wild are normal (p > 0.05) but GVZ and VA are not.
# TZ n=9, wild n=45, GVZ n=41, VA n=15
# TZ may have come out as normal just because of very small sample size 

#try transformation
full_21$log.SMI <- log(full_21$SMI)
plot(full_21$log.SMI)
hist(full_21$log.SMI) #looks slightly better
qqnorm(full_21$log.SMI) #better but still not great
shapiro.test(subset(full_21, population == "GVZ")$log.SMI) # p-value: 0.001568
shapiro.test(subset(full_21, population == "TZ")$log.SMI) # p-value: 0.942
shapiro.test(subset(full_21, population == "VA")$log.SMI) # p-value: 0.0185
shapiro.test(subset(full_21, population == "wild")$log.SMI) # p-value: 0.2525
###  GVZoo and VA are still not normal - but much closer than they were. 
# GVZoo has highest sample size of other zoos
# Will still use log transformed data - with higher sample size GVZoo can be considered normal ???

#####################################################################################################

# Test assumption of homogeneity of variance on transformed data
# Null hypothesis of bartlett: variance is equal

plot(log.SMI~population, data = full_21)
bartlett.test(log.SMI~population, data = full_21)
# p-value: 0.01035 - barely reject the null hypothesis
# the variance is not homogenous so we cannot proceed to ANOVA
# Alternative tests are Welch-ANOVA or Kruskal-Wallis
# Kruskal-Wallis does not require normality or homoscedasticity of variances
# Let's use our untransformed data in Kruskal-Wallis test because KW does not require normality

#####################################################################################################
# Kruskal-Wallis
kruskal.test(SMI~population, data = full_21)
# chi-squared = 54.8578, p-value = 7.3628e-12
# p < 0.05 so we reject the null hypotheses of all mean SMI being equal between populations
# but we still do not know WHICH population means differ - just know at least one is different

# need to do post-hoc test to determine which populations differ
# Use Dunn Test (though could also use pairwise.wilcox.test())
library(FSA) #see citation('FSA') if used in publication
dunnTest(SMI~population, data = full_21, method = "holm")
# if last column (adjusted p-value) is < 0.05 then the indicated means differ significantly
# GVZ - TZ: significant
# GVZ - VA: significant 
# TZ - VA: NOT significant = means of VA and TZ do not differ significantly (though it is close)
# GVZ - wild: significant (meaning GVZ is significantly lower)
# TZ - wild: significant (barely though - much closer than others)
# VA - wild: significant

###### we see that all means differ significantly EXCEPT VA from TZ - but ALL zoos differ significantly from wild

# Plot statistical results to show significant differences
library(ggstatsplot)

ggbetweenstats(
  data = full_21,
  x = population,
  y = SMI,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

######### Case 2: SMI ~ wild pop #######
######## do each of the wild populations (in 2021) significantly differ in SMI? #####
# Continuous dependent variable: SMI
# Categorical independent variables: (pop) CH, MS, MV (also GVZ, TZ, VA)
#Unbalanced design, unequal sample sizes

##################################################################
# Test assumptions
hist(full_21$SMI) #long right tail
qqnorm(full_21$SMI) #quite sloped

#Shapiro Wilks test: null hypothesis = variance of data is normally distributed
shapiro.test(subset(full_21, pop == "CH")$SMI) # p-value: 0.1874 (n=7)
shapiro.test(subset(full_21, pop == "MS")$SMI) # p-value: 0.9055 (n=25)
shapiro.test(subset(full_21, pop == "MV")$SMI) # p-value: 0.07678 (n=13)
shapiro.test(subset(full_21, pop == "GVZ")$SMI) # p-value: 2.918e-05 (n=41)
shapiro.test(subset(full_21, pop == "TZ")$SMI) # p-value: 0.769 (n=9)
shapiro.test(subset(full_21, pop == "VA")$SMI) # p-value: 0.001311 (n=15)
## GVZ and VA are not normal - others are.
# try with transformed data
shapiro.test(subset(full_21, pop == "CH")$log.SMI) # p-value: 0.1224
shapiro.test(subset(full_21, pop == "MS")$log.SMI) # p-value: 0.78999
shapiro.test(subset(full_21, pop == "MV")$log.SMI) # p-value: 0.4164
shapiro.test(subset(full_21, pop == "GVZ")$log.SMI) # p-value: 0.001568
shapiro.test(subset(full_21, pop == "TZ")$log.SMI) # p-value: 0.942
shapiro.test(subset(full_21, pop == "VA")$log.SMI) # p-value: 0.0185
## GVZ and VA are closer but still not normal

# Test variance. Null hypothesis of bartlett: variance is equal
plot(log.SMI~pop, data = full_21) #a few potential outliers but sample sizes too low to remove
bartlett.test(log.SMI~pop, data = full_21) 
## p-value = 0.00102 means reject the null so variance is not homogenous. Cannot use ANOVA

kruskal.test(SMI~pop, data = full_21)
# chi-squared = 59.0722, p-value = 1.8897e-11 (reject the null - mean SMI are not equal)

# Post hoc: Dunn Test 
dunnTest(SMI~pop, data = full_21, method = "holm")
# specifically interested in significant differences between wild populations (CH, MS, MV)
## CH-MS not significant (p: 0.986), CH-MV not significant (p: 1.00), MS-MV not significant (p: 0.315)
## interestingly we also see few sig. diffs between wild pops and GVZ and TZ. 
### is this different when using full (2012-2022) wild dataset - compared to just 2021 captures? 
ggbetweenstats(
  data = full_21,
  x = pop,
  y = SMI,
  type = "nonparametric", # Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)


######## Case 2b: look at SMI ~ pop for full wild dataset (2012-2022) #######
############ use wild.move dataframe to compare only wild pops  #############
# Test assumptions
hist(wild.move$SMI) #bit of a right tail
qqnorm(wild.move$SMI) # looks pretty good actually
wild.move$pop <- as.factor(wild.move$pop)
summary(wild.move$pop)

shapiro.test(subset(wild.move, pop == "CH")$SMI) # p-value: 0.1874 (n=7)
shapiro.test(subset(wild.move, pop == "MS")$SMI) # p-value: 0.01266 (n=161)
shapiro.test(subset(wild.move, pop == "MT")$SMI) # insufficient sample size (n=1)
shapiro.test(subset(wild.move, pop == "MV")$SMI) # p-value: 5.452e-12 (n=330)
shapiro.test(subset(wild.move, pop == "ST")$SMI) # p-value: 0.3217 (n=3)
## MV and MS are not normal, and they have largest sample size 
# try with transformed data
wild.move$log.SMI <- log(wild.move$SMI)
shapiro.test(subset(wild.move, pop == "CH")$log.SMI) # p-value: 0.1224
shapiro.test(subset(wild.move, pop == "MS")$log.SMI) # p-value: 0.5132 ## now normal
shapiro.test(subset(wild.move, pop == "MV")$log.SMI) # p-value: 5.993e-05 ### still not normal
shapiro.test(subset(wild.move, pop == "ST")$log.SMI) # p-value: 0.29422 ### this is actually worse ??
## Bartlett test
plot(log.SMI~pop, data = wild.move) # quite a few potential outliers in MS and MV especially - look into this
wild.noMT <- wild.move %>% 
  filter(pop != "MT")
bartlett.test(log.SMI~pop, data = wild.noMT) 
## p-value = 2.4011e-09 variance is not homogenous. Kruskal-Wallis needed
kruskal.test(log.SMI~pop, data = wild.noMT)
# chi-squared = 6.1874, p-value = 0.1028
# Interpretation: p-value is > 0.05 we accept the null hypothesis
##### There are NO significant differences in SMI mean between our wild populations ###
########## Interesting contrast to using only 2021 data for wild pops - where we saw sig differences ########
#############################################################################################################


######### Case 3: SMI ~ population + season ####################################
### SMI of all zoos and combined wild (since no sig diffs were found above) ####
# Need to subset full_pp dataset by "when" to compare within correct season ####
full_pp$log.SMI <- log(full_pp$SMI)
full_f <- subset(full_pp, when == "fall") ### only has SMI from VA (n=80) and GVZ (n=82)
full_f$population <- factor(full_f$population) #should drop off the extra levels (i.e. TZ and wild)
levels(full_f$population) #now only have GVZ and VA
full_s <- subset(full_pp, when == "spring") #GVZ (n=40), TZ (n=9), VA (n=15), wild (n=502)
hist(full_f$SMI) # looks quite normal 
hist(full_s$SMI) # right tail
qqnorm(full_f$SMI) # pretty straight!
qqnorm(full_s$SMI) # quite sloped

#Shapiro Wilks test: null hypothesis = variance of data is normally distributed
shapiro.test(subset(full_f, population == "GVZ")$SMI) # p-value: 0.0172 (n=84) (not normal)
shapiro.test(subset(full_f, population == "VA")$SMI) # p-value: 6.156e-05 (n=90) (not normal)

shapiro.test(subset(full_s, population == "GVZ")$SMI) # p-value: 4.499e-05 (n=42) (not normal)
shapiro.test(subset(full_s, population == "VA")$SMI) # p-value: 0.00266 (n=45) (not normal)
shapiro.test(subset(full_s, population == "TZ")$SMI) # p-value: 0.769 (n=9) (normal)
shapiro.test(subset(full_s, population == "wild")$SMI) # p-value: 5.782e-10 (n=502) (not normal)
# try transformed data
shapiro.test(subset(full_f, population == "GVZ")$log.SMI) # p-value: 0.2964
shapiro.test(subset(full_s, population == "VA")$log.SMI) # p-value: 0.0234

shapiro.test(subset(full_s, population == "GVZ")$log.SMI) # p-value: 0.002292
shapiro.test(subset(full_s, population == "VA")$log.SMI) # p-value: 0.0234
shapiro.test(subset(full_s, population == "TZ")$log.SMI) # p-value: 0.942
shapiro.test(subset(full_s, population == "wild")$log.SMI) # p-value: 0.00363
## closer but still not normal for most

plot(log.SMI~population, data = full_f) #a few potential outliers but not bad
plot(log.SMI~population, data = full_s) # a couple outliers in VA but less sample size here
bartlett.test(log.SMI~population, data = full_f) #p-value: 0.91111 (= homogenous variance)
bartlett.test(log.SMI~population, data = full_s) #p-value: 0.4789 (= homogenous)
## even though variance was homogenous, most of the data was not normal so let's stick with KW test
kruskal.test(SMI~population, data = full_f)
# chi-squared = 0.0755, p-value = 0.7835 so we accept the null
kruskal.test(SMI~population, data = full_s)
# chi-squared = 65.301, p-value = 4.325e-14 so we reject the null

##compare with anova
aov.21f <- aov(SMI~population, data = full_f) 
summary(aov.21f)# p-value: 0.8321
aov.21s <- aov(SMI~population, data = full_s)
summary(aov.21s) # p-value: 2.22e-16
######### accept the null for fall - meaning there are no significant differences

library(rstatix) #had this for previous tests that were removed... still need? 

##### for full_s group we have >2 groups so we can do regular Dunn Test
dunnTest(SMI~population, data = full_s, method = "holm")
# all significant except TZ-VA, and GVZ-wild

ggbetweenstats(
  data = full_s,
  x = population,
  y = SMI,
  type = "nonparametric", #  Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant", #not including centrality plotting means we see the n=, but also medians
  bf.message = FALSE
)

######### Case 3: SMI ~ population + EB ####################################
### use combine_prepo dataset to compare EB vs other in VA and GVZoo in 2021 data ####
# Need to subset by "when" to compare within correct season ####
EB21_f <- subset(combine_prepo, when == "fall") # GVZ (n=82), VA (n=80); EB (27), OK (135)
EB21_s <- subset(combine_prepo, when == "spring") # GVZ (n=40), VA (n=15), EB (9), OK (46)
hist(EB21_f$SMI) # looks quite normal 
hist(EB21_s$SMI) # right tail
qqnorm(EB21_f$SMI) # pretty straight!
qqnorm(EB21_s$SMI) # quite sloped

#Shapiro Wilks test: we know previously data was not normal but was better log transformed 
EB21_f$log.SMI <- log(EB21_f$SMI)
EB21_s$log.SMI <- log(EB21_s$SMI)

#Bartlett by EB status (could use an F test for homogeneity... is this better for only 2 samples?)
bartlett.test(SMI~EB, data = EB21_f) #p = 0.5192 (homogenous)
bartlett.test(SMI~EB, data = EB21_s) #p = 0.5911 (homogenous)
bartlett.test(SMI~pop, data = EB21_f) #p = 0.7579 (homogenous)
bartlett.test(SMI~pop, data = EB21_s) #p = 0.0168 (NOT homogenous...)
### spring variances are not homogenous among population groups... shouldn't really use ANOVA for this? 

one.way.f <- aov(SMI~EB, data = EB21_f)
summary(one.way.f) #p = 0.727 ### accept the null. not significantly different
one.way.s <- aov(SMI~EB, data = EB21_s)
summary(one.way.s) #p= 0.135 ### accept the null. not significantly different

## try one-way using pop as independent
one.f.pop <- aov(SMI ~ pop, data = EB21_f)
one.s.pop <- aov(SMI ~ pop, data = EB21_s)

## Test some two-way ANOVAs to see which is the best fit (highest F?)
### two-way ANOVA without interaction
two.way.f <- aov(SMI ~ EB + pop, data = EB21_f)
two.way.s <- aov(SMI ~ EB + pop, data = EB21_s)
### with interaction
interaction.f <- aov(SMI ~ EB * pop, data = EB21_f)
interaction.s <- aov(SMI ~ EB * pop, data = EB21_s)

### with blocking interaction not applicable?

## compare fall models using AIC (Akaike information criterion) model selection
library(AICcmodavg)
model.set.f <- list(one.way.f, one.f.pop, two.way.f, interaction.f)
model.names <- c("one way EB", "one way pop", "two way", "interaction")
aictab(model.set.f, modnames = model.names)
# lists the model with the best fit first - here it is: interaction
# look at summary of best model (interaction)
summary(interaction.f)

##### can see only EB:pop explains a significant amount of the variation (*)
##### need to do post hoc to see which levels are different from one another 
############# Tukey's Honestly-Significant-Difference (TukeyHSD) test ##############
TukeyHSD(interaction.f)
# look for p < 0.05 to indicate groups with significant differences
## We see NO significant differences, even in the EB:pop - though this is supposed to explain a sig part of variance
##### and, there are no within-pop significant differences. Aka egg bound and other within a pop are not significantly different
## So does EB:pop really even explain variance? Does it just explain more than either pop or EB alone?

#plot the interaction model's post hoc results
tukey.plot.aov<-aov(SMI ~ EB:pop, data=EB21_f)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1) #can add error bars and mean values in different plot eventually

## let's fit a model to the spring data
model.set.s <- list(one.way.s, one.s.pop, two.way.s, interaction.s)
model.names <- c("one way EB", "one way pop", "two way", "interaction")
aictab(model.set.s, modnames = model.names)
# best fit model: one-way pop
summary(one.s.pop)
# p-value: 1.923e-10 (***) no need for post-hoc 

#################################### Future thought: look at each pop and use season as independent variable #


######### Case 4: SMI ~ age + EB #######
#####SMI of VanAqua frogs by age from 2011-2022#####
# Continuous dependent variable: SMI
# Categorical independent variable: years (0-9)
#Unbalanced design, unequal sample sizes
VA_age$age <- factor(VA_age$age, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ordered = TRUE))

hist(VA_age$SMI)
qqnorm(VA_age$SMI) #a little curved
VA_age$log.SMI <- log(VA_age$SMI)
hist(VA_age$log.SMI) #a nice bell curve!
qqnorm(VA_age$log.SMI)

shapiro.test(subset(VA_age, age == "0")$SMI) # p-value: 5.864e-15
shapiro.test(subset(VA_age, age == "1")$SMI) # p-value: 1.528e-14
shapiro.test(subset(VA_age, age == "2")$SMI) # p-value: 0.00826
shapiro.test(subset(VA_age, age == "3")$SMI) # p-value: 6.781e-07
shapiro.test(subset(VA_age, age == "4")$SMI) # p-value: 0.004899
shapiro.test(subset(VA_age, age == "5")$SMI) # p-value: 0.3218 (normal)
shapiro.test(subset(VA_age, age == "6")$SMI) # p-value: 0.4897 (normal)
shapiro.test(subset(VA_age, age == "7")$SMI) # p-value: 0.0290
shapiro.test(subset(VA_age, age == "8")$SMI) # p-value: 0.9422 (normal)
shapiro.test(subset(VA_age, age == "9")$SMI) # p-value: 0.2777 (normal)
#### about half are normal and half are not. older age groups seem to be normal. 
#### smaller age groups have bigger sample size though
shapiro.test(subset(VA_age, EB == "0")$SMI) # p-value: 7.414e-10
shapiro.test(subset(VA_age, EB == "1")$SMI) # p-value: <2.2e-16
## neither of these are anywhere near normal

plot(SMI ~ age, data = VA_age) ##some major outliers here. particularly in the younger age groups. 
plot(SMI ~ EB, data = VA_age) ## lots of spread
VA_outliers <- VA_age %>% 
  group_by(age) %>% 
  identify_outliers("SMI")
view(VA_outliers) # 8 extreme outliers
### these outliers are mostly young ages... could be mix-up of year recorded or mix up btw frogs... 
EB_outliers <- VA_age %>% 
  group_by(EB) %>% 
  identify_outliers("SMI")
view(EB_outliers) # 8 "extreme outliers" - mostly not egg bound (6:2)
### test out with and without outliers 

bartlett.test(SMI ~ age, data = VA_age) #p-value: 3.5903e-14
bartlett.test(SMI ~ EB, data = VA_age) #p-value: 3.3311e-08
### variance is not equal. So we cannot use ANOVA. Friedman test is the nonparametric alternative
## but our blocks are incomplete... use Skillings-Mack test instead
# for un-balanced two-way block designs
# each block must get at least two treatments
# H0: treatments are equal
# 
#########################################################
# library(Skillings.Mack)
# skmk <- VA_age %>% 
#   Ski.Mack(SMI, groups=EB, blocks=frog_id) 
#########################################################
# cannot figure above homie out. Not a lot of references on this. 
#### Alternative plan: subset dataframe to a few frogs with the greatest number of repeated measures
##################### create a complete block
#########################################################
# data is not super complete... could have a decent subset for ages 1-3, but older ages have few and incomplete

view(VA_byage_wide)
block_13 <- VA_byage_wide %>%
  select(1:8, 11:16) %>% 
  filter_at(vars(mass1, sul1, mass2, sul2, mass3, sul3), all_vars(!is.na(.)))
view(block_13)
# all_vars ages 1-3 leaves n=43

# try with different age brackets
block_24 <- VA_byage_wide %>%
  select(1:8, 13:18) %>% 
  filter_at(vars(mass2, sul2, mass3, sul3, mass4, sul4), all_vars(!is.na(.)))
view(block_24)
# all_vars ages 1-4 leaves n=9
# all_vars ages 2-4 leaves n=14
block_67 <- VA_byage_wide %>%
  select(1:8, 21:24) %>% 
  filter_at(vars(mass6, sul6, mass7, sul7), all_vars(!is.na(.)))
view(block_67)
# ages 6-7 leaves n=7

#transform these datasets - the best blocks 
block_13 <- melt(setDT(block_13), id=1:8, measure=patterns("^mass", "^sul"),
               value.name=c("mass", "sul"), variable.name="age", na.rm = TRUE) %>% 
  relocate(age, .after = EB)  
view(block_13)
block_13$frog_id <- as.factor(block_13$frog_id) #factor 43 levels (individual frogs)
block_13$EB <- as.factor(block_13$EB) #as factor. Two levels: 0=EB 1=not (OK)
block_13$pop <- as.factor(block_13$pop) #factor with 1 level (all from VA)
block_13$source <- as.factor(block_13$source) #factor with 10 levels
block_13$birth_yr <- as.factor(block_13$birth_yr) #5 levels
str(block_13)

block_24 <- melt(setDT(block_24), id=1:8, measure=patterns("^mass", "^sul"),
                   value.name=c("mass", "sul"), variable.name="ageof", na.rm = TRUE)
block_24 <- block_24 %>% 
  mutate(age = case_when(ageof == "1" ~ "2", ageof == "2" ~ "3", ageof == "3" ~ "4")) %>%
  relocate(age, .after = EB)
view(block_24)
block_24$ageof <- NULL #no longer need this column
block_24$frog_id <- as.factor(block_24$frog_id) #factor 14 levels (individual frogs)
block_24$EB <- as.factor(block_24$EB) #as factor. Two levels: 0=EB 1=not (OK)
block_24$pop <- as.factor(block_24$pop) #factor with 1 level (all from VA)
block_24$source <- as.factor(block_24$source) #factor with 7 levels
block_24$birth_yr <- as.factor(block_24$birth_yr) #4 levels
str(block_24)

block_67 <- melt(setDT(block_67), id=1:8, measure=patterns("^mass", "^sul"),
                   value.name=c("mass", "sul"), variable.name="ageof", na.rm = TRUE)
block_67 <- block_67 %>% 
  mutate(age = case_when(ageof == "1" ~ "6", ageof == "2" ~ "7")) %>%
  relocate(age, .after = EB)
view(block_67)
block_67$ageof <- NULL 
block_67$frog_id <- as.factor(block_67$frog_id) #factor 43 levels (individual frogs)
block_67$EB <- as.factor(block_67$EB) #as factor. Two levels: 0=EB 1=not (OK)
block_67$pop <- as.factor(block_67$pop) #factor with 1 level (all from VA)
block_67$source <- as.factor(block_67$source) #factor with 10 levels
block_67$birth_yr <- as.factor(block_67$birth_yr) #5 levels
str(block_67)

##### need to convert sul to svl, and calculate SMI for each block
block_13$svl <- SULconvert(block_13$sul)
block_24$svl <- SULconvert(block_24$sul)
block_67$svl <- SULconvert(block_67$sul)
block_13$SMI <- ScaledMassIndex.VA(block_13$svl, block_13$mass)
block_24$SMI <- ScaledMassIndex.VA(block_24$svl, block_24$mass)
block_67$SMI <- ScaledMassIndex.VA(block_67$svl, block_67$mass)

######### Remember, we are blocking by age but we are actually interested in EB vs OK ############
# let's take a look at these data sets now
hist(block_13$SMI) #looks normal! 
hist(block_24$SMI) #right tail
hist(block_67$SMI) #bimodal
ggqqplot(block_13$SMI)
ggqqplot(block_24$SMI)
ggqqplot(block_67$SMI) #all three of these look quite good

# test normality
shapiro.test(block_13$SMI) # p-value: 0.01024 (not normal) - surprising because this looked normal
shapiro.test(block_24$SMI) # p-value: 0.000131 (not normal)
shapiro.test(block_67$SMI) # p-value: 0.0597 (barely normal)
shapiro.test(subset(block_13, EB = "0")$SMI) #subset by EB...
shapiro.test(subset(block_13, EB = "1")$SMI) #both p = 0.01 = not normal

####################################################################################
#### data mostly not normal, mostly not homogenous. not repeated measures but between-subject factors?
#### Mixed ANOVA? 
####################################################################################
#### don't think I can do the below kruskal wallis... only two groups. and subjects are not independent.
#### FIND correct statistical test #####################################################################
#homogeneity 
plot(SMI ~ EB, data = block_13)
bartlett.test(SMI ~ EB, data = block_13) # p-value: 0.0279 - not homogenous
plot(SMI ~ EB, data = block_24)
bartlett.test(SMI ~ EB, block_24) # p-value: 0.02811 - not homogenous 
plot(SMI ~ EB, data = block_67)
bartlett.test(SMI ~ EB, block_67) # p-value: 0.6669 - variance is homogenous
# let's do Kruskal-Wallis to be consistent across
kruskal.test(SMI ~ EB, data = block_13) #chi-squared: 1.6400, p-value: 0.2003
kruskal.test(SMI ~ EB, data = block_24) #chi-squared: 0.0284, p-value: 0.8661
kruskal.test(SMI ~ EB, data = block_67) #chi-squared: 6.6667, p-value: 0.0098
### we ACCEPT the null (no sig diffs) for block_13 and block_24 but NOT for block_67
dunnTest(SMI ~ EB, data = block_67)

#########################################################################
# select for frogs with most years of data
# 2009 GRP2-02, 2011 CBGHA#7-03, 2012 WCEM09-02, 2014 WCEM02-02 (maybe) - the oldest frogs ? lived the longest
#########################################################################
#########################################################################



#########################################################################
## Quick regression by age for wild frogs - only have age for a few #####
#########################################################################
wild_age <- read.csv("wild_withage.csv", header = TRUE)
view(wild_age)
wild_age$SMI <- ScaledMassIndex.w(wild_age$svl, wild_age$mass)
str(wild_age)
## trying to predict age based on other measurements that we have
plot(wild_age$mass, wild_age$age, pch = 16, cex = 1.3,
     xlab = "Mass (g)", ylab = "Age (yrs)")
abline(lm(age ~ mass, data = wild_age)) 
lm(age~mass, data = wild_age)
summary(lm(age~mass, data = wild_age)) #R-squared 0.4725, p-value: 0.0168

#try transformation
wild_age$log.mass <- log(wild_age$mass)
wild_age$log.svl <- log(wild_age$svl)
wild_age$log.shank <- log(wild_age$shank)
wild_age$log.SMI <- log(wild_age$SMI)
wild_age$log.age <- log(wild_age$age)

plot(wild_age$log.mass, wild_age$age, pch = 16, cex = 1.3,
     xlab = "Log Mass (g)", ylab = "Age (yrs)")
abline(lm(age ~ log.mass, data = wild_age)) 
lm(age~log.mass, data = wild_age)
summary(lm(age~log.mass, data = wild_age)) #R-squared: 0.3892, p-value: 0.0319

plot(wild_age$svl, wild_age$age, pch = 16, cex = 1.3,
     xlab = "Snout-vent length (mm)", ylab = "Age (yrs)")
abline(lm(age~svl, data = wild_age))
lm(age~svl, data = wild_age) # looks quite good
summary(lm(age~svl, data = wild_age)) # R-squared: 0.8156, p-value 0.0002117
#try transformation
plot(wild_age$log.svl, wild_age$age, pch = 16, cex = 1.3,
     xlab = "Log Snout-vent length (mm)", ylab = "Age (yrs)")
abline(lm(age~log.svl, data = wild_age))
lm(age~log.svl, data = wild_age) # looks quite good
summary(lm(age~log.svl, data = wild_age)) #R-squared is better without transformation

plot(wild_age$log.svl, wild_age$log.age, pch = 16, cex = 1.3,
     xlab = "Log Snout-vent length (mm)", ylab = "Age (yrs)")
abline(lm(log.age~log.svl, data = wild_age))
lm(log.age~log.svl, data = wild_age) 
summary(lm(log.age~log.svl, data = wild_age)) ## R-squared is only getting lower

plot(wild_age$shank, wild_age$age, pch = 16, cex = 1.3,
     xlab = "Shank length (mm)", ylab = "Age (yrs)")
abline(lm(age~shank, data = wild_age))
lm(age~shank, data = wild_age)
summary(lm(age~shank, data = wild_age)) #R-squared: 0.3283, p-value: 0.0487

plot(wild_age$SMI, wild_age$age, pch = 16, cex = 1.3,
     xlab = "Scaled Mass Index (g)", ylab = "Age (yrs)")
abline(lm(age ~ SMI, data = wild_age)) #see opposite relationship to all other regressions
lm(age ~ SMI, data = wild_age)
summary(lm(age ~ SMI, data = wild_age)) # R-squared: 0.7554, p-value: 0.000672

#### The highest R-squared thus far is from svl~age, (0.8156) ####
######## The equation of this line is: ###########################
#######################  (age)  y = 0.2049x - 10.922 #############
plot(wild_age$svl, wild_age$age, pch = 16, cex = 1.3,
     xlab = "Snout-vent length (mm)", ylab = "Age (yrs)")
abline(lm(age~svl, data = wild_age))
summary(lm(age~svl, data = wild_age))
legend("topleft", title = "Regression line", c("n=10, Adjusted R-squared: 0.8156, p-value: 0.00021"),
       fill = 'black', cex = 0.8)
legend("bottomright", title = "Equation of line", c("y = 0.205x - 10.922"),
       fill = 'blue', cex = 0.9)
title("Age of VIE-marked wild frogs")

# use this equation to predict age of wild frogs based on their measured SVL
AGE <- 
  function(x) {
    0.20487957 * x - 10.92177008
  }

wild_SMI$pred_age <- AGE(wild_SMI$svl) #temporary problem: replaced old age column with predicted age
view(wild_SMI)
plot(pred_age~svl, wild_SMI)
