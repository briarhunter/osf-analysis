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

SMI_zoos <- read.csv("2021SMI.csv", header = TRUE)
head(SMI_zoos)
view(SMI_zoos)
dim(SMI_zoos) #110 rows, 8 columns
str(SMI_zoos) #set EB as factor ##set birth_yr as date/numeric?
SMI_zoos$EB <- as.factor(SMI_zoos$EB) #as factor. Two levels: 0=EB 1=not (OK)
#yr_birth <- as.Date(SMI_zoos$birth_yr, format = "%Y") # make new column or set birth_yr as numeric?
summary(SMI_zoos)

mass_zoos <- read.csv("2021mass.csv", header = TRUE)
head(mass_zoos)
view(mass_zoos)
dim(mass_zoos) #64 rows, 6 columns
str(mass_zoos) #set EB as factor
mass_zoos$EB <- as.factor(mass_zoos$EB) #as factor. Two levels: 0=EB 1=not (OK)
summary(mass_zoos)

wild_SMI <- read.csv("2021wildSMI.csv", header = TRUE)
head(wild_SMI)
view(wild_SMI)
dim(wild_SMI) #502 rows, 5 columns
str(wild_SMI)
summary(wild_SMI)

prepo <-read.csv("prepost_SMI.csv", header = TRUE) #have now added calculated svl (using linear regressions below)
head(prepo)
view(prepo)
dim(prepo) #87 rows, 12 columns
prepo$pop <- as.factor(prepo$pop)
prepo$EB <- as.factor(prepo$EB)
str(prepo)
summary(prepo)

#visualize the data. looking for normalcy. 
hist(SMI_zoos$age_yrs, main = "OSF age")
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

plot(SVL ~ mass, data = SMI_zoos, main = "SVL by mass")
plot(svl_pre20 ~ mass_pre20, data = prepo, main = "svl by mass pre-brumation 2020")
plot(svl_post21 ~ mass_post21, data = prepo, main = "svl by mass post-brumation 2021")
plot(svl_pre21 ~ mass_pre_21, data = prepo, main = "svl by mass pre-brumation 2021")

ggplot(data = mass_zoos, aes(pop, mass_2021, #total hours in amplexus with status (ER, EB, OK)
                    colour = EB))+ #how to make EB a factor for this plot and/or permanently
  geom_bar()+ ##started trying to make box plot - fix
  labs(title = "OSF mass in 2021 by population") #should be a box-whiskers or box plot
###not working - also want to turn into box plot? 

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
shapiro.test(prepo$length_post21)
shapiro.test(prepo$svl_post21)
shapiro.test(prepo$mass_pre_21) #normal
shapiro.test(prepo$svl_pre21) #nearly normal but not quite

#done visualizing. get your bootstraps fastened
#boostrap ci's for non-normal variables (aka all of them)

my_mean <- function(data, indices) { #create function to calculate mean
  return( mean(data[indices]))
}

set.seed(100) #reproducibility

#filter out NAs 
age_na <- SMI_zoos %>% 
  filter(!is.na(age_yrs))

age_boot <- boot(age_na$age_yrs, my_mean, 1000) #bootstrap the mean for age
age_boot 
plot(age_boot)
boot.ci(boot_out <- age_boot, #confidence intervals for age
        type = c("norm", "basic", "perc", "bca"))

mass_boot <- boot(SMI_zoos$mass, my_mean, 1000) 
mass_boot
plot(mass_boot)
boot.ci(boot_out <- mass_boot,
        type = c("norm", "basic", "perc", "bca"))

SVL_boot <- boot(SMI_zoos$SVL, my_mean, 1000)
SVL_boot 
plot(SVL_boot)
boot.ci(boot_out <- SVL_boot,
        type = c("norm", "basic", "perc", "bca"))

mass22_boot <- boot(mass_zoos$mass_2022, my_mean, 1000) 
mass22_boot
plot(mass22_boot)
boot.ci(boot_out <- mass22_boot,
        type = c("norm", "basic", "perc", "bca"))

mass21_boot <- boot(mass_zoos$mass_2021, my_mean, 1000) 
mass21_boot
plot(mass21_boot)
boot.ci(boot_out <- mass21_boot,
        type = c("norm", "basic", "perc", "bca"))

Wmass_boot <- boot(wild_SMI$avg_mass, my_mean, 1000) 
Wmass_boot
plot(Wmass_boot)
boot.ci(boot_out <- Wmass_boot,
        type = c("norm", "basic", "perc", "bca"))

Wsvl_boot <- boot(wild_SMI$avg_SVL, my_mean, 1000) 
Wsvl_boot
plot(Wsvl_boot)
boot.ci(boot_out <- Wsvl_boot,
        type = c("norm", "basic", "perc", "bca"))

prelen_boot <- boot(prepo$svl_pre20, my_mean, 1000) 
prelen_boot
plot(prelen_boot)
boot.ci(boot_out <- prelen_boot,
        type = c("norm", "basic", "perc", "bca"))

#filter out NAs
post_2021 <- prepo %>% 
  filter(!is.na(mass_post21))

pomass_boot <- boot(post_2021$mass_post21, my_mean, 1000) #filter out NA first
pomass_boot
plot(pomass_boot)
boot.ci(boot_out <- pomass_boot,
        type = c("norm", "basic", "perc", "bca"))

#filter out NAs
po_2021 <- prepo %>% 
  filter(!is.na(svl_post21))

posvl_boot <- boot(po_2021$svl_post21, my_mean, 1000) 
posvl_boot
plot(posvl_boot)
boot.ci(boot_out <- posvl_boot,
        type = c("norm", "basic", "perc", "bca"))

#filter out NAs
pre_2021 <- prepo %>% 
  filter(!is.na(svl_pre21))

pre_boot <- boot(pre_2021$svl_pre21, my_mean, 1000) 
pre_boot
plot(pre_boot)
boot.ci(boot_out <- pre_boot,
        type = c("norm", "basic", "perc", "bca"))

#for normal distribution (pre-brum mass 2020, pre-brum mass 2021)
prebru_m_mean <- mean(prepo$mass_pre20)
prebru_m_mean #45.89575
prebru_m_error <- qnorm(0.975)*sd(prepo$mass_pre20)/sqrt(length(prepo$mass_pre20))
prebru_m_error
left <- prebru_m_mean - prebru_m_error
right <- prebru_m_mean + prebru_m_error
print(c(left,right)) #(41.74126 50.05024)

prebru_mean <- mean(pre_2021$mass_pre_21) #use same NA filtered out dataset as svl_pre21
prebru_mean #52.5812
prebru_error <- qnorm(0.975)*sd(pre_2021$mass_pre_21)/sqrt(length(pre_2021$mass_pre_21))
prebru_error
left <- prebru_mean - prebru_error
right <- prebru_mean + prebru_error
print(c(left,right)) #(49.10911 56.05329)

####linear regressions to predict SVL for VanAqua SUL values
SULreg <- read.csv("SULregression.csv", header = TRUE) #only VanAqua frogs
head(SULreg)
view(SULreg)
dim(SULreg) #30 rows, 5 columns
str(SULreg)
summary(SULreg)
#eventually hopefully replace these SUL/SVL measurements with ones taken on same day from same individuals in triplicate
#for now using what we have (as of Sept 2022)

#linear regression of SVL by SUL - to predict SVL from SUL
lmSVL1 <- lm(SVL_spr2021 ~ SUL_2020, data = SULreg)
summary(lmSVL1)

lmSVL2 <- lm(SVL_Nov2021 ~ SUL_2020, data = SULreg)
summary(lmSVL2)
#plot variables again but add regression line
plot(SVL_spr2021 ~ SUL_2020 , data = SULreg, main = "2021 spring SVL by 2020 SUL") # n=30
abline(lmSVL1)
plot(lmSVL1$residuals, pch = 16, col = "red", main = "sprSVL ~ SUL residuals") #plot the residuals - do they look random?
plot(cooks.distance(lmSVL1), pch = 16, col = "blue", main = "sprSVL ~ SUL cooks distance") #plot cooks distance - see at least one point clearly not following pattern
#2017 WCEM04-05 has a spring SVL of 57.86mm which is much lower than most others - might be outlier? Check if recorded correctly

plot(SVL_Nov2021 ~ SUL_2020, data = SULreg, main = "2021 fall SVL by 2020 SUL") # n=20
abline(lmSVL2)
plot(lmSVL2$residuals, pch = 16, col = "red") #these residuals look random

#SVL_spr2021 on SUL is the better model, but still not great. 
#let's see if any transformations improve the model fit
ggplot(SULreg, aes(SUL_2020, SVL_spr2021))+
  geom_point()+
  scale_x_log10() + scale_y_log10() #doesn't look any different

#try a log transformation on the regression
lmSVL3 <- lm(log(SVL_spr2021) ~ log(SUL_2020), data = SULreg)
summary(lmSVL3) #improves the model but only slightly
plot(log(SVL_spr2021) ~ log(SUL_2020), data = SULreg, main = "2021 spring SVL by 2020 SUL") 
abline(lmSVL3)
plot(lmSVL3$residuals, pch = 16, col = "red", main = "sprSVL ~ SUL residuals") 
plot(cooks.distance(lmSVL3), pch = 16, col = "blue", main = "sprSVL ~ SUL cooks distance") #now we see two influential points...

#using non-transformed regression our equation is: SVL = 0.889x + 5.1677, R2 = 0.4573
#create predicted SVL column using equation and SUL values from SULreg
lmSVL1 <- lm(SVL_spr2021 ~ SUL_2020, data = SULreg)
summary(lmSVL1) #lmSVL1 is our model
predSVL <- data.frame(SUL_2020 = prepo[43:87, 6]) #turn SUL values into data frame so we can predict their SVL
predicted <- predict(lmSVL1, newdata = predSVL)
predSVL <- predSVL %>% 
  mutate(SUL_SVL = predicted) #predict SVL from previously defined SUL values

# prepo <- prepo %>% 
#   mutate(SUL_SVL = rbind(prepo[1:42, 6], predicted) %>% 
#   relocate(SUL_SVL, .after = 6)

####trying to bind GVZ SVL values (rows 1:42) to predicted VA values for 43:87 but not working
##may need to just edit data in excel and sort out reproducible method later

predSVL2 <- data.frame(SUL_2020 = prepo[43:87, 8]) #turn SUL 2021 values into data frame so we can predict their SVL
predicted21 <- predict(lmSVL1, newdata = predSVL2)
predSVL <- predSVL %>% 
  mutate(SUL21_SVL = predicted21) #added column of predicted 2021 SUL values from VanAqua
view(predSVL)


####calculate some SMI values - ref. Chen-Pan Liao R code inspired by Peig & Green (2009)

#linear regressions of mass by length

plot(mass_pre20 ~ svl_pre20, data = prepo, main = "mass by svl pre-brumation 2020")

ggplot(data = prepo, aes(svl_pre20, mass_pre20, colour = pop))+
  geom_point(size = 2)+ #no obvious outliers
  geom_line(aes(svl_pre20, mass_pre20), data = lmSVL1)
  labs(title = "pre-brumation 2020 mass by svl", x = "svl (mm)", y = "mass (g)") #add a line of best fit?



plot(mass_post21 ~ svl_post21, data = prepo, main = "mass by svl post-brumation 2021")
plot(mass_pre_21 ~ svl_pre21, data = prepo, main = "mass by svl pre-brumation 2021")

###########################
#Use steph's code 
###########################

str(prepo)
#subset prepo data frame by treatment
prepo.GVZ <- subset(prepo, pop == "GVZ")
str(prepo.GVZ)
prepo.VA <- subset(prepo, pop == "VA")
str(prepo.VA)

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
legend("topright", title = "Regression lines", c("n=55, Adjusted R-squared: 0.04873, p-value: 0.0576",
                                                "n=40, Adjusted R-squared: -0.02128, p-value: 0.6675",
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
legend("topright", title = "Regression lines", c("n=55, Adjusted R-squared: 0.08053, p-value: 0.02025",
                                                    "n=40, Adjusted R-squared: -0.02392, p-value: 0.7672",
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

################################################################################
# Calculate scaled mass index (SMI)
################################################################################

# Scaled mass index, Pieg & Green 2009

# Scaled mass index (SMI): ^Mi = Mi (Lo/Li)^bsma
# where: ^Mi is scaled mass: predicted body mass for individual i when the linear body measure is standardized to L0.
#         Mi is mass of individual i
#         Lo is an arbitraty value of length (e.g. arithmetric mean of L for the population of study)
#         Li is the linear body measurement of individual i
#         bsma is the scaling exponent: calculated indirectly by dividing the slope from an OLS regression (bOLS) by the Pearson’s correlation coefficient r (LaBarbera 1989), or directly using online software (Bohonak and van der Linde 2004).

# calculate bsma manually (described in Pieg & Green 2009 p.1886)