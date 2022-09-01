getwd()

library(assertr)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(boot)
library(mvnormtest)

theme_set(theme_bw())

VA_data <- read.csv("OSF_VA22females.csv", header = TRUE)
#stacked version of duration data better for repeated measures statistics
head(VA_data)
view(VA_data)
dim(VA_data) #35 rows, 24 columns
str(VA_data) 
summary(VA_data)
attach(VA_data)

#visualize the data. looking for normalcy. 
hist(VA_data$mass_fall, main = "Fall Mass (g)")
ggdensity(mass_fall, main = "normality test for Fall Mass", xlab = "mass (g)")

hist(VA_data$mass_fall[VA_data$treatment == "C"], main = "Fall Mass in treatment C")
hist(VA_data$mass_fall[VA_data$treatment == "H"], main = "Fall Mass in treatment H")
hist(VA_data$mass_fall[VA_data$treatment == "V"], main = "Fall Mass in treatment V")
hist(VA_data$mass_fall[VA_data$treatment == "HV"], main = "Fall Mass in treatment HV")
hist(VA_data$mass_fall[VA_data$treatment == "PHV"], main = "Fall Mass in treatment PHV")

#I cannot get this multivariate shapiro test to work: mshapiro.test(VA_data[, 8:10])
shapiro.test(VA_data$age)
shapiro.test(VA_data$mass_fall)#can add [VA_data$treatment == "C"] for treatment p-value
shapiro.test(VA_data$mass_spr)
shapiro.test(VA_data$eggs_laid)
shapiro.test(VA_data$to_first)
shapiro.test(VA_data$hrs_apx)
shapiro.test(VA_data$avg_tch_dur)
shapiro.test(VA_data$avg_atmp_dur)
shapiro.test(VA_data$avg_suc_dur)
shapiro.test(VA_data$tot_num_amplex)
shapiro.test(VA_data$tot_num_suc) 
#messy messy - clean up later

#run shapiro tests by treatment
age_SW <- aggregate(formula = age ~ treatment,
          data = VA_data,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
view(age_SW)

fall_SW <- aggregate(formula = mass_fall ~ treatment,
          data = VA_data,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
view(fall_SW)

spr_SW <- aggregate(formula = mass_spr ~ treatment,
          data = VA_data,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
view(spr_SW)

hrs_SW <- aggregate(formula = hrs_apx ~ treatment,
          data = VA_data,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
view(hrs_SW)

num_amp_SW <- aggregate(formula = tot_num_amplex ~ treatment,
          data = VA_data,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
view(num_amp_SW)

num_suc_SW <- aggregate(formula = tot_num_suc ~ treatment,
          data = VA_data,
          FUN = function(x) {y <- shapiro.test(x); c(y$statistic, y$p.value)})
#check by individual treatments - says "all 'x' values are identical"

#only mass_fall and mass_spr are really normal - double check this on qq plot
ggqqplot(VA_data$mass_fall)
ggqqplot(VA_data$mass_spr) #looks good. these two can be treated as normal distribution

####calculate confidence intervals

#should I use t-distribution or normal distribution? t-test assumes working with sample stddev instead of exact
#qt(0.975,df=length(VA_data$mass_fall))*sd ##alternative t-test method
mass_fall_mean <- mean(VA_data$mass_fall)
mass_fall_mean #51.21686
mass_fall_error <- qnorm(0.975)*sd(VA_data$mass_fall)/sqrt(length(VA_data$mass_fall))
mass_fall_error
left <- mass_fall_mean - mass_fall_error
right <- mass_fall_mean + mass_fall_error
print(c(left,right)) # (44.99385, 57.43986)

mass_spr_mean <- mean(VA_data$mass_spr)
mass_spr_mean #47.53971
mass_spr_error <- qnorm(0.975)*sd(VA_data$mass_spr)/sqrt(length(VA_data$mass_spr))
mass_spr_error
left <- mass_spr_mean - mass_spr_error
right <- mass_spr_mean + mass_spr_error
print(c(left,right)) # (41.56360, 53.51583)

#boostrap ci's for non-normal variables: age, eggs_laid, to_first, hrs_apx, tot_num_amplex, tot_num_suc

#use boot function to find the R bootstrap of the mean
my_mean <- function(data, indices) { #create function to calculate mean
  return( mean(data[indices]))
}

set.seed(100) #reproducibility
age_boot <- boot(VA_data$age, my_mean, 1000) #bootstrap the mean for age
age_boot
plot(age_boot)
boot.ci(boot_out <- age_boot, #confidence intervals for age
        type = c("norm", "basic", "perc", "bca"))

eggs_laid_boot <- boot(VA_data$eggs_laid, my_mean, 1000) #eggs_laid
eggs_laid_boot
plot(eggs_laid_boot)
boot.ci(boot_out <- eggs_laid_boot,
        type = c("norm", "basic", "perc", "bca"))


to_first_boot <- boot(VA_data$to_first, my_mean, 1000) #need to omit PHV - all -1 values
to_first_boot
plot(to_first_boot)
boot.ci(boot_out <- to_first_boot, 
        type = c("norm", "basic", "perc", "bca"))
###error. did not work. come back

hrs_apx_boot <- boot(VA_data$hrs_apx, my_mean, 1000) #not working???
hrs_apx_boot
plot(hrs_apx_boot)
boot.ci(boot_out <- hrs_apx_boot, #confidence intervals for touch phase
        type = c("norm", "basic", "perc", "bca"))
### error for plot and boot.ci

tot_num_amplex_boot <- boot(VA_data$tot_num_amplex, my_mean, 1000) #total number of amplectants
tot_num_amplex_boot
plot(tot_num_amplex_boot)
boot.ci(boot_out <- tot_num_amplex_boot, #error: all values of t1* are NA
        type = c("norm", "basic", "perc", "bca"))
###error

tot_num_suc_boot <- boot(VA_data$tot_num_suc, my_mean, 1000) #bootstrap the mean for touch phase
tot_num_suc_boot
plot(tot_num_suc_boot)
boot.ci(boot_out <- tot_num_suc_boot, #confidence intervals for touch phase
        type = c("norm", "basic", "perc", "bca"))
###error

#bunch of errors in calculating CI's that I need to come back to
#maybe need to separate into treatments to calculate CI? 
#why so many 'equal' values or NA?





##try some scatterplots to visualize
ggplot(VA_data, aes(mass_spr, to_first,
                       colour = treatment))+
  geom_point(size = 2)+
  labs(title = "time to first amplexus by weight")

ggplot(VA_data, aes(treatment, to_first,
                    colour = treatment))+
  geom_point(size = VA_data$age)+ #size of point is associated with age of frog. cannot differentiate overlaid points
  labs(title = "time to first amplexus by treatment")

ggplot(VA_data, aes(treatment, to_first, #time to first by status (ER, EB, OK)
                    colour = status))+
  geom_point(size = VA_data$age)+ 
  labs(title = "time to first amplexus by treatment")

ggplot(VA_data, aes(treatment, hrs_apx, #total hours in amplexus with status (ER, EB, OK)
                    colour = status))+
  geom_point(size = VA_data$age)+ 
  labs(title = "hours in amplexus by treatment")

ggplot(VA_data, aes(treatment, hrs_apx, #total hours in amplexus with status (ER, EB, OK)
                    colour = status))+
  geom_point(size = VA_data$tot_num_amplex)+ 
  labs(title = "Hours in amplexus by treatment", subtitle = "size as number of amplectants")
#size as the total number of amplexus phases - doesn't quite work. Some have too many

ggplot(VA_data, aes(tot_num_amplex, hrs_apx, #total hours in amplexus on number of amplectants
                    colour = status))+
  geom_point(size = VA_data$age)+ 
  labs(title = "Hours in amplexus by total number of amplectants") #maybe this should be bar chart or box?

