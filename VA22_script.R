getwd()

library(assertr)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(boot)
library(mvnormtest)
library(crunch)

theme_set(theme_bw())

VA_data <- read.csv("OSF_VA22females.csv", header = TRUE)
#stacked version of duration data better for repeated measures statistics
head(VA_data)
view(VA_data)
dim(VA_data) #35 rows, 24 columns
VA_data$grade_01 <- factor(VA_data$grade_01, levels = c("0", "1", "2", "3", "4"), ordered = TRUE)
VA_data$frog_code <- as.factor(VA_data$frog_code)
VA_data$frog_id <- as.factor(VA_data$frog_id)
VA_data$treatment <- factor(VA_data$treatment, levels = c(
  "C", "H", "V", "HV", "PHV"), ordered = FALSE)
VA_data$source <- as.factor(VA_data$source)
VA_data$birth_type <- as.factor(VA_data$birth_type)
VA_data$status <- as.factor(VA_data$status)
str(VA_data)
# percent_hrs is currently character - should be numeric but leave for now 
summary(VA_data)
attach(VA_data)

#visualize the data. looking for normalcy. 
hist(VA_data$mass_fall, main = "Fall Mass (g)")
ggdensity(mass_fall, main = "normality test for Fall Mass", xlab = "mass (g)") #looks bimodal
ggdensity(mass_spr, main = "normality test for Spring Mass", xlab = "mass (g)") #looks fairly normal

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

plot(mass_fall ~ mass_spr, main = "Fall Mass by Spring Mass")

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

#subset my to_first to exclude PHV since they are all -1
first_noPHV <- VA_data %>% 
  filter(to_first > 0) #has filtered out PHV and F treatments 
view(first_noPHV)

to_first_boot <- boot(first_noPHV$to_first, my_mean, 1000) #PHV and F excluded
to_first_boot
plot(to_first_boot)
boot.ci(boot_out <- to_first_boot, 
        type = c("norm", "basic", "perc", "bca"))

#filter out NA values for hrs_apx
hrs_na <- VA_data %>% 
  filter (!is.na(hrs_apx))

hrs_apx_boot <- boot(hrs_na$hrs_apx, my_mean, 1000) #no F treatment
hrs_apx_boot
plot(hrs_apx_boot)
boot.ci(boot_out <- hrs_apx_boot, 
        type = c("norm", "basic", "perc", "bca"))

#filter out NA values
num_amplex <- VA_data %>% 
  filter (!is.na(tot_num_amplex))

tot_num_amplex_boot <- boot(num_amplex$tot_num_amplex, my_mean, 1000) #total number of amplectants
tot_num_amplex_boot #no F treatment
plot(tot_num_amplex_boot)
boot.ci(boot_out <- tot_num_amplex_boot,
        type = c("norm", "basic", "perc", "bca"))

#filter out NA values
suc_amplex <- VA_data %>% 
  filter(!is.na(tot_num_suc))

tot_num_suc_boot <- boot(suc_amplex$tot_num_suc, my_mean, 1000) #no F
tot_num_suc_boot
plot(tot_num_suc_boot)
boot.ci(boot_out <- tot_num_suc_boot,
        type = c("norm", "basic", "perc", "bca"))


###do I want confidence intervals by treatment? Or leave as full pop CI's ?


##try some scatterplots to visualize
ggplot(VA_data, aes(mass_spr, to_first,
                       colour = treatment))+
  geom_point(size = 2)+
  labs(title = "time to first amplexus by spring weight", x = "mass (g)", y = "time (hrs)")

ggplot(VA_data, aes(treatment, to_first,
                    colour = treatment))+
  geom_point(size = VA_data$age)+ #size of point is associated with age of frog. cannot differentiate overlaid points
  labs(title = "time to first amplexus by treatment", caption = "size = age of frog", y = "time (hrs)")

ggplot(VA_data, aes(treatment, to_first, #time to first by status (ER, EB, OK)
                    colour = status))+
  geom_point(size = VA_data$age)+ 
  labs(title = "time to first amplexus by treatment",
       caption = "size = age of frog", y = "time (hrs)")

ggplot(VA_data, aes(treatment, hrs_apx, #total hours in amplexus with status (ER, EB, OK)
                    colour = status))+
  geom_point(size = VA_data$age)+ 
  labs(title = "total hours in amplexus by treatment",
       y = "time (hrs)", caption = "size = age of frog")

ggplot(VA_data, aes(status, hrs_apx, #total hours in amplexus with status (ER, EB, OK)
                    colour = treatment))+
  geom_point(size = VA_data$age)+ 
  labs(title = "Total hours in amplexus by status", caption = "size = age of frog",
       y = "time (hrs)")

ggplot(VA_data, aes(tot_num_amplex, hrs_apx, #total hours in amplexus on number of amplectants
                    colour = status))+
  geom_point(size = (VA_data$tot_num_suc)+2)+ 
  labs(title = "Total hours in amplexus by total number of amplectants",
       x = "number of amplectants", y = "time (hrs)",
       caption = "size = number of successful amplexes + 2 (0 is min, 2 is max)") #maybe this should be bar chart or box?

ggplot(VA_data, aes(treatment, to_first, #time to first by status (ER, EB, OK)
                    colour = status))+
  geom_point(size = VA_data$grade_01 +1)+ 
  labs(title = "time to first amplexus by treatment",
       caption = "size = (ultrasound grade 01) + 1", y = "time (hrs)")

ggplot(VA_data, aes(treatment, hrs_apx, 
                    colour = status))+
  geom_point(size = VA_data$grade_01 +1)+ 
  labs(title = "Total hours in amplexus by treatment",
       caption = "size = (ultrasound grade 01) + 1", y = "time (hrs)")

ggplot(VA_data, aes(status, hrs_apx, 
                    colour = treatment))+
  geom_point(size = VA_data$grade_01 +1)+ 
  labs(title = "Total hours in amplexus by status",
       caption = "size = (ultrasound grade 01) + 1", y = "time (hrs)")

##########################################################################
### nicer graphs of amplexus behaviour (multiple metrics) by treatment ###
##########################################################################
#filter out F treatment - not applicable for amplexus behaviours
VA_noF <- VA_data %>% 
  filter(treatment != "F")
view(VA_noF) # n=30
### time to first ###
ggplot(data = VA_noF, aes(x = treatment, y = to_first, fill=factor(treatment, 
                                                                    labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Time (hrs)", title = "Time until first contact")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC", "#006633", "#006699", "#FF9933", "#33CC99"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#facet by treatment
ggboxplot(VA_noF, x = "frog_code", y = "to_first", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Time (hrs)", title = "Time to first contact")
#### hrs apx ###
ggplot(data = VA_noF, aes(x = treatment, y = hrs_apx, fill=factor(treatment, 
                                                                   labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Time (hrs)", title = "Total time in amplexus", subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC", "#006633", "#006699", "#FF9933", "#33CC99"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#facet by treatment
ggboxplot(VA_noF, x = "frog_code", y = "hrs_apx", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Time (hrs)", title = "Total time in amplexus")
#### average touch duration ###
ggplot(data = VA_noF, aes(x = treatment, y = avg_tch_dur, fill=factor(treatment, 
                                                                  labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Time (hrs)", title = "Average Touch Phase Duration", subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC", "#006633", "#006699", "#FF9933", "#33CC99"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
ggboxplot(VA_noF, x = "frog_code", y = "avg_tch_dur", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Time (hrs)", title = "Average Touch Phase Duration")
### average attempt ###
ggplot(data = VA_noF, aes(x = treatment, y = avg_atmp_dur, fill=factor(treatment, 
                                                                      labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Time (hrs)", title = "Average Attempt Phase Duration", subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC", "#006633", "#006699", "#FF9933", "#33CC99"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
ggboxplot(VA_noF, x = "frog_code", y = "avg_atmp_dur", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Time (hrs)", title = "Average Attempt Phase Duration")
### average successful attempts ###
ggplot(data = VA_noF, aes(x = treatment, y = avg_suc_dur, fill=factor(treatment, 
                                                                      labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Time (hrs)", title = "Average Successful Phase Duration", subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC", "#006699", "#FF9933", "#33CC99"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
ggboxplot(VA_noF, x = "frog_code", y = "avg_suc_dur", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Time (hrs)", title = "Average Successful Phase Duration")
### total number of amplectants ###
ggplot(data = VA_noF, aes(x = treatment, y = tot_num_amplex, fill=factor(treatment, 
                                                                      labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Frequency", title = "Total number of amplectants", subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC","#006633", "#006699", "#FF9933", "#33CC99"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
ggboxplot(VA_noF, x = "frog_code", y = "tot_num_amplex", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Frequency", title = "Total number of amplectants")
### total number of successful amplex ###
ggplot(data = VA_noF, aes(x = treatment, y = tot_num_suc, fill=factor(treatment, 
                                                                         labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Frequency", title = "Total number of successful amplex",
       subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC","#006633", "#006699", "#FF9933", "#33CC99"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
ggboxplot(VA_noF, x = "frog_code", y = "tot_num_suc", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Frequency", title = "Total number of successful amplex")
