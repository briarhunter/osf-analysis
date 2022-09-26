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
  "C", "H", "V", "HV", "PHV", "F"), ordered = FALSE)
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
### time to first #####################################################
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
#### hrs apx #####################################################
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
#### average touch duration #####################################################
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
### average attempt ######################################################
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
### average successful attempts #####################################################
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
### total number of amplectants #####################################################
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
########################################
####### Ultrasound grade by treatment ##
########################################
ggplot(data = VA_data, aes(x = treatment, y = grade_01, colour = treatment))+
  geom_point(size = 3, show.legend = FALSE)+
  labs(x = "Treatment", y = "Grade", title = "Day 00 Ultrasound Grade")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12)) ###########scatter plot - not best visual
####let's try a box plot - where grade is numeric
ggplot(data = VA_data, aes(x = treatment, y = as.numeric(as.character(grade_01)), fill = treatment))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Grade", title = "Day 00 Ultrasound Grade")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC","#006633", "#006699", "#FF9933", "#33CC99", "#CC9933"))+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV", "F" = "Pourya's"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

ggdotplot(VA_data, x = "frog_code", y = "grade_01", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Grade", title = "Ultrasound Grade on Day00")

################# Ultrasound grade by age ##############
ggplot(data = VA_data, aes(x = age, y = grade_01, colour = treatment, scales = "free_x"))+
  geom_jitter(size = 2)+
  facet_grid(~birth_type)+
  labs(x = "Age (yrs)", y = "Follicular Grade", title = "Day 00 Ultrasound Grade")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC","#006633", "#006699", "#FF9933", "#33CC99", "#CC9933"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

ggplot(data = VA_data, aes(x = age, y = grade_01, colour = birth_type))+
  geom_jitter(size = 2, alpha = 0.7)+
  labs(x = "Age (yrs)", y = "Follicular Grade", title = "Day 00 Ultrasound Grade of VA OSF in 2022")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC","#006633", "#006699", "#FF9933", "#33CC99", "#CC9933"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

ggdotplot(VA_data, x = "age", y = "grade_01", facet.by = "treatment")+
  labs(x = "Age (yrs)", y = "Follicular grade", title = "Ultrasound Grade on Day00")

###################################################################################
############### Amplexus behaviours by ultrasound grade1 (day 00)##################
###################################################################################
#using VA_noF dataframe because F group does not have measured amplexus behaviours
#####start with time to first contact
ggplot(data = VA_noF, aes(x = grade_01, y = to_first, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Time until first contact")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#facet grid by treatment
ggplot(data = VA_noF, aes(x = grade_01, y = to_first, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~treatment)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Time until first contact")+
  scale_x_discrete()
#### hrs apx #####################################################
ggplot(data = VA_noF, aes(x = grade_01, y = hrs_apx, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Total time in amplexus")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#
ggplot(data = VA_noF, aes(x = grade_01, y = hrs_apx, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~treatment)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Total time in amplexus")+
  scale_x_discrete()
#### touch duration #####################################################
ggplot(data = VA_noF, aes(x = grade_01, y = avg_tch_dur, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Touch Phase Duration")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#
ggplot(data = VA_noF, aes(x = grade_01, y = avg_tch_dur, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~treatment)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Touch Phase Duration")+
  scale_x_discrete()
#### attempt duration #####################################################
ggplot(data = VA_noF, aes(x = grade_01, y = avg_atmp_dur, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Attempt Phase Duration")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#
ggplot(data = VA_noF, aes(x = grade_01, y = avg_atmp_dur, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~treatment)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Attempt Phase Duration")+
  scale_x_discrete()
#### successful duration #####################################################
ggplot(data = VA_noF, aes(x = grade_01, y = avg_suc_dur, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Successful Phase Duration")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#
ggplot(data = VA_noF, aes(x = grade_01, y = avg_suc_dur, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~treatment)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Successful Phase Duration")+
  scale_x_discrete()
#### total number of amplectants #####################################################
ggplot(data = VA_noF, aes(x = grade_01, y = tot_num_amplex, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Total number of amplectants")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#
ggplot(data = VA_noF, aes(x = grade_01, y = tot_num_amplex, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~treatment)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Total number of amplectants")+
  scale_x_discrete()
#### number of successful amplectants #####################################################
ggplot(data = VA_noF, aes(x = grade_01, y = tot_num_suc, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Number of successful amplectants")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#
ggplot(data = VA_noF, aes(x = grade_01, y = tot_num_suc, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  facet_grid(~treatment)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Number of successful amplectants")+
  scale_x_discrete()


#### ultrasound grade01 by status #####################################################
ggplot(data = VA_noF, aes(x = status, y = as.numeric(as.character(grade_01)), fill=status))+
  geom_boxplot(alpha = 0.5, show.legend = FALSE)+
  labs(x = "Final Status", y = "Follicular Grade", title = "Day 01 follicular grade to final status")+
  theme_classic()+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("OK" = "Other (no mature eggs)", "ER" = "Egg retention",
                            "EB" = "Egg bound"))+
  stat_summary(fun=mean, geom="point", shape=20, size=10, color="black", fill="black")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

######not sure how I feel about below graph - showing follicular grade by status on bubble
ggplot(data = VA_data, aes(x = status, y = treatment, colour = treatment))+
  geom_count(aes(size =..prop..), colour = "lightgrey")+
  geom_count(alpha = 0.7, aes(size = ..prop.., group = treatment))+
  labs(x = "Status", y = "Treatment", title = "Final Status per Treatment")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC","#006633", "#006699", "#FF9933", "#33CC99", "#CC9933"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#try box plot instead
ggplot(VA_data, aes(x = treatment, fill = status))+
  geom_bar(alpha = 0.5)+
  labs(x = "Treatment", title = "Number of egg bound per treatment")+
  theme_classic()+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV", "F" = "Pourya's"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

####### bar chart of status by age ##########
ggplot(VA_data, aes(x = age, fill = status))+
  geom_bar(alpha = 0.5)+
  labs(x = "Age (yrs)", title = "Impact of age on final status")+
  theme_classic()+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

ggplot(VA_data, aes(x = age, fill = status))+
  geom_bar(alpha = 0.5)+
  facet_grid(~birth_type)+
  labs(x = "Age (yrs)", title = "Number of egg bound by birth type")+
  theme_classic()+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

### age ~ mass relationship - "you must be this mass to breed" 
ggplot(data = VA_noF, aes(x = as.factor(age), y = mass_spr, fill=birth_type))+
  geom_boxplot(alpha = 0.5)+
  labs(x = "Age (yrs)", y = "Mass (g)", title = "Spring mass of VanAqua females")+
  theme_classic()+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

########### grade01 ~ status relationship with or without facet by birth type ######
###### does particular day 00 grade predict egg binding down the road? #############

################### Set nice colour scheme for ppt plots #############
library(devtools)
devtools::install_github('Mikata-Project/ggthemr')
library(ggthemr)
ggthemr('dust')

lighten_swatch(0.2)
ggplot(VA_data, aes(x = grade_01, fill = status))+
  geom_bar(alpha = 0.5)+
  labs(x = "Follicular Grade", title = "Does Follicular Grade Predict Final Status?")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
    scale_fill_manual(values = c("#CC0000", "#FF9933", "#663300", "#CC9999"),
                      labels = c("EB" = "egg bound", "ER" = "egg retention", "OK" = "other"))+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

darken_swatch(0.1)
ggplot(VA_data, aes(x = age, fill = grade_01))+
  geom_bar(alpha = 0.9)+
  facet_grid(~status)+
  labs(x = "Age (yrs)")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
  scale_fill_discrete(name = "Follicular Grade")+
  scale_fill_manual(values = c("#CC0000", "#FF9933", "#663300", "#CC9999"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))


###### does heavier weight = higher follicular grade on day 00?
##is size/weight a proxy for males to determine which females are ready to oviposit?
ggplot(VA_data, aes(x = mass_spr, y = grade_01, colour = status))+
  geom_point(size = 2) ####no clear pattern but maybe flip axes

ggplot(VA_data, aes(x = grade_01, y = mass_spr, colour = status, fill = NA))+
  geom_boxplot(alpha = 0.7)+
  labs(x = "Follicular Grade", y = "Mass (g)", title = "Interaction of mass and follicular grade")+
  theme_classic()+
  scale_x_discrete(labels=c("0" = "0: not gravid", "1" = "1: eary gravid", "2" = "2: mid-gravid",
                            "3" = "3: late gravid"))+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#### expected higher grade to imply higher mass but we do not see this
  ## HOWEVER ##
##could be that most of the grade 3's were 3yo females and this is only a mass metric, not SMI
  # these 3yo females at grade 3 might be much heavier than other 3yo females at lower grade...
##subset by age
VA_3 <- subset(VA_data, age == "3")

ggplot(VA_3, aes(x = grade_01, y = mass_spr, colour = status))+
  geom_boxplot(alpha = 0.7)+
  labs(x = "Treatment", title = "Interaction of mass and follicular grade")+
  theme_classic()+
  scale_x_discrete(labels=c("0" = "0: not gravid", "1" = "1: eary gravid", "2" = "2: mid-gravid",
                            "3" = "3: late gravid"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
###add a facet by birth type
ggplot(VA_3, aes(x = grade_01, y = mass_spr, colour = status))+
  geom_boxplot(alpha = 0.7)+
  facet_grid(~birth_type)+
  labs(x = "Follicular Grade", y = "Mass (g)", title = "Interaction of mass and follicular grade")+
  theme_classic()+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
####we see our three wild caught 2019 females do not differ too much from our captive ones
# the two WC 2019's who developed mature eggs (grade 2 or 3) 
# weigh about the same as the grade 2 & 3 captive born frogs
#### BUT - the wild caught frogs do not become egg bound. They end with egg retention - but do not die
ggplot(VA_data, aes(x = grade_01, y = mass_spr, colour = status))+
  geom_boxplot(alpha = 0.7)+
  facet_grid(~prev_breed)+
  labs(x = "Follicular Grade", y = "Mass (g)", title = "Previous breeding impact on status")+
  theme_classic()+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#### see the only two females with previous breeding experience who became EGG BOUND were very high mass
## were they also old? 
lighten_swatch(0.3)
ggplot(VA_data, aes(x = as.factor(age), y = mass_spr, colour = status, fill = NA))+
  geom_boxplot(alpha = 0.6, size = 0.5)+
  labs(x = "Age (yrs)", y = "Mass (g)", title = "Is age indicative of mass?")+
  theme_classic()+
  scale_fill_discrete(labels = c("EB" = "egg bound", "ER" = "egg retention", "OK" = "other"))+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#### we actually see that for the 3yo frogs most of the egg bound ones were lower weight than the ER ones
##but in older ages we see those two obese frogs who become egg bound
#subset into 3 age groups (young - adult - old) by adding "maturity" column
VA_data <- VA_data %>% 
  mutate(maturity = case_when(age == "3" ~ "young", age == "4" ~"young",
                              age == "5" ~ "adult", age == "6" ~ "adult", 
                              age == "7" ~ "adult", age == "8" ~ "old", 
                              age == "9" ~ "old", age == "10" ~ "old", 
                              age == "11" ~ "old", age == "12" ~ "old", age == "13" ~ "old")) %>% 
  relocate(maturity, .after = age)
VA_data$maturity <- factor(VA_data$maturity, levels = c("young", "adult", "old"), ordered = TRUE)
str(VA_data)
view(VA_data)
# recreate previous graph but now view by maturity instead of age
ggplot(VA_data, aes(x = maturity, y = mass_spr, colour = status))+
  geom_boxplot(alpha = 0.7)+
  labs(x = "Age (yrs)", y = "Mass (g)", title = "Is age indicative of mass?")+
  theme_classic()+
  scale_x_discrete(labels=c("young" = "young (3-4)", "adult" = "adult (5-7)", "old" = "old (8+)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
##### see the trend of low EB mass in young frogs, high mass in older EB frogs more clearly


ggboxplot(VA_noF, x = "frog_code", y = "to_first", facet.by = "treatment", scales = "free_x")+
  labs(x = "Individual females", y = "Time (hrs)", title = "Time to first contact")


ggplot(VA_data, aes(x = grade_01, y = mass_spr, colour = status))+
  geom_point(size = VA_data$age)+
  labs(x = "Treatment", title = "Number of egg bound per treatment")+
  theme_classic()+
  scale_x_discrete(labels=c("0" = "not gravid", "1" = "eary gravid", "2" = "mid-gravid",
                            "3" = "late gravid"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

##### is age more connected to time to first than ultrasound grade ? 
ggplot(VA_data, aes(x = age, y = to_first, colour = treatment))+
  geom_jitter(alpha = 0.5)
##something to above graph but don't know what to make of it yet 