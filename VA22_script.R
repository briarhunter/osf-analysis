getwd()

library(assertr)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(boot)
library(mvnormtest)
library(crunch)
library(EnvStats)

theme_set(theme_bw())

VA_data <- read.csv("OSF_VA22females.csv", header = TRUE)
#stacked version of duration data better for repeated measures statistics
head(VA_data)
view(VA_data)
dim(VA_data) #35 rows, 24 columns
VA_data$grade_01 <- factor(VA_data$grade_01, levels = c("0", "1", "2", "3", "4"), ordered = TRUE)
VA_data$grade_02 <- factor(VA_data$grade_02, levels = c("0", "1", "2", "3", "4"), ordered = TRUE)
VA_data$grade_03 <- factor(VA_data$grade_03, levels = c("0", "1", "2", "3", "4"), ordered = TRUE)
VA_data$grade_04 <- factor(VA_data$grade_04, levels = c("0", "1", "2", "3", "4"), ordered = TRUE)
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

########################## calculate confidence intervals ###########################################
#####################################################################################################
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

################### boostrap ci's for non-normal variables: ############################
################### age, eggs_laid, to_first, hrs_apx, tot_num_amplex, tot_num_suc #####

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

######################################################################################
######### Visualize on some scatterplots #############################################
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

##########################################################################
### nicer graphs of amplexus behaviour (multiple metrics) by treatment ###
##########################################################################
#filter out F treatment - not applicable for amplexus behaviours
VA_noF <- VA_data %>% 
  filter(treatment != "F")
view(VA_noF) # n=30
########################### Amplexus behaviours by treatment ###############################
### time to first contact ##################################################################
ggplot(data = VA_noF, aes(x = treatment, y = to_first, fill=factor(treatment, 
                                                                   labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Time (hrs)", title = "Time until first contact")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#### total hours in amplexus #####################################################
ggplot(data = VA_noF, aes(x = treatment, y = hrs_apx, fill=factor(treatment, 
                                                                   labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Time (hrs)", title = "Total time in amplexus", subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

##### if we are interested in amplexus phases by treatment... could do simple distribution ############
#### would need to convert amplexus data to long form #################################################
# but maybe we are no longer interested in answering this question? ## 

### total number of amplectants #####################################################
ggplot(data = VA_noF, aes(x = treatment, y = tot_num_amplex, fill=factor(treatment, 
                                                                      labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Frequency", title = "Total number of amplectants", subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
### total number of successful amplex #################################################
ggplot(data = VA_noF, aes(x = treatment, y = tot_num_suc, fill=factor(treatment, 
                                                                         labels = c("C", "H", "V", "HV", "PHV"))))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Frequency", title = "Total number of successful amplex",
       subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

########################################
####### Ultrasound grade by treatment ##
######################################## Unclear if we actually want this comparison anymore. maybe irrelevant question/answer 
ggplot(data = VA_data, aes(x = treatment, y = as.numeric(as.character(grade_01)), fill = treatment))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Treatment", y = "Grade", title = "Day 00 Ultrasound Grade")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV", "F" = "Pourya's"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
### if you cycle through the above graphic but swap out the ultrasound grades (01 for 02, 03, 04) 
### you see HV lags behind - the other treatments have evened out by grade 02-03 but HV is still much lower
### just happens to have all the immature frogs or caused all the immature frogs? 

################# Ultrasound grade by age - faceted by Captive vs Wild born ##############
ggplot(data = VA_data, aes(x = age, y = grade_01, colour = treatment, scales = "free_x"))+
  geom_jitter(size = 2)+
  facet_grid(~birth_type)+
  labs(x = "Age (yrs)", y = "Follicular Grade", title = "Day 00 Ultrasound Grade")+
  theme_classic()+
  scale_fill_brewer(palette = "Dark2")+ #### the proper colour palette is not working - not applying to the points ? 
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

ggplot(data = VA_data, aes(x = age, y = grade_01, colour = birth_type))+
  geom_jitter(size = 2, alpha = 0.7)+
  labs(x = "Age (yrs)", y = "Follicular Grade", title = "Day 00 Ultrasound Grade of VA OSF in 2022")+
  theme_classic()+
  scale_fill_manual(values = c("#9933CC","#006633", "#006699", "#FF9933", "#33CC99", "#CC9933"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

###################################################################################
############### Amplexus behaviours by ultrasound grade1 (day 00)##################
###################################################################################
#using VA_noF dataframe because F group does not have measured amplexus behaviours



#####start with time to first contact
ggplot(data = VA_noF, aes(x = grade_01, y = to_first, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Time until first contact")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))+
  stat_n_text()

#### hrs apx #####################################################
## how long did females of each grade 01 spend in amplexus over first 59 hours? ## 

# if substituting below grade_01 with other grades, use: data = subset(VA_noF, !is.na(grade_02))
# this will remove NA values for that particular grade
ggplot(data = VA_noF, aes(x = grade_01, y = hrs_apx, fill=grade_01))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Time (hrs)", title = "Total time spent in amplexus by grade 01",
       subtitle = "over first 59 hours")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#### convert amplexus data to long form and look at distributions ############################
#### does grade 01 predict whether female spends more time in touch or attempt or successful phase? #### 

#### total number of amplectants #####################################################
ggplot(data = subset(VA_noF, !is.na(grade_04)), aes(x = grade_04, y = tot_num_amplex, fill=grade_04))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Count", title = "Total number of amplectants")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#### number of successful amplectants #####################################################
ggplot(data = subset(VA_noF, !is.na(grade_04)), aes(x = grade_04, y = tot_num_suc, fill=grade_04))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Follicular Grade", y = "Count", title = "Number of successful amplectants")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("0" = "0 (non-gravid)", "1" = "1 (early-gravid)",
                            "2" = "2 (mid-gravid)", "3" = "3 (late-gravid)"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
##### run through grades 01-04, what we see: 
## grade 01: 0-1 seem to be sig lower # than 2-3
## grade 02: maybe none sig? weird spread
## 03: 4 is higher but maybe not sig cause of spread? 04: " 

#######################################################################################
#### ultrasound grade01 by status #####################################################
ggplot(data = VA_noF, aes(x = status, y = as.numeric(as.character(grade_01)), fill=status))+
  geom_boxplot(alpha = 0.5, show.legend = FALSE)+
  labs(x = "Final Status", y = "Follicular Grade", title = "Day 01 follicular grade to final status")+
  theme_classic()+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.title = element_blank())+
  scale_x_discrete(labels=c("OK" = "Other (no mature eggs)", "ER" = "Egg retention",
                            "EB" = "Egg bound"))+
  stat_summary(fun=mean, geom="point", shape=20, size=10, color="black", fill="black")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
### maybe not the best graph - shows the trend very clearly but doesn't tell the full story

### final status distribution by treatment 
## Q: does overwintering treatment group (i.e. exposure to males) impact the final status of females? 
## does greater exposure to males lead to more egg binding? Or does hormonal exposure (via water) lead to more egg binding than just visual?
ggplot(VA_data, aes(x = treatment, fill = status))+
  geom_bar(alpha = 0.5)+
  labs(x = "Treatment", title = "Number of egg bound per treatment")+
  theme_classic()+
  scale_fill_brewer(palette = "Dark2")+ ###could maybe flip colours so EB is red... or set different colours per status
  scale_x_discrete(labels=c("C" = "Control", "V" = "Visual", "H" = "Hormone",
                            "HV" = "H + V", "PHV" = "Physical + HV", "F" = "Pourya's"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#### we see no frogs in the Visual treatment group became egg bound. All had egg retention but no egg binding. Hormone group had lots of egg binding
## but Control group also had lots of egg binding... the most interesting trend seems to be the lack of EB in the Visual group (only one without) 

####### bar chart of status by age ##########
ggplot(VA_data, aes(x = age, fill = status))+
  geom_bar(alpha = 0.5)+
  labs(x = "Age (yrs)", title = "Impact of age on final status")+
  theme_classic()+
  scale_fill_brewer(palette = "Dark2")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
#### see most EB were 3yo frogs - in their first breeding season. 

## is there a difference in Captive born vs Wild caught frogs? 
ggplot(VA_data, aes(x = age, fill = status))+
  geom_bar(alpha = 0.5)+
  facet_grid(~birth_type)+
  labs(x = "Age (yrs)", title = "Number of egg bound by birth type")+
  theme_classic()+
  scale_fill_brewer(palette = "Dark2")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
##### we see most EB frogs were captive born... but what does this really tell us? 
## most wild caught frogs still grew up in captivity and probably still had similar mass as captive born frogs... 

###########################################################################################
######### age ~ mass relationship - "you must be this mass to breed" ######################
## faceting by origins of Captive vs Wild born 
ggplot(data = VA_noF, aes(x = as.factor(age), y = mass_spr, colour = birth_type))+
  geom_boxplot(alpha = 0.5)+
  labs(x = "Age (yrs)", y = "Mass (g)", title = "Spring mass of VanAqua females")+
  theme_classic()+
  scale_fill_brewer(palette = "Dark2")+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

########### grade01 ~ status relationship with or without facet by birth type ######
###### does particular day 00 grade predict egg binding down the road? #############

################### Set nice colour scheme for ppt plots #############
# library(devtools)
# devtools::install_github('Mikata-Project/ggthemr')
# library(ggthemr)
# ggthemr('dust')
###### I do like this colour scheme (matches my earthy palatte) but not always the easiest contrasts to visualize 

##### plotting with earthy ('dust') colour scheme and background to match typical ppt background
ggplot(VA_data, aes(x = grade_01, fill = status))+
  geom_bar(alpha = 0.5)+
  labs(x = "Follicular Grade", title = "Does Follicular Grade Predict Final Status?")+
  theme_classic()+
  scale_fill_brewer(palette = "Pastel2")+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
    scale_fill_manual(values = c("#CC0000", "#FF9933", "#663300", "#CC9999"),
                      labels = c("EB" = "egg bound", "ER" = "egg retention", "OK" = "other"))+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#### same plot as above but with brewer palate of "Dark2" - still has background colour from ppt
ggplot(VA_data, aes(x = grade_01, fill = status))+
  geom_bar(alpha = 0.5)+
  labs(x = "Follicular Grade", title = "Does Initial Follicular Grade Predict Final Status?")+
  theme_classic()+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
  scale_fill_manual(labels = c("EB" = "egg bound", "ER" = "egg retention", "OK" = "other"))+
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

### Same plot as above but subbing in grades 02-04 without NA values included
# this is because NA represents females who have already died - should not be included
ggplot(data = subset(VA_noF, !is.na(grade_02)), aes(x = grade_02, fill = status))+
  geom_bar(alpha = 0.5)+
  labs(x = "Follicular Grade", title = "Does Follicular Grade 02 Predict Final Status?")+
  theme_classic()+
  scale_fill_manual(labels = c("EB" = "egg bound", "ER" = "egg retention", "OK" = "other"))+
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(breaks = c(2, 4, 6, 8, 10, 12))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
########################################################
#### we have designated the final u/s image for females who become egg bound or ER as a grade 4 ....
#### but this means we don't see what the grade actually was leading up to egg binding... i.e. if it was 2 or 3
### may need to just add a column for 'designated final grade' etc to indicate grade 4s, but still have them count as 2 or 3
########################################################
# do some more complicated model (GLM?) to see which u/s grade (01-04) best predicts final status 

### what is the distribution of grade 01 across age and status? 
## i.e. do age and grade predict final status? ##
ggplot(VA_data, aes(x = age, fill = grade_01))+
  geom_bar(alpha = 0.9)+
  facet_grid(~status)+
  labs(x = "Age (yrs)")+
  theme_classic()+
  scale_fill_discrete(name = "Follicular Grade")+
  scale_fill_brewer(palette = "Pastel2")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))


###### does heavier weight = higher follicular grade on day 00?
##is size/weight a proxy for males to determine which females are ready to oviposit?
ggplot(VA_data, aes(x = grade_01, y = mass_spr, colour = status, fill = NA))+
  geom_boxplot(alpha = 0.7)+
  labs(x = "Follicular Grade", y = "Mass (g)", title = "Interaction of mass and follicular grade")+
  theme_classic()+
  scale_x_discrete(labels=c("0" = "0: not gravid", "1" = "1: eary gravid", "2" = "2: mid-gravid",
                            "3" = "3: late gravid"))+
  theme(panel.background = element_rect(fill = "#F6F0ED"), ### background fill to match ppt 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
  scale_fill_brewer(palette = "Dark2")+ #don't think this fill is actually applied here
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))

#### expected higher grade to imply higher mass but we do not see this
  ## HOWEVER ##
##could be that most of the grade 3's were 3yo females and this is only a mass metric, not SMI
  # these 3yo females at grade 3 might be much heavier than other 3yo females at lower grade...
ggplot(VA_data, aes(x = as.factor(age), y = mass_spr, colour = status, fill = NA))+
  geom_boxplot(alpha = 0.6, size = 0.5)+
  labs(x = "Age (yrs)", y = "Mass (g)", title = "Is age indicative of mass?")+
  theme_classic()+
  scale_fill_discrete(labels = c("EB" = "egg bound", "ER" = "egg retention", "OK" = "other"))+
  theme(panel.background = element_rect(fill = "#F6F0ED"), 
        plot.background = element_rect(fill = "#F6F0ED"), 
        legend.background = element_rect(fill = "#F6F0ED"),
        legend.box.background = element_rect(fill = "#F6F0ED"))+
  scale_fill_brewer(palette = "Dark2")+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))+
  stat_n_text()
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
  theme(text = element_text(size = 12))+
  stat_n_text()
##### see the trend of low EB mass in young frogs, high mass in older EB frogs more clearly

### can we visualize the interact of mass and grade 01 on what the final status is? 
# previously predicted higher grade (i.e. 3) would correlate to higher mass (vice-versa) but did not see this in the data. 
# look at it a different way than previous
ggplot(VA_data, aes(x = grade_01, y = mass_spr, colour = status))+
  geom_point(size = VA_data$age)+
  labs(x = "Follicular Grade", y = "Mass (g)", title = "Mass and grade01 interaction on status",
       caption = "size = age")+
  theme_classic()+
  scale_x_discrete(labels=c("0" = "0: not gravid", "1" = "1: eary gravid", "2" = "2: mid-gravid",
                            "3" = "3: late gravid"))+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))
# we see, again, no real correlation between higher grade = higher mass
# but see our two high-mass frogs who were EB clearly above the others .. but one of them tarted at a grade 1... 
# already very high mass at grade 1 would suggest it was not the reproductive load but just excess fat ... also one of our oldest frogs

##### is age more connected to time to first than ultrasound grade ? 
ggplot(VA_data, aes(x = age, y = to_first, colour = treatment))+
  geom_jitter(alpha = 0.5)
##something to above graph but don't know what to make of it yet 


################################################################################################
####### Let's do some ANOVAs and/or Kruskal-Wallis tests to see if differences (in mass) are significant
################################################################################################
####### Compare mass across age groups, grade01s, and final status 

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

######### Case 1: spring mass ~ age #######
# Continuous dependent variable: mass
# Categorical independent variables: age 3:13
#Unbalanced design, unequal sample sizes

##################################################################
# Test assumptions
hist(VA_data$mass_spr) #looking for bell shape
# looks almost normal, slight divet in the top makes it almost bimodal 

qqnorm(VA_data$mass_spr) #looking for straight line
#looks quite good

#Shapiro Wilks test: null hypothesis = variance of data is normally distributed
#check by each age group
shapiro.test(subset(VA_data, age == "3")$mass_spr) # p-value: 0.9207
shapiro.test(subset(VA_data, age == "4")$mass_spr) # p-value: 0.1619
shapiro.test(subset(VA_data, age == "5")$mass_spr) # p-value: 0.6329
shapiro.test(subset(VA_data, age == "6")$mass_spr) # p-value: 0.2081
# shapiro.test(subset(VA_data, age == "7")$mass_spr) # not a large enough sample size for the other age groups (n has to be >= 3)

### all of these are normally distributed, BUT sample sizes are so small, these results might be biased. Can we still do ANOVA? 
# cannot do bartlett or ANOVA (?) on sample sizes less than 2 
# need to group ages by maturity level

##################################################################
####### let's try Case 2: spring mass ~ maturity #################
# does spring mass differ by maturity (age group)?
# Continuous dependent variable: mass
# Categorical independent variables: maturity = young, adult, old
#Unbalanced design, unequal sample sizes

##################################################################
# Test assumptions = the same as above

# check Shapiro test by maturity groups, not individual ages
shapiro.test(subset(VA_data, maturity == "young")$mass_spr) # p-value: 0.5375
shapiro.test(subset(VA_data, maturity == "adult")$mass_spr) # p-value: 0.9028
shapiro.test(subset(VA_data, maturity == "old")$mass_spr) # p-value: 0.5645
# all p-values are >0.05 so all variances are normally distributed, but again, still fairly small sample sizes
# young n=19, adult n=12, old n=4

##################################################################
# Test assumption of homogeneity of variance
# Null hypothesis of Bartlett: variance is equal

plot(mass_spr~maturity, data = VA_data) #visual check - do not see any significant outliers
bartlett.test(mass_spr~maturity, data = VA_data)
# p-value: 0.6778 - accept the null hypothesis
# the variance is homogenous so we can proceed to the ANOVA

### one-way ANOVA
oneway.test(mass_spr~maturity, data = VA_data, var.equal = TRUE)
# F = 17.412, num df = 2, denom df = 32, p-value = 7.646e-06
# easy method to switch from ANOVA (where variances are equal) to Welch ANOVA (where var.equal = FALSE)

### second method of ANOVA using aov (also a one-way test)
mat_aov <- aov(mass_spr~maturity, data = VA_data)
summary(mat_aov)
## this method prints the full ANOVA table and results can be saved for use later in post-hoc tests
# all the variation not explained by the independent variables (maturity) is called residual variance

# Diagnostic plots - evaluate normality assumption
plot(mat_aov, which = 2, add.smooth = FALSE) #Q-Q plot, look for straight line
plot(mat_aov, which = 3, add.smooth = TRUE) #constant variance assumption with scale-location plot - look for no patterns
plot(mat_aov, which = 4, add.smooth = FALSE) # evaluate if influential points with Cook's distance plot - look for outliers
###### not entirely sure of interpretation of these diagnostic plots... might need to look into and come back to ############

# Interpretation of results
# since p-value is < 0.05, we reject the null hypothesis so we conclude at least one age group's mass is different than the others

# Tukey HSD test - post-hoc test
library(multcomp)
post_mat_aov <- glht(mat_aov, linfct = mcp(maturity = "Tukey"))
summary(post_mat_aov)
## adult and young significantly differ (***), old - young (**), old-adult not significantly different (p-value: 0.99898)
par(mar = c(3, 8, 3, 3))
plot(post_mat_aov)

# alternative method using: TukeyHSD(mat_aov) - not as much information provided

### plot the ANOVA and post-hoc tests
library(ggstatsplot)
### significance of mass differences between maturity groups (young, adult, old)
ggbetweenstats(
  data = VA_data,
  x = maturity,
  y = mass_spr,
  type = "parametric", # ANOVA or Kruskal-Wallis
  var.equal = TRUE, # ANOVA or Welch ANOVA
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

# we want to see if differences in mass by status (EB, ER, OK) are significant - but need to subset by age groups first
## subset by maturity
VA_young <- subset(VA_data, maturity == "young")
view(VA_young) # n=19
VA_adult <- subset(VA_data, maturity == "adult")
view(VA_adult) # n=12
VA_old <- subset(VA_data, maturity == "old")
view(VA_old) # n=4

#try subsets of 3yo against all others
VA_3 <- subset(VA_data, age == 3)
view(VA_3) #n=15
VA_no3 <- VA_data %>% 
  filter(age != 3)
view(VA_no3) # n=20

#alternatively just make new column to separate groups - so same data frame
VA_data <- VA_data %>% 
  mutate(age_grp = case_when(age == "3" ~ "3", age == "4" ~"older",
                             age == "5" ~ "older", age == "6" ~ "older", 
                             age == "7" ~ "older", age == "8" ~ "older", 
                             age == "9" ~ "older", age == "10" ~ "older", 
                             age == "11" ~ "older", age == "12" ~ "older", age == "13" ~ "older")) %>% 
  relocate(age_grp, .after = maturity)
view(VA_data)

# check Shapiro test by maturity groups, not individual ages
shapiro.test(subset(VA_3, status == "EB")$mass_spr) # p-value: 0.1912
shapiro.test(subset(VA_3, status == "ER")$mass_spr) # p-value: 0.284
shapiro.test(subset(VA_3, status == "OK")$mass_spr) # p-value: 0.7721
# all p-values are >0.05 so all variances are normally distributed, but small sample sizes
# ER n=6, ER n=6, OK n=3

##################################################################
# Test assumption of homogeneity of variance
# Null hypothesis of Bartlett: variance is equal

plot(mass_spr~status, data = VA_3) #see one significant outlier in the EB group - keep because of small sample sizes?
bartlett.test(mass_spr~status, data = VA_3)
# p-value: 0.1686 - accept the null hypothesis
# the variance is homogenous so we can proceed to the ANOVA

### ANOVA using aov (one-way test)
aov_3 <- aov(mass_spr~status, data = VA_3)
summary(aov_3)
# F = 662.1, num df = 2, p-value = 0.0128 *

# Interpretation of results
# since p-value is < 0.05 (just barely), we reject the null hypothesis
# we conclude at least one status group's mass is different than the others - do post-hoc test

# Tukey HSD test - post-hoc test
library(multcomp)
post_aov_3 <- glht(aov_3, linfct = mcp(status = "Tukey"))
summary(post_aov_3)
## ER - EB are not significantly different but OK - EB and OK - ER are slightly significantly different (*)
par(mar = c(3, 8, 3, 3))
plot(post_aov_3)

### plot the ANOVA and post-hoc tests
library(ggstatsplot)
### significance of mass differences between status groups (EB, ER, OK) for 3 yo 
ggbetweenstats(
  data = VA_3,
  x = status,
  y = mass_spr,
  type = "parametric", # ANOVA 
  var.equal = TRUE, 
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)
#########################################################################################
#### this graph is automatically displaying results of a Fisher test - not ANOVA. 
#### Should we have been doing a Fisher and not an ANOVA? But Fisher seems to be for categorical data not nominal...
#########################################################################################

########### Do same thing as above but with "older" group
##### Here is where we expect to see the significant differences by status grouping
shapiro.test(subset(VA_no3, status == "EB")$mass_spr) # insufficient sample size for comparison (n=2)
shapiro.test(subset(VA_no3, status == "ER")$mass_spr) # p-value: 0.7311
shapiro.test(subset(VA_no3, status == "OK")$mass_spr) # p-value: 0.927
# both p-values are >0.05 so all variances are normally distributed, but small sample sizes. EB group is only n=2
# ER n=2, ER n=15, OK n=3

##################################################################
# Test assumption of homogeneity of variance
# Null hypothesis of Bartlett: variance is equal

plot(mass_spr~status, data = VA_no3) #see no significant outliers 
bartlett.test(mass_spr~status, data = VA_no3)
# p-value: 0.6257 - accept the null hypothesis - proceed to ANOVA

aov_no3 <- aov(mass_spr~status, data = VA_no3)
summary(aov_no3)
# F = 6.85, num df = 2, p-value = 0.00658 **

# Interpretation of results
# since p-value is < 0.05 (just barely), we reject the null hypothesis
# we conclude at least one status group's mass is different than the others - do post-hoc test

# Tukey HSD test - post-hoc test
post_aov_no3 <- glht(aov_no3, linfct = mcp(status = "Tukey"))
summary(post_aov_no3)
## ER - EB are  significantly different (**), OK - EB also significantly different (**)
## but OK - ER are not significantly different

### plot the ANOVA and post-hoc tests
ggbetweenstats(
  data = VA_no3,
  x = status,
  y = mass_spr,
  type = "parametric", # should be ANOVA but it's not??? - automatically selected
  var.equal = TRUE, 
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)


##### can we do a Two-Way ANOVA to see differences in mass by status on one graph?

### ANOVA by aov, adding status for two-way ANOVA: modeling mass as a function of maturity and status
matst_aov <- aov(mass_spr ~ maturity + status, data = VA_data)
summary(matst_aov)
# seems to have improved the model. The residual sum of squares has gone down and both maturity and status are significant (p < 0.05)

# test whether the two variables have an interaction effect in ANOVA
interaction <- aov(mass_spr ~ maturity*status, data = VA_data)
summary(interaction)

####################### Whenever I get around to it... adding GVZoo and VA ultrasounds... ###########
## need to subset VA dataset with frog info and grades 01 and 02 and combine with GVZoo separate dataset

##################################################################
####### Case 4(?): amplexus behaviour(s) ~ ultrasound #################
# 1. does time to first contact (to_first) differ by grade 01? 
# 2. does total hours in amplexus (hrs_apx) differ by ultrasound grade (01-04)?
# 3. does number of amplectants (tot_num_amplex) differ by ultrasound grade (01-04)?

# 1: Continuous dependent variable: time to first (to_first)
# Categorical independent variables: ultrasound grade
#Unbalanced design, unequal sample sizes
##################################################################
shapiro.test(subset(VA_data, grade_01 == "0")$to_first) # p-value: 0.1069
shapiro.test(subset(VA_data, grade_01 == "1")$to_first) # p-value: 0.0414
shapiro.test(subset(VA_data, grade_01 == "2")$to_first) # p-value: 0.01232
shapiro.test(subset(VA_data, grade_01 == "3")$to_first) # p-value: 0.2805
# p-values for grades 0 and 3 are >0.05 (normal) but 1 and 2 are not normal 
# small sample sizes all around but particularly 0 and 3. 0: n=4, 1: n=12, 2: n=12, 3: n=7

##################################################################
# Test assumption of homogeneity of variance
# Null hypothesis of Bartlett: variance is equal

plot(to_first~grade_01, data = VA_data) #visual check - grade 2 seems to have a couple significant outliers
bartlett.test(to_first~grade_01, data = VA_data)
# p-value: 0.5403 - accept the null hypothesis
# the variance is homogenous but data was mostly not normal so let's use a Kruskal-Wallis

### Kruskal-Wallis test
kruskal.test(to_first~grade_01, data = VA_data)
# chi-squared = 7.122, p-value = 0.06811
# p > 0.05 (barely) so we ACCEPT the null hypothesis of all times to first contact being equal

### on the previous plot it did appear there were differences between times
### perhaps the large variation in times is what makes them not significant... what if the two outliers in grade 2 were removed? 

### can we really assess this question (to_first) if we are including the PHV treatment? Big bias no? 
# use dataset first_noPHV
plot(to_first~grade_01, data = first_noPHV) #gr 2's times jump up and suddenly include way more spread
kruskal.test(to_first~grade_01, data = first_noPHV)
# chi-squared = 5.836, p-value = 0.1199
##################### we see the significance go way down here (p > 0.05) and visually boxes look much closer
##################### confirms previous test: no differences between times by ultrasound grade 01

## One-way ANOVA can tell me if treatment impacts time spent in amplexus
## and Two-way ANOVA can tell me if treatment and ultrasound grade impact time spent in amplexus
## (time spent in amplexus probably only relates to grade 01...)

##### is there a way to use Two-way ANOVA to see if ultrasound grade and day (of ultrasound grade - i.e. day 01, or ....)
## impacts time spent in amplexus or number of amplectants ? 
## how do I designate "day" without moving ultrasound grades into long form ... how to do long form without replicating other data?

# 2: Continuous dependent variable:total hours in amplexus (hrs_apx)
# Categorical independent variables: ultrasound grade (and treatment?)
#Unbalanced design, unequal sample sizes
##################################################################
shapiro.test(subset(VA_data, grade_01 == "0")$hrs_apx) # p-value: 0.004132 = significant (n=4) = not-normal
shapiro.test(subset(VA_data, grade_01 == "1")$hrs_apx) # p-value: 0.001602 = significant (n=12)
shapiro.test(subset(VA_data, grade_01 == "2")$hrs_apx) # p-value: 0.007329 = significant (n=12)
shapiro.test(subset(VA_data, grade_01 == "3")$hrs_apx) # p-value: 0.4475 = not-significant (n=7) = normal

shapiro.test(subset(VA_data, grade_02 == "0")$hrs_apx) # p-value: 0.0002105 (n=5)
shapiro.test(subset(VA_data, grade_02 == "1")$hrs_apx) # p-value: cannot calculate (n=3)
shapiro.test(subset(VA_data, grade_02 == "2")$hrs_apx) # p-value: 0.002582
shapiro.test(subset(VA_data, grade_02 == "3")$hrs_apx) # p-value: 0.003731 (n=8)
shapiro.test(subset(VA_data, grade_02 == "4")$hrs_apx) # p-value: cannot calculate (n=2)

qqnorm(VA_data$hrs_apx) #does not look normal - not straight at all
hist(VA_data$hrs_apx) #large right tail
# can I do a transformation of the data to make it more normal? log doesn't seem to work

plot(hrs_apx~grade_01, data = VA_data)
plot(hrs_apx~grade_02, data = VA_data)
plot(hrs_apx~grade_03, data = VA_data)
plot(hrs_apx~grade_04, data = VA_data)

# check homogeneity of variances with Bartlett test
bartlett.test(hrs_apx~grade_01, data = VA_data) #p = 0.0001404 (reject null)
bartlett.test(hrs_apx~grade_02, data = VA_data) # p cannot be calculated (too small sample sizes)
bartlett.test(hrs_apx~grade_03, data = VA_data) #p = 0.8502 (accept null = homogenous)
bartlett.test(hrs_apx~grade_04, data = VA_data) #p = 0.05544 (barely accept null = homogenous)
# use Kruskal-Wallis since most data is not normal and some don't have homogenous variance

kruskal.test(hrs_apx~grade_01, data = VA_data) ## p-value = 0.04075 #### reject null = some differ
kruskal.test(hrs_apx~grade_02, data = VA_data) ## p-value = 0.5236 ### accept null. none differ
kruskal.test(hrs_apx~grade_03, data = VA_data) ## p-value = 0.6787 ### accept null. none differ
kruskal.test(hrs_apx~grade_04, data = VA_data) ## p-value = 0.2937 ### accept null. none differ

### only need to do post-hoc test on grade_01 since it was the only one with significant K-W test
library(FSA)
dunnTest(hrs_apx~grade_01, data = VA_data, method = "holm")
#### the only means that differ significantly are 0 to 2 with a p.adj = 0.04755612
#### this is barely significant

## plot statistical results to show minor significant differences
ggbetweenstats(
  data = VA_data,
  x = grade_01,
  y = hrs_apx,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

############## Same thing for total number of amplectants ################
# 3: Continuous dependent variable: total number of amplectants (tot_num_amplex)
# Categorical independent variables: ultrasound grade
#Unbalanced design, unequal sample sizes
##################################################################
shapiro.test(subset(VA_data, grade_01 == "0")$tot_num_amplex) # p-value: 0.00717
shapiro.test(subset(VA_data, grade_01 == "1")$tot_num_amplex) # p-value: 0.03373
shapiro.test(subset(VA_data, grade_01 == "2")$tot_num_amplex) # p-value: 0.2073
shapiro.test(subset(VA_data, grade_01 == "3")$tot_num_amplex) # p-value: 0.1275
# p-values for grades 0 and 3 are >0.05 (normal) but 1 and 2 are not normal 

shapiro.test(subset(VA_data, grade_02 == "0")$tot_num_amplex) # p-value: 0.04007 = significant = not-normal
shapiro.test(subset(VA_data, grade_02 == "1")$tot_num_amplex) # small sample size
shapiro.test(subset(VA_data, grade_02 == "2")$tot_num_amplex) # p-value: 0.03789 = significant 
shapiro.test(subset(VA_data, grade_02 == "3")$tot_num_amplex) # p-value: 0.05132 = not-significant = normal
shapiro.test(subset(VA_data, grade_02 == "4")$tot_num_amplex) # sample size too small

shapiro.test(subset(VA_data, grade_03 == "0")$tot_num_amplex) # p-value: 0.02778
shapiro.test(subset(VA_data, grade_03 == "1")$tot_num_amplex) # p-value: 0.1736
shapiro.test(subset(VA_data, grade_03 == "2")$tot_num_amplex) # p-value: 0.07599
shapiro.test(subset(VA_data, grade_03 == "3")$tot_num_amplex) # p-value: 0.3864
shapiro.test(subset(VA_data, grade_03 == "4")$tot_num_amplex) # p-value: 0.4638
## grade_03 are all normal - but may be biased by very small sample sizes

qqnorm(VA_data$tot_num_amplex) #not very straight
hist(VA_data$tot_num_amplex) #large right tail
# can I do a transformation of the data to make it more normal? log doesn't seem to work
library(car)
qqPlot(VA_data$tot_num_amplex, id = FALSE) #a very wide confidence band but most points are within it

plot(tot_num_amplex~grade_01, data = VA_data)
plot(tot_num_amplex~grade_02, data = VA_data)
plot(tot_num_amplex~grade_03, data = VA_data)
plot(tot_num_amplex~grade_04, data = VA_data) #see one potentially significant outlier in grade 4 here

# check homogeneity of variances with Bartlett test
bartlett.test(tot_num_amplex~grade_01, data = VA_data) #p = 0.1965 (accept null = homogenous)
bartlett.test(tot_num_amplex~grade_02, data = VA_data) # p cannot be calculated (too small sample sizes)
bartlett.test(tot_num_amplex~grade_03, data = VA_data) #p = 0.6622 (accept null = homogenous)
bartlett.test(tot_num_amplex~grade_04, data = VA_data) #p = 0.6718 ( accept null = homogenous)
# we can use ANOVA because data seems to be normal and variance is homogenous (we don't know this for grade_02 for sure)

### ANOVA using aov (one-way test)
aov_num01 <- aov(tot_num_amplex~grade_01, data = VA_data)
summary(aov_num01)
# F = 1.43 , num df = 3, p-value = 0.257

# Interpretation of results
# since p-value is > 0.05, we accept the null hypothesis
#### There are no significant differences between means
aov_num02 <- aov(tot_num_amplex~grade_02, data = VA_data)
summary(aov_num02) # F = 0.543 , num df = 4, p-value = 0.706 #### accept null
aov_num03 <- aov(tot_num_amplex~grade_03, data = VA_data)
summary(aov_num03) # F = 0.645 , num df = 4, p-value = 0.645 #### accept null
aov_num04 <- aov(tot_num_amplex~grade_04, data = VA_data)
summary(aov_num04) # F = 0.6 , num df = 4, p-value = 0.558 #### accept null

##### Conclusion: there are no significant differences between means in any grade_ 
##### No post-hoc test necessary

##### can we do a Two-Way ANOVA to see differences in number of amplectants by grade and treatment ?
### ANOVA by aov, adding treatment for two-way ANOVA: modeling tot_num as a function of grade and treatment
numt_aov <- aov(tot_num_amplex ~ grade_01 + treatment, data = VA_data)
summary(numt_aov)
# seems to have improved the model. The residual sum of squares has gone down
# only treatment is significant (p < 0.05)

# test whether the two variables have an interaction effect in ANOVA
interaction <- aov(tot_num_amplex ~ grade_01*treatment, data = VA_data)
summary(interaction)
# interaction is not significant so - no? 
# try two-way ANOVA with other grade_
numt_aov <- aov(tot_num_amplex ~ grade_02 + treatment, data = VA_data)
summary(numt_aov)
numt_aov <- aov(tot_num_amplex ~ grade_03 + treatment, data = VA_data)
summary(numt_aov)
numt_aov <- aov(tot_num_amplex ~ grade_04 + treatment, data = VA_data)
summary(numt_aov)
### all seem to improve the models, but treatment is still only one significant
### which makes sense because the one-way ANOVA said grade was not significant

################################################################################
### Let's look at amplexus behaviours ~ treatment for significance #############
######### Case 1: to_first ~ treatment 
hist(VA_data$to_first) #looks pretty bimodal
qqnorm(VA_data$to_first) #not straight
shapiro.test(subset(VA_data, treatment == "C")$to_first) #p-value: 0.02249 (not normal)
shapiro.test(subset(VA_data, treatment == "H")$to_first) #p-value: 0.2719 (not normal)
shapiro.test(subset(VA_data, treatment == "V")$to_first) #p-value: 0.1582 (not normal)
shapiro.test(subset(VA_data, treatment == "HV")$to_first) #p-value: 0.09262 (not normal)
shapiro.test(subset(VA_data, treatment == "PHV")$to_first) #cannot calculate (all values same)
#### could I do a transformation? couldn't get it to work previously
plot(to_first~treatment, data = VA_data) #see quite a few potentially significant outliers. keep because small sample sizes
bartlett.test(to_first~treatment, data = VA_data)
# p-value: 2.2e-16 = reject the null hypothesis so variance is not homogenous. Use Kruskal-Wallis test
kruskal.test(to_first~treatment, data = VA_data)
# chi-squared = 19.603, df = 4, p-value = 0.0005981
# p < 0.05 so we reject null of all means being equal. need to do post-hoc test (Dunn)
dunnTest(to_first~treatment, data = VA_data, method = "holm")
## see only  slightly significant values (p < 0.05) for H-PHV, HV-PHV, and PHV-V
## so really, PHV differs significantly from all treatments except Control

ggbetweenstats(
  data = VA_data,
  x = treatment,
  y = to_first,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

######### Case 2: hrs_apx ~ treatment 
hist(VA_data$hrs_apx) #not curved
qqnorm(VA_data$hrs_apx) #not straight
shapiro.test(subset(VA_data, treatment == "C")$hrs_apx) #p-value: 0.02372 (not normal)
shapiro.test(subset(VA_data, treatment == "H")$hrs_apx) #p-value: 0.6921 (not normal)
shapiro.test(subset(VA_data, treatment == "V")$hrs_apx) #p-value: 0.3695 (not normal)
shapiro.test(subset(VA_data, treatment == "HV")$hrs_apx) #p-value: 3.733e-05 (normal)
shapiro.test(subset(VA_data, treatment == "PHV")$hrs_apx) #p-value: 0.0003212 (normal)
#### could I do a transformation? couldn't get it to work previously
plot(hrs_apx~treatment, data = VA_data) #see  a couple potentially significant outliers. keep because small sample sizes
bartlett.test(hrs_apx~treatment, data = VA_data)
# p-value: 3.84e-08 = reject the null hypothesis so variance is not homogenous. Use Kruskal-Wallis test
kruskal.test(hrs_apx~treatment, data = VA_data)
# chi-squared = 13.877, df = 4, p-value = 0.007699
# p < 0.05 so we reject null of all means being equal. need to do post-hoc test (Dunn)
dunnTest(hrs_apx~treatment, data = VA_data, method = "holm")
## only see significant difference between H-PHV

ggbetweenstats(
  data = VA_data,
  x = treatment,
  y = hrs_apx,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

######### Case 3: tot_num_amplex ~ treatment 
hist(VA_data$tot_num_amplex) #not curved
qqnorm(VA_data$tot_num_amplex) #not straight
shapiro.test(subset(VA_data, treatment == "C")$tot_num_amplex) #p-value: 0.4195 (not normal)
shapiro.test(subset(VA_data, treatment == "H")$tot_num_amplex) #p-value: 0.1494 (not normal)
shapiro.test(subset(VA_data, treatment == "V")$tot_num_amplex) #p-value: 0.6656 (not normal)
shapiro.test(subset(VA_data, treatment == "HV")$tot_num_amplex) #p-value: 01638 (normal)
shapiro.test(subset(VA_data, treatment == "PHV")$tot_num_amplex) #p-value: 0.1012 (not normal)
#### could I do a transformation? couldn't get it to work previously
plot(tot_num_amplex~treatment, data = VA_data) #no obvious outliers
bartlett.test(tot_num_amplex~treatment, data = VA_data)
# p-value: 7.088e-05 = reject the null hypothesis so variance is not homogenous. Use Kruskal-Wallis test
kruskal.test(tot_num_amplex~treatment, data = VA_data)
# chi-squared = 12.458, df = 4, p-value = 0.01425
# p < 0.05 so we reject null of all means being equal. need to do post-hoc test (Dunn)
dunnTest(tot_num_amplex~treatment, data = VA_data, method = "holm")
## only see significant difference between PHV-V

ggbetweenstats(
  data = VA_data,
  x = treatment,
  y = tot_num_amplex,
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE
)

