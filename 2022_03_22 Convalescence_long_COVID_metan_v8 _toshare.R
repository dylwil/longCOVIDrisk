
######################################################################################################################## 

# Script for meta-analyses and forest plotting of long COVID risk factors analyses from 10 Longitudinal Study samples and EHR data from OpenSAFELY

# Dylan Williams (dylan.williams@ucl.ac.uk) on behalf of the CONVALESCENCE consortium (iterations of meta-analyses conducted June - December 2021)

########################################################################################################################

# Packages:

my_packages = c("dplyr", "readr", "meta", "forcats")   # Packages used in following script
lapply(my_packages, require, character.only = TRUE)    # Load packages if present


#  Data loading and some reformatting steps prior to analysis: ------------

    # Edit the following file path to import data from where these are deposited:

File_path = "C:/Users/rmgpdmw/OneDrive - University College London/Projects/COVID/LHW programme in the NCS/ARQ2 - long COVID/Data"

res = read_delim(file.path(File_path, "/2022_03_22_model2_regressions_reformatted_toshare.txt"), "\t", escape_double = FALSE, trim_ws = TRUE)

res = rename(res, Study = cohort, se = coef_se, lci = lower_ci, uci = upper_ci, Exposure = exposure, Outcome = outcome) # colnames renamed to match previous script

# Adding log-odds to data:

res = res %>% mutate(logodds = log(coef))


# Changing BMI coefficients to be per 5kg

res = res %>% 
    mutate(logodds = ifelse(Grouping == "BMI", logodds*5, logodds)) %>%
    mutate(coef = ifelse(Grouping == "BMI", exp(logodds), coef)) %>%
    mutate(lci = ifelse(Grouping == "BMI", exp(logodds - (1.96*5*se)), lci)) %>%
    mutate(uci = ifelse(Grouping == "BMI", exp(logodds + (1.96*5*se)), uci))


# Removing education results for MCS (omitted because not all students are old enough to have completed degree):

res= subset(res, !(Study=="MCS" & Grouping=="Education"))

# Removing IMD result for Scotland because this is not directly comparable to English IMD values

res= subset(res, !(Study=="GS" & Exposure=="SIMD"))


# Renaming Next Steps as NS:
res = res %>% 
  mutate(Study = ifelse(Study == "Next Steps", "NS", Study))


# Fixing study order to match legend order:
table(res$Study)
res$study_order = NA
res = mutate(res,
                        study_order = case_when(
                          Study == "ALSPAC G0" ~ 1,
                          Study == "ALSPAC G1" ~ 2,
                          Study == "BCS70" ~ 3,
                          Study == "GS" ~ 4,
                          Study == "MCS" ~ 5,
                          Study == "NS" ~ 6,
                          Study == "NCDS" ~ 7,
                          Study == "OpenSAFELY" ~ 8,
                          Study == "TwinsUK" ~ 9,
                          Study == "USoc" ~ 10))

# Subsets of res for 12+ and 4+ week outcomes:

res_12w = filter(res, Outcome=="12+ weeks")

res_4w = filter(res, Outcome!="12+ weeks")


# OpenSAFELY res alone
res_open = filter(res, Study=="OpenSAFELY") 

    # Reordering results:


res_open = res_open %>%
  mutate(Exposure = fct_reorder(Exposure, Order)) 



######## Meta-analyses ########

# 1) Plotting 4+ week outcome results for cohorts -----------------------


  # removing results that won't be meta-analysed:
res_4w_cohorts = res_4w %>% 
  filter(Grouping!="Age") %>% 
  filter(Exposure!="pre-pandemic mh score")  %>% 
  filter(Exposure!="prepandemic mh score NO SOMATIC SUBSCALE")  %>% 
  filter(Study!="OpenSAFELY") %>%
  filter(Grouping!="Age") %>%
  filter(Grouping!="Social_class")

# Labelling grouping variables:

table(res_4w_cohorts$Grouping)

res_4w_cohorts$Grouping_f = factor(res_4w_cohorts$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking', 'Gen_health', 
                                                                     'M_Health', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 
                                                                     'Asthma'))



  # Splitting results to plot demographic and health results separately:

res_4w_cohorts_demog = filter(res_4w_cohorts, as.numeric(Grouping_f)<6)

res_4w_cohorts_demog$Group_labs = factor(res_4w_cohorts_demog$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking'))

levels(res_4w_cohorts_demog$Group_labs) = c("Female sex (ref. male)", "Other ethnicity (ref. white)", "No higher education (ref. degree attained)", "Index of multiple deprivation (per one decile)", "Current smoking")


res_4w_cohorts_health = filter(res_4w_cohorts, as.numeric(Grouping_f)>=6)

res_4w_cohorts_health$Group_labs = factor(res_4w_cohorts_health$Grouping, levels = c('Gen_health', 'M_Health', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 'Asthma'))

levels(res_4w_cohorts_health$Group_labs) = c("Poor overall health", "Psychological distress", "BMI (per 5 kg/m^2)", "Overweight or obese", "Diabetes", "Hypertension", "High cholesterol", "Asthma")

  # M-A of 4+ week demography res from cohorts:


setwd("C:/Users/rmgpdmw/OneDrive - University College London/Projects/COVID/LHW programme in the NCS/ARQ2 - long COVID/Output/Plots/Test of final script to share")

attach(res_4w_cohorts_demog)

meta_4w_cohorts_demog <-metagen(log(coef), 
                                lower = log(lci),
                                upper = log(uci),
                  studlab=paste(Study),
                  subgroup = Group_labs, 
                  sm = "OR",
                  title = n)

summary(meta_4w_cohorts_demog)

tiff("v8_4plus_model2_demog.tiff", width = 300, height = 380, units = "mm", res = 100)
forest(meta_4w_cohorts_demog, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_4w_cohorts_demog$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off() 

detach(res_4w_cohorts_demog)

    # creating table of fixed effects MA output for subgroup res:

str(meta_4w_cohorts_demog, list.len = 153) # see full number of items in list
attach(meta_4w_cohorts_demog)
fixed_4w_demog_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_4w_demog_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_4w_cohorts_demog)


# M-A of 4+ health res from cohorts:

attach(res_4w_cohorts_health)

meta_4w_cohorts_health <-metagen(log(coef),
                                lower = log(lci),
                                upper = log(uci),
                                level.ci = 0.95,
                                studlab=paste(Study),
                                subgroup = Group_labs, 
                                sm = "OR",
                                title = n)


    # Plot
tiff("v8_4plus_model2_health.tiff", width = 300, height = 600, units = "mm", res = 150)
forest(meta_4w_cohorts_health, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       sortvar = w.fixed,
       bylab = meta_4w_cohorts_health$Group_labs,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       test.overall.fixed = FALSE,
       test.overall.random = FALSE,
       ff.lr = "bold", 
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off() 

detach(res_4w_cohorts_health)


attach(meta_4w_cohorts_health)
fixed_4w_health_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_4w_health_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_4w_cohorts_health)

  # Combining demog and health res together for export:
fixed_4w_res = rbind(fixed_4w_demog_res, fixed_4w_health_res)
write.csv(fixed_4w_res, "2021_12_06_4plus_model2_fixedmeta_res_v8.csv", row.names=F, quote=F)




# 2) Plotting 12+ week results for cohorts -------------------------------

res_12w_cohorts = res_12w %>% 
  filter(Grouping!="Age") %>% 
  filter(Exposure!="pre-pandemic mh score")  %>% 
  filter(Exposure!="prepandemic mh score NO SOMATIC SUBSCALE")  %>% 
  filter(Study!="OpenSAFELY") %>%
  filter(Grouping!="Age") %>%
  filter(Grouping!="Social_class")


# Labelling grouping variables:

res_12w_cohorts$Grouping_f = factor(res_12w_cohorts$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking', 'Gen_health', 
                                                                         'M_Health', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 
                                                                         'Asthma'))
# Splitting results to plot demographic and health results separately:

res_12w_cohorts_demog = filter(res_12w_cohorts, as.numeric(Grouping_f)<6)


res_12w_cohorts_demog$Group_labs = factor(res_12w_cohorts_demog$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking'))

levels(res_12w_cohorts_demog$Group_labs) = c("Female sex (ref. male)", "Other ethnicity (ref. white)", "No higher education (ref. degree attained)", "Index of multiple deprivation (per one decile)", "Current smoker")


res_12w_cohorts_health = filter(res_12w_cohorts, as.numeric(Grouping_f)>=6)
table(res_12w_cohorts_health$Grouping)

res_12w_cohorts_health$Group_labs = factor(res_12w_cohorts_health$Grouping, levels = c('Gen_health', 'M_Health', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 'Asthma'))

levels(res_12w_cohorts_health$Group_labs) = c("Poor overall health", "Psychological distress", "BMI (per 5 /m^2)", "Overweight or obese", "Diabetes", "Hypertension", "High cholesterol", "Asthma")


########## M-A of 12+ week demography res from cohorts:

attach(res_12w_cohorts_demog)

meta_12w_cohorts_demog <-metagen(log(coef), 
                                lower = log(lci),
                                upper = log(uci),
                                studlab=paste(Study),
                                subgroup = Group_labs, 
                                sm = "OR",
                                title = n)

summary(meta_12w_cohorts_demog)

tiff("v8_12plus_model2_demog.tiff", width = 300, height = 360, units = "mm", res = 100)
forest(meta_12w_cohorts_demog, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_12w_cohorts_demog$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 12+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off()

detach(res_12w_cohorts_demog)


attach(meta_12w_cohorts_demog)
fixed_12w_demog_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_12w_demog_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_12w_cohorts_demog)




# M-A of 12+ health res from cohorts:

attach(res_12w_cohorts_health)

meta_12w_cohorts_health <- metagen(log(coef),
                                   lower = log(lci),
                                   upper = log(uci),
                                   level.ci = 0.95,
                                   studlab=paste(Study),
                                   subgroup = Group_labs, 
                                   sm = "OR",
                                   title = n)



summary(meta_12w_cohorts_health)

tiff("v8_12plus_model2_health.tiff", width = 300, height = 600, units = "mm", res = 150)
forest(meta_12w_cohorts_health, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_12w_cohorts_health$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 12+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off()

detach(res_12w_cohorts_health)

attach(meta_12w_cohorts_health)
fixed_12w_health_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_12w_health_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_12w_cohorts_health)


# Combining demog and health res together for export:
fixed_12w_res = rbind(fixed_12w_demog_res, fixed_12w_health_res)
write.csv(fixed_12w_res, "2021_12_06_12plus_model2_fixedmeta_res_v8.csv", row.names=F, quote=F)


# 3) Forest plots for age results only ------------------------------------

res_age = res %>% 
  filter(Grouping=="Age")

res_age_12w = filter(res_age, Outcome=="12+ weeks")

res_age_4w = filter(res_age, Outcome!="12+ weeks")


res_age_4w = filter(res_age_4w, Study!="MCS")
res_age_4w = filter(res_age_4w, Study!="NS")



      # Plotting 4 week results inc OpenSAFELY:

# Reassigning study order:

res_age_4w$study_order = NA
res_age_4w$study_order[res_age_4w$Study == "USoc"] = 1
res_age_4w$study_order[res_age_4w$Study == "TwinsUK"] = 2
res_age_4w$study_order[res_age_4w$Study == "GS"] = 3
res_age_4w$study_order[res_age_4w$Study == "ALSPAC G0"] = 4
res_age_4w$study_order[res_age_4w$Study == "OpenSAFELY"] = 5

table(res_age_4w$Exposure)
res_age_4w$Exposure[res_age_4w$Exposure == "age(70+vs18-44)"] = "70+"
res_age_4w$Exposure[res_age_4w$Exposure == "age(45-69+vs18-44)"] = "45-69"
res_age_4w$Exposure[res_age_4w$Exposure == "age(linear)"] = "Per year increase"
res_age_4w$Exposure[res_age_4w$Exposure == "age (cont)"] = "Per year increase"
res_age_4w$Exposure[res_age_4w$Exposure == "age (60-76 vs 45-54)"] = "60-76"
res_age_4w$Exposure[res_age_4w$Exposure == "age (55-59 vs 45-54)"] = "55-59"
res_age_4w$Exposure[res_age_4w$Exposure == "age(45+vs25-34)"] = "45+"
res_age_4w$Exposure[res_age_4w$Exposure == "age(35-44+vs25-34)"] = "35-44"
res_age_4w$Exposure[res_age_4w$Exposure == "25-34 (all vs 18-24)"] = "25-34"
res_age_4w$Exposure[res_age_4w$Exposure == "age (70+vs18-44)"] = "70+"
res_age_4w$Exposure[res_age_4w$Exposure == "age (45-69+vs18-44)"] = "45-69"


res_age_4w = res_age_4w %>%
  arrange(study_order)

res_age_4w$Group_labs = factor(res_age_4w$Study, levels = c('USoc', 'TwinsUK', 'GS', 'ALSPAC G0', 'OpenSAFELY'))

levels(res_age_4w$Group_labs) = c("USoc (ref. 18-44y)", "TwinsUK (ref. 18-44y)", "GS (ref. 18-44y)", "ALSPAC G0 (ref. 45-54y)", "OpenSAFELY (ref. 18-24y)")


attach(res_age_4w)

meta_4w_age <- metagen(log(coef),
                                   lower = log(lci),
                                   upper = log(uci),
                                   level.ci = 0.95,
                                   studlab=paste(Exposure),
                                   subgroup = Group_labs, 
                                   sm = "OR")

tiff("v8_4plus_ageonly.tiff", width = 280, height = 250, units = "mm", res = 100)
forest(meta_4w_age, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab"),
       leftlabs=c("Study"),
       bylab = res_age_4w$Group_labs,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       weight.study = "same",
       squaresize = 0,
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold", 
       just = "right",
       subgroup = FALSE,
       print.byvar =FALSE,
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       hetstat = FALSE)
dev.off()

detach(res_age_4w)

 

      # Plotting 12 week results without OpenSAFELY:

# Reassigning study order:


res_age_12w = filter(res_age_12w, Study!="MCS")
res_age_12w = filter(res_age_12w, Study!="NS")

res_age_12w$study_order = NA
res_age_12w$study_order[res_age_12w$Study == "USoc"] = 1
res_age_12w$study_order[res_age_12w$Study == "TwinsUK"] = 2
res_age_12w$study_order[res_age_12w$Study == "GS"] = 3
res_age_12w$study_order[res_age_12w$Study == "ALSPAC G0"] = 4

table(res_age_12w$Exposure)
res_age_12w$Exposure[res_age_12w$Exposure == "age(70+vs18-44)"] = "70+"
res_age_12w$Exposure[res_age_12w$Exposure == "age(45-69+vs18-44)"] = "45-69"
res_age_12w$Exposure[res_age_12w$Exposure == "age(linear)"] = "Per year increase"
res_age_12w$Exposure[res_age_12w$Exposure == "age (cont)"] = "Per year increase"
res_age_12w$Exposure[res_age_12w$Exposure == "age (60-76 vs 45-54)"] = "60-76"
res_age_12w$Exposure[res_age_12w$Exposure == "age (55-59 vs 45-54)"] = "55-59"
res_age_12w$Exposure[res_age_12w$Exposure == "age(45+vs25-34)"] = "45+"
res_age_12w$Exposure[res_age_12w$Exposure == "age(35-44+vs25-34)"] = "35-44"
res_age_12w$Exposure[res_age_12w$Exposure == "25-34 (all vs 18-24)"] = "25-34"
res_age_12w$Exposure[res_age_12w$Exposure == "age (70+vs18-44)"] = "70+"
res_age_12w$Exposure[res_age_12w$Exposure == "age (45-69+vs18-44)"] = "45-69"


res_age_12w = res_age_12w %>%
  arrange(study_order)

res_age_12w$Group_labs = factor(res_age_12w$Study, levels = c('USoc', 'TwinsUK', 'GS', 'ALSPAC G0'))

levels(res_age_12w$Group_labs) = c("USoc (ref. 18-44y)", "TwinsUK (ref. 18-44y)", "GS (ref. 18-44y)", "ALSPAC G0 (ref. 45-54y)")



attach(res_age_12w)

meta_12w_age <- metagen(log(coef),
                       lower = log(lci),
                       upper = log(uci),
                       level.ci = 0.95,
                       studlab=paste(Exposure),
                       subgroup = Group_labs, 
                       sm = "OR")

tiff("v8_12plus_ageonly.tiff", width = 280, height = 250, units = "mm", res = 100)
forest(meta_12w_age, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab"),
       leftlabs=c("Study"),
       bylab = res_age_12w$Group_labs,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 12+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       weight.study = "same",
       squaresize = 0,
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold", 
       just = "right",
       subgroup = FALSE,
       print.byvar =FALSE,
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       hetstat = FALSE)
dev.off()

detach(res_age_12w)




# 4) Meta-analysis of factors associated with C19 risk --------
# (to inform IPW strategy for index event bias)

res = read_delim(file.path(File_path,"/2021_06_16_C19status_regressions_reformatted.txt"), "\t", escape_double = FALSE, trim_ws = TRUE)

res = rename(res, Study = cohort, se = coef_se, lci = lower_ci, uci = upper_ci, Exposure = exposure, Outcome = outcome) # colnames renamed to match previous script

# Adding log-odds to data:

res = res %>% mutate(logodds = log(coef))


# Removing education results for MCS (probably not meaningful because not all students are old enough to have completed degree):

res= subset(res, !(Study=="MCS" & Grouping=="Education"))


# Renaming Next Steps as NS:
res = res %>% 
  mutate(Study = ifelse(Study == "Next Steps", "NS", Study))




# Fixing study order to match legend order:

res$study_order = NA
res <- mutate(res,
              study_order = case_when(
                Study == "ALSPAC G0" ~ 1,
                Study == "ALSPAC G1" ~ 2,
                Study == "BCS70" ~ 3,
                Study == "GS" ~ 4,
                Study == "MCS" ~ 5,
                Study == "NS" ~ 6,
                Study == "NCDS" ~ 7,
                Study == "TwinsUK" ~ 8,
                Study == "USoc" ~ 9))

  # Removing age before M-As:
res_c19 = res %>% 
  filter(Grouping!="Age") 


table(res_c19$Grouping)
res_c19$Grouping_f = factor(res_c19$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Social_class', 'Smoking',  'Gen_health', 
                                                                       'M_Health_cont', 'M_Health_binary', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 
                                                                       'Asthma'))

# Splitting results to plot demographic and health results separately:

res_c19_demog = filter(res_c19, as.numeric(Grouping_f)<7)

res_c19_demog$Group_labs = factor(res_c19_demog$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Social_class', 'Smoking'))

levels(res_c19_demog$Group_labs) = c("Female sex (ref. male)", "Other ethnicity (ref. white)", "No higher education (ref. degree attained)", "Index of multiple deprivation (per one decile)", "Low social class", "Current smoking")



res_c19_health = filter(res_c19, as.numeric(Grouping_f)>=7)

res_c19_health$Group_labs = factor(res_c19_health$Grouping, levels = c('Gen_health', 'M_Health_cont', 'M_Health_binary', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 'Asthma'))

levels(res_c19_health$Group_labs) = c("Poor overall health", "Mental health score", "Psychological distress", "BMI (per 5 kg/m^2)", "Overweight or obese", "Diabetes", "Hypertension", "High cholesterol", "Asthma")


  # M-A for demographic characteristics:

attach(res_c19_demog)

meta_c19_demog <-metagen(log(coef), 
                                lower = log(lci),
                                upper = log(uci),
                                studlab=paste(Study),
                                subgroup = Group_labs, 
                                sm = "OR",
                                title = n)
summary(meta_c19_demog)

tiff("v8_C19_model3_demog.tiff", width = 300, height = 450, units = "mm", res = 100)
forest(meta_c19_demog, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_c19_demog$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for COVID-19",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher COVID-19 risk",
       label.left="Lower COVID-19 risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off() 

detach(res_c19_demog)


attach(meta_c19_demog)
fixed_c19_demog_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_c19_demog_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_c19_demog)




# M-A for health characteristics:

attach(res_c19_health)

meta_c19_health <-metagen(log(coef), 
                         lower = log(lci),
                         upper = log(uci),
                         studlab=paste(Study),
                         subgroup = Group_labs, 
                         sm = "OR",
                         title = n)

summary(meta_c19_health)

tiff("v8_C19_model3_health.tiff", width = 300, height = 600, units = "mm", res = 100)
forest(meta_c19_health, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       sortvar = w.fixed,
       bylab = meta_c19_health$Group_labs,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for COVID-19",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher COVID-19 risk",
       label.left="Lower COVID-19 risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off() 

detach(res_c19_health)

attach(meta_c19_health)
fixed_c19_health_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_c19_health_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_c19_health)


# Combining demog and health res together for export:
fixed_c19_res = rbind(fixed_c19_demog_res, fixed_c19_health_res)
write.csv(fixed_c19_res, "2021_12_06_c19_model3_fixedmeta_res_v8.csv", row.names=F, quote=F)




# 6) combining results from OpenSAFELY alongside findings from cohort meta-analyses  --------


res_all <- read_delim(file.path(File_path,"2021_12_06_model2_combining_all_OS_MA_res.txt"),
                      "\t", escape_double = FALSE, trim_ws = TRUE)

res_all = res_all %>% 
  mutate(Group = ifelse(Group == "Overweight or obesity", "Overweight and obesity", Group))

attach(res_all)

meta_all <- metagen(log(or),
                    lower = log(lci),
                    upper = log(uci),
                    level.ci = 0.95,
                    studlab=paste(Measure),
                    subgroup = Group, 
                    sm = "OR", 
                    title = k)
summary(meta_all)

tiff("v8_model2_overall_forest.tiff", width = 300, height = 380, units = "mm", res = 100)
forest(meta_all, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "# Studies"),
       bylab = meta_all$Group,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "12cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       weight.study = "same",
       squaresize = 0,
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold", 
       just = "right",
       subgroup = FALSE,
       hetstat = FALSE,
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       at=c(0.5,0.75,1,1.5,2))
dev.off()


detach(res_all)



# 7) Meta-analyses of long COVID results generated with IPW for COVID status, --------


# Loading in new res with IPW:

res <- read_delim(file.path(File_path, "2021_06_20_IPW_model2_regressions_reformatted.txt"),
                  "\t", escape_double = FALSE, trim_ws = TRUE)

res = rename(res, Study = cohort, se = coef_se, lci = lower_ci, uci = upper_ci, Exposure = exposure, Outcome = outcome) # colnames renamed to match previous script

# Adding log-odds to data:

res = res %>% mutate(logodds = log(coef))

# Changing BMI coefficients to be per 5kg

res = res %>% 
  mutate(logodds = ifelse(Grouping == "BMI", logodds*5, logodds)) %>%
  mutate(coef = ifelse(Grouping == "BMI", exp(logodds), coef)) %>%
  mutate(lci = ifelse(Grouping == "BMI", exp(logodds - (1.96*5*se)), lci)) %>%
  mutate(uci = ifelse(Grouping == "BMI", exp(logodds + (1.96*5*se)), uci))




# Removing education results for MCS (probably not meaningful because not all students are old enough to have completed degree):

res= subset(res, !(Study=="MCS" & Grouping=="Education"))


# Removing Gen Scot from IMD (scottish IMD not comparable)

res= subset(res, !(Study=="GS" & Exposure=="SIMD"))




# Renaming Next Steps as NS:
res = res %>% 
  mutate(Study = ifelse(Study == "Next Steps", "NS", Study))


# Subsets of res for 12+ and 4+ week outcomes:

res_12w = filter(res, Outcome=="12+ weeks")

res_4w = filter(res, Outcome!="12+ weeks")



# Plotting 4+ week results for cohorts:

table(res_4w$Exposure)

# removing results that won't be meta-analysed:
res_4w_cohorts = res_4w %>% 
  filter(Grouping!="Age") %>% 
  filter(Exposure!="pre-pandemic mh score")  %>% 
  filter(Exposure!="prepandemic mh score NO SOMATIC SUBSCALE")  %>% 
  filter(Study!="OpenSAFELY") %>%
  filter(Grouping!="Age") %>%
  filter(Grouping!="Social_class")

table(res_4w_cohorts$Grouping)
table(res_4w_cohorts$Exposure)

# Labelling grouping variables:

table(res_4w_cohorts$Grouping)


res_4w_cohorts$Grouping_f = factor(res_4w_cohorts$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking', 'Gen_health', 
                                                                       'M_Health_binary', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 
                                                                       'Asthma'))



# Splitting results to plot demographic and health results separately:

res_4w_cohorts_demog = filter(res_4w_cohorts, as.numeric(Grouping_f)<6)

res_4w_cohorts_demog$Group_labs = factor(res_4w_cohorts_demog$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking'))

levels(res_4w_cohorts_demog$Group_labs) = c("Female sex (ref. male)", "Other ethnicity (ref. white)", "No higher education (ref. degree attained)", "Index of multiple deprivation (per one decile)", "Current smoking")


res_4w_cohorts_health = filter(res_4w_cohorts, as.numeric(Grouping_f)>=6)
table(res_4w_cohorts_health$Grouping)

res_4w_cohorts_health$Group_labs = factor(res_4w_cohorts_health$Grouping, levels = c('Gen_health', 'M_Health_binary', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 'Asthma'))

levels(res_4w_cohorts_health$Group_labs) = c("Poor overall health", "Psychological distress", "BMI (per 5 kg/m^2)", "Overweight or obese", "Diabetes", "Hypertension", "High cholesterol", "Asthma")




# M-A of 4+ week demography res from cohorts:

attach(res_4w_cohorts_demog)

meta_4w_cohorts_demog <-metagen(log(coef), 
                                lower = log(lci),
                                upper = log(uci),
                                studlab=paste(Study),
                                subgroup = Group_labs, 
                                sm = "OR",
                                title = n)

summary(meta_4w_cohorts_demog)

tiff("v8_IPW_4plus_model2_demog.tiff", width = 300, height = 380, units = "mm", res = 100)
forest(meta_4w_cohorts_demog, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_4w_cohorts_demog$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off() 

detach(res_4w_cohorts_demog)

# creating table of fixed effects MA output for subgroup res:

str(meta_4w_cohorts_demog, list.len = 153) # see full number of items in list
attach(meta_4w_cohorts_demog)
fixed_4w_demog_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_4w_demog_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_4w_cohorts_demog)


# M-A of 4+ health res from cohorts:

attach(res_4w_cohorts_health)

meta_4w_cohorts_health <-metagen(log(coef),
                                 lower = log(lci),
                                 upper = log(uci),
                                 level.ci = 0.95,
                                 studlab=paste(Study),
                                 subgroup = Group_labs, 
                                 sm = "OR",
                                 title = n)


summary(meta_4w_cohorts_health)


# Plot
tiff("v8_IPW_4plus_model2_health.tiff", width = 300, height = 600, units = "mm", res = 100)
forest(meta_4w_cohorts_health, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       sortvar = w.fixed,
       bylab = meta_4w_cohorts_health$Group_labs,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off() 

detach(res_4w_cohorts_health)


attach(meta_4w_cohorts_health)
fixed_4w_health_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_4w_health_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_4w_cohorts_health)

# Combining demog and health res together for export:
fixed_4w_res = rbind(fixed_4w_demog_res, fixed_4w_health_res)
write.csv(fixed_4w_res, "2021_12_06_IPW_4plus_model2_fixedmeta_res_v8.csv", row.names=F, quote=F)



# Plotting 12+ week results for cohorts:

table(res_12w$Exposure)

res_12w_cohorts = res_12w %>% 
  filter(Grouping!="Age") %>% 
  filter(Exposure!="pre-pandemic mh score")  %>% 
  filter(Exposure!="prepandemic mh score NO SOMATIC SUBSCALE")  %>% 
  filter(Study!="OpenSAFELY") %>%
  filter(Grouping!="Age") %>%
  filter(Grouping!="Social_class")


# Labelling grouping variables:

table(res_12w_cohorts$Grouping)

res_12w_cohorts$Grouping_f = factor(res_12w_cohorts$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking', 'Gen_health', 
                                                                         'M_Health_binary', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 
                                                                         'Asthma'))
# Splitting results to plot demographic and health results separately:

res_12w_cohorts_demog = filter(res_12w_cohorts, as.numeric(Grouping_f)<6)


res_12w_cohorts_demog$Group_labs = factor(res_12w_cohorts_demog$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking'))

levels(res_12w_cohorts_demog$Group_labs) = c("Female sex (ref. male)", "Other ethnicity (ref. white)", "No higher education (ref. degree attained)", "Index of multiple deprivation (per one decile)", "Current smoker")


res_12w_cohorts_health = filter(res_12w_cohorts, as.numeric(Grouping_f)>=6)
table(res_12w_cohorts_health$Grouping)

res_12w_cohorts_health$Group_labs = factor(res_12w_cohorts_health$Grouping, levels = c('Gen_health', 'M_Health_binary', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 'Asthma'))

levels(res_12w_cohorts_health$Group_labs) = c("Poor overall health", "Psychological distress", "BMI (per 5 kg/m^2)", "Overweight or obese", "Diabetes", "Hypertension", "High cholesterol", "Asthma")


########## M-A of 12+ week demography res from cohorts:

attach(res_12w_cohorts_demog)

meta_12w_cohorts_demog <-metagen(log(coef), 
                                 lower = log(lci),
                                 upper = log(uci),
                                 studlab=paste(Study),
                                 subgroup = Group_labs, 
                                 sm = "OR",
                                 title = n)

summary(meta_12w_cohorts_demog)

tiff("v8_IPW_12plus_model2_demog.tiff", width = 300, height = 360, units = "mm", res = 100)
forest(meta_12w_cohorts_demog, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_12w_cohorts_demog$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 12+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold", 
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off()

detach(res_12w_cohorts_demog)


attach(meta_12w_cohorts_demog)
fixed_12w_demog_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_12w_demog_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_12w_cohorts_demog)




# M-A of 12+ health res from cohorts:

attach(res_12w_cohorts_health)

meta_12w_cohorts_health <- metagen(log(coef),
                                   lower = log(lci),
                                   upper = log(uci),
                                   level.ci = 0.95,
                                   studlab=paste(Study),
                                   subgroup = Group_labs, 
                                   sm = "OR",
                                   title = n)



summary(meta_12w_cohorts_health)

tiff("v8_IPW_12plus_model2_health.tiff", width = 300, height = 600, units = "mm", res = 100)
forest(meta_12w_cohorts_health, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_12w_cohorts_health$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 12+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off()

detach(res_12w_cohorts_health)

attach(meta_12w_cohorts_health)
fixed_12w_health_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_12w_health_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_12w_cohorts_health)


# Combining demog and health res together for export:
fixed_12w_res = rbind(fixed_12w_demog_res, fixed_12w_health_res)
write.csv(fixed_12w_res, "2021_12_06_IPW_12plus_model2_fixedmeta_res_v8.csv", row.names=F, quote=F)



# Forest plots for age results only with IPW added 

res_age = res %>% 
  filter(Grouping=="Age")

res_age_12w = filter(res_age, Outcome=="12+ weeks")

res_age_4w = filter(res_age, Outcome!="12+ weeks")



# Plotting 4 week results inc OpenSAFELY:

# Reassigning study order:

res_age_4w$study_order = NA
res_age_4w$study_order[res_age_4w$Study == "MCS"] = 1
res_age_4w$study_order[res_age_4w$Study == "NS"] = 2
res_age_4w$study_order[res_age_4w$Study == "USoc"] = 3
res_age_4w$study_order[res_age_4w$Study == "TwinsUK"] = 4
res_age_4w$study_order[res_age_4w$Study == "GS"] = 5
res_age_4w$study_order[res_age_4w$Study == "ALSPAC G0"] = 6
res_age_4w$study_order[res_age_4w$Study == "OpenSAFELY"] = 7

table(res_age_4w$Exposure)
res_age_4w$Exposure[res_age_4w$Exposure == "age(70+vs18-44)"] = "70+"
res_age_4w$Exposure[res_age_4w$Exposure == "age(45-69+vs18-44)"] = "45-69"
res_age_4w$Exposure[res_age_4w$Exposure == "age(linear)"] = "Per year increase"
res_age_4w$Exposure[res_age_4w$Exposure == "age (cont)"] = "Per year increase"
res_age_4w$Exposure[res_age_4w$Exposure == "age (60-76 vs 45-54)"] = "60-76"
res_age_4w$Exposure[res_age_4w$Exposure == "age (55-59 vs 45-54)"] = "55-59"
res_age_4w$Exposure[res_age_4w$Exposure == "age(45+vs25-34)"] = "45+"
res_age_4w$Exposure[res_age_4w$Exposure == "age(35-44+vs25-34)"] = "35-44"
res_age_4w$Exposure[res_age_4w$Exposure == "25-34 (all vs 18-24)"] = "25-34"
res_age_4w$Exposure[res_age_4w$Exposure == "age (70+vs18-44)"] = "70+"
res_age_4w$Exposure[res_age_4w$Exposure == "age (45-69+vs18-44)"] = "45-69"


res_age_4w = res_age_4w %>%
  arrange(study_order)

res_age_4w$Group_labs = factor(res_age_4w$Study, levels = c('MCS', 'NS', 'USoc', 'TwinsUK', 'GS', 'ALSPAC G0', 'OpenSAFELY'))

levels(res_age_4w$Group_labs) = c("MCS (linear)", "NS (linear)", "USoc (ref. 18-44y)", "TwinsUK (ref. 18-44y)", "GS (ref. 18-44y)", "ALSPAC G0 (ref. 45-54y)", "OpenSAFELY (ref. 18-24y)")


attach(res_age_4w)

meta_4w_age <- metagen(log(coef),
                       lower = log(lci),
                       upper = log(uci),
                       level.ci = 0.95,
                       studlab=paste(Exposure),
                       subgroup = Group_labs, 
                       sm = "OR")

tiff("v8_IPW_4plus_ageonly.tiff", width = 280, height = 250, units = "mm", res = 100)
forest(meta_4w_age, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab"),
       leftlabs=c("Study"),
       bylab = res_age_4w$Group_labs,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       weight.study = "same",
       squaresize = 0,
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold", 
       just = "right",
       subgroup = FALSE,
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       hetstat = FALSE)
dev.off()

detach(res_age_4w)



# Plotting 12 week results without OpenSAFELY:

# Reassigning study order:

res_age_12w$study_order = NA
res_age_12w$study_order[res_age_12w$Study == "MCS"] = 1
res_age_12w$study_order[res_age_12w$Study == "NS"] = 2
res_age_12w$study_order[res_age_12w$Study == "USoc"] = 3
res_age_12w$study_order[res_age_12w$Study == "TwinsUK"] = 4
res_age_12w$study_order[res_age_12w$Study == "GS"] = 5
res_age_12w$study_order[res_age_12w$Study == "ALSPAC G0"] = 6

table(res_age_12w$Exposure)
res_age_12w$Exposure[res_age_12w$Exposure == "age(70+vs18-44)"] = "70+"
res_age_12w$Exposure[res_age_12w$Exposure == "age(45-69+vs18-44)"] = "45-69"
res_age_12w$Exposure[res_age_12w$Exposure == "age(linear)"] = "Per year increase"
res_age_12w$Exposure[res_age_12w$Exposure == "age (cont)"] = "Per year increase"
res_age_12w$Exposure[res_age_12w$Exposure == "age (60-76 vs 45-54)"] = "60-76"
res_age_12w$Exposure[res_age_12w$Exposure == "age (55-59 vs 45-54)"] = "55-59"
res_age_12w$Exposure[res_age_12w$Exposure == "age(45+vs25-34)"] = "45+"
res_age_12w$Exposure[res_age_12w$Exposure == "age(35-44+vs25-34)"] = "35-44"
res_age_12w$Exposure[res_age_12w$Exposure == "25-34 (all vs 18-24)"] = "25-34"
res_age_12w$Exposure[res_age_12w$Exposure == "age (70+vs18-44)"] = "70+"
res_age_12w$Exposure[res_age_12w$Exposure == "age (45-69+vs18-44)"] = "45-69"


res_age_12w = res_age_12w %>%
  arrange(study_order)

res_age_12w$Group_labs = factor(res_age_12w$Study, levels = c('MCS', 'NS', 'USoc', 'TwinsUK', 'GS', 'ALSPAC G0'))

levels(res_age_12w$Group_labs) = c("MCS (linear)", "NS (linear)", "USoc (ref. 18-44y)", "TwinsUK (ref. 18-44y)", "GS (ref. 18-44y)", "ALSPAC G0 (ref. 45-54y)")



attach(res_age_12w)

meta_12w_age <- metagen(log(coef),
                        lower = log(lci),
                        upper = log(uci),
                        level.ci = 0.95,
                        studlab=paste(Exposure),
                        subgroup = Group_labs, 
                        sm = "OR")

tiff("v8_IPW_12plus_ageonly.tiff", width = 280, height = 250, units = "mm", res = 100)
forest(meta_12w_age, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab"),
       leftlabs=c("Study"),
       bylab = res_age_12w$Group_labs,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 12+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       weight.study = "same",
       squaresize = 0,
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold", 
       just = "right",
       subgroup = FALSE,
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       hetstat = FALSE)
dev.off()

detach(res_age_12w)



# 8) Plotting of seropositive / PCR+ subgroup analysis --------



res_pos <- read_delim(file.path(File_path,"/2021_12_06_seroPCR_positive_model2_unweighted_regressions_reformatted.txt"),
                  "\t", escape_double = FALSE, trim_ws = TRUE)


res_pos = rename(res_pos, Study = cohort, se = coef_se, lci = lower_ci, uci = upper_ci, Exposure = exposure, Outcome = outcome) # colnames renamed to match previous script


# Adding log-odds to data:

res_pos = res_pos %>% mutate(logodds = log(coef))


res_pos = res_pos %>% 
  mutate(logodds = ifelse(Grouping == "BMI", logodds*5, logodds)) %>%
  mutate(coef = ifelse(Grouping == "BMI", exp(logodds), coef)) %>%
  mutate(lci = ifelse(Grouping == "BMI", exp(logodds - (1.96*5*se)), lci)) %>%
  mutate(uci = ifelse(Grouping == "BMI", exp(logodds + (1.96*5*se)), uci))




res_pos <- mutate(res_pos,
              study_order = case_when(
                Study == "ALSPAC G0" ~ 1,
                Study == "ALSPAC G1" ~ 2,
                Study == "TwinsUK" ~ 3))



# Subsets of res for 12+ and 4+ week outcomes:

res_12w = filter(res_pos, Outcome=="12+ weeks")

res_4w = filter(res_pos, Outcome!="12+ weeks")




# Plotting 4+ week results for sero/PCR+ subgroup:

# removing results that won't be meta-analysed:
res_4w_cohorts = res_4w %>% 
  filter(Grouping!="Age") %>% 
  filter(Exposure!="pre-pandemic mh score")  %>% 
  filter(Exposure!="prepandemic mh score NO SOMATIC SUBSCALE")  %>% 
  filter(Study!="OpenSAFELY") %>%
  filter(Grouping!="Social_class")

table(res_4w_cohorts$Grouping)
table(res_4w_cohorts$Exposure)

# Labelling grouping variables:

table(res_4w_cohorts$Grouping)


res_4w_cohorts$Grouping_f = factor(res_4w_cohorts$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking', 'Gen_health', 
                                                                       'M_Health', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 
                                                                       'Asthma'))



# Splitting results to plot demographic and health results separately:

res_4w_cohorts_demog = filter(res_4w_cohorts, as.numeric(Grouping_f)<6)

res_4w_cohorts_demog$Group_labs = factor(res_4w_cohorts_demog$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking'))

levels(res_4w_cohorts_demog$Group_labs) = c("Female sex (ref. male)", "Other ethnicity (ref. white)", "No higher education (ref. degree attained)", "Index of multiple deprivation (per one decile)", "Current smoking")


res_4w_cohorts_health = filter(res_4w_cohorts, as.numeric(Grouping_f)>=6)
table(res_4w_cohorts_health$Grouping)

res_4w_cohorts_health$Group_labs = factor(res_4w_cohorts_health$Grouping, levels = c('Gen_health', 'M_Health', 'BMI', 'Overweight_obese', 'Diabetes', 'Hypertension', 'Cholesterol', 'Asthma'))

levels(res_4w_cohorts_health$Group_labs) = c("Poor overall health", "Psychological distress", "BMI (per 5 kg/m^2)", "Overweight or obese", "Diabetes", "Hypertension", "High cholesterol", "Asthma")

# M-A of 4+ week demography res from cohorts:

attach(res_4w_cohorts_demog)

meta_4w_cohorts_demog <-metagen(log(coef), 
                                lower = log(lci),
                                upper = log(uci),
                                studlab=paste(Study),
                                subgroup = Group_labs, 
                                sm = "OR",
                                title = n)

summary(meta_4w_cohorts_demog)

tiff("v8_4plus_seroPCRpos_subgroup_model3_demog.tiff", width = 300, height = 240, units = "mm", res = 100)
forest(meta_4w_cohorts_demog, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_4w_cohorts_demog$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold",
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off() 

detach(res_4w_cohorts_demog)

# creating table of fixed effects MA output for subgroup res:

str(meta_4w_cohorts_demog, list.len = 153) # see full number of items in list
attach(meta_4w_cohorts_demog)
fixed_4w_demog_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_4w_demog_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_4w_cohorts_demog)


# M-A of 4+ health res from cohorts:

attach(res_4w_cohorts_health)

meta_4w_cohorts_health <-metagen(log(coef),
                                 lower = log(lci),
                                 upper = log(uci),
                                 level.ci = 0.95,
                                 studlab=paste(Study),
                                 subgroup = Group_labs, 
                                 sm = "OR",
                                 title = n)


summary(meta_4w_cohorts_health)


# Plot
tiff("v8_4plus_seroPCRpos_subgroup_model3_health.tiff", width = 300, height = 300, units = "mm", res = 100)
forest(meta_4w_cohorts_health, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       sortvar = w.fixed,
       bylab = meta_4w_cohorts_health$Group_labs,
       colgap = "4mm",
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 4+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       test.overall.fixed = FALSE,
       test.overall.random = FALSE,
       ff.lr = "bold", 
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off() 

detach(res_4w_cohorts_health)


attach(meta_4w_cohorts_health)
fixed_4w_health_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_4w_health_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_4w_cohorts_health)

# Combining demog and health res together for export:
fixed_4w_res = rbind(fixed_4w_demog_res, fixed_4w_health_res)
write.csv(fixed_4w_res, "2021_12_06_4plus_seroPCRpos_subgroup_model3_fixedmeta_res_v8.csv", row.names=F, quote=F)


# Plotting 12+ week results for sero/PCR+ subgroup analysis:

res_12w_cohorts = res_12w %>% 
  filter(Grouping!="Age") %>% 
  filter(Exposure!="pre-pandemic mh score")  %>% 
  filter(Exposure!="prepandemic mh score NO SOMATIC SUBSCALE")  %>% 
  filter(Study!="OpenSAFELY") %>%
  filter(Grouping!="Social_class")


# Labelling grouping variables:

res_12w_cohorts$Grouping_f = factor(res_12w_cohorts$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking', 'Gen_health', 
                                                                         'M_Health', 'BMI', 'Overweight_obese', 'Hypertension', 'Cholesterol', 
                                                                         'Asthma'))
# Splitting results to plot demographic and health results separately:

res_12w_cohorts_demog = filter(res_12w_cohorts, as.numeric(Grouping_f)<6)


res_12w_cohorts_demog$Group_labs = factor(res_12w_cohorts_demog$Grouping, levels = c('Sex', 'Ethnicity', 'Education', 'SEP', 'Smoking'))

levels(res_12w_cohorts_demog$Group_labs) = c("Female sex (ref. male)", "Other ethnicity (ref. white)", "No higher education (ref. degree attained)", "Index of multiple deprivation (per one decile)", "Current smoker")


res_12w_cohorts_health = filter(res_12w_cohorts, as.numeric(Grouping_f)>=6)
table(res_12w_cohorts_health$Grouping)

res_12w_cohorts_health$Group_labs = factor(res_12w_cohorts_health$Grouping, levels = c('Gen_health', 'M_Health', 'BMI', 'Overweight_obese', 'Hypertension', 'Cholesterol', 'Asthma'))

levels(res_12w_cohorts_health$Group_labs) = c("Poor overall health", "Psychological distress", "BMI (per 5 /m^2)", "Overweight or obese", "Hypertension", "High cholesterol", "Asthma")


########## M-A of 12+ week demography res from cohorts:

attach(res_12w_cohorts_demog)

meta_12w_cohorts_demog <-metagen(log(coef), 
                                 lower = log(lci),
                                 upper = log(uci),
                                 studlab=paste(Study),
                                 subgroup = Group_labs, 
                                 sm = "OR",
                                 title = n)

summary(meta_12w_cohorts_demog)

tiff("v8_12plus_seroPCRpos_subgroup_model3_demog.tiff", width = 300, height = 240, units = "mm", res = 100)
forest(meta_12w_cohorts_demog, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_12w_cohorts_demog$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 12+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold", 
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off()

detach(res_12w_cohorts_demog)


attach(meta_12w_cohorts_demog)
fixed_12w_demog_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_12w_demog_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_12w_cohorts_demog)




# M-A of 12+ health res from cohorts:

attach(res_12w_cohorts_health)

meta_12w_cohorts_health <- metagen(log(coef),
                                   lower = log(lci),
                                   upper = log(uci),
                                   level.ci = 0.95,
                                   studlab=paste(Study),
                                   subgroup = Group_labs, 
                                   sm = "OR",
                                   title = n)



summary(meta_12w_cohorts_health)

tiff("v8_12plus_seroPCRpos_subgroup_model3_health.tiff", width = 300, height = 300, units = "mm", res = 100)
forest(meta_12w_cohorts_health, 
       overall = FALSE,
       overall.hetstat = FALSE,
       leftcols=c("studlab", "title"),
       leftlabs=c("", "N"),
       bylab = meta_12w_cohorts_health$Group_labs,
       colgap = "4mm",
       sortvar = w.fixed,
       colgap.studlab = "20mm",
       plotwidth = "10cm",
       xlab = "Odds ratio for symptoms lasting 12+ weeks",
       ff.test.subgroup = "italic",
       col.by="dimgray",
       label.right="Higher risk",
       label.left="Lower risk",
       colgap.forest = "8mm",
       bottom.lr = FALSE,
       ff.lr = "bold", 
       test.subgroup.fixed = FALSE,
       test.subgroup.random = FALSE,
       label.test.subgroup.fixed = FALSE,
       label.test.subgroup.random = FALSE,
       print.subgroup.name = F, 
       just = "right")
dev.off()

detach(res_12w_cohorts_health)

attach(meta_12w_cohorts_health)
fixed_12w_health_res = as.data.frame(cbind(bylevs, k.w, exp(TE.fixed.w), exp(lower.fixed.w), exp(upper.fixed.w), pval.fixed.w, I2.w))
colnames(fixed_12w_health_res) = c("Trait", "# Studies", "OR", "lower CI", "upper CI", "P", "I2")
detach(meta_12w_cohorts_health)


# Combining demog and health res together for export:
fixed_12w_res = rbind(fixed_12w_demog_res, fixed_12w_health_res)
write.csv(fixed_12w_res, "2021_12_06_12plus_seroPCRpos_subgroup_model3_fixedmeta_res_v8.csv", row.names=F, quote=F)
