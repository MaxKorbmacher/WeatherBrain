#! /usr/bin/env Rscript
################################
# Data cleaning and assembling #
################################
# load packages
library(dplyr)
# load data
pheno = read.csv("/cluster/projects/p33/users/maxk/UKB/environment/data/demo/environment.csv")
dMRI_cross = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/cross_sectional.csv")
dMRI_T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/T1.csv")
dMRI_T2 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/T2.csv")

# rename phenotypes

# we use baseline income (that is not assessed right at the imaging time point, but before)
# why? due to missingness in the other variables!
names(pheno)[names(pheno) == "X738.0.0"] <- "income"
pheno$income = ifelse(pheno$income == 1, "less_18k",pheno$income)
pheno$income = ifelse(pheno$income == 2, "18k_31k",pheno$income)
pheno$income = ifelse(pheno$income == 3, "31k_52k",pheno$income)
pheno$income = ifelse(pheno$income == 4, "52k_100k",pheno$income)
pheno$income = ifelse(pheno$income == 5, "more_100k",pheno$income)
pheno$income = ifelse(pheno$income == -1, "idk",pheno$income)
pheno$income = ifelse(pheno$income == -3, "no_answer",pheno$income)

# names(pheno)[names(pheno) == "X738.4.0"] <- "income2" # does not exist
# pheno$income = ifelse(pheno$income == 1, "less_18k",pheno$income)
# pheno$income = ifelse(pheno$income == 2, "18k_31k",pheno$income)
# pheno$income = ifelse(pheno$income == 3, "31k_52k",pheno$income)
# pheno$income = ifelse(pheno$income == 4, "52k_100k",pheno$income)
# pheno$income = ifelse(pheno$income == 5, "more_100k",pheno$income)
# pheno$income = ifelse(pheno$income == -1, "idk",pheno$income)
# pheno$income = ifelse(pheno$income == -3, "no_answer",pheno$income)


#names(pheno)[names(pheno) == "X31.0.0"] <- "sex"
names(pheno)[names(pheno) == "X33.0.0"] <- "birthday"
names(pheno)[names(pheno) == "X34.0.0"] <- "birthyear"
names(pheno)[names(pheno) == "X34.0.0"] <- "birthyear"
names(pheno)[names(pheno) == "X52.0.0"] <- "birthmonth"

#names(pheno)[names(pheno) == "X130.0.0"] <- "pob_UK"
#names(pheno)[names(pheno) == "X129.0.0"] <- "pob_non_UK"
#names(pheno)[names(pheno) == "X120115.0.0"] <- "country_of_birth"
names(pheno)[names(pheno) == "X1647.0.0"] <- "country_of_birth"

pheno$country_of_birth = ifelse(pheno$country_of_birth == 1, "England",pheno$country_of_birth)
pheno$country_of_birth = ifelse(pheno$country_of_birth == 2, "Wales",pheno$country_of_birth)
pheno$country_of_birth = ifelse(pheno$country_of_birth == 3, "Scotland",pheno$country_of_birth)
pheno$country_of_birth = ifelse(pheno$country_of_birth == 4, "Northern_Ireland",pheno$country_of_birth)
pheno$country_of_birth = ifelse(pheno$country_of_birth == 5, "Republic_Ireland",pheno$country_of_birth)
pheno$country_of_birth = ifelse(pheno$country_of_birth == 6, "Elsewhere",pheno$country_of_birth)
pheno$country_of_birth = ifelse(pheno$country_of_birth == -1, NA,pheno$country_of_birth)
pheno$country_of_birth = ifelse(pheno$country_of_birth == -3, NA,pheno$country_of_birth)
table(pheno$country_of_birth)


# merge
T1 = merge(pheno %>% select(eid,birthday, birthmonth, birthyear,country_of_birth, income),dMRI_T1 %>% select(eid, sex, age, site))
T2 = merge(pheno %>% select(eid,birthday, birthmonth, birthyear,country_of_birth, income),dMRI_T2 %>% select(eid, sex, age, site))
cross = merge(pheno %>% select(eid,birthday, birthmonth, birthyear,country_of_birth, income),dMRI_cross %>% select(eid, sex, age, site))

# rename income at tp2... if it would exist ...
# names(T2)[names(T2) == "income2"] <- "income"

# write tables
write.csv(T1, "/cluster/projects/p33/users/maxk/UKB/environment/data/demo/T1.csv")
write.csv(T2, "/cluster/projects/p33/users/maxk/UKB/environment/data/demo/T2.csv")
write.csv(cross, "/cluster/projects/p33/users/maxk/UKB/environment/data/demo/cross.csv")
