# data wrangling to connect environmental and MRI data
# 
# Max Korbmacher, July 2025
#
# ------------------------------------------------- #
# -----------------Contents------------------------ #
# ------------------------------------------------- #
# 1. Prep------------------------------------------ #
# 1.1 Packages------------------------------------- #
# 1.2 Load data------------------------------------ #
# 1.3 Label data----------------------------------- #
# 1.4 Merge data----------------------------------- #
# ------------------------------------------------- #
# ------------------------------------------------- #
#
#
#
# 1.Prep-------------------------------------------

# 1.1 Packages-------------------------------------

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, dplyr, ggplot2, reshape2,lme4,ggpubr,
               psych, lmerTest, performance, broom)

# 1.2 Load data------------------------------------

# Specify the data path
datapath = "/cluster/projects/p33/users/maxk/UKB/environment/data/"
# get mental health data
mh = read.csv("/cluster/projects/p33/users/maxk/UKB/mentalhealth/long_data.csv")
mh_cross = read.csv("/cluster/projects/p33/users/maxk/UKB/mentalhealth/pheno_cross.csv")
names(mh_cross)
# brain age data
ba1 = read.csv("/cluster/projects/p33/users/maxk/UKB/longBAG/data/dMRI_test1.csv")
ba2 = read.csv("/cluster/projects/p33/users/maxk/UKB/longBAG/data/dMRI_test2.csv")
ba_long = rbind(ba1,ba2)
ba_train = read.csv("/cluster/projects/p33/users/maxk/UKB/longBAG/data/dMRI_train.csv")
# poligenic risk scores
pgrs = read.csv("/cluster/projects/p33/users/maxk/UKB/genetics/PRS.csv")
# demo
demo = read.csv(paste(datapath,"demo/cross.csv",sep=""))
#Put all csv files from the env_data and demo directory into a list
extract_data = function(data_frame_location){
  filelist = list.files(pattern="*.csv$",
                      path = paste(datapath,data_frame_location,sep=""),
                      full.names = T)
  out_list <- lapply(filelist, fread)
  #get the filenames, remove path and .csv extension
  names(out_list) <- gsub(filelist, pattern=paste(datapath,data_frame_location,"/",sep=""), replacement="")
  names(out_list) <- gsub(names(out_list), pattern=".csv", replacement="")
  return(out_list)
}
env_list = extract_data("env_data/")
MRI_list = extract_data("demo/")
# load csv from env_data/two_timepoints/ 
# those are currently the cross-sectional data
c1 = read.csv(paste(datapath,"env_data/two_timepoints/climate_data_two_timepoints.csv",sep=""))
c2 = read.csv(paste(datapath,"env_data/two_timepoints/climate_data_two_timepoints_2.csv",sep=""))

# Now, load MRI data
## diffusion
dMRI_cross = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/cross_sectional.csv")
dMRI_T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/T1.csv")
dMRI_T2 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/T2.csv")
## T1w
T1w = read.delim("/cluster/projects/p33/users/maxk/UKB/data/T1w_50k/subcorticalstats.UKBB.txt",sep=" ")
T1w_long = read.csv("/cluster/projects/p33/groups/imaging/ukbio_R/recon/aseg_stats.txt",sep="\t")

## QC (T1w)
euler = read.delim("/cluster/projects/p33/users/maxk/UKB/data/T1w_50k/allEuler_UKB.csv",sep=" ")

## 4SD over
euler = na.omit(euler)

# 0 for keep, 1 throw out
a = ifelse(abs(euler$euler_lh-mean(euler$euler_lh)) < sd(euler$euler_lh)*3,0,1)
b = ifelse(abs(euler$euler_rh-mean(euler$euler_rh)) < sd(euler$euler_rh)*3,0,1)
euler$kick = a+b
kick=euler %>% filter(kick == 1) %>% pull(eid)
T1w = T1w[!T1w$eid %in% kick,]

colnames(T1w_long)[1] = "eid"
T1w_long$eid = gsub("_20216_3_0/","",T1w_long$eid)
T1w_long = T1w_long[!T1w_long$eid %in% kick,]
#dMRI_cross = merge(T1w,dMRI_cross,by="eid")
dMRI_T1 = merge(T1w,dMRI_T1,by="eid")
dMRI_T2 = merge(T1w_long,dMRI_T2,by="eid")

# 1.3 Label data-----------------------------------
#
# Data descriptor extracted_climate_data_T* files
# 
# avg_slhtf = Time-mean surface latent heat flux
# avg_snswrf = surface net short-wave radiation flux
# avg_sdlwrf = Time-mean surface downward long-wave radiation flux
# avg_lsprate = Time-mean large-scale precipitation rate
# uvb = Surface downward UV radiation
# avg_sdirswrf = Time-mean surface direct short-wave radiation flux
# avg_sduvrf = Time-mean surface downward UV radiation flux
# avg_snlwrf = Time-mean surface net long-wave radiation flux
# tp = Total precipitation
env_list["extracted_climate_data_T1"]
#
#
#
# Data descriptor extracted_climate_data_avg* files
#
# si10 = 10 metre wind speed
# d2m = 2 metre dewpoint temperature
# t2m = 2 metre temperature
# u10 = 10 metre U wind component
# v10 = 10 metre V wind component
# sp = surface pressure
#
env_list["extracted_climate_data_avgua"]
#
# 1.4 Merge data-----------------------------------
T1_env = merge(env_list["exctracted_climate_data_T1_avgua"]$exctracted_climate_data_T1_avgua
               %>% select(-birthyear,-birthmonth,-age),
               env_list["extracted_climate_data_T1"]$extracted_climate_data_T1, by = "eid")
T1_env = merge(dMRI_T1 %>% dplyr::select(eid,sex,site, contains("Mean"),CortexVol,CorticalWhiteMatterVol, SurfaceHoles, EstimatedTotalIntraCranialVol), T1_env, by = "eid")
T2_env = merge(env_list["exctracted_climate_data_T2_avgua"]$exctracted_climate_data_T2_avgua
               %>% select(-birthyear,-birthmonth,-age),
               env_list["extracted_climate_data_T2"]$extracted_climate_data_T2, by = "eid")
#T2_env = merge(MRI_list["T2"]$T2 %>% select(eid,sex,site, contains("Mean")), T2_env, by = "eid")
T2_env = merge(dMRI_T2 %>% select(eid,sex,site, contains("Mean"),CortexVol,CorticalWhiteMatterVol, SurfaceHoles, EstimatedTotalIntraCranialVol), T2_env, by = "eid")

cross_env = merge(data.frame(env_list["extracted_climate_data_avgua"]$extracted_climate_data_avgua),
                  dMRI_cross, by = "eid")

#env_list["extracted_climate_data_avgua"]$extracted_climate_data_avgua$eid %in% dMRI_cross$eid


# cross_env = merge(MRI_list["cross"]$cross%>% select(-birthyear,-birthmonth,-age), 
#                   env_list["extracted_climate_data_avgua"]$extracted_climate_data_avgua, by = "eid")
T1_env$session = 1
T2_env$session = 2
long_env = rbind(T1_env,T2_env, fill=T)

cross_env = merge(c1,c2, by = "eid")
cross_env = merge(cross_env,dMRI_cross, by = "eid")
cross_env = merge(cross_env,T1w, by = "eid")

# add covariates
# mental health
## N = N12 neuroticism score
mh = mh %>% select(eid,RDS,N,TP)
colnames(mh)[4] = "session"
mh$session = mh$session+1
long_env = merge(mh,long_env,by=c("eid","session"))
cross_env = merge(mh_cross,cross_env,by=c("eid"))

# polygenic risk scores
colnames(pgrs)[2] = "eid"
pgrs = pgrs %>% select(-X)
long_env = base::merge(pgrs,long_env,by.x="eid",by.y="eid")
cross_env = base::merge(pgrs,cross_env,by.x="eid",by.y="eid")

# brain age (only in longitudinal data!)
## train model
ba_train = ba_train %>% select(-c(sex, site, eid))
model = lm(age ~., data = ba_train)
## prep prediction data
ba_long$session = c(replicate(nrow(ba1), 0), replicate(nrow(ba2), 1))
## predict
ba_long$BrainAge = predict(model,ba_long)
ba_long = ba_long %>% select(eid, session, BrainAge)
ba_long$session = ba_long$session+1
long_env = merge(long_env,ba_long,by=c("eid","session"))

# cross_env$income
cross_env = merge(cross_env,demo%>%select(eid,income), by = "eid")
#merge(long_env,demo%>%select(eid,income), by = "eid")
#
write.csv(cross_env,paste(datapath,"final_cross.csv",sep=""))
write.csv(long_env,paste(datapath,"final_long.csv",sep=""))