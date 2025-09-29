# get T1-w stats
# process them so that eids match

library(dplyr)
library(stringr)
library("magrittr")


# LOAD AND FIX EID NAME FOR EULER NUMBERS
Euler = read.csv("/cluster/projects/p33/groups/imaging/ukbio/stats/batch_50k/allEuler_UKB.csv")
# euler nb
Euler = Euler %>% dplyr::rename("eid"= "subject") 
# remove "FS_" prefix
Euler$eid = substring(Euler[,1],4)
Euler = Euler %>% subset(select = c(eid,euler_lh, euler_rh))
write.table(Euler, file = "/cluster/projects/p33/users/maxk/UKB/data/T1w_50k/allEuler_UKB.txt", col.names = T, row.names = F, sep = " ")


# LOAD AND FIX EID NAME FOR FS OUTPUTS
files = c("lh.area.UKBB.txt",
"lh.thicknessstd.UKBB.txt",
"lh.thickness.UKBB.txt",
"lh.volume.UKBB.txt",
"rh.area.UKBB.txt",
"rh.thicknessstd.UKBB.txt",
"rh.thickness.UKBB.txt",
"rh.volume.UKBB.txt",
"subcorticalstats.UKBB.txt")
T1 = list()
for (i in files){
  T1[[i]] = read.delim(paste("/cluster/projects/p33/groups/imaging/ukbio/stats/batch_50k/regularFSstats/",i, sep=""))
  T1[[i]][,1] = substring(T1[[i]][,1],4)
  colnames(T1[[i]])[1] = "eid"
}
# save files
for (i in files){
  write.table(T1[[i]], file = paste("/cluster/projects/p33/users/maxk/UKB/data/T1w_50k/",i,sep=""), col.names = T, row.names = F, sep = " ")
}
