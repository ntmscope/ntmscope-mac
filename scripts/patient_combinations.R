#Script for possible combinations of patients

#Load libraries
library(reshape2)
library(lubridate)
library(DescTools)


#Possible combinations without repetition
patients_combs <- t(combn(single_patients,2))
patients_combs <- data.frame(pat1=patients_combs[,1],pat2=patients_combs[,2])
#patients_combs$transmission <- 1 #Cave put coding here for possibility of transmission
patients_combs
patients_combs_inverse <- patients_combs
patients_combs_inverse$pat1 <- patients_combs$pat2
patients_combs_inverse$pat2 <- patients_combs$pat1
patients_combs_double <- rbind(patients_combs,patients_combs_inverse)