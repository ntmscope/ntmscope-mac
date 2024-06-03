#Load libraries
library(lubridate)
library(DescTools)

#For loop for possible transmission occasion between all possible patients-------
setwd("YOUR_FOLDER")

#Create file for output
patients_combs$tm_no <- NULL #reset results
sink("ta_results.txt", split=TRUE)
#Start for loop
for (x in 1:nrow(patients_combs)){
  print(paste(patients_combs[x,1],"vs.",patients_combs[x,2]), quote=FALSE)
  print("---------------------------------", quote=FALSE)
  print("", quote=FALSE)
  #Object for the two patients of the comparison
  pat1_time <- sum_frame[sum_frame$newkey1==paste(patients_combs[x,1]),]
  pat2_time <- sum_frame[sum_frame$newkey1==paste(patients_combs[x,2]),]
  count <- 0 # Set transmission opportunity count to zero
  #Nested for loop for the two regarding patients
  for (i in 1:length(pat1_time$date1)){
    for (y in 1:length(pat2_time$date1)){
      result <- as.Date(c(pat1_time$date1[i],pat1_time$date2[i])) %overlaps% as.Date(c(pat2_time$date1[y],pat2_time$date2[y]))
      if(result){
        count <- count+1
        print(count)
        print(paste(pat1_time$date1[i]), quote=FALSE)
        print(paste(pat1_time$date2[i]), quote=FALSE)
        print(paste(pat1_time$ward[i]), quote=FALSE)
        print(pat2_time$date1[y], quote=FALSE)
        print(pat2_time$date2[y], quote=FALSE)
        print(paste(pat2_time$ward[y]), quote=FALSE)
        print(result, quote=FALSE)
        print("", quote=FALSE)
      }
    }
    patients_combs$tm_no[x] <- count
  }
}
sink()
View(patients_combs)
write.csv(file="patients_combs_result.csv",patients_combs)
