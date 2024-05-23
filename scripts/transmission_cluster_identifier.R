#Load libraries:
library(tidyverse)
library(reshape2)

#d12 clusters-----------------------------------------------------------------------
#Read in d12 groups
d12_groups <- read.delim("d12_groups.csv")

#Delete ungrouped isolates
d12_groups <- filter(d12_groups,group!="ungrouped")

#Join metadata to groups
d12_groups <- left_join(d12_groups,metadata_all,by=c("sequence"="nrz_id_1"))

#Identify groups containing patient IDs
d12_groups$group <- as.factor(d12_groups$group)
length(levels(d12_groups$group)) 

#For loop to identify clusters with more than one patient
transmission_frame <- d12_groups[1,]  
transmission_frame <- transmission_frame[-1,]  
for (i in levels(d12_groups$group)){
  print(i)
  group_frame <- filter(d12_groups,group==i)
  if (length(unique(group_frame$patient_id))>1){
    transmission_frame <- rbind(transmission_frame,group_frame)
  }
}

#How many clusters?
length(unique(transmission_frame$group))
length(unique(transmission_frame$patient_id))
length(unique(transmission_frame$site))
length(unique(transmission_frame$country))
write.csv(transmission_frame[,c("sequence","group")],file="transmission_frame.csv",row.names=FALSE)

#For loop to identify clusters with more than one site
transmission_frame_site <- d12_groups[1,]  
transmission_frame_site <- transmission_frame_site[-1,]  
for (i in levels(d12_groups$group)){
  print(i)
  group_frame <- filter(d12_groups,group==i)
  if (length(unique(group_frame$site))>1){
    transmission_frame_site <- rbind(transmission_frame_site,group_frame)
  }
}

#How many clusters?
length(unique(transmission_frame_site$group))
length(unique(transmission_frame_site$patient_id))
length(unique(transmission_frame_site$site))

write.csv(transmission_frame_site, file="transmission_frame_diff_site.csv",row.names=FALSE)

#d12 clusters by country------------------------------------------------------------
#For loop to identify clusters with more than one country
transmission_frame_country <- d12_groups[1,]  
transmission_frame_country <- transmission_frame_country[-1,]  
for (i in levels(d12_groups$group)){
  print(i)
  group_frame <- filter(d12_groups,group==i)
  if (length(unique(group_frame$country))>1){
    transmission_frame_country <- rbind(transmission_frame_country,group_frame)
  }
}

#How many clusters?
length(unique(transmission_frame_country$group))
length(unique(transmission_frame_country$patient_id))
length(unique(transmission_frame_country$country))
length(unique(transmission_frame_country$site))

#d12 clusters by continent----------------------------------------------------------
#For loop to identify clusters with more than one continent
transmission_frame_continent <- d12_groups[1,]  
transmission_frame_continent <- transmission_frame_continent[-1,]  
for (i in levels(d12_groups$group)){
  print(i)
  group_frame <- filter(d12_groups,group==i)
  if (length(unique(group_frame$continent))>1){
    transmission_frame_continent <- rbind(transmission_frame_continent,group_frame)
  }
}

#How many clusters?
length(unique(transmission_frame_continent$group))
length(unique(transmission_frame_continent$patient_id))
length(unique(transmission_frame_continent$continent))
length(unique(transmission_frame_continent$site))

write.csv(transmission_frame_continent, file="transmission_frame_continent.csv",row.names=FALSE)
