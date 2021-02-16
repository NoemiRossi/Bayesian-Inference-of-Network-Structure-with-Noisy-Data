#Import the first dataset
dataset1<-read.table("Datasets/High-School_data_2013.csv", quote="\"", comment.char="")

#Extract the nodes
nodes=c(dataset1$V2, dataset1$V3)
nodes<-unique(nodes)

#Extract the classes
classes=dataset1$V4

#Rescale the time so that the beginning is time 0.
dataset1[,1]=dataset1[,1]-min(dataset1[,1])

#Add the day of the observation in the dataset
for (i in 1:dim(dataset1)[1]){
  if(dataset1[i,1]<18000)
    dataset1[i,6]=1
  if(dataset1[i,1]>71000 & dataset1[i,1]<105000 )
    dataset1[i,6]=2
  if(dataset1[i,1]>158000 & dataset1[i,1]<191000 )
    dataset1[i,6]=3
  if(dataset1[i,1]>244000 & dataset1[i,1]<278000 )
    dataset1[i,6]=4
  if(dataset1[i,1]>331000)
    dataset1[i,6]=5
}

write.csv(dataset1,"Datasets/dataset1_withdays.csv", row.names=FALSE)

#Find all the possible couple of nodes
possible_couples<-as.data.frame(unique(cbind(dataset1[,2],dataset1[,3])))

#Create a dataset of the interaction of every couple in day1. 
#We consider an interaction if two people meet for 20 seconds, they stay together for any time, and if they meet again after at least 70 sec there is a new interaction.
primo_giorno<-dataset1[dataset1$V6==1,]
possible_couples1<-as.data.frame(unique(cbind(primo_giorno[,2],primo_giorno[,3])))


for(i in 1:dim(possible_couples1)[1]){
  prevoius=-20
  count=0
  couple<-primo_giorno[primo_giorno$V2==possible_couples1[i,1] & primo_giorno$V3==possible_couples1[i,2],]
  for(j in 1:dim(couple)[1]){
    if(couple[j,1]-prevoius>70)
      count=count+1
    prevoius=couple[j,1]
  }
  possible_couples1[i,3]=count
}

for(i in 1:dim(possible_couples1)[1]){
  if(possible_couples1[i,3]==0){
    possible_couples1[i,3]=1
  }
}

write.csv(possible_couples1,"Datasets/day1.csv", row.names=FALSE)


#Create a dataset of the interaction of every couple in day2. 
secondo_giorno<-dataset1[dataset1$V6==2,]
possible_couples2<-as.data.frame(unique(cbind(secondo_giorno[,2],secondo_giorno[,3])))


for(i in 1:dim(possible_couples2)[1]){
  prevoius=-20
  count=0
  couple<-secondo_giorno[secondo_giorno$V2==possible_couples2[i,1] & secondo_giorno$V3==possible_couples2[i,2],]
  for(j in 1:dim(couple)[1]){
    if(couple[j,1]-prevoius>70)
      count=count+1
    prevoius=couple[j,1]
  }
  possible_couples2[i,3]=count
}

write.csv(possible_couples2,"Datasets/day2.csv", row.names=FALSE)

#Create a dataset of the interaction of every couple in day3. 
terzo_giorno<-dataset1[dataset1$V6==3,]
possible_couples3<-as.data.frame(unique(cbind(terzo_giorno[,2],terzo_giorno[,3])))


for(i in 1:dim(possible_couples3)[1]){
  prevoius=-20
  count=0
  couple<-terzo_giorno[terzo_giorno$V2==possible_couples3[i,1] & terzo_giorno$V3==possible_couples3[i,2],]
  for(j in 1:dim(couple)[1]){
    if(couple[j,1]-prevoius>70)
      count=count+1
    prevoius=couple[j,1]
  }
  possible_couples3[i,3]=count
}

write.csv(possible_couples3,"Datasets/day3.csv", row.names=FALSE)

#Create a dataset of the interaction of every couple in day4. 
quarto_giorno<-dataset1[dataset1$V6==4,]
possible_couples4<-as.data.frame(unique(cbind(quarto_giorno[,2],quarto_giorno[,3])))


for(i in 1:dim(possible_couples4)[1]){
  prevoius=-20
  count=0
  couple<-quarto_giorno[quarto_giorno$V2==possible_couples4[i,1] & quarto_giorno$V3==possible_couples4[i,2],]
  for(j in 1:dim(couple)[1]){
    if(couple[j,1]-prevoius>70)
      count=count+1
    prevoius=couple[j,1]
  }
  possible_couples4[i,3]=count
}

write.csv(possible_couples4,"Datasets/day4.csv", row.names=FALSE)

#Create a dataset of the interaction of every couple in day5. 
quinto_giorno<-dataset1[dataset1$V6==5,]
possible_couples5<-as.data.frame(unique(cbind(quinto_giorno[,2],quinto_giorno[,3])))


for(i in 1:dim(possible_couples5)[1]){
  prevoius=-20
  count=0
  couple<-quinto_giorno[quinto_giorno$V2==possible_couples5[i,1] & quinto_giorno$V3==possible_couples5[i,2],]
  for(j in 1:dim(couple)[1]){
    if(couple[j,1]-prevoius>70)
      count=count+1
    prevoius=couple[j,1]
  }
  possible_couples5[i,3]=count
}

write.csv(possible_couples5,"Datasets/day5.csv", row.names=FALSE)

#Create a new dataset in which the interaction in every day is recorded
day1<-read.csv("Datasets/day1.csv")
day2<-read.csv("Datasets/day2.csv")
day3<-read.csv("Datasets/day3.csv")
day4<-read.csv("Datasets/day4.csv")
day5<-read.csv("Datasets/day5.csv")

possible_couples<-as.data.frame(unique(cbind(dataset1[,2],dataset1[,3])))

for(i in 1:dim(possible_couples)[1]){
  possible_couples[i,3]=0
  possible_couples[i,4]=0
  possible_couples[i,5]=0
  possible_couples[i,6]=0
  possible_couples[i,7]=0
  
}

for(i in 1:dim(possible_couples)[1]){
  if(length(possible_couples1[possible_couples1$V1==possible_couples[i,1] & possible_couples1$V2==possible_couples[i,2],3]))
    possible_couples[i,3]=possible_couples1[possible_couples1$V1==possible_couples[i,1] & possible_couples1$V2==possible_couples[i,2],3] 
  if(length(possible_couples2[possible_couples2$V1==possible_couples[i,1] & possible_couples2$V2==possible_couples[i,2],3]))  
    possible_couples[i,4]=possible_couples2[possible_couples2$V1==possible_couples[i,1] & possible_couples2$V2==possible_couples[i,2],3] 
  if(length(possible_couples3[possible_couples3$V1==possible_couples[i,1] & possible_couples3$V2==possible_couples[i,2],3]))  
    possible_couples[i,5]=possible_couples3[possible_couples3$V1==possible_couples[i,1] & possible_couples3$V2==possible_couples[i,2],3] 
  if(length(possible_couples4[possible_couples4$V1==possible_couples[i,1] & possible_couples4$V2==possible_couples[i,2],3] ))  
    possible_couples[i,6]=possible_couples4[possible_couples4$V1==possible_couples[i,1] & possible_couples4$V2==possible_couples[i,2],3] 
  if(length(possible_couples5[possible_couples5$V1==possible_couples[i,1] & possible_couples5$V2==possible_couples[i,2],3] ))  
    possible_couples[i,7]=possible_couples5[possible_couples5$V1==possible_couples[i,1] & possible_couples5$V2==possible_couples[i,2],3] 
    
}

write.csv(possible_couples,"Datasets/allday.csv", row.names=FALSE)

#Create a new dataset in which for every couple there is the sum of the interaction happened in the 5 days.
allday<-read.csv("Datasets/allday.csv")

final_data<-allday[,1:2]

for(i in 1:dim(allday)[1]){
  final_data[i,3]=sum(allday[i,3:7])
  
}

write.csv(final_data,"Datasets/new_dataset1.csv", row.names=FALSE)

new_dataset1=read.csv(("Datasets/new_dataset1.csv"))

#Associate the class to every couple
for(i in 1:dim(data_poisson)[1]) {
  new_dataset1[i,4]=dataset1[dataset1$V2==new_dataset1[i,1],4][1]
  new_dataset1[i,5]=dataset1[dataset1$V2==new_dataset1[i,2],5][1]
  
}  

new_dataset1<-new_dataset1[(new_dataset1$V4=='2BIO1' | new_dataset1$V4=='2BIO2'| new_dataset1$V4=='2BIO3'),]
new_dataset1<-new_dataset1[(new_dataset1$V5=='2BIO1' | new_dataset1$V5=='2BIO2'| new_dataset1$V5=='2BIO3'),]
new_dataset1<-new_dataset1[-which(is.na(new_dataset1$V4)),]
new_dataset1<-new_dataset1[-which(is.na(new_dataset1$V5)),]

write.csv(new_dataset1,"Datasets/new_dataset1_bio.csv", row.names=FALSE)
