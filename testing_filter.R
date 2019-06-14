#make fake data set that is structured similarly to what I think we'll have
library(dplyr)
library(tidyr)

subject_id<-c(1,2, 2, 1, 2, 3, 3)
choice<-c(1,1, 2, 2, 3, 1, 2)
classification<-c("moose", "coyote", "coyote", "moose", "wolf", "wtd", "wtd")

df<-as_tibble(cbind(subject_id, classification))
df1<-as_data_frame(cbind(subject_id, choice, classification))
df1$subject_id<-as.factor(df1$subject_id)
df1$choice<- as.factor(df1$choice)
					 #now I need a loop that will go through df1, group a subset of data according to sharing the same subject_id,
#and widen to a format where there is only 1 row per subject id with as many columns as there are classfications

#start by making second data frame to accept the data
numRows<-length(unique(subject_id))
df2<-tibble(subject_id = numeric(numRows),
			 sp1 = character(numRows),
			 sp2 = character(numRows),
			 sp3 = character(numRows))
#now build out the process, but not as a loop; then "loopify it"
library(dplyr)
#separate events
class_event<-filter(df1, subject_id == 1)
class_event2<-filter(df1, subject_id == 2)
class_event3<-filter(df1, subject_id ==3)
#go to first event and assign species
df2$subject_id[1]<-1
df2$sp1[1]<-class_event$classification[1]
df2$sp2[1]<- class_event$classification[2]

#now go to second event and assign species

df3<-df1
df3$subject_id<-as.factor(df3$subject_id)
library(tidyr)
df4<-spread(df1, choice, classification)


#so, the take home message is that if I can add a column that identifies each choice for each subject
#as a factor, I can spread the data with spread (tidyr)

#now let's think about how to add choice numbers to a df
#begin by sorting by subject id
df[order(df$subject_id),]
unique(df$subject_id)


library(dplyr)
df2<-df %>% group_by(subject_id) %>% summarize(
	numObs= length((classification))
)


ids<-unique(df$subject_id)


#hint from Ivan (SOLUTION from ivan.)

library(dplyr)
df<-df %>%
	group_by(subject_id) %>%
	mutate(choice = letters[1:length(subject_id)])

library(tidyr)
spread(df, choice, classification)
