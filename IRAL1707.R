#### Code for Dowding et al in prep, Palaeontology and earth science database diversity dynamics.



# Database diversity and range (based on activity, last update, or if not available digitally,last citation)
install.packages(c("ggplot2","tidyverse", "dplyr", "divDyn","zoo", "Cairo"))


library(ggplot2)
library(tidyverse)
library(dplyr)
library(divDyn)
library(zoo)
library(Cairo)



#Set up working environment
setwd("xxxxxxxxxxxxx")
#Load in data from wd 
dat<-read.csv("IRAL.csv")

#alternative, read from GIT
```{r load_data}
path_dat <- "~/IRAL/Data"
#dat <- read.csv("https://github.com/dowdingem/main/IRAL/IRAL.csv")
dat <- read.csv(paste0(path_dat, "/IRAL.csv"))

# set up dataframe
df <- data.frame(
  Title = dat$Title,
  Start = dat$Inception,
  End = dat$Last.Ref.UPDATE,
  Omit = dat$Omit
)

# Clean Data

unique(df$Start)
unique(df$End)

clean.df <- na.omit(df)

clean.df <- clean.df[clean.df$Start != "",]

clean.df <- clean.df[clean.df$End != "",]

clean.df <- clean.df[clean.df$Start != "2006-7",]

clean.df <- clean.df[clean.df$End != "No record" & clean.df$End != "Absent" & 
                       clean.df$End != "absent" & clean.df$End != "2023/4?",]

clean.df <- clean.df[clean.df$Title != "",]

unique(clean.df$Start)
unique(clean.df$End)

df <- clean.df


# Set years to numeric

df$Start <- as.numeric(df$Start)
df$End <- as.numeric(df$End)

# Remove singletons 

df <- df[df$Start!=df$End,]

# Organize DBs based on start date

df <- df[order(df$Start),]

dateTitle <- factor(df$Title, levels = unique(df$Title[order(df$Start)]))


#############################################################################################################################

# Database Range Plot

ggplot(df, aes(y = dateTitle)) +
  geom_segment(aes(x = Start, xend = End, yend = dateTitle), color = "black", size = 2) +
  geom_point(aes(x = Start), color = "grey", size = 1) +   # Start point
  geom_point(aes(x = End), color = "grey", size = 1) +   # End point
  labs(title = "Database Range Plot",
       x = "Year",
       y = "DB active phase") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_blank(),  # Remove the panel background
    axis.line = element_line(color = "black"),  # Add axis lines for clarity
    axis.text.y = element_blank()
  )



#########################################################################################################

# Range through diversity plot

# Create stgs 
min_start <- min(df$Start)
max_end <- max(df$End)
stg <- seq(from = min_start, to = max_end, by = 2)

# Dataframe with Title duplicated for each stage
stg_df <- data.frame(Title = character(), stg = numeric(), stringsAsFactors = FALSE)

for (i in 1:nrow(df)) {
  title <- df$Title[i]
  start <- df$Start[i]
  end <- df$End[i]
  
  # Find all stg intervals the title is present in
  relevant_stg <- stg[stg >= start & stg <= end]
  
  # Create a dataframe for the current title and its relevant stg intervals
  temp_df <- data.frame(Title = rep(title, length(relevant_stg)), stg = relevant_stg)
  
  # Append to the main dataframe
  stg_df <- rbind(stg_df, temp_df)
}

# Convert stg intervals to ranks
stg_df$rank <- as.numeric(factor(stg_df$stg, levels = sort(unique(stg_df$stg))))


# basic plot using divDyn. 

dd <-divDyn(stg_df, bin="rank", tax="Title")

ggplot(data=dd, aes(x=rank, y=divRT, group=1)) +
  geom_line()+
  geom_point()


# create a stages like object? use stg_df
rank <- seq_along(stg)
stages <- data.frame(stg, rank)

dd$stg <- stg

#tsplot(stages, ylab="Range Through Richness", bottom = stages$stg)

png("RTDiv.png")
# base R plot
plot(dd$divRT, xlim = range(min(dd$stg), max(dd$stg)), xlab="Year", ylab="Number of Databases", main="Range Through Richness")
points(dd$stg, dd$divRT, pch=16)
lines(dd$stg, dd$divRT)
dev.off()

####################################################################################
#Extinction with or without rolling mean

#NEW DF for rank by year
### PER YEAR

# need to iterate years so can do extinction by which year of life db died in
# need to do one year per rank

stg_df$Years <- NA

i=0 # title
j=1 # row number
k=1 # number of years

for (i in 1:length(unique(stg_df$Title))) {
  
  
  j = 1
  k = 1
  
  current_title <- unique(stg_df$Title)[i]
  
  i = i + 1
  
  for (j in 1:nrow(stg_df)) {
    
    if (stg_df$Title[j] == current_title) stg_df$Years[j] <- k else k=0
    
    j = j + 1
    k = k + 1
    
  }
}

stg_df$Years

dd2 <-divDyn(stg_df, bin="Years", tax="Title")

#####EXTINCTION PC = iterate through by changing the selected column in "DD2" 

plot(dd2$extPC, xlim = range(min(dd2$Years), max(dd2$Years)), xlab="Number of Years of Database Activity", ylab="Rate", main="PC Extinction")
points(dd2$Years, dd2$extPC, pch=16)
lines(dd2$Years, dd2$extPC)


## EXTINCTION WITH ROLLING MEAN

# Extinction - PC w/ rolling mean =  iterate through by changing the selected column in "DD2" 
# temporary object
temp.zoo <- zoo(dd2$extPC, dd2$Years)
# rolling mean
roll.mean.extPC<-rollmean(temp.zoo, 3, fill = list(NA, NULL, NA))
# add to dataframe
dd2$roll.mean.extPC=coredata(roll.mean.extPC)

# plot PC
plot(dd2$roll.mean.extPC, xlim = range(min(dd2$Years), max(dd2$Years)), xlab="Number of Years of Database Activity", ylab="Rate", main="Extinction (PC) with Rolling Mean")
points(dd2$Years, dd2$roll.mean.extPC, pch=16)
lines(dd2$Years, dd2$roll.mean.extPC)

#div by years lived
### DIVERSITY

# SIB Diversity by # years DB lived

plot(dd2$divSIB, xlim = range(min(dd2$Years), max(dd2$Years)), xlab="Number of Years of Database Activity", ylab="Number of Databases", main="Diversity (SIB)")
points(dd2$Years, dd2$divSIB, pch=16)
lines(dd2$Years, dd2$divSIB)
#dev.off()


########################################################################################

#Four panel plot of origination and diversity without or without rolling means.

par(mfrow=c(2,2))

# Origination - PC 
# scale
dd$scale.oriPC <- scale(dd$oriPC)
dd$scale.divRT <- scale(dd$divRT)
# plot
plot(dd$scale.oriPC, xlim = range(min(dd$stg), max(dd$stg)), xlab="Years", ylab="Rate")
points(dd$stg, dd$scale.oriPC, pch=16)
lines(dd$stg, dd$scale.oriPC)
points(dd$stg, dd$scale.divRT, pch=16, col="blue")
lines(dd$stg, dd$scale.divRT, col="blue")
legend("topleft", bg="white", legend=c("Scaled Origination Rates (PC)", "Scaled Diversity (RT)"), 
       col=c("black", "blue"), lwd=3, inset=c(0.01,0.01), cex=0.75) 


# Origination - t.Ori 
# scale
dd$scale.tOri <- scale(dd$tOri)
dd$scale.divRT <- scale(dd$divRT)
# plot
plot(dd$scale.tOri, xlim = range(min(dd$stg), max(dd$stg)), xlab="Years", ylab="Rate")
points(dd$stg, dd$scale.tOri, pch=16)
lines(dd$stg, dd$scale.tOri)
points(dd$stg, dd$scale.divRT, pch=16, col="blue")
lines(dd$stg, dd$scale.divRT, col="blue")
legend("topleft", bg="white", legend=c("Scaled Origination Rates (tOri)", "Scaled Diversity (RT)"), 
       col=c("black", "blue"), lwd=3, inset=c(0.01,0.01), cex=0.75) 

# Origination with Rolling Mean - PC 
# temporary object
temp.zoo <- zoo(dd$oriPC, dd$stg)
# rolling mean
roll.mean.oriPC<-rollmean(temp.zoo, 3, fill = list(NA, NULL, NA))
# add to dataframe
dd$roll.mean.oriPC=coredata(roll.mean.oriPC)
# scale
dd$roll.mean.scale.oriPC <- scale(dd$roll.mean.oriPC)
dd$scale.divRT <- scale(dd$divRT)
# plot
plot(dd$roll.mean.scale.oriPC, xlim = range(min(dd$stg), max(dd$stg)), xlab="Years", ylab="Rate")
points(dd$stg, dd$roll.mean.scale.oriPC, pch=16)
lines(dd$stg, dd$roll.mean.scale.oriPC)
points(dd$stg, dd$scale.divRT, pch=16, col="blue")
lines(dd$stg, dd$scale.divRT, col="blue")
legend("topleft", bg="white", legend=c("Scaled Origination Rates (PC) with Rolling Mean", "Scaled Diversity (RT)"), 
       col=c("black", "blue"), lwd=3, inset=c(0.01,0.01), cex=0.75) 

# Origination with Rolling Mean - tOri 
# temporary object
temp.zoo <- zoo(dd$tOri, dd$stg)
# rolling mean
roll.mean.tOri<-rollmean(temp.zoo, 3, fill = list(NA, NULL, NA))
# add to dataframe
dd$roll.mean.tOri=coredata(roll.mean.tOri)
# scale
dd$roll.mean.scale.tOri <- scale(dd$roll.mean.tOri)
dd$scale.divRT <- scale(dd$divRT)

# plot
plot(dd$roll.mean.scale.tOri, xlim = range(min(dd$stg), max(dd$stg)), xlab="Years", ylab="Rate")
points(dd$stg, dd$roll.mean.scale.tOri, pch=16)
lines(dd$stg, dd$roll.mean.scale.tOri)
points(dd$stg, dd$scale.divRT, pch=16, col="blue")
lines(dd$stg, dd$scale.divRT, col="blue")
legend("topleft", bg="white", legend=c("Scaled Origination Rates (tOri) with Rolling Mean", "Scaled Diversity (RT)"), 
       col=c("black", "blue"), lwd=3, inset=c(0.01,0.01), cex=0.75)






