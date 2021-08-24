library("dplyr")
library("tidyverse")
library(tidyr)
library(lubridate)
library(naniar)

#upload weather data
weather_2006 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2006_P1D.csv")
weather_2007 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2007_P1D.csv")
weather_2008 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2008_P1D.csv")
weather_2009 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2009_P1D.csv")
weather_2010 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2010_P1D.csv")
weather_2011 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2011_P1D.csv")
weather_2012 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2012_P1D.csv")
weather_2013 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2013_P1D.csv")
weather_2014 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2014_P1D.csv")
weather_2015 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2015_P1D.csv")
weather_2016 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2016_P1D.csv")
weather_2017 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2017_P1D.csv")
weather_2018 <- read.csv("../data/raw data/fr_climat_quotidiennes_QC_7022579_2018_P1D.csv")

#merge all weather dataframes
weather <- Reduce(function(x,y) merge(x, y, all=TRUE), list(weather_2006, weather_2007, weather_2008,
                                                            weather_2009, weather_2010, weather_2011,
                                                            weather_2012, weather_2013, weather_2014,
                                                            weather_2015, weather_2016, weather_2017,
                                                            weather_2018))
#replace missing values with NA
weather <-replace_na(weather)

#format as.Date
weather$Date.Heure <- as.Date(weather$Date.Heure)

#rename cols & keep only what is useful to me
weather <- weather %>% rename(Date = Date.Heure, Temp.moy = Temp.moy...C., Precip.tot = Précip..tot...mm.) %>%
                        select(Date, Temp.moy, Precip.tot) 

#convert commas to periods & keep values as.numeric
weather$Temp.moy <- as.numeric(gsub(",", ".", weather$Temp.moy))

weather$Precip.tot <- as.numeric(gsub(",", ".", weather$Precip.tot))

# write.csv(weather, "weather.csv")

##get mean temp of 7 days leading to date
#function to select range of t-7:t
get_date_range <- function(x){
  weather[weather$Date >= as.Date(x) - 7 & weather$Date <= as.Date(x),]
}
get_date_range(samples$Date[100])

#function to get mean temp from t-7:t
get_mean_temp <- function(x){
  y = get_date_range(x)
  return(mean(y$Temp.moy))
}

#function to get mean precipitation from t-7:t
get_mean_prec <- function(x){
  y = get_date_range(x)
  return(mean(y$Precip.tot))
}

#get mean temp for each sample
for (i in 1:length(samples$Date)){
  samples$Mean_temperature_t0_t7[i] <- get_mean_temp(samples$Date[i])
}

#get mean prec for each sample
for (i in 1:length(samples$Date)){
  samples$Cumulative_precipitation_t1_t7_mm[i] <- get_mean_prec(samples$Date[i])
}

##########
meta2 <- read.csv("data/metadata2.csv") #already edited
#View(meta2)

#replace missing values with NA
meta2 <- meta2 %>% mutate_all(~replace(., . ==0, NA))

#rename
meta2<- rename(meta2, "Date" = "Sample")

#remove everything before 4th period in sampleID (just to keep date)
meta2$Date <- gsub("*........(.*)\\-.*","\\1", meta2$SampleID)

#fix dates that didn't correctly format
meta2$Description <- gsub("-Aug-17", "-08-2017", meta2$Description)
meta2$Description <- gsub("-Aug-18", "-08-2018", meta2$Description)
meta2$Description <- gsub("-Sep-17", "-09-2017", meta2$Description)
meta2$Description <- gsub("-Sep-18", "-09-2018", meta2$Description)
meta2$Description <- gsub("-Oct-18", "-10-2018", meta2$Description)

meta2$Date[200:210] <- meta2$Description[200:210]
meta2$Date[26] <- "02-06-2007"
meta2$Date[33] <- "29-08-2007"
meta2$Date[43] <- "01-06-2008"
meta2$Date[44] <- "02-06-2008"
meta2$Date[47] <- "02-07-2008"
meta2$Date[50] <- "28-07-2008"
meta2$Date[54] <- "18-08-2008"
meta2$Date[70] <- "19-06-2009"
meta2$Date[71] <- "19-06-2009"
meta2$Date[79] <- "29-07-2009"
meta2$Date[80] <- "08-08-2009"
meta2$Date[110] <- "12-06-2011"
meta2$Date[143] <- "27-08-2013"
meta2$Date[158] <- "29-07-2015"
meta2$Date[162] <- "31-07-2015"
meta2$Date[165] <- "05-08-2015"
meta2$Date[167] <- "06-08-2015"
meta2$Date[180] <- "19-07-2016"
meta2$Date[181] <- "20-07-2016"
meta2$Date[182] <- "21-07-2016"
meta2$Date[183] <- "25-07-2016"
meta2$Date[184] <- "27-07-2016"
meta2$Date[185] <- "01-08-2016"
meta2$Date[186] <- "03-08-2016"
meta2$Date[187] <- "09-08-2016"
meta2$Date[188] <- "16-08-2016"
meta2$Date[189] <- "18-08-2016"
meta2$Date[190] <- "23-08-2016"
meta2$Date[191] <- "01-09-2016"
meta2$Date[192] <- "03-09-2016"
meta2$Date[193] <- "08-09-2016"
meta2$Date[194] <- "15-09-2016"
meta2$Date[195] <- "19-09-2016"
meta2$Date[196] <- "22-09-2016"
meta2$Date[197] <- "25-10-2016"
meta2$Date[198] <- "26-09-2016"
meta2$Date[199] <- "06-10-2016"

#extract proper year
meta2$Years <- gsub(".*-(.*)\\-.*", "\\1", meta2$SampleID, perl=T)

#extract proper month
#meta2$Months <- gsub(".*-.(.*)\\-.*-.*", "\\1", meta2$SampleID, perl=T) #"." before ( indicates to remove the first character between xx-Xx-xxxx
meta2$Months <- gsub(".*-(.*)\\-.*-.*", "\\1", meta2$SampleID, perl=T)

#change month from numeric to title
meta2$Months <- gsub("03", "March", meta2$Months)
meta2$Months <- gsub("04", "April", meta2$Months)
meta2$Months <- gsub("05", "May", meta2$Months)
meta2$Months <- gsub("06", "June", meta2$Months)
meta2$Months <- gsub("07", "July", meta2$Months)
meta2$Months <- gsub("08", "August", meta2$Months)
meta2$Months <- gsub("09", "September", meta2$Months)
meta2$Months <- gsub("10", "October", meta2$Months)

#assign season depending on date
getSeason <- function(DATES) {
  WS <- as.Date("15-12-2012", format = "%d-%m-%Y") # Winter Solstice
  SE <- as.Date("15-3-2012",  format = "%d-%m-%Y") # Spring Equinox
  SS <- as.Date("15-6-2012",  format = "%d-%m-%Y") # Summer Solstice
  FE <- as.Date("15-9-2012",  format = "%d-%m-%Y") # Fall Equinox
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Fall")))
}
#format as.Date
meta2$Date <- as.Date(meta2$Date, "%d-%m-%Y")

#assign season depending on date
meta2$Period <- getSeason(meta2$Date)

#get avg temp for each sample
for (i in 1:length(meta2$Date)){
  meta2$Mean_temperature_t0_t7[i] <- get_mean_temp(meta2$Date[i])
}

#get avg prec for each sample
for (i in 1:length(meta2$Date)){
  meta2$Cumulative_precipitation_t1_t7_mm[i] <- get_mean_prec(meta2$Date[i])
}

write.csv(meta2, "metadata2.csv")





##### merge metatable (from metadata_w_cmd.csv) -- add microcystin data ######
env_tab <- read.table("../data/raw data/File_S1_Environmental_Table.txt")
head(env_tab)
head(meta)

colnames(env_tab) <- c("SampleID", "Julia", "Week", "Months", "Years", "Site", "Season", "Bloom", "TP", "TN", "DP", 
                       "DN", "Precipitation", "Temperature", "Microcystin", "Description")

#match sample dates
meta2 <- meta

meta2$ID <- rownames(meta2)

table(rownames(meta2)) #check for duplicates -- won't be cause the sample ID is different even tho date is the same
rownames(meta2)[rownames(meta2) == "FLD0295_15_05_2011_1"] <- "FLD0295_15_05_2011_2" #dates were duplicated therefore need to correct

#remove sample ID at beginning
meta2$description <- sub("*._*._*._*._*._*._*._","", row.names(meta2))
#change "_" to "."
meta2$description <- gsub("_", ".", meta2$description)

#make sure meta matches new meta samples
nrow(meta2)
nrow(env_tab)

env_tab$Description %in% meta2$description
intersect(env_tab$Description, meta2$description) #which ones are the same
length(setdiff(env_tab$Description, meta2$description)) ##count how many are different
 
sum(!is.na(env_tab$Microcystin)) #how many non NA values in Microcystin

#merge env_tab Microcystin based on Description
metamerge <- merge(meta2, env_tab[,c("Description", "Microcystin")], by.x="description", by.y = "Description", all.x = T)
rownames(metamerge) <- metamerge$ID

#rm a cols
drop <- c("ID")
metamerge <- metamerge[, !(names(metamerge) %in% drop)]

head(metamerge, n=2)
metamerge <- metamerge %>% select(description, everything()) #move Description col to first position in df
write.csv(metamerge, "Metadata.csv")
