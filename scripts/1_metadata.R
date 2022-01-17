library("dplyr")
library("tidyverse")
library(tidyr)
library(lubridate)
library(naniar)

#upload weather data
weather_2006 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2006_P1D.csv")
weather_2007 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2007_P1D.csv")
weather_2008 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2008_P1D.csv")
weather_2009 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2009_P1D.csv")
weather_2010 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2010_P1D.csv")
weather_2011 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2011_P1D.csv")
weather_2012 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2012_P1D.csv")
weather_2013 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2013_P1D.csv")
weather_2014 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2014_P1D.csv")
weather_2015 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2015_P1D.csv")
weather_2016 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2016_P1D.csv")
weather_2017 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2017_P1D.csv")
weather_2018 <- read.csv("data/raw data/fr_climat_quotidiennes_QC_7022579_2018_P1D.csv")

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




#### import viral ASV table ####
samples <- read.table('data/ASVs_counts_copy.tsv', header = T, sep = '\t', row.names = 1)

#create function to isolate date from the sampleID
get_date <- function(samp){
  date <- sub("*._*._*._*._*._*._*._","", samp) #remove sample ID at beginning
  date <- sub('^([^_]+_[^_]+_[^_]+).*', '\\1', date) #keep only dates (rm everything after third _)
  date <- as.Date(format(dmy(date), '%Y-%m-%d')) #format to date
  return(date)
}

dates <- get_date(colnames(samples))
str(dates)
(dates_unique <- unique(dates)) #remove duplicated dates


#function to select range of t-7:t (7 days leading up to date of interest)
get_date_range <- function(x){
  weather[weather$Date >= as.Date(x) - 7 & weather$Date <= as.Date(x),]
}
get_date_range(dates_unique[100]) #7 days leading up to and incl 2016-09-22

#function to get mean temp of 7 days leading to date
get_mean_temp <- function(x){
  y = get_date_range(x)
  return(mean(y$Temp.moy))
}
get_mean_temp(dates_unique[100])

#function to get cumulative precipitation from t-7:t
get_cumul_prec <- function(x){
  y = get_date_range(x)
  return(sum(y$Precip.tot))
}
get_cumul_prec(dates_unique[100])


#get mean temp for each sample date
# create new metadata 
meta <- as.data.frame(dates_unique)
colnames(meta) = c('Date') #change col name to Date
head(meta)

#add new columns to meta df 
for (i in 1:length(meta$Date)){
  meta$Mean_temp_t0_t7[i] <- get_mean_temp(meta$Date[i])
}

#get mean prec for each sample
for (i in 1:length(meta$Date)){
  meta$Cumul_precip_t1_t7_mm[i] <- get_cumul_prec(meta$Date[i])
}

head(meta)



# get all samples 
samp_meta <- as.data.frame(colnames(samples))
colnames(samp_meta) <- 'sampleID'
samp_meta$Date <- get_date(samp_meta$sampleID)

#create new columns
samp_meta$description <- sub("*._*._*._*._*._*._*._","", samp_meta$sampleID) #remove sample ID at beginning
samp_meta$description <- gsub("_", ".", samp_meta$description) #replace _ with .
samp_meta$Years <- format(as.Date(samp_meta$Date, format='%Y-%m-%d'), "%Y") #get year
samp_meta$Month <- format(as.Date(samp_meta$Date, format='%Y-%m-%d'), "%m") #get month

#change month from numeric to title
samp_meta$Month <- gsub("03", "March", samp_meta$Month)
samp_meta$Month <- gsub("04", "April",samp_meta$Month)
samp_meta$Month <- gsub("05", "May", samp_meta$Month)
samp_meta$Month <- gsub("06", "June", samp_meta$Month)
samp_meta$Month <- gsub("07", "July", samp_meta$Month)
samp_meta$Month <- gsub("08", "August", samp_meta$Month)
samp_meta$Month <- gsub("09", "September", samp_meta$Month)
samp_meta$Month <- gsub("10", "October", samp_meta$Month)

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

#assign season depending on date
samp_meta$Period <- getSeason(samp_meta$Date)

head(samp_meta)

env_data <- read.table('data/raw data/File_S1_Environmental_Table.txt')
colnames(env_data) <- c('SampleID', 'Julian', 'Week', 'Month', 'Year', 'Site', 'Season', 'Bloom', 'Tot.P_ug', 'Tot.N_mg', "Dissolved.P",
                        'Dissolved.N', 'Precipitation', 'Temperature', 'Microcystin', 'Description')
head(env_data)
env_data_filtered <- env_data[colnames(env_data) %in% c('Julian', 'Week', 'Site', 'Bloom', 'Tot.P_ug', 'Tot.N_mg', "Dissolved.P",
                                                        'Dissolved.N', 'Microcystin', 'Description')]


#explore the data
#make sure meta matches new meta samples
nrow(samp_meta)
nrow(env_data_filtered)

# env_data_filtered[which(env_data_filtered$Description %in% samp_meta$description), ] #see the rows from env_data_filtered that are in samp_meta
# intersect(env_data_filtered$Description, samp_meta$description) #which ones are the same
# length(setdiff(env_data_filtered$Description, samp_meta$description)) ##count how many are different
# sum(!is.na(env_data_filtered$Microcystin)) #how many non NA values in Microcystin


# compare metadata tables to see which information is accurate
weather_2006[which(weather_2006$Date == '2006-06-13'),] #ground truth weather information from specific date
env_data[which(env_data$SampleID == '13.06.2006.1c'), ]
meta[which(meta$Date == '2006-06-13'),]
get_date_range('2006-06-13') #7 days leading up to date of interest



#merge dataframes to get final metadata table
nrow(samp_meta)
nrow(meta)
meta_all <- merge(samp_meta, meta, by = 'Date', all.x=T)
nrow(meta_all)
nrow(env_data_filtered)
meta_all <- merge(meta_all, env_data_filtered, by.x = 'description', by.y = 'Description', all.x=T)
nrow(meta_all)

#clean up meta_all df
row.names(meta_all) <- meta_all$sampleID
meta_all <- meta_all[, !colnames(meta_all) == 'sampleID']

head(meta_all)
#write.csv(meta_all, 'data/metadata.csv')

