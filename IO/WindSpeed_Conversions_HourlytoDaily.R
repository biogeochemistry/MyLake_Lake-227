library(plyr)
Wind <- read.csv("WindSpeed_hourly_1969-2016.csv")
head(Wind)

winddaily <- ddply(Wind, .(measurement.date), summarize, dailyavgwind = mean(speed.at.10m))
winddaily$measurement.date <- as.Date(winddaily$measurement.date, format = "%m/%d/%y")
winddaily <- arrange(winddaily, measurement.date)

#this is the overland speed at 10 m height. We need to convert to the over-water speed at 10 m height. 
# The conversion can be done with a ratio, which for Lake 470 (similar size, Solinske 1982) is 0.527.

winddaily <- mutate(winddaily, windoverwater = dailyavgwind * 0.527)

# Convert to m s-1 (files are in km h-1)

winddaily$dailyavgwind <- winddaily$dailyavgwind * 1000 / 3600
winddaily$windoverwater <- winddaily$windoverwater * 1000 / 3600

winddaily$dailyavgwind <- round(winddaily$dailyavgwind, digits = 1)
winddaily$windoverwater <- round(winddaily$windoverwater, digits = 1)

write.csv(winddaily, file = "WindDaily_1969to2016.csv", row.names = F)

