
library(plyr)
humidity <- read.csv("Humidity_Hourly_2010to2016.csv")
head(humidity)

humiditydaily <- ddply(humidity, .(measurement.date), summarize, dailyavghumidity = mean(relative.humidity))
humiditydaily$measurement.date <- as.Date(humiditydaily$measurement.date, format = "%m/%d/%y")
humiditydaily <- arrange(humiditydaily, measurement.date)
humiditydaily$dailyavghumidity <- round(humiditydaily$dailyavghumidity, digits = 1)

write.csv(humiditydaily, file = "Humidity_Daily_2010to2016.csv", row.names = F)

