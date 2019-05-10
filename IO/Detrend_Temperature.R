
library(tidyverse)
library(lubridate)
#library(fpp)
#library(forecast)
Inputs <- read.table("./L227_input_interpolated_constantfertilization.txt", skip = 1, header = T)

Temp <- Inputs %>%
  select(AirTemperature) 

Temp_ts <- ts(Temp[[1]], frequency = 365)
Temp_Decomposed <- stl(Temp_ts, "periodic")
Temp_Seasonal <- as.data.frame(Temp_Decomposed$time.series[,1])
Temp_Trend <- as.data.frame(Temp_Decomposed$time.series[,2])
Temp_Random <- as.data.frame(Temp_Decomposed$time.series[,3])

plot(Temp_Decomposed)

Temp_Seasonal<- mutate(Temp_Seasonal, 
                       Date = as_date(paste(Inputs$Year, Inputs$Month, Inputs$Day, sep = "-")))
Temp_Trend <- mutate(Temp_Trend, 
                       Date = as_date(paste(Inputs$Year, Inputs$Month, Inputs$Day, sep = "-")))
Temp_Random <- mutate(Temp_Random, 
                       Date = as_date(paste(Inputs$Year, Inputs$Month, Inputs$Day, sep = "-")))
Temp_Observed <- mutate(Temp, 
                        Date = as_date(paste(Inputs$Year, Inputs$Month, Inputs$Day, sep = "-")))

colnames(Temp_Seasonal) <- c("Temp.Seasonal", "Date")
colnames(Temp_Trend) <- c("Temp.Trend", "Date")
colnames(Temp_Random) <- c("Temp.Random", "Date")
colnames(Temp_Observed) <- c("Temp.Observed", "Date")

Temp_All <- 
  left_join(Temp_Observed, Temp_Seasonal, by = "Date") %>%
  left_join(., Temp_Trend, by = "Date") %>%
  left_join(., Temp_Random, by = "Date") 

Temp_All <- mutate(Temp_All, Temp.Detrended = Temp.Observed - Temp.Trend)
Temp_All_Period3 <- Temp_All[7495:17355,]
summary(Temp_All_Period3$Temp.Trend)
sd(Temp_All_Period3$Temp.Trend)


ggplot(Temp_All) +
  geom_line(aes(y = Temp.Observed, x = Date), color = "#240c4cff", size = 0.25) +
  geom_line(aes(y = Temp.Detrended, x = Date), color = "#f99d15ff", size = 0.25, alpha = 0.6) +
  geom_line(aes(y = Temp.Trend, x = Date), color = "black") +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Temperature") +
  theme_classic()

Temp.Period3 <-
ggplot(Temp_All[7494:17355,]) +
  geom_line(aes(y = Temp.Observed, x = Date), color = "black", size = 0.3) +
  geom_line(aes(y = Temp.Detrended, x = Date), color = "#f99d15ff", size = 0.3, alpha = 0.6) +
  geom_line(aes(y = Temp.Trend, x = Date), color = "#d14a42ff") +
  geom_hline(yintercept = 0, lty = 2) +
  ylab(expression("Air Temperature " ( degree*C))) +
  theme_classic()
print(Temp.Period3)

ggsave("DecomposedTemp.Period3.pdf", Temp.Period3, dpi = 300, width = 6.5, height = 3, units = "in")


write.csv(file = "Temperature_Trends.csv", Temp_All, row.names = FALSE)
