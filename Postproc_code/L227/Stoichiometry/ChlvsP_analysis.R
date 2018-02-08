
setwd("/Users/krsalkgu/Documents/SourceTree/Lake227/Postproc_code/L227/Stoichiometry")

L227chemistry <- read.csv("./IISDELA_L227_chemistry_cleaned_1969-2016.csv")
L227chemistry$Date <- as.Date(L227chemistry$Date, format = "%Y-%m-%d")
L227chemistry <- mutate(L227chemistry, Month = format(L227chemistry$Date, "%m"))
L227chemistry$Month <- as.numeric(L227chemistry$Month)

L227chemistry_epi_icefree <- 
  L227chemistry %>%
  filter(Stratum == "EPI") %>%
  filter(Month > 4 & Month < 11) %>%
  filter(Susp.P !=-1111 & Susp.P != -200 & Susp.P <100) %>%
  filter(chla !=-1111 & chla != -200 & chla <100)

 
PPvschl <- lm(L227chemistry_epi_icefree$chla ~ L227chemistry_epi_icefree$Susp.P)  
summary(PPvschl)

PPvschlbymonth <- lm(L227chemistry_epi_icefree$chla ~ L227chemistry_epi_icefree$Susp.P + L227chemistry_epi_icefree$Month)  
summary(PPvschlbymonth)

library(ggplot2)
ggplot(L227chemistry_epi_icefree, aes(x = Susp.P, y = chla, color = Month)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", color = "black")
