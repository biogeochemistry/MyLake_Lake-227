
NtoPLake <- read.csv("NP_Stoichiometry_L227.csv")
NtoPLake$Date <- as.Date(NtoPLake$Date, "%m/%d/%y")
library(tidyverse)

Zoop <- read.csv("Zooplankton_Biomass.csv")
Zoop$Date <- as.Date(Zoop$Date, "%m/%d/%y")

ZoopandStoic <- merge(NtoPLake, Zoop, by = "Date")
ZoopandStoic <- select(ZoopandStoic, Date, PP, Chl:Total.biomass)
ZoopandStoic <- na.omit(ZoopandStoic)
ZoopandStoic <- mutate(ZoopandStoic, Year = format(ZoopandStoic$Date, "%Y"), Month = format(ZoopandStoic$Date, "%m"))
ZoopandStoic <- filter(ZoopandStoic, Month == "05" | Month == "06" | Month == "07" | Month == "08" | Month == "09" | Month == "10")

ZoopandStoicMeans <- aggregate(ZoopandStoic, list(Year = ZoopandStoic$Year), mean)
ZoopandStoicMeans <- ZoopandStoicMeans[c(1, 3, 4, 11)]
CumulativePP <- end.of.season.cumulative[c(3, 6)]

ZoopandStoicMeans2 <- merge(ZoopandStoicMeans, CumulativePP, by = "Year")


ggplot(ZoopandStoic) + 
  geom_point(aes(x = Total.biomass, y = PP))

ggplot(ZoopandStoic) + 
  geom_point(aes(x = Total.biomass, y = Chl))

ggplot(ZoopandStoic) + 
  geom_point(aes(x = log(Total.biomass), y = log(PP)))

ggplot(ZoopandStoic) + 
  geom_point(aes(x = log(Total.biomass), y = log(Chl)))

ggplot(ZoopandStoicMeans2) + 
  geom_point(aes(x = Total.biomass, y = PP))

ggplot(ZoopandStoicMeans2) + 
  geom_point(aes(x = Total.biomass, y = Chl))

ggplot(ZoopandStoicMeans2) + 
  geom_point(aes(x = Total.biomass, y = Cum.Sum.obs.PP))

ggplot(ZoopandStoicMeans2) + 
  geom_point(aes(x = log(Total.biomass), y = log(Cum.Sum.obs.PP)))

summary(lm(ZoopandStoic$PP ~ ZoopandStoic$Total.biomass))
summary(lm(ZoopandStoic$Chl ~ ZoopandStoic$Total.biomass))
summary(lm(log(ZoopandStoic$PP) ~ log(ZoopandStoic$Total.biomass)))
summary(lm(log(ZoopandStoic$Chl) ~ log(ZoopandStoic$Total.biomass)))

summary(lm(ZoopandStoicMeans2$PP ~ ZoopandStoicMeans2$Total.biomass))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      33.56675    2.55491  13.138 2.76e-16 ***
#   ZoopandStoicMeans2$Total.biomass -0.13478    0.08356  -1.613    0.114    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 11.31 on 41 degrees of freedom
# Multiple R-squared:  0.05967,	Adjusted R-squared:  0.03674 
# F-statistic: 2.602 on 1 and 41 DF,  p-value: 0.1144

summary(lm(ZoopandStoicMeans2$Chl ~ ZoopandStoicMeans2$Total.biomass))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                       29.1544     2.7948   10.43 4.18e-13 ***
#   ZoopandStoicMeans2$Total.biomass  -0.0850     0.0914   -0.93    0.358    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 12.37 on 41 degrees of freedom
# Multiple R-squared:  0.02066,	Adjusted R-squared:  -0.00323 
# F-statistic: 0.8648 on 1 and 41 DF,  p-value: 0.3579

summary(lm(ZoopandStoicMeans2$Cum.Sum.obs.PP ~ ZoopandStoicMeans2$Total.biomass))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                       5666.67     418.42  13.543   <2e-16 ***
#   ZoopandStoicMeans2$Total.biomass   -25.83      13.68  -1.888   0.0662 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1852 on 41 degrees of freedom
# Multiple R-squared:  0.07996,	Adjusted R-squared:  0.05752 
# F-statistic: 3.563 on 1 and 41 DF,  p-value: 0.06616

summary(lm(log(ZoopandStoicMeans2$PP) ~ log(ZoopandStoicMeans2$Total.biomass)))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            3.65227    0.19037  19.185   <2e-16 ***
#   log(ZoopandStoicMeans2$Total.biomass) -0.10154    0.06521  -1.557    0.127    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.307 on 41 degrees of freedom
# Multiple R-squared:  0.05585,	Adjusted R-squared:  0.03282 
# F-statistic: 2.425 on 1 and 41 DF,  p-value: 0.1271

summary(lm(log(ZoopandStoicMeans2$Chl) ~ log(ZoopandStoicMeans2$Total.biomass)))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            3.45147    0.26691  12.931 4.68e-16 ***
#   log(ZoopandStoicMeans2$Total.biomass) -0.08420    0.09142  -0.921    0.362    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4305 on 41 degrees of freedom
# Multiple R-squared:  0.02027,	Adjusted R-squared:  -0.003626 
# F-statistic: 0.8483 on 1 and 41 DF,  p-value: 0.3624

summary(lm(log(ZoopandStoicMeans2$Cum.Sum.obs.PP) ~ log(ZoopandStoicMeans2$Total.biomass)))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            8.84632    0.18606  47.545   <2e-16 ***
#   log(ZoopandStoicMeans2$Total.biomass) -0.12923    0.06373  -2.028   0.0491 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3001 on 41 degrees of freedom
# Multiple R-squared:  0.09115,	Adjusted R-squared:  0.06899 
# F-statistic: 4.112 on 1 and 41 DF,  p-value: 0.04911


