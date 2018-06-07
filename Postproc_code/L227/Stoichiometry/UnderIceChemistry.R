
library(tidyverse)

L227chemistry <- read.csv("./Stoichiometry/IISDELA_L227_chemistry_cleaned_1969-2016.csv")
L227chemistry <- select(L227chemistry, Date, Depth, Stratum:TDP)
L227chemistry$Date <- as.Date(L227chemistry$Date, "%Y-%m-%d")
L227chemistry <- L227chemistry %>%
  mutate(Month = as.numeric(format(L227chemistry$Date, "%m"))) %>%
  mutate(Year = as.numeric(format(L227chemistry$Date, "%Y"))) %>%
  filter(Month == 1 | Month == 2 | Month == 3) %>%
  filter(Depth == 0 | Depth == 1)

L227chemistry$Susp.N <- ifelse(L227chemistry$Susp.N == -1111, NA, L227chemistry$Susp.N)
L227chemistry$TDN <- ifelse(L227chemistry$TDN == -200, NA, L227chemistry$TDN)
L227chemistry$Susp.P <- ifelse(L227chemistry$Susp.P == -1111, NA, L227chemistry$Susp.P)
L227chemistry$TDP <- ifelse(L227chemistry$TDP == -200, NA, L227chemistry$TDP)

CumulativePP <- end.of.season.cumulative[c(3, 6)]
L227chemistrycomparisons <- merge(L227chemistry, CumulativePP, by = "Year")
L227chemistrycomparisons2 <- merge(L227chemistrycomparisons, ZoopandStoicMeans, by = "Year")

ggplot(L227chemistry) + 
  geom_point(aes(x = Date, y = NO3))

ggplot(L227chemistry) + 
  geom_point(aes(x = Date, y = NH4))

ggplot(L227chemistry) + 
  geom_point(aes(x = Date, y = Susp.N))

ggplot(L227chemistry) + 
  geom_point(aes(x = Date, y = TDN))

ggplot(L227chemistry) + 
  geom_point(aes(x = Date, y = Susp.P))

ggplot(L227chemistry) + 
  geom_point(aes(x = Date, y = TDP))

ggplot(L227chemistrycomparisons) +
  geom_point(aes(x = NO3, y = Cum.Sum.obs.PP))

ggplot(L227chemistrycomparisons) +
  geom_point(aes(x = NH4, y = Cum.Sum.obs.PP))

ggplot(L227chemistrycomparisons) +
  geom_point(aes(x = Susp.N, y = Cum.Sum.obs.PP))

ggplot(L227chemistrycomparisons) +
  geom_point(aes(x = TDN, y = Cum.Sum.obs.PP))

ggplot(L227chemistrycomparisons) +
  geom_point(aes(x = Susp.P, y = Cum.Sum.obs.PP))

ggplot(L227chemistrycomparisons) +
  geom_point(aes(x = TDP, y = Cum.Sum.obs.PP))

# only significant correlations are displayed
summary(lm(L227chemistrycomparisons$Cum.Sum.obs.PP ~ L227chemistrycomparisons$NO3))

summary(lm(L227chemistrycomparisons$Cum.Sum.obs.PP ~ L227chemistrycomparisons$NH4))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  4226.3680   550.2624   7.681 3.64e-11 ***
#   L227chemistrycomparisons$NH4    2.1204     0.7959   2.664  0.00935 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2113 on 79 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.08244,	Adjusted R-squared:  0.07082 
# F-statistic: 7.098 on 1 and 79 DF,  p-value: 0.009353

summary(lm(L227chemistrycomparisons$Cum.Sum.obs.PP ~ L227chemistrycomparisons$Susp.N))
summary(lm(L227chemistrycomparisons$Cum.Sum.obs.PP ~ L227chemistrycomparisons$TDN))
summary(lm(L227chemistrycomparisons$Cum.Sum.obs.PP ~ L227chemistrycomparisons$Susp.P))
summary(lm(L227chemistrycomparisons$Cum.Sum.obs.PP ~ L227chemistrycomparisons$TDP))


summary(lm(L227chemistrycomparisons2$PP ~ L227chemistrycomparisons2$NH4))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   24.658979   3.332655   7.399 1.65e-10 ***
#   L227chemistrycomparisons2$NH4  0.014739   0.004758   3.098  0.00274 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 12.55 on 75 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.1135,	Adjusted R-squared:  0.1016 
# F-statistic: 9.598 on 1 and 75 DF,  p-value: 0.002741

summary(lm(L227chemistrycomparisons2$Chl ~ L227chemistrycomparisons2$NH4))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   39.788833   3.547742  11.215  < 2e-16 ***
#   L227chemistrycomparisons2$NH4 -0.015216   0.005065  -3.004  0.00362 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 13.36 on 75 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.1074,	Adjusted R-squared:  0.09552 
# F-statistic: 9.026 on 1 and 75 DF,  p-value: 0.003617

summary(lm(L227chemistrycomparisons2$PP ~ L227chemistrycomparisons2$NO3))
summary(lm(L227chemistrycomparisons2$Chl ~ L227chemistrycomparisons2$NO3))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   27.28578    2.07840  13.128   <2e-16 ***
#   L227chemistrycomparisons2$NO3  0.02729    0.01323   2.063   0.0425 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 13.67 on 76 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.05304,	Adjusted R-squared:  0.04058 
# F-statistic: 4.257 on 1 and 76 DF,  p-value: 0.04251

summary(lm(L227chemistrycomparisons2$PP ~ L227chemistrycomparisons2$Susp.N))
summary(lm(L227chemistrycomparisons2$Chl ~ L227chemistrycomparisons2$Susp.N))
summary(lm(L227chemistrycomparisons2$PP ~ L227chemistrycomparisons2$TDN))
summary(lm(L227chemistrycomparisons2$Chl ~ L227chemistrycomparisons2$TDN))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   47.001861   5.152658   9.122  9.8e-14 ***
#   L227chemistrycomparisons2$TDN -0.013695   0.004023  -3.404  0.00108 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 13.23 on 74 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.1354,	Adjusted R-squared:  0.1237 
# F-statistic: 11.59 on 1 and 74 DF,  p-value: 0.001076

summary(lm(L227chemistrycomparisons2$PP ~ L227chemistrycomparisons2$Susp.P))
summary(lm(L227chemistrycomparisons2$Chl ~ L227chemistrycomparisons2$Susp.P))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                       37.5163     3.0249  12.402  < 2e-16 ***
#   L227chemistrycomparisons2$Susp.P  -0.3040     0.1085  -2.802  0.00649 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 13.52 on 74 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.09589,	Adjusted R-squared:  0.08367 
# F-statistic: 7.849 on 1 and 74 DF,  p-value: 0.006487

summary(lm(L227chemistrycomparisons2$PP ~ L227chemistrycomparisons2$TDP))
summary(lm(L227chemistrycomparisons2$Chl ~ L227chemistrycomparisons2$TDP))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    38.9198     4.1751   9.322 4.11e-14 ***
#   L227chemistrycomparisons2$TDP  -1.0570     0.4708  -2.245   0.0277 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 13.76 on 74 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.06378,	Adjusted R-squared:  0.05113 
# F-statistic: 5.042 on 1 and 74 DF,  p-value: 0.02773
