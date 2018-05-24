setwd("/Users/krsalkgu/Documents/SourceTree/Lake227/Postproc_code/L227/Stoichiometry")
library(tidyverse)

NtoPLake <- read.csv("NP_Stoichiometry_L227.csv")
attach(NtoPLake)
head(NtoPLake)
NtoPInflow <- read.csv("NP_Stoichiometry_Inflow.csv")
attach(NtoPInflow)
head(NtoPInflow)
CyanoPercent <- read.csv("../Phytoplankton/L227_cyanopercent.csv")

#Change concentrations from mass (ug/L or mg/d) to molar (umol/L or mmol/d)
Fert_TP_molar <- Fert_TP_mg.d/30.97
Inflow_TP_molar <- Inflow_TP_mg.d/30.97
Fert_TN_molar <- Fert_TN_mg.d/14.01
Inflow_TN_molar <- Inflow_TN_mg.d/14.01
NO3_molar <- NO3/14.01
NH4_molar <- NH4/14.01
DIN_molar <- NO3_molar + NH4_molar
PN_molar <- PN/14.01
TDN_molar <- TDN/14.01
PP_molar <- PP/30.97
TDP_molar <- TDP/30.97
TN_molar <- TDN_molar + PN_molar
TP_molar <- TDP_molar + PP_molar

#Combine inflows + fertilization
Input_TN_molar <- Fert_TN_molar + Inflow_TN_molar
Input_TP_molar <- Fert_TP_molar + Inflow_TP_molar

#Create N:P stoichiometric ratios in lake and inputs
Fert_NtoP <- Fert_TN_molar/Fert_TP_molar
Inflow_NtoP <- Inflow_TN_molar/Inflow_TP_molar
Input_NtoP <- Input_TN_molar/Input_TP_molar
DINtoTDP <- DIN_molar/TDP_molar
TDNtoTDP <- TDN_molar/TDP_molar
PNtoPP <- PN_molar/PP_molar
TNtoTP <- TN_molar/TP_molar

#Convert dates from factor to date format, add month 
Datelake <- as.Date(NtoPLake$Date, "%m/%d/%y")
Dateinflow <- as.Date(NtoPInflow$Date, "%m/%d/%y")
Monthlake <- as.numeric(format(Datelake, "%m"))
Monthinflow <- as.numeric(format(Dateinflow, "%m"))
CyanoPercent$Date <- as.Date(CyanoPercent$Date, "%m/%d/%y")

#Create data frames for N:P stoichiometry in lake and inflows for May-October
NtoPinsitu <- data.frame(Datelake, Monthlake, TNtoTP, PNtoPP, DINtoTDP, TDNtoTDP, DIN_molar, PP)
NtoPinsitu <- filter(NtoPinsitu, Monthlake > 4 & Monthlake < 11)
NtoPinput <- data.frame(Dateinflow, Monthinflow, Fert_NtoP, Inflow_NtoP, Input_NtoP, Input_TN_molar)
NtoPinput <- filter(NtoPinput, Monthinflow > 4 & Monthinflow < 11)

# Set up analysis basics
library(ggplot2)
theme_std <- function (base_size = 16, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.ticks = element_line(colour = "black", size = 1), 
          legend.key = element_rect(colour = "white"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = NA),
          axis.line = element_line(size = 0.5, colour = "black"),
          panel.grid.major = element_line(NA), 
          panel.grid.minor = element_line(NA), 
          strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
          axis.text  = element_text(size=rel(0.9)),
          axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm"),size=rel(1)),
          axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"),size=rel(1), angle = 90),
          strip.text = element_text(size = rel(1.15), colour = "black", face = "bold"),
          plot.margin=unit(c(10,10,10,10),"pt"))}
theme_set(theme_std())

library(strucchange)

#### cyanopercent ----
shapiro.test(CyanoPercent$Nfixcyano.percent)
    # Shapiro-Wilk normality test
    # 
    # data:  CyanoPercent$Nfixcyano.percent
    # W = 0.88162, p-value < 2.2e-16

#### PP ----
shapiro.test(NtoPinsitu$PP)
    # Shapiro-Wilk normality test
    # 
    # data:  NtoPinsitu$PP_molar
    # W = 0.74035, p-value < 2.2e-16

#### TN:TP ----
TNtoTPbreak = breakpoints (TNtoTP ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
TNtoTPbreak = breakpoints (log(TNtoTP) ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
summary(TNtoTPbreak)
  # Breakpoints at observation number:
  # m = 1   470

TNtoTPbreak2 = breakpoints (TNtoTP ~ Datelake, h = 20, data = NtoPinsitu, breaks=2)
summary(TNtoTPbreak2)
  # Breakpoints at observation number:
  # m = 1       470
  # m = 2   112 132

shapiro.test(NtoPinsitu$TNtoTP)
    # Shapiro-Wilk normality test
    # 
    # data:  NtoPinsitu$TNtoTP
    # W = 0.83693, p-value < 2.2e-16
qqnorm(NtoPinsitu$TNtoTP)

TNtoTPbydate <- lm(NtoPinsitu$TNtoTP ~ NtoPinsitu$Datelake)
summary(TNtoTPbydate)
  # Coefficients:
  #   Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)         61.1493052  1.6310814  37.490  < 2e-16 ***
  #   NtoPinsitu$Datelake -0.0005242  0.0001748  -2.998  0.00281 ** 
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # Residual standard error: 23.58 on 750 degrees of freedom
  # (220 observations deleted due to missingness)
  # Multiple R-squared:  0.01184,	Adjusted R-squared:  0.01052 
  # F-statistic: 8.987 on 1 and 750 DF,  p-value: 0.002808


TNtoTPperiod1 <- TNtoTPbydate <- lm(NtoPinsitu$TNtoTP[1:469] ~ NtoPinsitu$Datelake[1:469])
summary(TNtoTPperiod1)
    # Coefficients:
    # Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                5.813e+01  2.908e+00  19.993   <2e-16 ***
    #   NtoPinsitu$Datelake[1:469] 1.846e-04  7.457e-04   0.248    0.805    
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 27.19 on 329 degrees of freedom
    # (138 observations deleted due to missingness)
    # Multiple R-squared:  0.0001862,	Adjusted R-squared:  -0.002853 
    # F-statistic: 0.06129 on 1 and 329 DF,  p-value: 0.8046

TNtoTPperiod2 <- TNtoTPbydate <- lm(NtoPinsitu$TNtoTP[470:972] ~ NtoPinsitu$Datelake[470:972])
summary(TNtoTPperiod2)
    # Coefficients:
    # Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                  68.6466462  3.6337370  18.891  < 2e-16 ***
    #   NtoPinsitu$Datelake[470:972] -0.0011304  0.0003034  -3.726 0.000221 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 20.21 on 419 degrees of freedom
    # (82 observations deleted due to missingness)
    # Multiple R-squared:  0.03208,	Adjusted R-squared:  0.02977 
    # F-statistic: 13.89 on 1 and 419 DF,  p-value: 0.0002209

TNtoTPplot <- 
  ggplot(NtoPinsitu, aes(x = Datelake, y = log(TNtoTP))) +
  geom_point(aes(y = TNtoTP), size = 0.5) +
  #ylim(0,200) +
  ylab(expression(TN:TP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) + 
  geom_smooth(method = 'lm', data = NtoPinsitu[1:469,], aes(x = Datelake, y = TNtoTP), se = FALSE, color = "black") + #non-significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[470:972,], aes(x = Datelake, y = TNtoTP), se = FALSE, color = "blue") + #significant slope
  geom_vline(xintercept = as.numeric(NtoPinsitu$Datelake[470]), lty = 5)

#### TDN:TDP ----
TDNtoTDPbreak = breakpoints (TDNtoTDP ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
TDNtoTDPbreak = breakpoints (log(TDNtoTDP) ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
summary(TDNtoTDPbreak)
    # Breakpoints at observation number:
    #   m = 1   558

shapiro.test(NtoPinsitu$TDNtoTDP)
    # Shapiro-Wilk normality test
    # 
    # data:  NtoPinsitu$TDNtoTDP
    # W = 0.21698, p-value < 2.2e-16
qqnorm(NtoPinsitu$TDNtoTDP)

TDNtoTDPbydate <- lm(NtoPinsitu$TDNtoTDP ~ NtoPinsitu$Datelake)
summary(TDNtoTDPbydate)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)         1.349e+02  1.421e+01   9.492   <2e-16 ***
    #   NtoPinsitu$Datelake 1.692e-03  1.532e-03   1.104     0.27    
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 211.2 on 779 degrees of freedom
    # (191 observations deleted due to missingness)
    # Multiple R-squared:  0.001563,	Adjusted R-squared:  0.0002811 
    # F-statistic: 1.219 on 1 and 779 DF,  p-value: 0.2698

TDNtoTDPperiod1 <- lm(NtoPinsitu$TDNtoTDP[1:557] ~ NtoPinsitu$Datelake[1:557])
summary(TDNtoTDPperiod1)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                1.146e+02  7.262e+00  15.783  < 2e-16 ***
    #   NtoPinsitu$Datelake[1:557] 6.165e-03  1.621e-03   3.803 0.000165 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 74.82 on 408 degrees of freedom
    # (147 observations deleted due to missingness)
    # Multiple R-squared:  0.03423,	Adjusted R-squared:  0.03187 
    # F-statistic: 14.46 on 1 and 408 DF,  p-value: 0.0001649

TDNtoTDPperiod2 <- lm(NtoPinsitu$TDNtoTDP[558:972] ~ NtoPinsitu$Datelake[558:972])
summary(TDNtoTDPperiod2)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)  
    # (Intercept)                  217.588484  67.213131   3.237  0.00132 **
    #   NtoPinsitu$Datelake[558:972]  -0.004776   0.005331  -0.896  0.37095   
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 295.8 on 369 degrees of freedom
    # (44 observations deleted due to missingness)
    # Multiple R-squared:  0.00217,	Adjusted R-squared:  -0.0005342 
    # F-statistic: 0.8025 on 1 and 369 DF,  p-value: 0.3709

TDNtoTDPplot <-
ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = TDNtoTDP), size = 0.5) +
  ylim(0,800) +
  ylab(expression(TDN:TDP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) +
  geom_smooth(method = 'lm', data = NtoPinsitu[1:557,], aes(x = Datelake, y = TDNtoTDP), se = FALSE, color = "blue") + #significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[558:972,], aes(x = Datelake, y = TDNtoTDP), se = FALSE, color = "black") + #non-significant slope
  geom_vline(xintercept = as.numeric(NtoPinsitu$Datelake[558]), lty = 5)

#### DIN:TDP (not used) ----
# DINtoTDPbreak = breakpoints (DINtoTDP ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
# summary(DINtoTDPbreak)
# #did not converge on a breakpoint (chose point 21 when h = 20, chose point 123 when h = 100)
# DINtoTDPbydate <- lm(DINtoTDP ~ Datelake)
# summary(DINtoTDPbydate)
# ggplot(NtoPinsitu, aes(x = Datelake)) +
#   geom_point(aes(y = DINtoTDP), size = 0.5) +
#   #ylim(0,500) +
#   ylab(expression(DIN:TDP)) +
#   xlab(" ") +
#   theme(legend.position = c(0.9,0.9))

#### PN:PP ----
PNtoPPbreak = breakpoints (PNtoPP ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
summary(PNtoPPbreak)
#did not converge on a breakpoint (chose point 20 when h = 20, chose point 100 when h = 100)

PNtoPPbydate <- lm(NtoPinsitu$PNtoPP ~ NtoPinsitu$Datelake)
summary(PNtoPPbydate)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)         36.3822046  1.8676529  19.480  < 2e-16 ***
    #   NtoPinsitu$Datelake -0.0005515  0.0002026  -2.722  0.00663 ** 
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 28.02 on 772 degrees of freedom
    # (198 observations deleted due to missingness)
    # Multiple R-squared:  0.009507,	Adjusted R-squared:  0.008224 
    # F-statistic:  7.41 on 1 and 772 DF,  p-value: 0.006633

PNtoPPplot <-
  ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = PNtoPP), size = 0.5) +
  ylim(0,100) +
  ylab(expression(PN:PP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) +
  geom_smooth(method = 'lm', data = NtoPinsitu, aes(x = Datelake, y = PNtoPP), se = FALSE, color = "blue") #significant slope

#### DIN ----
DINbreak = breakpoints (DIN_molar ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
summary(DINbreak)
    # Breakpoints at observation number:
    #   m = 1   176

DINbydate <- lm(NtoPinsitu$DIN_molar ~ NtoPinsitu$Datelake)
summary(DINbydate)
    # Coefficients:
    # Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)          8.076e+00  7.002e-01  11.534  < 2e-16 ***
    #   NtoPinsitu$Datelake -3.058e-04  7.651e-05  -3.997 7.02e-05 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 10.47 on 778 degrees of freedom
    # (192 observations deleted due to missingness)
    # Multiple R-squared:  0.02012,	Adjusted R-squared:  0.01886 
    # F-statistic: 15.98 on 1 and 778 DF,  p-value: 7.024e-05

DINperiod1 <- lm(NtoPinsitu$DIN_molar[1:176] ~ NtoPinsitu$Datelake[1:176])
summary(DINperiod1)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                12.329699   1.513276   8.148 1.48e-13 ***
    #   NtoPinsitu$Datelake[1:176] -0.003616   0.001040  -3.478 0.000664 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 10.59 on 147 degrees of freedom
    # (27 observations deleted due to missingness)
    # Multiple R-squared:  0.07604,	Adjusted R-squared:  0.06975 
    # F-statistic:  12.1 on 1 and 147 DF,  p-value: 0.0006643

DINperiod2 <- lm(NtoPinsitu$DIN_molar[176:972] ~ NtoPinsitu$Datelake[176:972])
summary(DINperiod2)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                   7.341e+00  1.012e+00   7.255 1.19e-12 ***
    #   NtoPinsitu$Datelake[176:972] -2.354e-04  9.975e-05  -2.360   0.0186 *  
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 10.36 on 630 degrees of freedom
    # (165 observations deleted due to missingness)
    # Multiple R-squared:  0.008762,	Adjusted R-squared:  0.007189 
    # F-statistic: 5.569 on 1 and 630 DF,  p-value: 0.01858

DIN_molarplot <-
ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = DIN_molar), size = 0.5) +
  ylim(0,125) +
  ylab(expression(DIN ~ (mu*M))) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) +
  geom_smooth(method = 'lm', data = NtoPinsitu[1:176,], aes(x = Datelake, y = DIN_molar), se = FALSE, color = "blue") + #significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[176:972,], aes(x = Datelake, y = DIN_molar), se = FALSE, color = "blue") + # significant slope
  geom_vline(xintercept = as.numeric(NtoPinsitu$Datelake[176]), lty = 5)

#### Fertilizer N:P ----
FertNtoPplot <-
ggplot(NtoPinput, aes(x = Dateinflow, y = Fert_NtoP)) +
  geom_point(size = 0.5) +
  ylim(0,50) +
  ylab(expression(FertN:FertP)) +
  xlab(" ") +
  theme(axis.text.x=element_blank()) 

#### Inflow N:P ----
Inflow_NtoPbreak = breakpoints (Inflow_NtoP ~ Dateinflow, h = 20, data = NtoPinput, breaks=1)
summary(Inflow_NtoPbreak)
    # Breakpoints at observation number:
    #   m = 1   884

Inflow_NtoPbydate <- lm(NtoPinput$Inflow_NtoP ~ NtoPinput$Dateinflow)
summary(Inflow_NtoPbydate)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept) 5.944e+00  3.253e-01  18.274  < 2e-16 ***
    #   Dateinflow  2.123e-04  3.789e-05   5.601  2.3e-08 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 9.387 on 3330 degrees of freedom
    # (4155 observations deleted due to missingness)
    # Multiple R-squared:  0.009334,	Adjusted R-squared:  0.009037 
    # F-statistic: 31.38 on 1 and 3330 DF,  p-value: 2.299e-08

Inflow_NtoPperiod1 <- lm(NtoPinput$Inflow_NtoP[1:883] ~ NtoPinput$Dateinflow[1:883])
summary(Inflow_NtoPperiod1)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)        6.5099257  0.5081759  12.810  < 2e-16 ***
    #   Dateinflow[1:883] -0.0017974  0.0006149  -2.923  0.00368 ** 
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 6.071 on 362 degrees of freedom
    # (519 observations deleted due to missingness)
    # Multiple R-squared:  0.02306,	Adjusted R-squared:  0.02036 
    # F-statistic: 8.546 on 1 and 362 DF,  p-value: 0.003681

Inflow_NtoPperiod2 <- lm(NtoPinput$Inflow_NtoP[884:7487] ~ NtoPinput$Dateinflow[884:7487])
summary(Inflow_NtoPperiod2)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)          6.389e+00  4.283e-01  14.917  < 2e-16 ***
    #   Dateinflow[884:7487] 1.692e-04  4.712e-05   3.591 0.000335 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 9.705 on 2966 degrees of freedom
    # (3636 observations deleted due to missingness)
    # Multiple R-squared:  0.004328,	Adjusted R-squared:  0.003993 
    # F-statistic: 12.89 on 1 and 2966 DF,  p-value: 0.000335

NtoPinputplot <-
ggplot(NtoPinput, aes(x = Dateinflow)) +
  geom_point(aes(y = Inflow_NtoP), size = 0.5) +
  ylab(expression(InflowN:InflowP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9))  #+
  #geom_smooth(method = 'lm', data = NtoPinput[1:883,], aes(x = Dateinflow, y = Inflow_NtoP), se = FALSE, color = "blue") + #significant slope
  #geom_smooth(method = 'lm', data = NtoPinput[884:7487,], aes(x = Dateinflow, y = Inflow_NtoP), se = FALSE, color = "blue") + #significant slope
  #geom_vline(xintercept = as.numeric(NtoPinput$Datelake[884]), lty = 5)

#### N:P fertilizer + inflow (not used) ----
ggplot(NtoPinput, aes(x = Dateinflow)) +
  geom_point(aes(y = Input_NtoP), size = 0.5) +
  ylim(0,50) +
  ylab(expression(InputN:InputP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) 

#### multiplot ----
library(gridExtra)
grid.arrange(TNtoTPplot, TDNtoTDPplot, PNtoPPplot, DIN_molarplot, FertNtoPplot, NtoPinputplot, ncol = 3)


#### non-parametric trend and changepoint testing for TN:TP, TDN:TDP, PP (non-normally distributed)
library(trend)
TNtoTPdataset <- NtoPinsitu %>%
  select(Datelake, TNtoTP) %>%
  drop_na()
rownames(TNtoTPdataset) <- NULL

TDNtoTDPdataset <- NtoPinsitu %>%
  select(Datelake, TDNtoTDP) %>%
  drop_na()
rownames(TDNtoTDPdataset) <- NULL

PPdataset <- NtoPinsitu %>%
  select(Datelake, PP) %>%
  drop_na()
rownames(PPdataset) <- NULL

# Pettitt's test: detects a single changepoint in hydrological/climate series with continuous data

TNtoTPvec <- as.vector(TNtoTPdataset$TNtoTP)
pettitt.test(TNtoTPvec)
    # data:  TNtoTPvec
    # U* = 30578, p-value = 3.796e-06
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 470 

TDNtoTDPvec <- as.vector(TDNtoTDPdataset$TDNtoTDP)
pettitt.test(TDNtoTDPvec)
    # data:  TDNtoTDPvec
    # U* = 27951, p-value = 0.0001079
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 120 

PPvec <- as.vector(PPdataset$PP)
pettitt.test(PPvec)
    # data:  PPvec
    # U* = 31464, p-value = 1.245e-05
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 179 

Cyanovec <- as.vector(CyanoPercent$Nfixcyano.percent)
pettitt.test(Cyanovec)
    # data:  Cyanovec
    # U* = 40536, p-value < 2.2e-16
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 87 

# Using shannon index requires running the phytoplankton portion of Unusual Year Analysis.
Shannonvec <- as.vector(Shannondataset$ShannonIndex)
pettitt.test(Shannonvec)
    # data:  Shannonvec
    # U* = 33428, p-value = 3.019e-15
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 232 

Shannonmeanvec <- as.vector(MeanShannon$shannon)
pettitt.test(Shannonmeanvec)
    # U* = 344, p-value = 0.0009788
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 16 

# Mann-Kendall test: detect monotonic trends in series of environmental/climate/hydrological data
mk.test(TNtoTPvec[1:469])
    # data:  TNtoTPvec[1:469]
    # z = 4.0215, n = 469, p-value = 5.782e-05
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 1.363800e+04 1.149894e+07 1.242727e-01 

mk.test(TNtoTPvec[470:752])
    # data:  TNtoTPvec[470:752]
    # z = 0.77053, n = 283, p-value = 0.441
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 1.227000e+03 2.531613e+06 3.075342e-02 

mk.test(TDNtoTDPvec[1:119])
    # data:  TDNtoTDPvec[1:119]
    # z = 1.4424, n = 119, p-value = 0.1492
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 6.290000e+02 1.895640e+05 8.960752e-02 

mk.test(TDNtoTDPvec[120:781])
    # data:  TDNtoTDPvec[120:781]
    # z = -3.5267, n = 662, p-value = 0.0004207
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -2.004700e+04  3.230782e+07 -9.167299e-02 

mk.test(PPvec[1:178])
    # data:  PPvec[1:178]
    # z = 1.3246, n = 178, p-value = 0.1853
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 1.053000e+03 6.307977e+05 6.781124e-02 

mk.test(PPvec[179:774])
    # data:  PPvec[179:774]
    # z = -6.4515, n = 596, p-value = 1.108e-10
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -3.131600e+04  2.356052e+07 -1.787924e-01 

mk.test(Cyanovec[1:86])
    # data:  Cyanovec[1:86]
    # z = 4.1489, n = 86, p-value = 3.34e-05
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 4.910000e+02 1.394833e+04 3.650356e-01 

mk.test(Cyanovec[87:611])
    # data:  Cyanovec[87:611]
    # z = -1.4936, n = 525, p-value = 0.1353
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -5.998000e+03  1.612026e+07 -4.370521e-02 

mk.test(Shannonvec[1:231])
    # data:  Shannonvec[1:231]
    # z = 4.5977, n = 231, p-value = 4.272e-06
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 5.399000e+03 1.378428e+06 2.032373e-01 

mk.test(Shannonvec[232:581])
    # data:  Shannonvec[232:581]
    # z = 3.5962, n = 350, p-value = 0.0003228
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 7.867000e+03 4.784208e+06 1.288088e-01 


# Sen's slope: computes slope and intercept, plus confidence limits
# sens.slope in trend does not provide intercept

library(colormap)
colormap(colormap = colormaps$viridis, nshades = 72, format = "hex",
         alpha = 1, reverse = FALSE)
scales::show_col(colormap(colormap = colormaps$magma, nshades=15))

TNtoTPearly <- TNtoTPdataset[1:469,]
TNtoTPlate <- TNtoTPdataset[470:752,]
TNtoTPplot <- 
  ggplot() +
  geom_point(data = TNtoTPearly, aes(x = Datelake, y = TNtoTP), size = 0.5, color = "#fc986cff") + #significant positive slope
  geom_point(data = TNtoTPlate, aes(x = Datelake, y = TNtoTP), size = 0.5, color = "gray40") + #non-significant slope
  ylim(0,200) +
  ylab(expression(TN:TP)) +
  xlab(" ") +
  theme(axis.text.x=element_blank()) + 
  geom_vline(xintercept = as.numeric(TNtoTPdataset$Datelake[470]), lty = 5)

TDNtoTDPearly <- TDNtoTDPdataset[1:119,]
TDNtoTDPlate <- TDNtoTDPdataset[120:781,]
TDNtoTDPplot <- 
  ggplot() +
  geom_point(data = TDNtoTDPearly, aes(x = Datelake, y = TDNtoTDP), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = TDNtoTDPlate, aes(x = Datelake, y = TDNtoTDP), size = 0.5, color = "#5d177dff") + #non-significant slope
  ylim(0,300) +
  ylab(expression(TDN:TDP)) +
  xlab(" ") +
  theme(axis.text.x=element_blank()) + 
  geom_vline(xintercept = as.numeric(TDNtoTDPdataset$Datelake[120]), lty = 5)

PPearly <- PPdataset[1:178,]
PPlate <- PPdataset[179:774,]
PPplot <- 
  ggplot() +
  geom_line(data = PPdataset, aes(x = Datelake, y = PP), color = "#b5367aff") +
  geom_point(data = PPearly, aes(x = Datelake, y = PP), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = PPlate, aes(x = Datelake, y = PP), size = 0.5, color = "#5d177dff") + #non-significant slope
  #ylim(0,200) +
  ylab(expression(PP)) +
  xlab(" ") +
  theme() + 
  geom_vline(xintercept = as.numeric(PPdataset$Datelake[179]), lty = 5)

cyanoplot <- 
  ggplot() +
  geom_point(data = CyanoPercent, aes(x = Date, y = Nfixcyano.percent), size = 0.5) +
  ylab(expression("Prop Cyano")) + 
  xlab("") + 
  theme(axis.text.x=element_blank())

Shannonearly <- Shannondataset[1:231,]
Shannonlate <- Shannondataset[232:581,]
Shannonearly$ShannonIndex <- round(Shannonearly$ShannonIndex, digits = 2)
Shannonlate$ShannonIndex <- round(Shannonlate$ShannonIndex, digits = 2)

Shannonplot <- 
  ggplot() +
  geom_point(data = Shannonearly, aes(x = Date, y = ShannonIndex), size = 0.5, color = "#fc986cff") + #significant positive slope
  geom_point(data = Shannonlate, aes(x = Date, y = ShannonIndex), size = 0.5, color = "#fc986cff") + #non-significant slope
  #ylim(0,200) +
  ylab(expression("Shannon Index")) +
  xlab(" ") +
  theme(axis.text.x=element_blank()) + 
  geom_vline(xintercept = as.numeric(Shannondataset$Date[232]), lty = 5)

Shannonmeanplot <- 
  ggplot(MeanShannon, aes(x = Year, y = shannon)) +
  geom_point(size = 0.5) + 
  geom_errorbar(aes(ymin = shannon - sdshannon, ymax = shannon + sdshannon)) + 
  ylab(expression("Shannon Index")) +
  xlab(" ")

var.test(Shannondataset$ShannonIndex[Shannondataset$Date < "1975-01-01"], 
         Shannondataset$ShannonIndex[Shannondataset$Date > "1989-12-31"])
var.test(Shannondataset$ShannonIndex[Shannondataset$Date < "1975-01-01"], 
         Shannondataset$ShannonIndex[Shannondataset$Date > "1975-01-01" & Shannondataset$Date < "1989-12-31"])
var.test(Shannondataset$ShannonIndex[Shannondataset$Date > "1975-01-01" & Shannondataset$Date < "1989-12-31"], 
         Shannondataset$ShannonIndex[Shannondataset$Date > "1989-12-31"])
var.test(Shannondataset$ShannonIndex[Shannondataset$Date < "1989-12-31"], 
         Shannondataset$ShannonIndex[Shannondataset$Date > "1989-12-31"])
var.test(PPdataset$PP[PPdataset$Datelake < "1989-12-31"], 
         PPdataset$PP[PPdataset$Datelake > "1989-12-31"])
var.test(PPdataset$PP[PPdataset$Datelake < "1975-01-01"], 
         PPdataset$PP[PPdataset$Datelake > "1989-12-31"])
var.test(PPdataset$PP[PPdataset$Datelake > "1975-01-01" & PPdataset$Datelake < "1989-12-31"], 
         PPdataset$PP[PPdataset$Datelake > "1989-12-31"])
var.test(PPdataset$PP[PPdataset$Datelake < "1975-01-01"], 
         PPdataset$PP[PPdataset$Datelake > "1975-01-01" & PPdataset$Datelake < "1989-12-31"])

library(gridExtra)
library(ggpubr)
library(gtable)
grid.arrange(FertNtoPplot, TNtoTPplot, cyanoplot, TDNtoTDPplot, Shannonplot, PPplot,  ncol = 2)

g = rbind(FertNtoPplot, TNtoTPplot, cyanoplot, TDNtoTDPplot, Shannonplot, PPplot, size = "first")
g$widths = unit.pmax(FertNtoPplot$widths, TNtoTPplot$widths, cyanoplot$widths, TDNtoTDPplot$widths, Shannonplot$widths, PPplot$widths)
