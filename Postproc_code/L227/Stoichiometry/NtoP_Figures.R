setwd("/Users/krsalkgu/Documents/SourceTree/Lake227/Postproc_code/L227/Stoichiometry")
library(tidyverse)

NtoPLake <- read.csv("NP_Stoichiometry_L227.csv")
attach(NtoPLake)
head(NtoPLake)
NtoPInflow <- read.csv("NP_Stoichiometry_Inflow.csv")
attach(NtoPInflow)
head(NtoPInflow)

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

#Create data frames for N:P stoichiometry in lake and inflows for May-October
NtoPinsitu <- data.frame(Datelake, Monthlake, TNtoTP, PNtoPP, DINtoTDP, TDNtoTDP, DIN_molar)
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

#### TN:TP ----
TNtoTPbreak = breakpoints (TNtoTP ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
summary(TNtoTPbreak)
    # Breakpoints at observation number:
    #   m = 1   470

TNtoTPbreak2 = breakpoints (TNtoTP ~ Datelake, h = 20, data = NtoPinsitu, breaks=2)
summary(TNtoTPbreak2)
    # Breakpoints at observation number:
    # m = 1       470
    # m = 2   112 132

TNtoTPbydate <- lm(NtoPinsitu$TNtoTP ~ NtoPinsitu$Datelake)
summary(TNtoTPbydate)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)         61.2304046  1.7452061  35.085   <2e-16 ***
    #   NtoPinsitu$Datelake -0.0005369  0.0002083  -2.577   0.0102 *  
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 23.92 on 686 degrees of freedom
    # (217 observations deleted due to missingness)
    # Multiple R-squared:  0.009588,	Adjusted R-squared:  0.008144 
    # F-statistic: 6.641 on 1 and 686 DF,  p-value: 0.01017

TNtoTPperiod1 <- TNtoTPbydate <- lm(NtoPinsitu$TNtoTP[1:469] ~ NtoPinsitu$Datelake[1:469])
summary(TNtoTPperiod1)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                5.813e+01  2.908e+00  19.993   <2e-16 ***
    #   NtoPinsitu$Datelake[1:469] 1.846e-04  7.457e-04   0.248    0.805    
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 27.19 on 329 degrees of freedom
    # (138 observations deleted due to missingness)
    # Multiple R-squared:  0.0001862,	Adjusted R-squared:  -0.002853 
    # F-statistic: 0.06129 on 1 and 329 DF,  p-value: 0.8046

TNtoTPperiod2 <- TNtoTPbydate <- lm(NtoPinsitu$TNtoTP[470:905] ~ NtoPinsitu$Datelake[470:905])
summary(TNtoTPperiod2)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                  72.3646322  4.3268408   16.73  < 2e-16 ***
    #   NtoPinsitu$Datelake[470:905] -0.0015216  0.0003931   -3.87 0.000129 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 20.27 on 355 degrees of freedom
    # (79 observations deleted due to missingness)
    # Multiple R-squared:  0.04049,	Adjusted R-squared:  0.03779 
    # F-statistic: 14.98 on 1 and 355 DF,  p-value: 0.0001293

TNtoTPplot <- 
  ggplot(NtoPinsitu, aes(x = Datelake, y = TNtoTP)) +
  geom_point(aes(y = TNtoTP), size = 0.5) +
  ylim(0,200) +
  ylab(expression(TN:TP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) + 
  geom_smooth(method = 'lm', data = NtoPinsitu[1:469,], aes(x = Datelake, y = TNtoTP), se = FALSE, color = "black") + #non-significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[470:905,], aes(x = Datelake, y = TNtoTP), se = FALSE, color = "blue") + #significant slope
  geom_vline(xintercept = as.numeric(NtoPinsitu$Datelake[470]), lty = 5)

#### TDN:TDP ----
TDNtoTDPbreak = breakpoints (TDNtoTDP ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
summary(TDNtoTDPbreak)
    # Breakpoints at observation number:
    #   m = 1   558

TDNtoTDPbydate <- lm(NtoPinsitu$TDNtoTDP ~ NtoPinsitu$Datelake)
summary(TDNtoTDPbydate)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)         1.251e+02  1.561e+01   8.013  4.6e-15 ***
    #   NtoPinsitu$Datelake 3.643e-03  1.881e-03   1.937   0.0531 .  
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 219.8 on 712 degrees of freedom
    # (191 observations deleted due to missingness)
    # Multiple R-squared:  0.005242,	Adjusted R-squared:  0.003845 
    # F-statistic: 3.752 on 1 and 712 DF,  p-value: 0.05313

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

TDNtoTDPperiod2 <- lm(NtoPinsitu$TDNtoTDP[558:905] ~ NtoPinsitu$Datelake[558:905])
summary(TDNtoTDPperiod2)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)  
    # (Intercept)                   1.778e+02  9.233e+01   1.926   0.0551 .
    # NtoPinsitu$Datelake[558:905] -9.062e-04  7.954e-03  -0.114   0.9094  
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 325.8 on 302 degrees of freedom
    # (44 observations deleted due to missingness)
    # Multiple R-squared:  4.298e-05,	Adjusted R-squared:  -0.003268 
    # F-statistic: 0.01298 on 1 and 302 DF,  p-value: 0.9094

TDNtoTDPplot <-
ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = TDNtoTDP), size = 0.5) +
  ylim(0,800) +
  ylab(expression(TDN:TDP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) +
  geom_smooth(method = 'lm', data = NtoPinsitu[1:557,], aes(x = Datelake, y = TDNtoTDP), se = FALSE, color = "blue") + #significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[558:905,], aes(x = Datelake, y = TDNtoTDP), se = FALSE, color = "black") + #non-significant slope
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
    # (Intercept)         37.3428413  2.0162596  18.521  < 2e-16 ***
    #   NtoPinsitu$Datelake -0.0007438  0.0002437  -3.053  0.00235 ** 
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 28.77 on 708 degrees of freedom
    # (195 observations deleted due to missingness)
    # Multiple R-squared:  0.01299,	Adjusted R-squared:  0.0116 
    # F-statistic: 9.319 on 1 and 708 DF,  p-value: 0.002352

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
    #   m = 1   185

DINbydate <- lm(NtoPinsitu$DIN_molar ~ NtoPinsitu$Datelake)
summary(DINbydate)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)          8.143e+00  7.607e-01  10.704  < 2e-16 ***
    #   NtoPinsitu$Datelake -3.200e-04  9.205e-05  -3.476  0.00054 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 10.82 on 718 degrees of freedom
    # (185 observations deleted due to missingness)
    # Multiple R-squared:  0.01655,	Adjusted R-squared:  0.01518 
    # F-statistic: 12.08 on 1 and 718 DF,  p-value: 0.0005398

DINperiod1 <- lm(NtoPinsitu$DIN_molar[1:184] ~ NtoPinsitu$Datelake[1:184])
summary(DINperiod1)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                12.2555964  1.4764695   8.301 5.07e-14 ***
    #   NtoPinsitu$Datelake[1:184] -0.0035193  0.0009763  -3.605 0.000422 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 10.44 on 153 degrees of freedom
    # (29 observations deleted due to missingness)
    # Multiple R-squared:  0.07828,	Adjusted R-squared:  0.07226 
    # F-statistic: 12.99 on 1 and 153 DF,  p-value: 0.0004222

DINperiod2 <- lm(NtoPinsitu$DIN_molar[185:905] ~ NtoPinsitu$Datelake[185:905])
summary(DINperiod2)
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)                   7.4890958  1.1929977   6.278 6.87e-10 ***
    #   NtoPinsitu$Datelake[185:905] -0.0002480  0.0001283  -1.932   0.0539 .  
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 10.84 on 563 degrees of freedom
    # (156 observations deleted due to missingness)
    # Multiple R-squared:  0.006586,	Adjusted R-squared:  0.004822 
    # F-statistic: 3.733 on 1 and 563 DF,  p-value: 0.05386

DIN_molarplot <-
ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = DIN_molar), size = 0.5) +
  ylim(0,125) +
  ylab(expression(DIN ~ (mu*M))) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) +
  geom_smooth(method = 'lm', data = NtoPinsitu[1:184,], aes(x = Datelake, y = DIN_molar), se = FALSE, color = "blue") + #significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[185:905,], aes(x = Datelake, y = DIN_molar), se = FALSE, color = "blue", lty = 5) + #marginally significant slope
  geom_vline(xintercept = as.numeric(NtoPinsitu$Datelake[185]), lty = 5)

#### Fertilizer N:P ----
FertNtoPplot <-
ggplot(NtoPinput, aes(x = Dateinflow, y = Fert_NtoP)) +
  geom_point(size = 0.5) +
  ylim(0,50) +
  ylab(expression(FertN:FertP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) 

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
