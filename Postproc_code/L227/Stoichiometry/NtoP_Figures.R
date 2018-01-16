setwd("/Users/krsalkgu/Documents/SourceTree/Lake227/Postproc_code/L227/Stoichiometry")

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

#Convert dates from factor to date format
Datelake <- as.Date(NtoPLake$Date, "%m/%d/%y")
Dateinflow <- as.Date(NtoPInflow$Date, "%m/%d/%y")

#Create data frames for N:P stoichiometry in lake and inflows
NtoPinsitu <- data.frame(Datelake, TNtoTP, PNtoPP, DINtoTDP, TDNtoTDP, DIN_molar)
NtoPinput <- data.frame(Dateinflow, Fert_NtoP, Inflow_NtoP, Input_NtoP, Input_TN_molar)

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
    #   m = 1   565

    #   m   0      1     
    # RSS 727962 703734
    # BIC   7654   7648
TNtoTPbreak2 = breakpoints (TNtoTP ~ Datelake, h = 20, data = NtoPinsitu, breaks=2)
summary(TNtoTPbreak2)
    # Breakpoints at observation number:
    # m = 1       565
    # m = 2   342 366

    #   m   0      1      2     
    # RSS 727962 703734 673869
    # BIC   7654   7648   7633
TNtoTPbydate <- lm(TNtoTP ~ Datelake)
summary(TNtoTPbydate)
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -61.957 -18.607  -6.346  10.281 250.469 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept) 72.4854883  2.0096124  36.069  < 2e-16 ***
    #   Datelake    -0.0012915  0.0002486  -5.195 2.61e-07 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 30.39 on 788 degrees of freedom
    # (225 observations deleted due to missingness)
    # Multiple R-squared:  0.03311,	Adjusted R-squared:  0.03189 
    # F-statistic: 26.99 on 1 and 788 DF,  p-value: 2.609e-07
TNtoTPperiod1 <- TNtoTPbydate <- lm(TNtoTP[1:564] ~ Datelake[1:564])
summary(TNtoTPperiod1)
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -58.920 -23.319 -10.876   9.573 255.278 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)     67.8801430  3.3360709  20.347   <2e-16 ***
    #   Datelake[1:564] -0.0001454  0.0008486  -0.171    0.864    
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 35.09 on 415 degrees of freedom
    # (147 observations deleted due to missingness)
    # Multiple R-squared:  7.071e-05,	Adjusted R-squared:  -0.002339 
    # F-statistic: 0.02935 on 1 and 415 DF,  p-value: 0.8641

TNtoTPperiod2 <- TNtoTPbydate <- lm(TNtoTP[565:1015] ~ Datelake[565:1015])
summary(TNtoTPperiod2)
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -63.020 -13.352  -3.307  10.433 123.269 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)        89.6086529  4.9549816  18.085  < 2e-16 ***
    #   Datelake[565:7487] -0.0028192  0.0004502  -6.262 1.05e-09 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 23.66 on 371 degrees of freedom
    # (6550 observations deleted due to missingness)
    # Multiple R-squared:  0.09559,	Adjusted R-squared:  0.09315 
    # F-statistic: 39.21 on 1 and 371 DF,  p-value: 1.052e-09

TNtoTPplot <- 
  ggplot(NtoPinsitu, aes(x = Datelake, y = TNtoTP)) +
  geom_point(aes(y = TNtoTP), size = 0.5) +
  ylim(0,200) +
  ylab(expression(TN:TP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) + 
  geom_smooth(method = 'lm', data = NtoPinsitu[1:564,], aes(x = Datelake, y = TNtoTP), se = FALSE, color = "black") + #non-significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[565:1015,], aes(x = Datelake, y = TNtoTP), se = FALSE, color = "blue") + #significant slope
  geom_vline(xintercept = as.numeric(NtoPinsitu$Datelake[565]), lty = 5)

#### TDN:TDP ----
TDNtoTDPbreak = breakpoints (TDNtoTDP ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
summary(TDNtoTDPbreak)
    # Breakpoints at observation number:
    #   m = 1   657
    #   m   0        1       
    # RSS 42649447 41521116
    # BIC    11226    11225
TDNtoTDPbydate <- lm(TDNtoTDP ~ Datelake)
summary(TDNtoTDPbydate)
    # Residuals:
    #   Min     1Q Median     3Q    Max 
    # -179.3  -77.0  -45.6    9.4 3799.2 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept) 1.728e+02  1.473e+01  11.730   <2e-16 ***
    #   Datelake    6.240e-04  1.838e-03   0.339    0.734    
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 228.6 on 816 degrees of freedom
    # (197 observations deleted due to missingness)
    # Multiple R-squared:  0.0001412,	Adjusted R-squared:  -0.001084 
    # F-statistic: 0.1152 on 1 and 816 DF,  p-value: 0.7344
TDNtoTDPperiod1 <- lm(TDNtoTDP[1:656] ~ Datelake[1:656])
summary(TDNtoTDPperiod1)
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -198.66  -76.53  -44.98   40.28 1232.76 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)     1.481e+02  1.219e+01  12.147  < 2e-16 ***
    #   Datelake[1:656] 7.783e-03  2.736e-03   2.845  0.00462 ** 
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 140.6 on 501 degrees of freedom
    # (153 observations deleted due to missingness)
    # Multiple R-squared:  0.0159,	Adjusted R-squared:  0.01394 
    # F-statistic: 8.094 on 1 and 501 DF,  p-value: 0.004622
TDNtoTDPperiod2 <- lm(TDNtoTDP[657:1015] ~ Datelake[657:1015])
summary(TDNtoTDPperiod2)
    # Residuals:
    #   Min     1Q Median     3Q    Max 
    # -174.6  -77.2  -46.9   -6.7 3802.5 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)  
    # (Intercept)        226.225795  89.333142   2.532   0.0118 *
    #   Datelake[657:1015]  -0.004439   0.007687  -0.578   0.5640  
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 322.5 on 313 degrees of freedom
    # (44 observations deleted due to missingness)
    # Multiple R-squared:  0.001064,	Adjusted R-squared:  -0.002127 
    # F-statistic: 0.3335 on 1 and 313 DF,  p-value: 0.564

TDNtoTDPplot <-
ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = TDNtoTDP), size = 0.5) +
  ylim(0,800) +
  ylab(expression(TDN:TDP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) +
  geom_smooth(method = 'lm', data = NtoPinsitu[1:656,], aes(x = Datelake, y = TDNtoTDP), se = FALSE, color = "blue") + #significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[657:1015,], aes(x = Datelake, y = TDNtoTDP), se = FALSE, color = "black") + #non-significant slope
  geom_vline(xintercept = as.numeric(NtoPinsitu$Datelake[657]), lty = 5)

#### DIN:TDP (not used) ----
# DINtoTDPbreak = breakpoints (DINtoTDP ~ Datelake, h = 20, data = NtoPinsitu, breaks=1)
# summary(DINtoTDPbreak)
# #did not converge on a breakpoint (chose point 20 when h = 20, chose point 100 when h = 100)
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

PNtoPPbydate <- lm(PNtoPP ~ Datelake)
summary(PNtoPPbydate)
    # Residuals:
    #   Min     1Q Median     3Q    Max 
    # -32.48 -10.05  -3.57   3.00 529.04 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept) 34.5587251  1.7431545  19.825   <2e-16 ***
    #   Datelake    -0.0005371  0.0002180  -2.463    0.014 *  
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 27.31 on 811 degrees of freedom
    # (202 observations deleted due to missingness)
    # Multiple R-squared:  0.007425,	Adjusted R-squared:  0.006201 
    # F-statistic: 6.067 on 1 and 811 DF,  p-value: 0.01398

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
    #   m = 1   537
    #   m   0      1     
    # RSS 290347 284694
    # BIC   7176   7180
DINbydate <- lm(DIN_molar ~ Datelake)
summary(DINbydate)
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -16.315 -10.718  -6.835   0.653 105.021 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept) 16.7277826  1.2075253  13.853  < 2e-16 ***
    #   Datelake    -0.0007954  0.0001511  -5.265 1.79e-07 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 18.82 on 820 degrees of freedom
    # (193 observations deleted due to missingness)
    # Multiple R-squared:  0.0327,	Adjusted R-squared:  0.03152 
    # F-statistic: 27.72 on 1 and 820 DF,  p-value: 1.79e-07
DINperiod1 <- lm(DIN_molar[1:536] ~ Datelake[1:536])
summary(DINperiod1)
    # Residuals:
    #   Min     1Q Median     3Q    Max 
    # -14.03 -13.16 -10.93   8.09  95.28 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)      1.426e+01  1.876e+00   7.603    2e-13 ***
    #   Datelake[1:536] -7.934e-05  5.121e-04  -0.155    0.877    
    # ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 20.25 on 409 degrees of freedom
    # (125 observations deleted due to missingness)
    # Multiple R-squared:  5.868e-05,	Adjusted R-squared:  -0.002386 
    # F-statistic: 0.024 on 1 and 409 DF,  p-value: 0.877

DINperiod2 <- lm(DIN_molar[537:1015] ~ Datelake[537:1015])
summary(DINperiod2)
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -13.369  -9.360  -5.143  -1.017 104.408 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)        22.8118412  3.2234141   7.077 6.45e-12 ***
    #   Datelake[537:1015] -0.0013616  0.0003014  -4.517 8.21e-06 ***
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 17.18 on 409 degrees of freedom
    # (68 observations deleted due to missingness)
    # Multiple R-squared:  0.04752,	Adjusted R-squared:  0.04519 
    # F-statistic:  20.4 on 1 and 409 DF,  p-value: 8.212e-06

DIN_molarplot <-
ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = DIN_molar), size = 0.5) +
  ylim(0,125) +
  ylab(expression(DIN ~ (mu*M))) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) +
  geom_smooth(method = 'lm', data = NtoPinsitu[1:536,], aes(x = Datelake, y = DIN_molar), se = FALSE, color = "black") + #non-significant slope
  geom_smooth(method = 'lm', data = NtoPinsitu[537:1015,], aes(x = Datelake, y = DIN_molar), se = FALSE, color = "blue") + #significant slope
  geom_vline(xintercept = as.numeric(NtoPinsitu$Datelake[537]), lty = 5)

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
    #   m   0      1     
    # RSS 293396 280929
    # BIC  24401  24280

Inflow_NtoPbydate <- lm(Inflow_NtoP ~ Dateinflow)
summary(Inflow_NtoPbydate)
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -6.984  -3.741  -1.979   0.453 150.860 
    # 
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

Inflow_NtoPperiod1 <- lm(Inflow_NtoP[1:883] ~ Dateinflow[1:883])
summary(Inflow_NtoPperiod1)
    # Residuals:
    #   Min     1Q Median     3Q    Max 
    # -5.936 -2.818 -0.851  1.342 81.907 
    # 
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

Inflow_NtoPperiod2 <- lm(Inflow_NtoP[884:7487] ~ Dateinflow[884:7487])
summary(Inflow_NtoPperiod2)
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -6.871  -3.926  -2.070   0.352 150.594 
    # 
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
