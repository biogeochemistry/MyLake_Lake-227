
#Interpolation script for L227 inflows
#All variables that have non-zero measurements need to be interpolated 
#(the model needs to receive inputs even if they weren't measured as part of the monitoring data)
#Fertilizer inputs (TP and NO3) do not need to be interpolated, 
#but inflows do. Thus, I will need to interpolate the inflow data, then
#put it back into the Fertilization Input Calculations.xlsx file to add
#the fertilizer back into the inflows. 

#For variables that have NaN as the first and last value in the dataset, 
#I added the average value of the dataset 
#1969-2016 averages
#TP: 18.37
#DOC: 28833.12
#NO3: 39.22
#NH4: 31.32
#SO4: 4112.95
#Fe2: 420.02
#Ca2: 2130.98
#pH: 4.80
#Fe3: 93.22
#SiO2: 4875.23

#Variables that have a random distribution should be aggregated 
#by the average value in the rest of the dataset. These include: 
#Cloud cover
#Humidity (there is a seasonal cycle, but 1973 has almost a whole year of no measurements)
#Air pressure
#Wind speed
#pH

#Variables that have a strong seasonal trend or shifts by year should be interpolated
#linearly. These include: 
#Air temperature
#Inflow temperature
#TP
#DOC
#NO3
#NH4
#SO4
#Fe2
#Ca2
#Fe3
#SiO2

dataset = read.csv("Inflow_for_interpolation.csv")
attach(dataset)
head(dataset)
library(zoo)

#### Variables aggregated by mean of remaining dataset ----
CloudI = na.aggregate(Cloud)
#Check
#Cloud2 = na.omit(Cloud)
#mean(Cloud2)
#CloudI[1372:1462]
HumidityI = na.aggregate(Humidity)
AirPressureI = na.aggregate(AirPressure)
WindSpeedI = na.aggregate(WindSpeed)
pHI = na.aggregate(pH)

#Variables interpolated linearly 
AirTempI = na.approx(AirTemp)
InflowTempI = na.approx(InflowTemp)
TPI = na.approx(TP)
TPI.constant = na.approx(TP, method = "constant")
DOCI = na.approx(DOC)
NO3I = na.approx(NO3)
NH4I = na.approx(NH4)
SO4I = na.approx(SO4)
Fe2I = na.approx(Fe2)
Ca2I = na.approx(Ca2)
Fe3I = na.approx(Fe3)
SiO2I = na.approx(SiO2)

#Write an output file
InterpolatedData = data.frame(Year, Month, Day, CloudI, AirTempI, HumidityI, AirPressureI, WindSpeedI, InflowTempI, TPI, DOCI, NO3I, NH4I, SO4I, Fe2I, Ca2I, pHI, Fe3I, SiO2I)
InterpolatedTP = data.frame(Year, Month, Day, TPI, TPI.constant)

#DOC was too high in the model inputs, which made the model crash. I will take all values > 50,000
#and replace them with 50,000

InterpolatedData$DOCI[InterpolatedData$DOCI > 50000] <- 50000

write.csv(InterpolatedData, "Inflow_Interpolated.csv", row.names = F)
write.csv(InterpolatedTP, "TP_Interpolated.csv", row.names = F)

