
datechar=datestr(MyLake_results.basin1.days);
dates=datevec(datechar);

Totalchl = (MyLake_results.basin1.concentrations.Chl + MyLake_results.basin1.concentrations.C)/0.42;
Chlintegratedepi = transpose(mean(Totalchl(1:8,1:3840)));
cyano = rdivide(MyLake_results.basin1.concentrations.C/0.42,Totalchl);
cyanointegratedepi = transpose(mean(cyano(1:8,1:3840)));

TDP = MyLake_results.basin1.concentrations.P + MyLake_results.basin1.concentrations.DOP;
TDPintegratedepi = transpose(mean(TDP(1:8,1:3840)));

TP = Totalchl + TDP + MyLake_results.basin1.concentrations.PP;
TPintegratedepi = transpose(mean(TP(1:8,1:3840)));

filename='Postproc_code/L227/Output_IntegratedEpi.csv';
M = [dates(:,1:3), Chlintegratedepi, cyanointegratedepi, TDPintegratedepi, TPintegratedepi];
csvwrite(filename,M);

Temp1m = transpose(MyLake_results.basin1.T(2,:));
Temp4m = transpose (MyLake_results.basin1.T(8,:));
Temp9m = transpose (MyLake_results.basin1.T(18,:));

Oxy2m = transpose(MyLake_results.basin1.concentrations.O2(4,:)/1000);
Oxy3m = transpose(MyLake_results.basin1.concentrations.O2(6,:)/1000);
Oxy4m = transpose(MyLake_results.basin1.concentrations.O2(8,:)/1000);
Oxy6m = transpose(MyLake_results.basin1.concentrations.O2(12,:)/1000);
Oxy8m = transpose(MyLake_results.basin1.concentrations.O2(16,:)/1000);
Oxy10m = transpose(MyLake_results.basin1.concentrations.O2(20,:)/1000);

TFe = MyLake_results.basin1.concentrations.Fe2 + MyLake_results.basin1.concentrations.Fe3;

filename='Postproc_code/L227/Output_Depths.csv';
M1 = [dates(:,1:3), Temp1m, Temp4m, Temp9m, Oxy2m, Oxy3m, Oxy4m, Oxy6m, Oxy8m, Oxy10m];
csvwrite(filename,M1);


