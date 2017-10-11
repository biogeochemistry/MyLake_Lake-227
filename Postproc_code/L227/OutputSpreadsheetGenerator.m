
datechar=datestr(MyLake_results.basin1.days);
dates=datevec(datechar);

Totalchl = (MyLake_results.basin1.concentrations.Chl + MyLake_results.basin1.concentrations.C)/0.42;
Chlintegratedepi = transpose(mean(Totalchl(1:4,1:3653)));
cyano = rdivide(MyLake_results.basin1.concentrations.C/0.42,Totalchl);
cyanointegratedepi = transpose(mean(cyano(1:4,1:3653)));

TDP = MyLake_results.basin1.concentrations.P + MyLake_results.basin1.concentrations.DOP;
TDPintegratedepi = transpose(mean(TDP(1:4,1:3653)));

TP = Totalchl + TDP + MyLake_results.basin1.concentrations.PP;
TPintegratedepi = transpose(mean(TP(1:4,1:3653)));

filename='Postproc_code/L227/Output_IntegratedEpi_old.csv';
M = [dates(:,1:3), Chlintegratedepi, cyanointegratedepi, TDPintegratedepi, TPintegratedepi];
csvwrite(filename,M);

Temp1m = transpose(MyLake_results.basin1.T(1,:));
Temp4m = transpose (MyLake_results.basin1.T(4,:));
Temp9m = transpose (MyLake_results.basin1.T(9,:));

Oxy1m = transpose(MyLake_results.basin1.concentrations.O2(1,:)/1000);
Oxy2m = transpose(MyLake_results.basin1.concentrations.O2(2,:)/1000);
Oxy3m = transpose(MyLake_results.basin1.concentrations.O2(3,:)/1000);
Oxy4m = transpose(MyLake_results.basin1.concentrations.O2(4,:)/1000);
Oxy5m = transpose(MyLake_results.basin1.concentrations.O2(5,:)/1000);
Oxy6m = transpose(MyLake_results.basin1.concentrations.O2(6,:)/1000);
Oxy7m = transpose(MyLake_results.basin1.concentrations.O2(7,:)/1000);
Oxy8m = transpose(MyLake_results.basin1.concentrations.O2(8,:)/1000);
Oxy9m = transpose(MyLake_results.basin1.concentrations.O2(9,:)/1000);
Oxy10m = transpose(MyLake_results.basin1.concentrations.O2(10,:)/1000);

filename='Postproc_code/L227/Output_Depths.csv';
M1 = [dates(:,1:3), Temp1m, Temp4m, Temp9m, Oxy1m, Oxy2m, Oxy3m, Oxy4m, Oxy5m, Oxy6m, Oxy7m, Oxy8m, Oxy9m, Oxy10m];
csvwrite(filename,M1);


