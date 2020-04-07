
datechar=datestr(MyLake_results.basin1.days);
dates=datevec(datechar);

epidepth = floor(MyLake_results.basin1.MixStat(12,:)*2)/2;
epidepth(isnan(epidepth)) = 10;
epidepthposition = 2*epidepth; %MyLake computes for every 0.5 m, so adding a factor of 2

%Totalchl = (MyLake_results.basin1.concentrations.Chl + MyLake_results.basin1.concentrations.C)/0.42;
%Chlintegratedepi = transpose(mean(Totalchl(1:10,:)));
%cyano = rdivide(MyLake_results.basin1.concentrations.C/0.42,Totalchl);
%cyanointegratedepi = transpose(mean(cyano(1:10,:)));

TDP = MyLake_results.basin1.concentrations.P + MyLake_results.basin1.concentrations.DOP;
TDPintegratedepi = zeros(1,length(TDP));
    for (i=1:length(TDP))
        TDPintegratedepi(i) = mean(TDP(1:epidepthposition(i), i));
    end %returns integrated epilimnion TPP measurement for each day (variable epilimnion depth)
TDPintegratedepi = transpose(TDPintegratedepi);

%TP = MyLake_results.basin1.concentrations.Chl + MyLake_results.basin1.concentrations.C + TDP + MyLake_results.basin1.concentrations.PP;
%TPintegratedepi = transpose(mean(TP(1:8,:)));

TPP = MyLake_results.basin1.concentrations.Chl + MyLake_results.basin1.concentrations.C + MyLake_results.basin1.concentrations.PP;

TPPintegratedepi = zeros(1,length(TPP));
    for (i=1:length(TPP))
        TPPintegratedepi(i) = mean(TPP(1:epidepthposition(i), i));
    end %returns integrated epilimnion TPP measurement for each day (variable epilimnion depth)
TPPintegratedepi = transpose(TPPintegratedepi);

TPPintegratedepi = zeros(1,length(TPP));
    for (i=1:length(TPP))
        TPPintegratedepi(i) = mean(TPP(1:epidepthposition(i), i));
    end %returns integrated epilimnion TPP measurement for each day (variable epilimnion depth)
TPPintegratedepi = transpose(TPPintegratedepi);

DOC = MyLake_results.basin1.concentrations.DOC;

DOCintegratedepi = zeros(1,length(DOC));
    for (i=1:length(DOC))
        DOCintegratedepi(i) = mean(DOC(1:epidepthposition(i), i));
    end %returns integrated epilimnion TPP measurement for each day (variable epilimnion depth)
DOCintegratedepi = transpose(DOCintegratedepi);


filename='Postproc_code/L227/Output_IntegratedEpi_NandP.csv';
M = [dates(:,1:3), TDPintegratedepi, TPPintegratedepi, DOCintegratedepi];
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
Fe4m = transpose(TFe(8,:)/1000);
Fe6m = transpose(TFe(12,:)/1000);
Fe8m = transpose(TFe(16,:)/1000);
Fe10m = transpose(TFe(20,:)/1000);

filename='Postproc_code/L227/Output_Depths_NandP.csv';
M1 = [dates(:,1:3), Temp1m, Temp4m, Temp9m, Oxy2m, Oxy3m, Oxy4m, Oxy6m, Oxy8m, Oxy10m, Fe4m, Fe6m, Fe8m, Fe10m];
csvwrite(filename,M1);



