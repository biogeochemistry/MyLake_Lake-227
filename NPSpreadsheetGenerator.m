
datechar=datestr(MyLake_results.basin1.days);
dates=datevec(datechar);

epidepth = floor(MyLake_results.basin1.MixStat(12,:)*2)/2;
epidepth(isnan(epidepth)) = 10;
epidepthposition = 2*epidepth; %MyLake computes for every 0.5 m, so adding a factor of 2

DIP = MyLake_results.basin1.concentrations.P;
DIPintegratedepi = zeros(1,length(DIP));
    for (i=1:length(DIP))
        DIPintegratedepi(i) = mean(DIP(1:epidepthposition(i), i));
    end %returns integrated epilimnion DIP measurement for each day (variable epilimnion depth)
DIPintegratedepi = transpose(DIPintegratedepi);


DOP = MyLake_results.basin1.concentrations.DOP;
DOPintegratedepi = zeros(1,length(DOP));
    for (i=1:length(DOP))
        DOPintegratedepi(i) = mean(DOP(1:epidepthposition(i), i));
    end %returns integrated epilimnion DOP measurement for each day (variable epilimnion depth)
DOPintegratedepi = transpose(DOPintegratedepi);

% TDP = MyLake_results.basin1.concentrations.P + MyLake_results.basin1.concentrations.DOP;
% TDPintegratedepi = zeros(1,length(TDP));
%     for (i=1:length(TDP))
%         TDPintegratedepi(i) = mean(TDP(1:epidepthposition(i), i));
%     end %returns integrated epilimnion TDP measurement for each day (variable epilimnion depth)
TDPintegratedepi = DIPintegratedepi + DOPintegratedepi;



TPP = MyLake_results.basin1.concentrations.Chl + MyLake_results.basin1.concentrations.C + MyLake_results.basin1.concentrations.PP;
TPPintegratedepi = zeros(1,length(TPP));
    for (i=1:length(TPP))
        TPPintegratedepi(i) = mean(TPP(1:epidepthposition(i), i));
    end %returns integrated epilimnion TPP measurement for each day (variable epilimnion depth)
TPPintegratedepi = transpose(TPPintegratedepi);

TPintegratedepi = TDPintegratedepi + TPPintegratedepi;


NO3 = MyLake_results.basin1.concentrations.NO3;
NO3integratedepi = zeros(1,length(NO3));
    for (i=1:length(NO3))
        NO3integratedepi(i) = mean(NO3(1:epidepthposition(i), i));
    end %returns integrated epilimnion NO3 measurement for each day (variable epilimnion depth)
NO3integratedepi = transpose(NO3integratedepi);

NH4 = MyLake_results.basin1.concentrations.NH4;
NH4integratedepi = zeros(1,length(NH4));
    for (i=1:length(NH4))
        NH4integratedepi(i) = mean(NH4(1:epidepthposition(i), i));
    end %returns integrated epilimnion NH4 measurement for each day (variable epilimnion depth)
NH4integratedepi = transpose(NH4integratedepi);

DINintegratedepi = NO3integratedepi + NH4integratedepi;


filename='Postproc_code/L227/Output_IntegratedEpi_Period3_NandP.csv';
M = [dates(:,1:3), DIPintegratedepi, DOPintegratedepi, TDPintegratedepi, TPPintegratedepi, TPintegratedepi, NO3integratedepi, NH4integratedepi, DINintegratedepi];
csvwrite(filename,M);





