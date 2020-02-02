
datechar=datestr(MyLake_results.basin1.days);
dates=datevec(datechar);

epidepth = floor(MyLake_results.basin1.MixStat(12,:)*2)/2;
epidepth(isnan(epidepth)) = 10;
epidepthposition = 2*epidepth; %MyLake computes for every 0.5 m, so adding a factor of 2

diazoPP = MyLake_results.basin1.concentrations.C;
nondiazoPP = MyLake_results.basin1.concentrations.Chl;

diazoPPintegratedepi = zeros(1,length(diazoPP));
    for (i=1:length(diazoPP))
        diazoPPintegratedepi(i) = mean(diazoPP(1:epidepthposition(i), i));
    end %returns integrated epilimnion diazoPP measurement for each day (variable epilimnion depth)
diazoPPintegratedepi = transpose(diazoPPintegratedepi);

nondiazoPPintegratedepi = zeros(1,length(nondiazoPP));
    for (i=1:length(nondiazoPP))
        nondiazoPPintegratedepi(i) = mean(nondiazoPP(1:epidepthposition(i), i));
    end %returns integrated epilimnion TPP measurement for each day (variable epilimnion depth)
nondiazoPPintegratedepi = transpose(nondiazoPPintegratedepi);


filename='Postproc_code/L227/Output_IntegratedEpi_Phyto.csv';
M = [dates(:,1:3), diazoPPintegratedepi, nondiazoPPintegratedepi];
csvwrite(filename,M);

