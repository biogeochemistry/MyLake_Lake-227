
%% Project specific setup script

% Population size for each generation of the genetic algorithm. If you use parallellization it will be
% more efficient if this is a multiple of the number of cores on the machine.
% It is generally more important to have a large
% population_size than it is to have a large max_generations since it is
% important to sample a large part of the parameter space before starting
% to sample around the best fits.
population_size = 36;  % 3 params 
% Max generations to run. The maximal amount of runs is population_size*max_generations.
max_generations = 8; % 
paralellize     = true; % Run many model evaluations in parallell (saves time if the computer has many cores).

m_start = [1971, 1, 1];
m_stop = [1975, 12, 31];

MyL_dates = datenum(m_start):datenum(m_stop);

% Loading observation data once so that it doesn't have to be loaded each time the model is evaluated.
Data = loadData(MyL_dates);

% Loading a priori K_lake and K_sedmient values for all parameters.
[K_lake, K_sediments] = load_params();

% varyindexes is the indexes of the parameters in K_lake and K_sediments we want to vary during
% optimization. The indexes start in K_lake and continue into K_sediments.
% Every value that has an index that belongs to the same column in
% varyindexes is given the same value. indexes in the second row or later
% can be NaN if you don't want the index in the first row of that column to
% covary with any in the next row.

% Example: (g_twty (parameter nr 50) is always given the same value as
% g_twty2 (parameter nr 58) and so on ..)
varyindexes = [25 71 73 74 75 76 77 78 85]; %I_scO, wc_factor k_chl, k_POP, k_POC, k_DOP, k_DOC, Km_O2, Kin_O2 

% Setting up the min and max boundaries for each covarying set of parameters.
minparam = [0.5, 0.05,  1, 1, 1, 1, 0.1, 0.003, 0.008];
maxparam = [1.2, 0.2, 1, 1, 1, 1, 0.6, 0.003, 0.008];

% The best initial guess for the values of each set of covarying parameters (can have
% multiple rows for multiple initial guesses. (up to population_size rows)
initial_guess = [1, 0.1, 1, 1, 1, 1, 0.04, 0.003, 0.008];

modeleval      = @MyLake_227_model_evaluation;
errfun         = @error_function_oxygen;
filenameprefix = 'O2'; % Prefix for the .mat file where the optimal parameters are saved in the end.

do_MyLake_optimization(m_start, m_stop, K_sediments, K_lake, Data, ...
    varyindexes, minparam, maxparam, initial_guess, modeleval, errfun,...
    population_size, max_generations, paralellize, filenameprefix);

%% Helper functions for project specific setup

function Data = loadData(MyL_dates)
    rawcsvdata = csvread('Postproc_code/L227/OutputForOptimization.csv', 1, 1);
    rawdatadates = datenum(rawcsvdata(:,9:11));
    withinmodelrange = (rawdatadates >= MyL_dates(1)) & (rawdatadates <= MyL_dates(end));
    rawO24m = rawcsvdata(:,4);
    rawO24m(rawO24m == 0) = NaN; %Unfortunately the readcsv function fills in 0s for missing fields. This fix only works if there are no legitimate 0s in the data. 
    Data.O24m = rawO24m(withinmodelrange);
    rawO26m = rawcsvdata(:,5);
    rawO26m(rawO26m == 0) = NaN; %Unfortunately the readcsv function fills in 0s for missing fields. This fix only works if there are no legitimate 0s in the data. 
    Data.O26m = rawO26m(withinmodelrange);
    rawO28m = rawcsvdata(:,6);
    rawO28m(rawO28m == 0) = NaN; %Unfortunately the readcsv function fills in 0s for missing fields. This fix only works if there are no legitimate 0s in the data. 
    Data.O28m = rawO28m(withinmodelrange);
    rawO210m = rawcsvdata(:,7);
    rawO210m(rawO210m == 0) = NaN; %Unfortunately the readcsv function fills in 0s for missing fields. This fix only works if there are no legitimate 0s in the data. 
    Data.O210m = rawO210m(withinmodelrange);
    dateswithinrange = rawdatadates(withinmodelrange);
    Data.date_mask = getmask(dateswithinrange, MyL_dates); % Used later to match model data to observed data by correct date.
end

function mask = getmask(subsetdates, alldates)
    [C,IA, mask] = intersect(subsetdates, alldates);
end

%% BEGIN project specific evaluation functions

% Model evaluation function. It should take the parameters
% (m_start, m_stop, K_sediment, K_lake) and return a ModelResult struct.
% The ModelResult struct should contain whatever the error function needs
% to compare the model result to data.
function ModelResult = MyLake_227_model_evaluation(m_start, m_stop, sediment_params, lake_params)
    run_INCA = 0; % 1- MyLake will run INCA, 0- No run
    use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files
    save_initial_conditions = false; % save final concentrations as initial for the next run
    run_ID = 0;
    clim_ID = 0;
    [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, save_initial_conditions); % runs the model and outputs obs and sim
    
    ModelResult.O24m = transpose(MyLake_results.basin1.concentrations.O2(8,:)/1000);
 %   ModelResult.O26m = transpose(MyLake_results.basin1.concentrations.O2(12,:)/1000);
 %   ModelResult.O28m = transpose(MyLake_results.basin1.concentrations.O2(16,:)/1000);
    ModelResult.O210m = transpose(MyLake_results.basin1.concentrations.O2(20,:)/1000);
end

% Error function. The error function takes a ModelResult
% struct and a Data struct. Remember that these two structs are user
% defined. They should contain whatever is needed for the error function to
% compare the model result to measured data. It has to return a positive
% number err, which is smaller the better fit the model is to the data.
function err = error_function_oxygen(ModelResult, Data)
    
    MatchedModelO24m = ModelResult.O24m(Data.date_mask);
   % MatchedModelO26m = ModelResult.O26m(Data.date_mask);
   % MatchedModelO28m = ModelResult.O28m(Data.date_mask);
    MatchedModelO210m = ModelResult.O210m(Data.date_mask);
%   err = nansum (((MatchedModelO24m - Data.O24m).^2) + ((MatchedModelO26m - Data.O26m).^2)+ ((MatchedModelO28m - Data.O28m).^2) + ((MatchedModelO210m - Data.O210m).^2));
   err = nansum (((MatchedModelO24m - Data.O24m).^2) + ((MatchedModelO210m - Data.O210m).^2)); % TOOK OUT 6 AND 8 M DEPTH
end

%% END project specific evaluation functions

%% The following two functions are general and should not have to be modified in each project

function do_MyLake_optimization(m_start, m_stop, K_values_sediment, K_values_lake, Data, varyindexes, minparam, maxparam, initial_guess, ...
    modeleval, errfun, max_generations, population_size, parallelize, filenameprefix)

    runfunc = @(varyparam)(MyLake_optimizer_single_run(m_start, m_stop, K_values_sediment, K_values_lake, Data, varyindexes, modeleval, errfun, varyparam));
    
    options = optimoptions('ga', 'MaxGenerations', max_generations, 'PopulationSize', population_size, 'UseParallel', parallelize, 'InitialPopulationMatrix', initial_guess);

    % Running the genetic optimizer algorithm
    [optimal_parameters, optimal_ss, exitflag, output, pop, scores] = ga(runfunc, size(varyindexes,2),...
        [], [], [], [], minparam, maxparam,...
        [], [], options);

    cl = clock;
    filename = sprintf('%s_optimized_parameters_%d_%d_%d.mat', filenameprefix, cl(3), cl(2), cl(1));
    save(filename, 'optimal_parameters', 'varyindexes');
end

function err = MyLake_optimizer_single_run(m_start, m_stop, K_sediments, K_lake, Data, varyindexes, modeleval, errfun, varyparam)

    % Inserting the varying parameters into K_lake and K_sediments
    for ii = 1:size(varyindexes,2)
        if varyindexes(1,ii) <= length(K_lake)
            for jj = 1:size(varyindexes, 1)
                if ~isnan(varyindexes(jj, ii))
                    K_lake{varyindexes(jj,ii), 1} = varyparam(ii);
                end
            end
        else
            for jj = 1:size(varyindexes, 1)
                if ~isnan(varyindexes(jj, ii))
                    K_sediments{varyindexes(jj,ii)-length(K_lake), 1} = varyparam(ii);
                end
            end
        end
    end

    % Running the model
    ModelResult = modeleval(m_start, m_stop, K_sediments, K_lake);

    % Evaluating the error
    err = errfun(ModelResult, Data);
    
    % Debug output.
    nf = java.text.DecimalFormat;
    ssstr = char(nf.format(err));
    fprintf(1, '\n');
    fprintf(1, '*******************************************************************************************\n');
    fprintf(1, '                Single model run finished. Error: %s\n', ssstr);
    fprintf(1, 'Parameters in this run:');
    for ii = 1:length(varyparam)
    fprintf(1, ' %.3g', varyparam(ii));
    end
    fprintf('\n');
    fprintf(1, '*******************************************************************************************\n');
    fprintf(1, '\n');
end