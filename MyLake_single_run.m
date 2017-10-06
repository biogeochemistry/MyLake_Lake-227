% for i=1:1000
tic
disp('Started at:')
disp(datetime('now'));


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

is_metrics = true; % print metrics in the end

m_start=[2000, 1, 1]; %
m_stop=[2000, 12, 31]; %

save_initial_conditions = false; % save final concentrations as initial for the next run
file_name = 'IO/test_C.mat'

[lake_params, sediment_params] = load_params();



% test C parameters are the same as for Chl
lake_params{56} = lake_params{47}; % 204.8121e-003  % 56    Settling velocity for Chl2 a (m day-1)
lake_params{57} = lake_params{49}; % 167.6746e-003   % 57    Loss rate (1/day) at 20 deg C
lake_params{58} = lake_params{50}; % 1.0985e+000   % 58    Specific growth rate (1/day) at 20 deg C
lake_params{59} = lake_params{53}; % 1.5525e+000   % 59    Half saturation growth P level (mg/m3)
lake_params{54} = lake_params{10};
lake_params{55} = lake_params{12};




run_ID = 0;
clim_ID = 0;
[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, save_initial_conditions); % runs the model and outputs obs and sim


disp('Saving results...')
save(file_name, 'MyLake_results', 'Sediment_results')
disp('Finished at:')
disp(datetime('now'));


toc
% end
%
