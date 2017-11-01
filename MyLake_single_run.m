tic
disp('Started at:')
disp(datetime('now'));


run_INCA = 0; % 1- MyLake will run INCA, 0- No run
use_INCA = 0; % 1- MyLake will take written INCA input, either written just now or saved before, and prepare inputs from them. 0- MyLake uses hand-made input files

m_start=[1969, 6, 27]; %
m_stop=[2009, 12, 31]; %

save_initial_conditions = false; % save final concentrations as initial for the next run

[lake_params, sediment_params] = load_params();


run_ID = 0;
clim_ID = 0;
[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params, lake_params, use_INCA, run_INCA, run_ID, clim_ID, save_initial_conditions); % runs the model and outputs obs and sim


disp('Saving results...')
save('IO/MyLakeResults_sediment.mat', 'MyLake_results', 'Sediment_results')
disp('Finished at:')
disp(datetime('now'));


toc
