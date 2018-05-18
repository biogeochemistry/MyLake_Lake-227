%% load_params: load parameters for MyLake and Sediment
function [lake_params, sediment_params] = load_params()

lake_params = {
    % PhysPar
    0.5, 'dz',                 % 1
    NaN, 'Kz_K1',     % 2     open water diffusion parameter (-)
    0.000898, 'Kz_K1_ice',     % 3     under ice diffusion parameter (-)
    7E-05, 'Kz_N0',            % 4     min. stability frequency (s-2)
    NaN, 'C_shelter',         % 5     wind shelter parameter (-)
    49.661944, 'lat',              % 6     latitude (decimal degrees)
    93.721944, 'lon',              % 7     longitude (decimal degrees)
    0.3, 'alb_melt_ice',       % 8     albedo of melting ice (-)
    0.4, 'alb_melt_snow',     % 9     albedo of melting snow (-)
    0.000228, 'PAR_sat',       % 10    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    0.45, 'f_par',             % 11    Fraction of PAR in incoming solar radiation (-)
    0.015, 'beta_chl',         % 12    Optical cross_section of chlorophyll (m2 mg-1)
    5, 'lambda_i',             % 13    PAR light attenuation coefficient for ice (m-1)
    15, 'lambda_s',            % 14    PAR light attenuation coefficient for snow (m-1)
    0.05, 'F_sed_sld',         % 15    volume fraction of solids in sediment (= 1-porosity)
    0.7359, 'I_scV',                % 16    scaling factor for inflow volume (-)
    0, 'I_scT',                % 17    adjusting delta for inflow temperature (-)
    1, 'I_scC',                % 18    scaling factor for inflow concentration of C (-)
    1, 'I_scPOC',              % 19    scaling factor for inflow concentration of POC (-)
    1, 'I_scTP',               % 20    scaling factor for inflow concentration of total P (-)
    1, 'I_scDOP',              % 21    scaling factor for inflow concentration of diss. organic P (-)
    1, 'I_scChl',              % 22    scaling factor for inflow concentration of Chl a (-)
    1, 'I_scDOC',              % 23    scaling factor for inflow concentration of DOC  (-)
    1, 'I_scPOP',              % 24    scaling factor for inflow concentration of POP  (-)
    1, 'I_scO',                % 25    Scaling factor for inflow concentration of O2 (-)
    1 , 'I_scDIC',             % 26    Scaling factor for inflow concentration of DOC  (-)
    1,  'I_scNO3',             % 27    Scaling factor for inflow concentration of NO3 (-)
    1,  'I_scNH4',             % 28    Scaling factor for inflow concentration of NH4 (-)
    1,  'I_scSO4',             % 29    Scaling factor for inflow concentration of SO4 (-)
    1,  'I_scFe2',             % 30    Scaling factor for inflow concentration of Fe2 (-)
    1,  'I_scCa2',             % 31    Scaling factor for inflow concentration of Ca2 (-)
    1,  'I_scpH',              % 32    Scaling factor for inflow concentration of pH (-)
    1,  'I_scCH4',             % 33    Scaling factor for inflow concentration of CH4 (-)
    1,  'I_scFe3',             % 34    Scaling factor for inflow concentration of Fe3 (-)
    1,  'I_scAl3',             % 35    Scaling factor for inflow concentration of Al3 (-)
    1,  'I_scFeS',             % 36    Scaling factor for inflow concentration of FeS (-)
    1,  'I_scCaCO3',           % 37    Scaling factor for inflow concentration of CaCO3 (-)
    1,  'I_scCH4g',            % 38    Scaling factor for inflow concentration of CH4g (-)
    2.5, 'swa_b0',             % 39     non-PAR light attenuation coeff. (m-1)
    1.05, 'swa_b1',            % 40     PAR light attenuation coeff. (m-1)
    3.30E-07, 'S_res_epi',     % 41     Particle resuspension mass transfer coefficient, epilimnion (m day-1, dry)
    3.30E-08, 'S_res_hypo',    % 42     Particle resuspension mass transfer coefficient, hypolimnion (m day-1, dry)
    0.03, 'H_sed',             % 43     height of active sediment layer (m, wet mass)
    2500, 'Psat_L',              % 44     NOTE: NOT USED: Half saturation parameter for Langmuir isotherm
    8000, 'Fmax_L',              % 45     NOTE: NOT USED: Scaling parameter for Langmuir isotherm !!!!!!!!!!!!
    0.1, 'w_s',                % 46     settling velocity for S (m day-1)
    0.0817, 'w_chl',             % 47     settling velocity for Chl a (m day-1)
    1, 'Y_cp',                 % 48     NOTE: NOT USED:  yield coefficient (chlorophyll to carbon) * (carbon to phosphorus) ratio (-)   1/55*112/1 = 1
    0.0375, 'm_twty',             % 49    loss rate (1/day) at 20 deg C
    1.85, 'g_twty',             % 50    specific growth rate (1/day) at 20 deg C
    2.00E-04, 'k_twty',        % 51    NOTE: NOT USED: specific Chl a to P transformation rate (1/day) at 20 deg C
    0, 'dop_twty',             % 52    NOTE: NOT USED: specific DOP to P transformation rate (day-1) at 20 deg C
    0.598, 'P_half',             % 53    Half saturation growth P level (mg/m3)
    0.000228, 'PAR_sat_2',     % 54    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    0.015, 'beta_chl_2',       % 55    Optical cross_section of chlorophyll (m2 mg-1)
    0.0817, 'w_chl_2',           % 56    Settling velocity for Chl a (m day-1)
    0.0375, 'm_twty_2',           % 57    Loss rate (1/day) at 20 deg C
    1.85, 'g_twty_2',           % 58    Specific growth rate (1/day) at 20 deg C
    0.598, 'P_half_2',         % 59    Half saturation growth P level (mg/m3)
    0.01, 'oc_DOC',            % 60    Optical cross-section of DOC (m2/mg DOC)
    0.1, 'qy_DOC',             % 61    Quantum yield (mg DOC degraded/mol quanta)
    0.01, 'k_BOD',              % 62    NOTE: NOT USED: Organic decomposition rate (1/d)
    5, 'w_CH4',                % 63    Methane gas rising velocity (m/d)
    1.047, 'theta_bod',        % 64    NOTE: NOT USED: Temperature adjustment coefficient for BOD, T ? 10 °C
    1.13, 'theta_bod_ice',     % 65    NOTE: NOT USED: Temperature adjustment coefficient for BOD, T < 10 °C
    1, 'theta_sod',            % 66    NOTE: NOT USED: Temperature adjustment coefficient for SOD, T ? 10 °C
    1, 'theta_sod_ice',        % 67    NOTE: NOT USED: Temperature adjustment coefficient for SOD, T < 10 °C
    4, 'BOD_temp_switch',      % 68    NOTE: NOT USED: Threshold for bod or bod_ice °C
    7.5, 'pH',                 % 69    Lake water pH
    2, 'Q10_wc',               % 70    Q10 for reactions of respiration
    1, 'wc_factor',            % 71    Scaling factor for rates in WC
    4.8497, 'T_ref_wc'};       % 72    Reference Temperature for rates


sediment_params = {
    0.4,   'k_Chl',                 % 1       % 1
    0.4,  'k_POP',                 % 2       % 1
    1.75, 'k_POC',                  % 3       % 0.01
    1.04,  'k_DOP',                 % 4       % 1
    0.02, 'k_DOC',                  % 5       % 1
    0.031, 'Km_O2',                 % 6       % Canavan, R. W (2006) rho=2.5
    0.01,  'Km_NO3',                % 7       % Canavan, R. W (2006) rho=2.5
    3.92,  'Km_Fe(OH)3',         % 8       % Canavan, R. W (2006) rho=2.5
    2415,  'Km_FeOOH',           % 9       % Canavan, R. W (2006) rho=2.5
    0.0293,  'Km_SO4',                 % 10       % Canavan, R. W (2006 rho=2.5
    0.001,'Km_oxao',                % 11       % the same as Km rho=2.5
    0.1,  'Km_amao',                % 12       % the same as Km rho=2.5
    0.329, 'Kin_O2',                % 13       % the same as Km rho=2.5
    0.01,  'Kin_NO3',               % 14       % the same as Km rho=2.5
    0.1,   'Kin_FeOH3',         % 15       % the same as Km rho=2.5
    0.1,   'Kin_FeOOH',         % 16       % the same as Km rho=2.5
    20000,    'k_amox',                % 17       % Canavan, R. W (2006)
    87000,  'k_Feox',                % 18       % Canavan, R. W (2006)
    0.1,   'k_Sdis',                % 19       %
    2500,  'k_Spre',                % 20       %
    3.17,   'k_FeS2pre',             % 21       % Canavan (2006)
    0.1,   'k_alum',                % 22
    1.35,     'k_pdesorb_a',         % 23
    1.35,     'k_pdesorb_b',         % 24
    20000,  'k_fesox',              % 25        % R23 %Canava
    1000,   'k_fes2ox',             % 26        % R23 % Katsev (2013)
    8,      'k_tS_Fe',              % 27      % Cappellen (1996) in Canavan, R. W (2006) the reaction is different
    2510,  'Ks_FeS',                % 28      % Canavan, R. W (2006)
    0.001, 'k_Fe_dis',              % 29      % Canavan, R. W (2006), Katsev, R. W (2013)
    21.3,'k_Fe_pre',             % 30         % Katsev, R. W (2013)
    0.37e-3,  'k_apa_pre',          % 31
    0.37,     'k_apa_dis',          % 32
    10^(-4.2249),  'K_apa',         % 33      % linl.dat PHREEQC
    0.1/2.5,  'k_CaCO3_pre',        % 34      % Katsev (2013)
    0.05,  'k_CaCO3_dis',           % 35      % Katsev (2013)
    5e-9,  'K_CaCO3',               % 36      %
    450/2.5,  'k_FeCO3_pre',        % 37      % Cappellen (1996)
    0.25,     'k_FeCO3_dis',        % 38      % Cappellen (1996)
    10^(-8.4),  'K_FeCO3',          % 39      % Cappellen (1996)
    0.37e-3,  'k_viv_pre',          % 40
    0.37,  'k_viv_dis',             % 41
    10^(-4.7237), 'K_viv',          % 42     % linl.dat PHREEQC
    1e-6,  'k_oms',                 % 43
    1e4,   'k_tsox',                % 44     % Canavan, R. W (2006)
    0.001, 'k_FeSpre',            % 45     % from "Non-steady state diagenesis of organic and inorganic sulfur in lake sediments Raoul-Marie Couture, Rachele Fischer b, Philippe Van Cappellen b, Charles Gobeil c
    1e7,   'k_ch4_o2',              % 46     % Canavan, R. W (2006)
    1e-1,  'k_ch4_so4',             % 47     % Canavan, R. W (2006)
    0.0015,  'Kh_CH4',              % 48     % Henry cobstant M/atm
    1e3,   'k_ch4_dis',             % 49
    5,     'w_CH4g',                % 50     % Rising velocity of methane
    0.034,  'Kh_CO2',               % 51     % Henry cobstant M/atm
    21.5,   'accel',                 % 52
    1e-6,   'f_pfe',                % 53
    1.35,   'k_pdesorb_c',          % 54
    0.95,   'fi_in',                % 55
    0.85,   'fi_f',                 % 56
    0.5,    'X_b',                  % 57
    1,      'tortuosity',           % 58
    0.1,    'w',                    % 59
    300,    'n',                    % 60
    30,     'depth',                % 61
    14.4/2,   'alfa0',                % 62
    236,    'Cx1',                  % 63           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    20,     'Ny1',                  % 64           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Pz1',                  % 65           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    236,    'Cx2',                  % 66           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    20,     'Ny2',                  % 67           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Pz2',                  % 68           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Cx3',                  % 69           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    0.1,    'Ny3',                  % 70           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    0,      'Pz3',                  % 71           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    9,     'effective_depth',     % 72           % depth below which the lake is affected by sediments, [m], if -1 (experimental) , then sediments below pycnocline
    96,     'n_ts',                 % 73           % (48 is the minimum, 24 for calibration) number of time steps during 1 day (fixed time step of MyLake) for chemical and sediment module (the modules should be in sync)
    0,      'pH algorithm',         % 74           % 0. Disabled  % 1. Phreeqc  % 2. Electro-neutrality Equation
    0, 'SO4 fx',
    };
