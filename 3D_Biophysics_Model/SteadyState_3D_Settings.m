function settings=SteadyState_3D_Settings
% Default model and solver settings for 3D steady state model of amino acid exchange
%% Settings: model parameters 
%
% AA1=pro, AA2=trp
% Type0 = DTry,  Type1 = DPro
% type 0 produces pro, growth limited by try
% type 1 produces try, growth limited by pro
% Boundary conditions: no flux on left, right, and top wall, [E]=0 on bottom wall
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020

lastUpdate='2020.02.06';
fprintf('Settings based on Dal Co et al, Nat Eco Evo 2020\n');

settings=struct;

%td_ouble

%non-scaled rates in sec
tD_s = 2790; %s (46.5min, mu=1.29)
ru_s = 7.0; % 1/s
rl_s = 3.1E-6; % 1/s
D0_s = 7.61E2; %um^2/s

%rescaled rates
ru = ru_s * tD_s;
rl = rl_s * tD_s;
D0 = D0_s * tD_s;

%% SET Initial grid settings
settings.initFracType1          = 0.5; %intial fraction of type 1 in random arrangement
settings.numReplacementInitGrid = 10;
settings.randNumFeedInitGrid    = 2018;


%% SET  model parameters
settings.Iconst1                = 20; %Concentration of produced amino-acid 
settings.Iconst2                = 20; %Concentration of produced amino-acid

settings.rho                    = 0.65; %Density
settings.D0                     = D0; %Diffusion constant in empty space
settings.cellSpacing            = 1.5; %Average cell spacing, used to rescale D0

settings.ru                     = ru;
settings.rl                     = rl; 

settings.delta_u                = 11.8; %ru2/ru1 = delta_u
settings.delta_l                = 1/26.3; %rl2/rl1 = delta_u
settings.delta_D                = 0.75; %D2/D1 = delta_u

settings.tD_s                   = tD_s;
settings.mu_wt                  = 3600 / tD_s;

fprintf('...rho=%#.2g ru=%#.2g rl=%#.2g d_D=%#.2g d_ru=%#.2g d_rl=%#.2g \n',...
    settings.rho,settings.ru,settings.rl,settings.delta_D,settings.delta_u,settings.delta_l);


%% SET grid properties
settings.gridSize               = [40, 40, 1]; %[Width (East-West), Depth (North-South), Height (Top-Bottom)]
settings.boundaryType           = [0, 0, 1, 0, 0, 0]; %Set boundary conditions
%boundaryType: 
% 0 is closed boundary (zero-flux)
% 1 is open boundary (zero concentration)
% order is: North, East, South, West, Top, Bottom

%settings.gridSizeCells          = 40; %size of grid in trems of cells (1x1um blocks)
settings.gridScalingBase        = 1; %scale grid by this factor to make diffusion process converge
settings.gridScaleFactor        = 2; %keep at 2 (set factor at which grid is refined)
settings.maxGridScaling         = 32;
settings.gridScaling            = settings.gridScalingBase;

%% SET properties of run
settings.DEBUG                  = 0;    % plot result at end of iterations
settings.DEBUG_Extensive        = 0;    % plot results during iterations
settings.timeOutTime            = 120;  % max time allowed to try to calc interaction range

%% SET properties SOR process
settings.tolerance              = 1e-3; %the max relative difference between previous and current solution should be less than tolerance for each grid site
settings.muTolerance            = 1e-2; %grid spacing is decreased untill mu profile changes by less than muTolerance (absolute difference)
settings.omega                  = 1.75; %hand picked omega
%NOTE: omega (between 0-2) sets relaxation of solution
%Values omega>1 are used to speed up convergence of a slow-converging process
%Values omega<1 are  used to help establish convergence of a diverging iterative process or speed up the convergence of an overshooting process.
%Use a value as large as possible without blowing up solution
%omega=2/(1+sin(pi/numGridPoint1D^2)); %optimal omega for Poisson process
