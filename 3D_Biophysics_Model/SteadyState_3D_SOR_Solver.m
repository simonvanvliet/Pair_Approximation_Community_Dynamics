function [output, gridE1, gridE2] = SteadyState_3D_SOR_Solver(settings, input, gridE1, gridE2, gridCellTypeScaled)
% Core solver, called from all other codes that solve 2D steady state problem
% Calculates steady state solution for 3D model of amino acid exchange
% This is the main handler, the actual computation is done by SteadyState_2D_SOR_Solver_2018_CORE
% This last function also has a c++ implementation for increased performance
%% settings: model parameters
%% input: initial guess of solution
%% output: steady state solution,
%
% two cell types:
%   type 0 produces AA1, growth limited by AA2
%   type 1 produces AA2, growth limited by AA1
%   type 0 produces pro, growth limited by try: Delta try on left
%   type 1 produces try, growth limited by pro: Delta pro on right
% internal concentration of produced AA is set to Iconst
% AA uptake is linear with [E], leakage is passive diffusion
% uptake, leakage, and diffusion rates can be assymetric between AA
% effective diffusion constant  =  (1-p)/(1+p/2), where p is density
% time scaled with growth rate (1h doubling time)
% space scaled with cell size (settings.cellSpacing um)
% concentrations scaled with Km of growth for each AA
%
% Boundary conditions: no flux on left, right, and top wall, [E] = 0 on bottom wall
%
% Solved with iterative finite difference scheme using successive over-relaxation (SOR) method
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020


%% get parameters model
Iconst1 = settings.Iconst1; %COncentration of produced amino-acid
Iconst2 = settings.Iconst2; %COncentration of produced amino-acid

D = (settings.D0/(settings.cellSpacing)^2)*(1-settings.rho)./(1+settings.rho/2); %effective diffusion constant
alpha = settings.rho/(1-settings.rho);

delta_u = settings.delta_u; %ru2/ru1  =  delta_u
delta_l = settings.delta_l; %rl2/rl1  =  delta_u
delta_D = settings.delta_D; %D2/D1  =  delta_u

ru = settings.ru; %uptake rate
rl = settings.rl; %leakage rate

%% get grid properties

%Grid size [Width (East-West), Depth (North-South), Height (Top-Bottom)]
    
%boundaryType: 
% 0 is closed boundary (zero-flux)
% 1 os open boundary (zero concentration)
% order is: North, East, South, West, Top, Bottom

boundaryType            = settings.boundaryType;  
gridSize                = settings.gridSize; %[W,D,H];
gridScaling             = settings.gridScaling; %scale grid by this factor to make diffusion process converge

%% get properties SOR process
tolerance = settings.tolerance; %the max relative difference between previous and current solution should be less than tolerance for each grid site
omega = settings.omega; 


%% calculate parameter relations
ru1 = ru / sqrt(delta_u);
ru2 = ru * sqrt(delta_u);
rl1 = rl / sqrt(delta_l);
rl2 = rl * sqrt(delta_l);
D1  = D / sqrt(delta_D);
D2  = D * sqrt(delta_D);


%% Define production functions
% *(1-type): selects type 1, kills type 2
% *(type):   selects type 2, kills type 1
IofE = @(E,ru,rl,D) (E*(ru+rl)   - rl + sqrt( (ru+rl)^2*E.^2 + 2*(rl+2)*(ru+rl)*E + rl^2 )) / (2*(1+rl));
muOfE = @(E1,E2,type) (1-type).* IofE(E2,ru2,rl2,D2)./(1+IofE(E2,ru2,rl2,D2)) + (type).*IofE(E1,ru1,rl1,D1)./(1+IofE(E1,ru1,rl1,D1));


%% get input grid
gridCellType = input.gridCellType;

%%Calc grid size
if gridSize(3)==1 %2D
    gridScaleVec = [gridScaling, gridScaling, 1];
    gridExtendVec = [2, 2, 0];
else %3D
    gridScaleVec = gridScaling*[1, 1, 1];
    gridExtendVec = [2, 2, 2];
end
    
gridSizeScaled          = gridSize .* gridScaleVec;
gridSizeExtended        = gridSizeScaled + gridExtendVec;
dx                      = 1/gridScaling; 


% %% check input grid
% if sum(size(gridE1) == gridSizeExpect) ~= expectedSize, error('E1 grid has wrong size'), end
% if sum(size(gridE2) == gridSizeExpect) ~= expectedSize, error('E2 grid has wrong size'), end
% if sum(size(gridCellTypeScaled) == gridSizeExpect) ~= expectedSize, disp(size(gridCellTypeScaled)), error('Cell Type Scaled grid has wrong size'),  end
% if sum(size(gridCellType) == gridSize) ~= expectedSize, error('Cell type grid has wrong size'), end


%extend grid to allow for boundary conditions
exGridE1 = zeros(gridSizeExtended);
exGridE2 = zeros(gridSizeExtended);
extendedGridType = zeros(gridSizeExtended);


if gridSize(3)==1 %2D
    exGridE1(2:end-1, 2:end-1) = gridE1;
    exGridE2(2:end-1, 2:end-1) = gridE2;
    extendedGridType(2:end-1, 2:end-1) = gridCellTypeScaled;
else %3D
    exGridE1(2:end-1, 2:end-1, 2:end-1) = gridE1;
    exGridE2(2:end-1, 2:end-1, 2:end-1) = gridE2;
    extendedGridType(2:end-1, 2:end-1, 2:end-1) = gridCellTypeScaled; 
end
    



%RUN SOR
%%
[exGridE1,exGridE2,numIt,errorCode] = SteadyState_3D_SOR_Solver_CoreC(...
    exGridE1,exGridE2, extendedGridType,...
    int32(gridSizeScaled), dx, boundaryType,...
    ru1 ,ru2, rl1, rl2, D1, D2, Iconst1, Iconst2, alpha,...
    tolerance, omega, settings.timeOutTime);

if numIt > 1000
    fprintf('num iterations = %i\n',numIt)
end
if errorCode == 1
    warning('Solution did not converge, skipping this grid\n'); 
    output = [];
else
    %remove padding of grid
    if gridSize(3)==1 %2D
        gridE1 = exGridE1(2:end-1, 2:end-1);
        gridE2 = exGridE2(2:end-1, 2:end-1);
    else %3D
        gridE1 = exGridE1(2:end-1, 2:end-1, 2:end-1);
        gridE2 = exGridE2(2:end-1, 2:end-1, 2:end-1);  
    end

    %calculate output
    if gridSize(3)==1 %2D
        gridE1_DownSampled = imresize(gridE1, 1/gridScaling, 'nearest');
        gridE2_DownSampled = imresize(gridE2, 1/gridScaling, 'nearest');
    else %3D
        gridE1_DownSampled = imresize3(gridE1, 1/gridScaling, 'nearest');
        gridE2_DownSampled = imresize3(gridE2, 1/gridScaling, 'nearest');
    end

    mu = muOfE(gridE1_DownSampled, gridE2_DownSampled, gridCellType);
    I1 = IofE(gridE1_DownSampled,ru1,rl1,D1);
    I2 = IofE(gridE2_DownSampled,ru2,rl2,D2);

    output = struct();
    output.gridE1_DownSampled = gridE1_DownSampled;
    output.gridE2_DownSampled = gridE2_DownSampled;
    output.I1 = I1;
    output.I2 = I2;
    output.mu = mu;
    output.gridCellType = gridCellType;
    output.settings = settings;
    
end


