function [output] = SteadyState_3D_InitGrid(settings, initialGridSetting, varargin)
% Initializes grid for 3D steady state model of amino acid exchange
%% Output: initialized grid
%% Settings: model parameters
%% initialGridSetting: type of initial grid, Options are:
%    SplitWorld: left half = type1, right half = type2
%    LineWorld1: all type 0, except for center line which is type 1
%    LineWorld0: all type 1, except for center line which is type 0
%    Type0Island: all type 1, except for center point which is type 1
%    Type1Island: all type 0, except for center point which is type 1
%    Random: random assignment with initFracType1 assigned to type 1
%    RandomClustered: creates random clusters: start with random and create clusters by iteratively flipping neighboring
%    sides to same state.
%    ProvidedGrid: provide external grid in varargin
%
% type 0 produces pro, growth limited by try
% type 1 produces try, growth limited by pro
% Boundary conditions: no flux on left, right, and top wall, [E] = 0 on bottom wall
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020

%% get model parameters
Iconst1 = settings.Iconst1; %COncentration of produced amino-acid
Iconst2 = settings.Iconst2; %COncentration of produced amino-acid

delta_u = settings.delta_u; %ru2/ru1  =  delta_u
delta_l = settings.delta_l; %rl2/rl1  =  delta_u

ru = settings.ru; %uptake rate
rl = settings.rl; %leakage rate
ru_1 = ru / sqrt(delta_u);
ru_2 = ru * sqrt(delta_u);
rl_1 = rl / sqrt(delta_l);
rl_2 = rl * sqrt(delta_l);



%% get grid properties
gridSize = settings.gridSize; %size of grid in terms of cells (1x1 blocks)
gridScaling = settings.gridScaling; %scale grid by this factor to make diffusion process converge

randNumFeed = settings.randNumFeedInitGrid;
initFracType1 = settings.initFracType1;

%% init grid
gridCellType = zeros(gridSize); %type: 0 produces AA1, 1 produces AA2

%% setup initial grid
switch initialGridSetting
    case 'SplitWorld'
        midPoint = floor(gridSize(1) / 2);
        gridCellType(:, midPoint+1:end, :) = 1;
    case 'Type0Island'
        gridCellType(:) = 1;
        midPoint = num2cell(ceil(gridSize / 2));
        gridCellType(midPoint{:}) = 0;        
    case 'Type1Island'
        gridCellType(:) = 0;
        midPoint = num2cell(ceil(gridSize / 2));
        gridCellType(midPoint{:}) = 1; 
    case 'Random'
        rng(randNumFeed)
        randMatrix = rand(gridSize);
        gridCellType(randMatrix < initFracType1) = 0;
        gridCellType(randMatrix >= initFracType1) = 1;
    case 'RandomClustered'
        gridCellType = SteadyState_3D_CreateClusteredGrid(settings);
    case 'ProvidedGrid2D'
        %check for input grid validity
        if isempty(varargin), error('no input grid provided'); 
        else, gridIn = varargin{1}; end
        if ~ismatrix(gridIn), error('input is not 2D matrix'); end
        if sum(size(gridIn) == [gridSizeCells gridSizeCells]) ~= 2, error('input is of wrong size'); end
        if (sum(gridIn(:) == 0)+sum(gridIn(:) == 1)) ~= gridSizeCells^2, error('input contains invalid numbers'); end
        gridCellType = gridIn;
    case 'ProvidedGrid'
        %check for input grid validity
        if isempty(varargin), error('no input grid provided'); 
        else, gridIn = varargin{1}; end
        if ~isnumeric(gridIn), error('input is not an array'); end
        if mean(size(gridIn) == gridSize(1:length(size(gridIn)))) ~= 1, error('input is of wrong size'); end
        if (sum(gridIn(:) == 0) + sum(gridIn(:) == 1)) ~= prod(gridSize), error('input contains invalid numbers'); end
        gridCellType = gridIn;
end

%% Guess initial solution
%scale type grid up to full grid
if settings.gridSize(3) == 1 %2D
    gridCellTypeScaled = round(imresize(gridCellType, gridScaling, 'nearest'));
else %3D
    gridCellTypeScaled = round(imresize3(gridCellType, gridScaling, 'nearest'));
end


%set initial E values
%in producing regions, set [E] to expected homogenous steady state value Iconst*rl_1/ru_1
%in non-producing regions, set [E] to expected homogenous steady state value devided by scaleFactorE
scaleFactorE = 100; %guess for difference in [E] between producing and non producing regions
E1Init = (1-gridCellTypeScaled)*Iconst1*rl_1/ru_1 +     gridCellTypeScaled*Iconst1*rl_1/(ru_1*scaleFactorE);
E2Init =     gridCellTypeScaled*Iconst2*rl_2/ru_2 + (1-gridCellTypeScaled)*Iconst2*rl_2/(ru_2*scaleFactorE);


%store output
output.gridE1 = E1Init;
output.gridE2 = E2Init;
output.gridCellType = gridCellType;
output.gridCellTypeScaled = gridCellTypeScaled;
output.settings = settings;

