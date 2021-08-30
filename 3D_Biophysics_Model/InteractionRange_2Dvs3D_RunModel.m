% Run 3D Model store to disk 
% Process output with InteractionRange_2Dvs3D_ProcessModel
%
% Written by Simon van Vliet & Alma Dal Co
% Eawag & ETH Zurich
%
% Initial development: 12.6.2018
% Last update: 2021.01.30 

% load default settings
settings = SteadyState_3D_Settings;
pathName = './ClusteredGridScan_20210202/';

% create 3D Grid 
%Set boundary conditions [North, East, South, West, Top, Bottom]
% 0 is closed boundary (zero-flux) / % 1 is open boundary (zero concentration)
settings.boundaryType            = [0, 0, 0, 0, 0, 0]; 
settings.omega                   = 1.6;

%set symmetric community
settings.ru = 1E4;
settings.rl = 1E-2;
settings.delta_u = 1;
settings.delta_l = 1;
settings.delta_D = 1;
settings.D0 = 2E6;

settings.randNumFeedInitGrid = 20210202;

settings2D = settings;
settings3D = settings;

settings2D.maxGridScaling          = 16; 
settings2D.timeOutTime             = 180;
settings3D.maxGridScaling          = 8; 
settings3D.timeOutTime             = 900;

settings2D.gridSize              = [40, 40, 1]; %[Width (East-West), Depth (North-South), Height (Top-Bottom)]
settings3D.gridSize              = [40, 40, 40]; %[Width (East-West), Depth (North-South), Height (Top-Bottom)]

%Set scan range
numChambers = 100; 
initFractions = linspace(0.2, 0.8, numChambers);
initFractions = initFractions(randperm(numChambers));

clusterSizes = 0:5:80;
numCluster = length(clusterSizes);

save([pathName, 'metaData.mat'],'settings2D','settings3D','initFractions','clusterSizes');


fprintf('Settings: rho=%#.2g ru=%#.2g rl=%#.2g d_D=%#.2g d_ru=%#.2g d_rl=%#.2g \n',...
    settings.rho,settings.ru,settings.rl,settings.delta_D,settings.delta_u,settings.delta_l);


parfor cc = 1:numChambers
    modelSolutions = cell(numCluster, 2);
    initFrac = initFractions(cc);
    
    for rr = 1:numCluster
        clusterR = clusterSizes(rr);
        %RUN 3D
        locSet3D = settings3D;
        locSet3D.numReplacementInitGrid = clusterR;
    
        %make grid
        grid3D = SteadyState_3D_ClusteredGrid(locSet3D, clusterR, initFrac);
        %reset grid scaling
        locSet3D.gridScaling = locSet3D.gridScalingBase;
        %setup grid
        initialCondition = SteadyState_3D_InitGrid(locSet3D, 'ProvidedGrid', grid3D);
        %run simulation
        modelSolutions{rr, 1} = SteadyState_3D_RefineGrid(locSet3D, initialCondition);
    
        %RUN 2D
        locSet2D = settings2D;
        locSet2D.numReplacementInitGrid = clusterR;
    
        %make grid
        grid2D = SteadyState_3D_ClusteredGrid(locSet2D, clusterR, initFrac);
        %reset grid scaling
        locSet2D.gridScaling = locSet2D.gridScalingBase;
        %setup grid
        initialCondition = SteadyState_3D_InitGrid(locSet2D, 'ProvidedGrid', grid2D);
        %run simulation
        modelSolutions{rr, 2} = SteadyState_3D_RefineGrid(locSet2D, initialCondition);
       
    end
    saveName = sprintf('%sGridDataChamber%03i.mat', pathName, cc);
    saveOutput(saveName, modelSolutions);
end
    
%java.lang.System.gc()



function saveOutput(saveName, modelSolutions)
    save(saveName, 'modelSolutions','-v7.3')
end


