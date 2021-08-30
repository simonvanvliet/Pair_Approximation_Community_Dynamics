% Plot grid arrangement 3D Model
%
% Written by Simon van Vliet & Alma Dal Co
% Eawag & ETH Zurich
%
% Initial development: 12.6.2018
% Last update: 2021.01.30 

% load default settings
settings = SteadyState_3D_Settings;

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

settings3D = settings;
settings3D.gridSize              = [40, 40, 40]; %[Width (East-West), Depth (North-South), Height (Top-Bottom)]
slice = 20;

%Set scan range
numChambers = 100; 
initFractions = [0.3 0.5 0.7];
numFrac = length(initFractions);

clusterSizes = 0:20:80;
numCluster = length(clusterSizes);


cY = [183 184 54]/256;
cP = [125 126 187]/256;

cmapYP = [cY;cP];

plotData = cell(numFrac, numCluster);

for ff = 1:numFrac
    for rr = 1:numCluster
        clusterR = clusterSizes(rr);
        %RUN 3D
        locSet3D = settings3D;
        locSet3D.numReplacementInitGrid = clusterR;

        %make grid
        grid3D = SteadyState_3D_ClusteredGrid(locSet3D, clusterR, initFractions(ff));
        
        gridSlice = squeeze(grid3D(:,:,slice));
        plotData{ff,rr} = gridSlice;
    end
end

%%

figure(101)
for ff = 1:numFrac
    for rr = 1:numCluster
        data = plotData{ff,rr};
        idx = (ff-1) * numCluster + rr;
        subplot(numFrac, numCluster, idx)
        imagesc(data)
        colormap(cmapYP)
        axis off
    end
end
     
