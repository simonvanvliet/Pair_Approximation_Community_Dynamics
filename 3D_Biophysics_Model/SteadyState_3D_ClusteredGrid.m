function typeGrid = SteadyState_3D_ClusteredGrid(settings, clusterSize, initFrac)
%creates 4D array with stack of artificial cell type arrangements
%INPUT: settings, model settings
%INPUT: numChambers, number of chambers to create
%INPUT: clusterSize, size of clusters to create
%
%OUTPUT: 4D array with stack of artificial cell type arrangements
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020

%SET GRID SETTINGS
settings.initialGridSetting = 'RandomClustered';
settings.numReplacementInitGrid = clusterSize;

%Generate grids
goodGrid = false;
settings.initFracType1 = initFrac;
settings.randNumFeedInitGrid = now; %give each chamber different rand feed

while ~goodGrid
    typeGrid = SteadyState_3D_CreateClusteredGrid(settings);
    if sum(typeGrid(:)==0)>0 && sum(typeGrid(:)==1)>0 %check that both cell types are there
        goodGrid=true;
    end
end



        
        
        