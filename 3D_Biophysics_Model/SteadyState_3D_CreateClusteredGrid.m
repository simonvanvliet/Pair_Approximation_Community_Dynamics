function gridCellType = SteadyState_3D_CreateClusteredGrid(settings)
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020

%% get grid properties
gridSize = settings.gridSize; %size of grid in terms of cells (1x1 blocks)

randNumFeed = settings.randNumFeedInitGrid;
initFracType1 = settings.initFracType1;

%% init grid
gridCellType = zeros(gridSize); %type: 0 produces AA1, 1 produces AA2

%% seed randomly
rng(randNumFeed)
randMatrix = rand(gridSize);
gridCellType(randMatrix < initFracType1) = 0;
gridCellType(randMatrix >= initFracType1) = 1;

%% cerate clusters
if settings.numReplacementInitGrid>0
    %perfrom settings.numReplacementInitGrid rounds of replacement where each grid site is visisted
    numIt = ceil(settings.numReplacementInitGrid * prod(gridSize));
    
    %get rand numbers
    pickRand = 1;
    while pickRand %check if valid random number, exclude 0 and 1
        randNumVec = rand(numIt,4);
        if sum(randNumVec(:) == 0) == 0 && sum(randNumVec(:) == 1) == 0
            pickRand = 0;
        end
    end
    
    if gridSize(3)==1 %2D
        tresholdVec = (1:6)/4; %excludes movements up and down
    else %3D
        tresholdVec = (1:6)/6;
    end
    movementVec = ['l','r','t','b','u','d'];
    
    %select random focal cell, select random neighbor and replace neighbor with type of focal cell
    %bias can be used to bias formation of clusters in direction towards chamber opening to create stripes
    for nn = 1:numIt
        %pick random cell
        rowCellToDevide = ceil(randNumVec(nn,1)*gridSize(1));
        colCellToDevide = ceil(randNumVec(nn,2)*gridSize(2));
        layCellToDevide = ceil(randNumVec(nn,3)*gridSize(3));
        
        %pick random direction of movement
        idxVec = (randNumVec(nn,4) <= tresholdVec) & (randNumVec(nn,4) > [0 tresholdVec(1:end-1)]);
        offspringSite = movementVec(idxVec);
        
        %move offspring
        offspringRow = rowCellToDevide;
        offspringCol = colCellToDevide;
        offspringLay = layCellToDevide;
        if offspringSite == 't', offspringRow = min(offspringRow+1, gridSize(1)); end
        if offspringSite == 'b', offspringRow = max(offspringRow-1, 1); end
        if offspringSite == 'r', offspringCol = min(offspringCol+1, gridSize(2)); end
        if offspringSite == 'l', offspringCol = max(offspringCol-1, 1); end
        if offspringSite == 'u', offspringCol = min(offspringLay+1, gridSize(3)); end
        if offspringSite == 'd', offspringCol = max(offspringLay-1, 1); end
        
        gridCellType(offspringRow, offspringCol, offspringLay) = ...
            gridCellType(rowCellToDevide, colCellToDevide, layCellToDevide);
    end
end