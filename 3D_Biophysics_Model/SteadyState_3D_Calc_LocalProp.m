function [muDTrp, fracDPForDTAtRadius, muDPro, fracDTForDPAtRadius] = ...
    SteadyState_3D_Calc_LocalProp(muGrid, typeGrid, gridSize, edge, radiusVec)%#codegen
%Extracts for all cells growth rate and spatial arrangement
%INPUT: muGrid: grid of growth rates
%INPUT: typeGrid: grid of cell types
%INPUT: edge: number of pixels on edge of chamber where cells are excluded from analyis
%INPUT: radiusVec: vector of radia (in pixels) for which spatial arrangement is quantified
%
%Output: muDTrp/muDPro: vector of growth rates for each DTrp or Dpro cell
%Output: fracDPForDTAtRadius/fracDTForDPAtRadius: matrix #cell*#radia
%... each row gives vector of fraction of other type within neighborhood as function of radius 
%
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020

%SET DEBUG

DTrpIndex = 0; %type 0 is Dtry
DProIndex = 1; %type 1 is Dpro

[xGrid, yGrid, zGrid] = meshgrid(1:gridSize(2), 1:gridSize(1), 1:gridSize(3));

if gridSize(3) == 1 %2D
    numZ = 1;
    zStart = 1;
    zEnd = 1;
else
    numZ = gridSize(3) - 2*edge;
    zStart = edge+1;
    zEnd = gridSize(3)-edge;
end

numInteriorPoints = (gridSize(1) - 2*edge) ...
    * (gridSize(2) - 2*edge) * numZ;

numRadia = length(radiusVec);

%init output
muVec = nan(numInteriorPoints, 1);
typeVec = nan(numInteriorPoints, 1);
DProWithinRadius = nan(numInteriorPoints, numRadia);
DTrpWithinRadius = nan(numInteriorPoints, numRadia);

%loop over all interior grid points
currIndex = 0;
for row = edge+1:gridSize(1)-edge
    for col = edge+1:gridSize(2)-edge
        for lay = zStart:zEnd
            currIndex = currIndex+1;

            %extract current type and growth rate
            currMu = muGrid(row, col, lay);
            currType = typeGrid(row, col, lay);

            muVec(currIndex) = currMu;
            typeVec(currIndex) = currType;

            %calculate distance to all other grid points
            currX = xGrid(row, col, lay);
            currY = yGrid(row, col, lay);
            currZ = zGrid(row, col, lay);
            
            rMat = sqrt((double(xGrid - currX)).^2 + ...
                (double(yGrid - currY)).^2 +...
                (double(zGrid - currZ)).^2);

            %extract distances to cells of certain type
            rMatDTrp = rMat(typeGrid == DTrpIndex);
            rMatDPro = rMat(typeGrid == DProIndex);

            %scan radia and store data
            for rr = 1:numRadia
                upperR = radiusVec(rr);
                DProWithinRadius(currIndex,rr) = sum(rMatDPro <= upperR);
                DTrpWithinRadius(currIndex,rr) = sum(rMatDTrp <= upperR);
            end
        end
    end
end

%split data into cell types and calculate fraction for each radius
muDTrp = muVec(typeVec == DTrpIndex);
fracDPForDTAtRadius = DProWithinRadius(typeVec == DTrpIndex,:) ./ ...
    (DProWithinRadius(typeVec == DTrpIndex,:) + DTrpWithinRadius(typeVec == DTrpIndex,:));

muDPro = muVec(typeVec == DProIndex);
fracDTForDPAtRadius = DTrpWithinRadius(typeVec == DProIndex,:) ./ ...
    (DProWithinRadius(typeVec == DProIndex,:) + DTrpWithinRadius(typeVec == DProIndex,:));


    
    