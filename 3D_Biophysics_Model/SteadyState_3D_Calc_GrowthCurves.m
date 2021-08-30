function [muDTrp,fracDPForDTAtRadius,muDPro,fracDTForDPAtRadius] = ...
    SteadyState_3D_Calc_GrowthCurves(muGrid, typeGrid, edge, radiusVec, gridSize)
%
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
DEBUG = 0;

DTrpIndex = 0; %type 0 is Dtry
DProIndex = 1; %type 1 is Dpro

colormap('parula')

[muVec, typeVec, DProWithinRadius, DTrpWithinRadius] = SteadyState_3D_Calc_LocalProp(muGrid, typeGrid, gridSize);


%split data into cell types and calculate fraction for each radius
muDTrp = muVec(typeVec == DTrpIndex);
fracDPForDTAtRadius = DProWithinRadius(typeVec == DTrpIndex,:) ./ ...
    (DProWithinRadius(typeVec == DTrpIndex,:) + DTrpWithinRadius(typeVec == DTrpIndex,:));

muDPro = muVec(typeVec == DProIndex);
fracDTForDPAtRadius = DTrpWithinRadius(typeVec == DProIndex,:) ./ ...
    (DProWithinRadius(typeVec == DProIndex,:) + DTrpWithinRadius(typeVec == DProIndex,:));

%%
if DEBUG
    figure(501);
    ax1 = subplot(1,2,1);
    imagesc(typeGrid)
    colormap(ax1,[1 0 0;0 1 0])
    axis image
    subplot(1,2,2);
    imagesc(muGrid)
    axis image
end

   