function functionOutput = SteadyState_3D_GetInteractionRange(mu1,mu2,fo1,fo2,radiusVec)
%Compares growth range to interaction range for single parameter set
%INPUT: settings: model settings
%INPUT: gridCellTypeStack: 4D matrix of grids to analyze
%INPUT: edgeToIgnore: number of pixels on edge of chamber where cells are excluded from analyis
%INPUT: radiusVec: vector of radia (in pixels) for which spatial arrangement is quantified
%
%functionOutput: structure with fields:
% -rho2AtRadius: spearmann correlation^2 of mu of 2 vs 1/(2+1) for different radia
% -rho1AtRadius: spearmann correlation^2 of mu of 1 vs 2/(2+1) for different radia
% -maxRho2: max correlation^2 2 vs 1/(2+1)
% -maxRho1: max correlation^2 1 vs 2/(2+1)
% -interactionRange2: interaction range 2 (radius of max correlation^2)
% -interactionRange1: interaction range 1 (radius of max correlation^2)
% -fractionOther2: 1/(1+2) within neighborhood of size that maximizes correlation
% -mu2: corresponding single cell growth rates
% -fractionOther1: 2/(1+2) within neighborhood of size that maximizes correlation
% -mu1: corresponding single cell growth rates
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 2021.02.02

%init output
functionOutput = struct;

DEBUG_GROWTHCURVE = 0;
DEBUG_CORRELATIONCURVE = 0;

correlationType  =  'Spearman';

%% calculate correlation between SA and mu for different radia
%init output
numRadia = length(radiusVec);
rho1AtRadius = nan(numRadia,1);
rho2AtRadius = nan(numRadia,1);

for rr = 1:numRadia %loop radia
    rho1 = corr(fo1(:,rr),mu1,'type',correlationType);
    rho1AtRadius(rr) = rho1^2;
    
    rho2 = corr(fo2(:,rr),mu2,'type',correlationType);
    rho2AtRadius(rr) = rho2^2;
end

functionOutput.rho1AtRadius = rho1AtRadius;
functionOutput.rho2AtRadius = rho2AtRadius;

%% find maxima of correlation
%first interpolate data to get better estimate
radiusVecFine = linspace(min(radiusVec), max(radiusVec), 1000);
RhoSq1Intrep = interp1(radiusVec, rho1AtRadius, radiusVecFine);
RhoSq2Intrep = interp1(radiusVec, rho2AtRadius, radiusVecFine);

%find interaction range as point of maximum correlation
%substract 1 as we need to measure distance between perimeters of cells
%(default distances are calculated center-center which adds 2* cell wi1h/2 = 1 to all distances
[maxRho1,idx] = max(RhoSq1Intrep);
interactionRange1 = radiusVecFine(idx) - 1;
[~, idx1Course] = min(abs(radiusVec - radiusVecFine(idx)));

[maxRho2,idx] = max(RhoSq2Intrep);
interactionRange2 = radiusVecFine(idx) - 1;
[~, idx2Course] = min(abs(radiusVec - radiusVecFine(idx)));

functionOutput.maxRho1 = maxRho1;
functionOutput.interactionRange1 = interactionRange1;

functionOutput.maxRho2 = maxRho2;
functionOutput.interactionRange2 = interactionRange2;

functionOutput.fractionOther2 = fo2(:,idx2Course);
functionOutput.mu2 = mu2;
functionOutput.fractionOther1 = fo1(:,idx1Course);
functionOutput.mu1 = mu1;


%%
if DEBUG_GROWTHCURVE
    figure(1012);
    clear clf
    plotGrowthCurve(radiusVec,fo1,fo2,mu1,mu2);
end

if DEBUG_CORRELATIONCURVE
    figure(1013)
    clear clf
    
    h = plot(radiusVec,rho1AtRadius,radiusVec,rho2AtRadius);
    set(h,'LineWi1h',2)
    set(h(1),'Color',redC)
    set(h(2),'Color',greenC)
    
    xlabel('radius (cells)')
    ylabel('rho')
    ylim([0 1])
    hold on
    k1 = line([muTryHM muTryHM],[0 1]);
    k2 = line([muProHM muProHM],[0 1]);
    
    set(k1,'Color',redC,'LineStyle','-','LineWi1h',1)
    set(k2,'Color',greenC,'LineStyle','-','LineWi1h',1)
    hold off
end



function plotGrowthCurve(radiusVec,fo1,fo2,mu1,mu2)
%plot options
greenC = [0,136,55]/256;
redC = [202,0,32]/256;

numRadia = length(radiusVec);

numRow = ceil(sqrt(numRadia)*0.4)*2;
numCol = ceil(numRadia/(numRow/2));

for rr = 1:numRadia
    %bin growth rate data
    fracBin = 0:0.05:1;
    fracBinCenter = (fracBin(2:end)+fracBin(1:end-1))/2;
    numBin = length(fracBinCenter);
    
    n1Bin = nan(numBin,1);
    n2Bin = nan(numBin,1);
    
    for bb = 1:numBin
        lTr = fracBin(bb);
        hTr = fracBin(bb+1);
        
        inBinR = fo1(:,rr)>= lTr & fo1(:,rr)<hTr;
        inBinG = fo2(:,rr)>= lTr & fo2(:,rr)<hTr;
        
        n1Bin(bb) = mean(mu1(inBinR));
        n2Bin(bb) = mean(mu2(inBinG));
    end
    
    %plot 1 data
    subplot(numRow,numCol,rr)
    plot(fo1(:,rr),mu1,'.b','MarkerSize',2)
    hold on
    k = plot(fracBinCenter,n1Bin,'o');
    set(k,'MarkerSize',5,'MarkerFaceColor',redC);
    hold off
    %ylabel('mu 1ry')
    %xlabel('2/(1+2)')
    title(sprintf('1 r = %.1g',radiusVec(rr)));
    %ylim([0 1])
    
    %plot 2 data
    subplot(numRow,numCol,round(numRow/2*numCol)+rr)
    plot(fo2(:,rr),mu2,'.b','MarkerSize',2)
    hold on
    k = plot(fracBinCenter,n2Bin,'o');
    set(k,'MarkerSize',5,'MarkerFaceColor',greenC);
    hold off
    %ylabel('mu 2')
    %xlabel('1/(1+2)')
    title(sprintf('2 r = %.1g',radiusVec(rr)));
    %ylim([0 1])
    
end




