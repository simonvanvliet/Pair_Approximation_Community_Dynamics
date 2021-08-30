function functionOutput = SteadyState_3D_Calc_InteractionRange(settings, gridCellTypeStack, ...
    edgeToIgnore, radiusVec, DEBUG, varargin)
%Compares growth range to interaction range for single parameter set
%INPUT: settings: model settings
%INPUT: gridCellTypeStack: 4D matrix of grids to analyze
%INPUT: edgeToIgnore: number of pixels on edge of chamber where cells are excluded from analyis
%INPUT: radiusVec: vector of radia (in pixels) for which spatial arrangement is quantified
%
%functionOutput: structure with fields:
% -rhoDPAtRadius: spearmann correlation^2 of mu of DP vs DT/(DP+DT) for different radia
% -rhoDTAtRadius: spearmann correlation^2 of mu of DT vs DP/(DP+DT) for different radia
% -maxRhoDP: max correlation^2 DP vs DT/(DP+DT)
% -maxRhoDT: max correlation^2 DT vs DP/(DP+DT)
% -interactionRangeDP: interaction range DP (radius of max correlation^2)
% -interactionRangeDT: interaction range DT (radius of max correlation^2)
% -fractionOtherDPro: DTrp/(DTrp+DPro) within neighborhood of size that maximizes correlation
% -muDPro: corresponding single cell growth rates
% -fractionOtherDTrp: DPro/(DTrp+DPro) within neighborhood of size that maximizes correlation
% -muDTrp: corresponding single cell growth rates
%
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020

if nargin  == 4; DEBUG = 0; end


%SET DEBUG options
%sets what is plotted if DEBUG = 1
DEBUG_GROWTHCURVE = 0; %plot growth curves
DEBUG_CORRELATIONCURVE = 1; %plot correlation curves

numChambers = size(gridCellTypeStack,4);
%init output
functionOutput = struct;

%plot options
greenC = [0,136,55]/256;
redC = [202,0,32]/256;

%init output
muDTrp = [];
fracDPforDTatRadius = [];
muDPro = [];
fracDTforDPatRadius = [];

%% First solve for steady state and calc mu and SA for all cells
%loop over chambers

%First loop creates job list for parallel pool
for nn = 1:numChambers
    %get grid
    currChamber = squeeze(gridCellTypeStack(:, :, :, nn));
    %send to parallel pool
    parOutput(nn)  =  parfeval(@parCompute, 1, currChamber, settings);
end

%SET timepout for task (in s)
maxTimeOut  =  settings.timeOutTime;

%Second loops gets results from parallel pool when they are done
for nn = 1:numChambers
    %get result, the loop waits untill next result comes in, then in goes though one step
    %this repeats until loop is done. When t>maxTimeOut, output is empty
    
%     %get grid
%     currChamber = squeeze(gridCellTypeStack(:, :, :, nn));
%     
%     %reset grid scaling
%     settings.gridScaling = settings.gridScalingBase;
%     %setup grid
%     [initialCondition] = SteadyState_3D_InitGrid(settings, 'ProvidedGrid', currChamber);
%     %run simulation
%     [output, ~, ~, ~] = SteadyState_3D_RefineGrid(settings, initialCondition);
%     
%     
    [~, output]  =  fetchNext(parOutput, maxTimeOut);
    
    if ~isempty(output)
        muGrid = output.mu;
        typeGrid = output.gridCellType;
        if sum(imag(muGrid(:))) == 0
            %compute growth curve (mu vs other/total) for current chamber
            [muDTrpCurr,fracDPForDTAtRadius,muDProCurr,fracDTForDPAtRadius] = ...
                SteadyState_3D_Calc_LocalProp(muGrid, typeGrid, settings.gridSize, edgeToIgnore, radiusVec);

            %combine data with other chambers
            muDTrp = cat(1, muDTrp, muDTrpCurr);
            fracDPforDTatRadius = cat(1, fracDPforDTatRadius, fracDPForDTAtRadius);
            muDPro = cat(1, muDPro, muDProCurr);
            fracDTforDPatRadius = cat(1, fracDTforDPatRadius, fracDTForDPAtRadius);

            clear output typeGrid muGrid initialCondition currChamber
            %run garabage collection
            java.lang.System.gc()
        else
            warning('Skipped chamber')
        end
    else
        warning('Skipped chamber')
    end
    
end

if isempty(muDTrp)
    functionOutput.rhoDTAtRadius = nan;
    functionOutput.rhoDPAtRadius = nan;
    functionOutput.maxRhoDT = nan;
    functionOutput.interactionRangeDT = nan;
    functionOutput.maxRhoDP = nan;
    functionOutput.interactionRangeDP = nan;
    functionOutput.fractionOtherDPro = nan;
    functionOutput.muDPro = nan;
    functionOutput.fractionOtherDTrp = nan;
    functionOutput.muDTrp = nan;
else
    
    %% calculate correlation between SA and mu for different radia
    %init output
    numRadia = length(radiusVec);
    rhoDTAtRadius = nan(numRadia,1);
    rhoDPAtRadius = nan(numRadia,1);
    
    if ~isempty(varargin)
        correlationType  =  varargin{1};
        fprintf('using %s correlation', correlationType);
    else
        correlationType  =  'Spearman';
    end
    
    for rr = 1:numRadia %loop radia
        rhoDT = corr(fracDPforDTatRadius(:,rr),muDTrp,'type',correlationType);
        rhoDTAtRadius(rr) = rhoDT^2;
        
        rhoDP = corr(fracDTforDPatRadius(:,rr),muDPro,'type',correlationType);
        rhoDPAtRadius(rr) = rhoDP^2;
    end
        
    functionOutput.rhoDTAtRadius = rhoDTAtRadius;
    functionOutput.rhoDPAtRadius = rhoDPAtRadius;
    
    %% find maxima of correlation
    %first interpolate data to get better estimate
    radiusVecFine = linspace(min(radiusVec), max(radiusVec), 1000);
    RhoSqTrpIntrep = interp1(radiusVec, rhoDTAtRadius, radiusVecFine);
    RhoSqProIntrep = interp1(radiusVec, rhoDPAtRadius, radiusVecFine);
    
    %find interaction range as point of maximum correlation
    %substract 1 as we need to measure distance between perimeters of cells
    %(default distances are calculated center-center which adds 2* cell width/2 = 1 to all distances
    [maxRhoDT,idx] = max(RhoSqTrpIntrep);
    interactionRangeDT = radiusVecFine(idx) - 1;
    [~, idxDTCourse] = min(abs(radiusVec - radiusVecFine(idx)));
    
    [maxRhoDP,idx] = max(RhoSqProIntrep);
    interactionRangeDP = radiusVecFine(idx) - 1;
    [~, idxDPCourse] = min(abs(radiusVec - radiusVecFine(idx)));
    
    functionOutput.maxRhoDT = maxRhoDT;
    functionOutput.interactionRangeDT = interactionRangeDT;
    
    functionOutput.maxRhoDP = maxRhoDP;
    functionOutput.interactionRangeDP = interactionRangeDP;
    
    functionOutput.fractionOtherDPro = fracDTforDPatRadius(:,idxDPCourse);
    functionOutput.muDPro = muDPro;
    functionOutput.fractionOtherDTrp = fracDPforDTatRadius(:,idxDTCourse);
    functionOutput.muDTrp = muDTrp;
    
    
    %%
    if DEBUG_GROWTHCURVE && DEBUG
        figure(1012);
        clear clf
        plotGrowthCurve(radiusVec,fracDPforDTatRadius,fracDTforDPatRadius,muDTrp,muDPro);
    end
    
    if DEBUG_CORRELATIONCURVE && DEBUG
        figure(1013)
        clear clf
        
        h = plot(radiusVec,rhoDTAtRadius,radiusVec,rhoDPAtRadius);
        set(h,'LineWidth',2)
        set(h(1),'Color',redC)
        set(h(2),'Color',greenC)
        
        xlabel('radius (cells)')
        ylabel('rho')
        ylim([0 1])
        hold on
        k1 = line([muTryHM muTryHM],[0 1]);
        k2 = line([muProHM muProHM],[0 1]);
        
        set(k1,'Color',redC,'LineStyle','-','LineWidth',1)
        set(k2,'Color',greenC,'LineStyle','-','LineWidth',1)
        hold off
        
    end
        
end
end

% 
%parallel execution function
function output  =  parCompute(currChamber, settings)
    %reset grid scaling
    settings.gridScaling = settings.gridScalingBase;
    %setup grid
    [initialCondition] = SteadyState_3D_InitGrid(settings, 'ProvidedGrid', currChamber);
    %run simulation
    [output, ~, ~, ~] = SteadyState_3D_RefineGrid(settings, initialCondition);
end


function plotGrowthCurve(radiusVec,fracDPforDTatRadius,fracDTforDPatRadius,muDTrp,muDPro)
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
    
    DTBin = nan(numBin,1);
    DPBin = nan(numBin,1);
    
    for bb = 1:numBin
        lTr = fracBin(bb);
        hTr = fracBin(bb+1);
        
        inBinR = fracDPforDTatRadius(:,rr)>= lTr & fracDPforDTatRadius(:,rr)<hTr;
        inBinG = fracDTforDPatRadius(:,rr)>= lTr & fracDTforDPatRadius(:,rr)<hTr;
        
        DTBin(bb) = mean(muDTrp(inBinR));
        DPBin(bb) = mean(muDPro(inBinG));
    end
    
    %plot DT data
    subplot(numRow,numCol,rr)
    plot(fracDPforDTatRadius(:,rr),muDTrp,'.b','MarkerSize',2)
    hold on
    k = plot(fracBinCenter,DTBin,'o');
    set(k,'MarkerSize',5,'MarkerFaceColor',redC);
    hold off
    %ylabel('mu Dtry')
    %xlabel('DP/(DT+DP)')
    title(sprintf('DTrp r = %.1g',radiusVec(rr)));
    %ylim([0 1])
    
    %plot DP data
    subplot(numRow,numCol,round(numRow/2*numCol)+rr)
    plot(fracDTforDPatRadius(:,rr),muDPro,'.b','MarkerSize',2)
    hold on
    k = plot(fracBinCenter,DPBin,'o');
    set(k,'MarkerSize',5,'MarkerFaceColor',greenC);
    hold off
    %ylabel('mu Dpro')
    %xlabel('DT/(DT+DP)')
    title(sprintf('DPro r = %.1g',radiusVec(rr)));
    %ylim([0 1])
    
end
end



