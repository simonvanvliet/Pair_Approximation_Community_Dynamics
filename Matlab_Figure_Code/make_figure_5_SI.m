%% LOAD ANALYTICAL FUNCTIONS
analytical_functions
colors = make_figure_plot_colors();

%% FIGURE SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.filename               = 'Fig5_SI';
settings.saveFig                = 1;
settings.saveFolder             = '~/Dropbox/Pair_Approximation_Manuscript/MainText/PairApproximationLatex/MatlabFigures/';

settings.MapFontColor           = [0 0 0];
settings.MapFontSize            = 7;
settings.AxisFontSize           = 7;
settings.LabelFontSize          = 7;
settings.LegendFontSize         = 6;

settings.alphaConfInt           = 0.3;
settings.legendLocation         = 'upperright';
settings.legendLineWidth        = 0.1;

maxRange                        = 20;
maxRelMu                        = 10;

nanColor = 0.45*[1 1 1];

%% SET margins in cm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 column wide (20.5 picas / 3.42? / 8.7 cm) 
%1.5 columns wide (27 picas / 4.5? / 11.4 cm) 
%2 columns wide (42.125 picas / 7? / 17.8 cm)

%set size of figure
totHeight    = 10;
totWidth     = 16.5;
numCol       = 3;
numRow       = 2;

%set margins
hMargin     = 2.5;
hMargin2    = 1.2;
vMargin     = 0.2;
lMargin     = 0.8;
tMargin     = 0.5;
bMargin     = 0.8;

cbMarg      = 0.2;
cbWid       = 0.3;
cbMarg2     = 0.9;
cbHeight    = 0.7;

rMargin = cbMarg+ cbWid+cbMarg2;

%calc panel size
marginsW  = (numCol-1)*hMargin + lMargin + rMargin;
marginsH  = (numRow-1)*vMargin + bMargin + tMargin;
plotH = (totHeight - marginsH) / numRow;
plotW = (totWidth - marginsW) / numCol;




%% calc rest of dimensions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotPos      = cell(numCol*numRow,1);
dd=1;
for rr = 1:numRow
    for cc = 1:numCol
        leftPos = (cc-1)*(plotW + hMargin) + lMargin;   
        botPos  = (numRow-rr)*(plotH + vMargin) + bMargin;
        plotPos{dd} = [leftPos botPos plotW  plotH];
        dd = dd+1;
    end
end

%% COLOR INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% PLOT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig                = figure('MenuBar','none', 'Name', settings.filename, 'NumberTitle', 'off');
hFig.Units          = 'centimeter';

RdBu = csvread('RdBu.csv');


%load data
CAData = load('simulation_CA_parameterScan.mat');
CAData = CAData.data;

%Setup data vectors
rangeA = CAData.r0;
muA = CAData.mumax;
rangeB = rangeA*CAData.neighborhoodRatio;
[nbBMat, muAMat] =  meshgrid(rangeB, muA);
XTickDiscrete = 1:length(rangeB);

%PA functions
freqDep = @(a,b,k,l) max(0, min(1, (a./k - b./l + a.*(k - 2)./k) ./ (b.*(l-2)./l + a.*(k-2)./k)));
clustering = @(r) (r - 2) ./ (r - 1);



%% ROW 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fPA = freqDep(1, muAMat, rangeA, nbBMat);
fCA = CAData.P1_end;

error = fPA - fCA;
error(~(fPA>0 & fPA<1 & fCA>0 & fCA<1)) = nan;

maxErr = ceil(max(abs(error(:))) / 0.05)*0.05;

dataSet = {fPA, fCA, error};
dataRange = {[0 1], [0 1], [-maxErr maxErr]};
alphaData = {(fPA>0 & fPA<1), (fCA>0 & fCA<1), (fPA>0 & fPA<1 & fCA>0 & fCA<1)};    

%plot frequencies
clear ax        
for idx = 1:2
    
    plotIdx = (idx-1)*3 + 1;
    
    ax(idx)      =  axes(hFig);

    imH = imagesc(XTickDiscrete, log10(muA), dataSet{idx}, dataRange{idx});
    colormap(ax(idx), flip(RdBu))
        
    set(imH, 'AlphaData', alphaData{idx})

    axis square
    axis xy

    make_figure_makeAxisNice(ax(idx), plotPos{plotIdx}, settings)

    xTick = CAData.neighborhoodRatio;
    yTick = {'1/8','1/4','1/2',1,2,4,8};
    yTickLoc = [1/8 1/4 1/2 1 2 4 8];
    
    ax(idx).XTick           = XTickDiscrete;    
    ax(idx).XTickLabel      = xTick;    
    ax(idx).YTick           = log10(yTickLoc);
    ax(idx).YTickLabel      = yTick;
    ax(idx).YLabel.String   = 'growth ratio         ';
    ax(idx).XLabel.String   = 'neighborhood size ratio        ';
    
    set(gca, 'Color', nanColor)
    
    %color bar data
    cbPos = ax(idx).Position;
    cbH  = cbPos(4)*cbHeight;
    dH = (cbPos(4)-cbH)/2;
    cbPos(1) = cbPos(1)+cbPos(3)+cbMarg;
    cbPos(2) = cbPos(2)+dH;
    cbPos(3) = cbWid;
    cbPos(4) = cbH;

    cb                      = colorbar(ax(idx));
    cb.Units                = 'centimeter';
    cb.Position             = cbPos;
    cb.Label.String         = 'P(A)';
    cb.Limits               = [0 1];
    cb.Ticks                = [0 0.5 1];

end



%color bar error 
% cbPos = ax(3).Position;
% cbH  = cbPos(4)*cbHeight;
% dH = (cbPos(4)-cbH)/2;
% cbPos(1) = cbPos(1)+cbPos(3)+cbMarg;
% cbPos(2) = cbPos(2)+dH;
% cbPos(3) = cbWid;
% cbPos(4) = cbH;
% 
% cb2                      = colorbar(ax(3));
% cb2.Units                = 'centimeter';
% cb2.Position             = cbPos;
% cb2.Label.String         = 'Difference PA-CA';
% cb2.Limits               = [-maxErr maxErr];
% cb2.Ticks                = [-maxErr 0 maxErr];



%% ROW 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clusteringLim = [0.8 1.2];
clusteringTick = [0.8 0.9 1 1.1 1.2];
clusteringTickLable = {'<0.8','0.9','1','1.1','>1.2'};

fPA = freqDep(1, muAMat, rangeA, nbBMat);
fCA = CAData.P1_end;

rangeAMat = rangeA * ones(size(nbBMat));


clPA = clustering(rangeAMat);
clCA = CAData.clustering1';

error = clPA - clCA;
error(~(fPA>0 & fPA<1 & fCA>0 & fCA<1) | isinf(error)) = nan;

maxErr = ceil(max(abs(error(:))) / 0.05)*0.05;

dataSet = {clPA, clCA, error};
dataRange = {clusteringLim, clusteringLim, [-maxErr maxErr]};
alphaData = {(fPA>0 & fPA<1), (fCA>0 & fCA<1), (fPA>0 & fPA<1 & fCA>0 & fCA<1)};   

%plot frequencies
clear ax        
for idx = 1:2
    plotIdx = (idx-1)*3 + 2;

    ax(idx)      =  axes(hFig);

    imH = imagesc(XTickDiscrete, log10(muA), dataSet{idx}, dataRange{idx});
    colormap(ax(idx), flip(RdBu))
        
    set(imH, 'AlphaData', alphaData{idx})

    axis square
    axis xy

    make_figure_makeAxisNice(ax(idx), plotPos{plotIdx}, settings)

    xTick = CAData.neighborhoodRatio;
    yTick = {'1/8','1/4','1/2',1,2,4,8};
    yTickLoc = [1/8 1/4 1/2 1 2 4 8];
    
    ax(idx).XTick           = XTickDiscrete;    
    ax(idx).XTickLabel      = xTick;    
    ax(idx).YTick           = log10(yTickLoc);
    ax(idx).YTickLabel      = yTick;
    ax(idx).YLabel.String   = 'growth ratio         ';
    ax(idx).XLabel.String   = 'neighborhood size ratio        ';
    
    set(gca, 'Color', nanColor)
    
    %color bar data
    cbPos = ax(idx).Position;
    cbH  = cbPos(4)*cbHeight;
    dH = (cbPos(4)-cbH)/2;
    cbPos(1) = cbPos(1)+cbPos(3)+cbMarg;
    cbPos(2) = cbPos(2)+dH;
    cbPos(3) = cbWid;
    cbPos(4) = cbH;

    cb                      = colorbar(ax(idx));
    cb.Units                = 'centimeter';
    cb.Position             = cbPos;
    cb.Label.String         = 'local/global frequency';
    cb.Limits               = clusteringLim;
    cb.Ticks                = clusteringTick;
    cb.TickLabels           = clusteringTickLable;
        
end



%color bar error 
% cbPos = ax(3).Position;
% cbH  = cbPos(4)*cbHeight;
% dH = (cbPos(4)-cbH)/2;
% cbPos(1) = cbPos(1)+cbPos(3)+cbMarg;
% cbPos(2) = cbPos(2)+dH;
% cbPos(3) = cbWid;
% cbPos(4) = cbH;
% 
% cb2                      = colorbar(ax(3));
% cb2.Units                = 'centimeter';
% cb2.Position             = cbPos;
% cb2.Label.String         = 'Difference PA-CA';
% cb2.Limits               = [-maxErr maxErr];
% cb2.Ticks                = [-maxErr 0 maxErr];

%% ROW 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clusteringLim = [0.8 1.2];
clusteringTick = [0.8 0.9 1 1.1 1.2];
clusteringTickLable = {'<0.8','0.9','1','1.1','>1.2'};

fPA = freqDep(1, muAMat, rangeA, nbBMat);
fCA = CAData.P1_end;

rangeAMat = rangeA * ones(size(nbBMat));


clPA = clustering(nbBMat);
clCA = CAData.clustering0';


error = clPA - clCA;
error(~(fPA>0 & fPA<1 & fCA>0 & fCA<1) | isinf(error)) = nan;

maxErr = ceil(max(abs(error(:))) / 0.05)*0.05;

dataSet = {clPA, clCA, error};
dataRange = {clusteringLim, clusteringLim, [-maxErr maxErr]};
alphaData = {(fPA>0 & fPA<1), (fCA>0 & fCA<1), (fPA>0 & fPA<1 & fCA>0 & fCA<1)};    

%plot frequencies
clear ax        
for idx = 1:2
    plotIdx = (idx-1)*3 + 3;

    ax(idx)      =  axes(hFig);

    imH = imagesc(XTickDiscrete, log10(muA), dataSet{idx}, dataRange{idx});
    colormap(ax(idx), flip(RdBu))
        
    set(imH, 'AlphaData', alphaData{idx})

    axis square
    axis xy

    make_figure_makeAxisNice(ax(idx), plotPos{plotIdx}, settings)

    xTick = CAData.neighborhoodRatio;
    yTick = {'1/8','1/4','1/2',1,2,4,8};
    yTickLoc = [1/8 1/4 1/2 1 2 4 8];
    
    ax(idx).XTick           = XTickDiscrete;    
    ax(idx).XTickLabel      = xTick;    
    ax(idx).YTick           = log10(yTickLoc);
    ax(idx).YTickLabel      = yTick;
    ax(idx).YLabel.String   = 'growth ratio         ';
    ax(idx).XLabel.String   = 'neighborhood size ratio        ';
    
    set(gca, 'Color', nanColor)
    
    %color bar data
    cbPos = ax(idx).Position;
    cbH  = cbPos(4)*cbHeight;
    dH = (cbPos(4)-cbH)/2;
    cbPos(1) = cbPos(1)+cbPos(3)+cbMarg;
    cbPos(2) = cbPos(2)+dH;
    cbPos(3) = cbWid;
    cbPos(4) = cbH;

    cb                      = colorbar(ax(idx));
    cb.Units                = 'centimeter';
    cb.Position             = cbPos;
    cb.Label.String         = 'local/global frequency';
    cb.Limits               = clusteringLim;
    cb.Ticks                = clusteringTick;
        cb.TickLabels           = clusteringTickLable;

        
end



%color bar error 
% cbPos = ax(3).Position;
% cbH  = cbPos(4)*cbHeight;
% dH = (cbPos(4)-cbH)/2;
% cbPos(1) = cbPos(1)+cbPos(3)+cbMarg;
% cbPos(2) = cbPos(2)+dH;
% cbPos(3) = cbWid;
% cbPos(4) = cbH;
% 
% cb2                      = colorbar(ax(3));
% cb2.Units                = 'centimeter';
% cb2.Position             = cbPos;
% cb2.Label.String         = 'Difference PA-CA';
% cb2.Limits               = [-maxErr maxErr];
% cb2.Ticks                = [-maxErr 0 maxErr];


%% ADJUST FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fg.figureWidth  = totWidth;
fg.figureHeight = totHeight;
hFig.Position   = [1 1 totWidth totHeight];

hFig.PaperUnits = 'centimeters';
hFig.PaperSize  = [totWidth totHeight];

% hFig.Color = 'none';
% hFig.InvertHardcopy = 'off';

set(gcf,'InvertHardCopy','Off');
set(gcf,'Color',[1 1 1]);


%% SAVE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.saveFig
    saveas( hFig, [settings.saveFolder settings.filename '.pdf']);
    disp([' * Saved figure in ' settings.filename '.pdf']);
end

%hFig.Color = [0.9 0.9 0.9];



function makeTextNice(tObj, settings)
    tObj.FontName             = 'Arial';
    tObj.FontWeight           = 'normal';
    tObj.FontSize             = settings.MapFontSize;
    tObj.HorizontalAlignment  = settings.HorizontalAlignment;
    tObj.VerticalAlignment   = 'top';
    tObj.Color = settings.MapFontColor;
end



