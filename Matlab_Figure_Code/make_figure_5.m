%% LOAD ANALYTICAL FUNCTIONS
analytical_functions
colors = make_figure_plot_colors();

%% FIGURE SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.filename               = 'Fig5';
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

%% SET margins in cm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 column wide (20.5 picas / 3.42? / 8.7 cm) 
%1.5 columns wide (27 picas / 4.5? / 11.4 cm) 
%2 columns wide (42.125 picas / 7? / 17.8 cm)

%set size of figure
totHeight    = 5;
totWidth     = 8.7;
numCol       = 2;
numRow       = 1;

%set margins
hMargin     = 1;
vMargin     = 0.8;
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


%% Panel AB  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RdBu = csvread('RdBu.csv');

%setup range vectors

rangeA = 10;
muA = logspace(log10(1/8),log10(8),21);
rangeB = 10:70;
[nbBMat, muAMat] =  meshgrid(rangeB, muA);

%calcualte steady state
freqDep = @(a,b,k,l) max(0, min(1, (a./k - b./l + a.*(k - 2)./k) ./ (b.*(l-2)./l + a.*(k-2)./k)));
densDep = @(a,b,k,l) max(0, min(1, (a - b + a.*(k - 2)) ./ (b.*(l-2) + a.*(k-2))));

dataSet =  {freqDep(muAMat, 1, rangeA, nbBMat), ...
            densDep(muAMat, 1, rangeA, nbBMat)};
clear ax        
        
for idx = 1:2
    ax(idx)      =  axes(hFig);

    imagesc(rangeB/rangeA, log10(muA), dataSet{idx}, [0 1])
    colormap(ax(idx), flip(RdBu))

    axis square
    axis xy

    make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

    xTick = [10, 30, 50, 70]/rangeA;
    yTick = {'1/8','1/4','1/2',1,2,4,8};
    yTickLoc = [1/8 1/4 1/2 1 2 4 8];
    
    ax(idx).XTick           = xTick;    
    ax(idx).XTickLabel      = xTick;    
    ax(idx).YTick           = log10(yTickLoc);
    ax(idx).YTickLabel      = yTick;
    ax(idx).YLabel.String   = 'growth ratio         ';
    ax(idx).XLabel.String   = 'neighborhood size ratio        ';

end

cbPos = ax(2).Position;
cbH  = cbPos(4)*cbHeight;
dH = (cbPos(4)-cbH)/2;
cbPos(1) = cbPos(1)+cbPos(3)+cbMarg;
cbPos(2) = cbPos(2)+dH;
cbPos(3) = cbWid;
cbPos(4) = cbH;

cb                      = colorbar();
cb.Units                = 'centimeter';
cb.Position             = cbPos;
cb.Label.String         = 'P(A)';
cb.Limits               = [0 1];
cb.Ticks                = [0 0.5 1];

%% ADJUST FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fg.figureWidth  = totWidth;
fg.figureHeight = totHeight;
hFig.Position   = [1 1 totWidth totHeight];

hFig.PaperUnits = 'centimeters';
hFig.PaperSize  = [totWidth totHeight];

% hFig.Color = 'none';
% hFig.InvertHardcopy = 'off';

%% SAVE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.saveFig
    saveas( hFig, [settings.saveFolder settings.filename '.pdf']);
    disp([' * Saved figure in ' settings.filename '.pdf']);
end

hFig.Color = [0.9 0.9 0.9];




function makeTextNice(tObj, settings)
    tObj.FontName             = 'Arial';
    tObj.FontWeight           = 'normal';
    tObj.FontSize             = settings.MapFontSize;
    tObj.HorizontalAlignment  = settings.HorizontalAlignment;
    tObj.VerticalAlignment   = 'top';
    tObj.Color = settings.MapFontColor;
end



