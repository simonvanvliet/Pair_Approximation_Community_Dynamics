analytical_functions
colors = make_figure_plot_colors();

%% FIGURE SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.filename               = 'Fig4_part2';
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

maxRange                        = 30;
maxRelMu                        = 10;

%% SET margins in cm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 column wide (20.5 picas / 3.42? / 8.7 cm) 
%1.5 columns wide (27 picas / 4.5? / 11.4 cm) 
%2 columns wide (42.125 picas / 7? / 17.8 cm)

%set size of figure
totHeight    = 4;
totWidth     = 7.5;
numCol       = 1;
numRow       = 1;

%set margins
hMargin     = 1;
vMargin     = 0.8;
lMargin     = 1;
rMargin     = 0.5;
tMargin     = 0.2;
bMargin     = 0.8;

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



%% Panel B  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear ax

rangeMuliplier = 1;
idx = 1;
ax(idx)      =  axes(hFig);

relMu = logspace(-1,1,100);

rangeB = rangeMuliplier * (12:-1:3);
rangeA = rangeMuliplier * 3 * ones(1,length(rangeB));

relRange = rangeB./rangeA;

[relMuMat, nbAMat] =  meshgrid(relMu, rangeA);
[~, nbBMat] =  meshgrid(relMu, rangeB);

fitness = Wrel(relMuMat, 1, nbAMat, nbBMat);
surf(log10(relMu), relRange, fitness,'LineStyle','none')
colormap(ax(idx), [0,0,0; parula(256)])
view(10, 35)
view(10, 35)
make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

xTick = [0.1 0.2 0.5 1 2 5 10];

ax(idx).XTick           = log10(xTick);    
ax(idx).XTickLabel      = xTick;    
ax(idx).YTick           = [1 2 3 4];
ax(idx).ZTick           = [0 0.5 1];
ax(idx).ZLim            = [0  1];
ax(idx).XLabel.String   = '';
ax(idx).YLabel.String   = '';
ax(idx).XTickLabel      = []; 
ax(idx).YTickLabel      = []; 
ax(idx).ZTickLabel      = []; 

%     cb                      = colorbar('NorthOutside');
%     cb.Label.String         = 'relative fitness';
%     cb.Limits               = [0 1];
%     cb.Ticks                = [0 0.5 1];
%

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



