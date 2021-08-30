%% LOAD ANALYTICAL FUNCTIONS
analytical_functions
colors = make_figure_plot_colors();

%% FIGURE SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.filename               = 'Fig1';
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

%% SET margins in cm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 column wide (20.5 picas / 3.42? / 8.7 cm) 
%1.5 columns wide (27 picas / 4.5? / 11.4 cm) 
%2 columns wide (42.125 picas / 7? / 17.8 cm)

%set size of figure
totHeight    = 3;
totWidth     = 3.5;
numCol       = 1;
numRow       = 2;

%set margins
hMargin     = 1.2;
vMargin     = 0.7;
lMargin     = 0.8;
rMargin     = 0.1;
tMargin     = 0.2;
bMargin     = 0.7;

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


%% PLOT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig                = figure('MenuBar','none', 'Name', settings.filename, 'NumberTitle', 'off');
hFig.Units          = 'centimeter';

%% Panel AB  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
data = load('data_growth_vs_arrangement.mat');
XDataNames = 'fractionVec';
YDataNames ='growthVec';
     
data.Y   = data.(YDataNames);
data.X   = data.(XDataNames);

fg = struct();
fg.LineStyle                    = {'-','-'};
fg.LineWidth                    = {2, 2}; 
fg.Color                        = {colors.DP, colors.DT*0.8}; 
fg.Marker                       = {'o', 'o'}; 
fg.MarkerSize                   = {4, 4};
fg.MarkerFaceColor              = {'k', 'k'}; 
fg.MarkerFaceAlpha              = {0.3, 0.3}; 
fg.MarkerEdgeColor              = {'none','none'}; 
fg.XLabel                       = {'local frequency \DeltaT', 'local frequency \DeltaP'};

for tt=1:2
    ax(tt)      =  axes(hFig);
    hold on
    
    %get data points
    notNan   = ~isnan(data.X{tt}) & ~isnan(data.Y{tt});
    xData = data.X{tt}(notNan);
    yData = data.Y{tt}(notNan);
    
    %linear fit of data
    polynom = polyfit(xData, yData, 1);
    fracVec = [0 1];
    fittedGrowth = fracVec*polynom(1) + polynom(2);
 
    scatter(xData, yData, ...
        'Marker', fg.Marker{tt}, ...
        'SizeData', fg.MarkerSize{tt}, ...
        'MarkerFaceColor', fg.MarkerFaceColor{tt}, ...
        'MarkerFaceAlpha', fg.MarkerFaceAlpha{tt}, ...
        'MarkerEdgeColor', fg.MarkerEdgeColor{tt});
       
    line(fracVec, fittedGrowth, ...
        'LineStyle', fg.LineStyle{tt}, ...
        'LineWidth', fg.LineWidth{tt}, ...
        'Color', fg.Color{tt}, ...
        'Marker', 'none');

    make_figure_makeAxisNice(ax(tt), plotPos{tt}, settings)
    %axis square

    ax(tt).XLim            = [0 1];
    ax(tt).YLim            = [0 0.8];
    ax(tt).XTick           = [0 0.5 1];    
    ax(tt).YTick           = [0 0.4 0.8];
    ax(tt).XLabel.String   = fg.XLabel{tt};
    ax(tt).YLabel.String   = 'growth 1/h';
end



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


