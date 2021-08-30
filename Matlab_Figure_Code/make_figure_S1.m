%load PA predictions
parametersModel = load('parameters_DP_DT_community_from_literature.mat','parameters');
parametersModel = parametersModel.parameters;

%% LOAD ANALYTICAL FUNCTIONS
analytical_functions
colors = make_figure_plot_colors();

%% FIGURE SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.filename               = 'FigS1';
settings.saveFig                = 1;
settings.saveFolder             = '~/Dropbox/Pair_Approximation_Manuscript/MainText/PairApproximationLatex/MatlabFigures/';

settings.AxisFontColor          = [1 1 1]*0.15;
settings.AxisFontSize           = 7;
settings.LabelFontSize          = 7;

settings.alphaConfInt           = 0.3;
settings.legendLocation         = 'upperright';
settings.legendLineWidth        = 0.1;

LineWidth                       = {0.5, 2};

maxRange                        = 21;
rangeTick                       = [3:6:21];

%% SET margins in cm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 column wide (20.5 picas / 3.42? / 8.7 cm)
%1.5 columns wide (27 picas / 4.5? / 11.4 cm)
%2 columns wide (42.125 picas / 7? / 17.8 cm)

%set size of figure
totHeight    = 10;
totWidth     = 17;
numCol       = 2;
numRow       = 2;

%set margins
hMargin     = 1.2;
vMargin     = 1;
lMargin     = 0.8;
rMargin     = 0.3;
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


%% PLOT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hFig                = figure('MenuBar','none', 'Name', settings.filename, 'NumberTitle', 'off');
hFig.Units          = 'centimeter';

%% Panel A CA replacing cells adjacent  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 1;
ax(idx)      =  axes(hFig);

%load data
dataA = load('simulation_CA_dynamics_ReplaceAdjNeighbor.mat');
dataA = dataA.data;

dataA.Y   = {dataA.DTinTime_mumax_r0(1,1:end,1,2), dataA.DTinTime_mumax_r0(1,1:end,2,2), dataA.DTinTime_mumax_r0(1,1:end,3,2), dataA.DTinTime_mumax_r0(1,1:end,4,2)};
dataA.X   = (0:100:100*(size(dataA.Y{1},2))-1)'/300;

fg = struct();
fg.LineStyle                    = {'-','-','-','-'};
fg.LineWidth                    = 1.5;
fg.Color                        = {colors.green, colors.MagentaDark, colors.Blue, colors.red};

dim = [.1 .65 .3 .3];
str = 'replication neighborhood {\itr_{R} = 8}';
annotation('textbox',dim,'String',str,'FitBoxToText','on');

for i = 1:size(dataA.Y,2)
    line(   dataA.X', dataA.Y{i}, ...
        'LineStyle', fg.LineStyle{i}, ...
        'LineWidth', fg.LineWidth, ...
        'Color', fg.Color{i} ); ...
end

make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0 400];
ax(idx).YLim            = [0 1];
ax(idx).XTick           = [0 200 400];
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'time [a.u.]';
ax(idx).YLabel.String   = 'frequency of type A';


%% Panel B CA replacing cell within interaction range  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 2;
ax(idx)      =  axes(hFig);

%load data
dataB = load('simulation_CA_dynamics_ReplaceInIntRange.mat');
dataB = dataB.data;

dataB.Y   = {dataB.DTinTime_mumax_r0(1,1:end,1,1), dataB.DTinTime_mumax_r0(1,1:end,2,1), dataB.DTinTime_mumax_r0(1,1:end,3,1), dataB.DTinTime_mumax_r0(1,1:end,4,1)};
dataB.X   = (0:100:100*(size(dataB.Y{1},2))-1)'/300;

fg = struct();
fg.LineStyle                    = {'-','-','-','-'};
fg.LineWidth                    = 1.5;
fg.Color                        = {colors.green, colors.MagentaDark, colors.Blue, colors.red};

dim = [.6 .65 .3 .3];
str = 'replication neighborhood {\it r_{R} = r_{A}}';
annotation('textbox',dim,'String',str,'FitBoxToText','on');

for i = 1:size(dataB.Y,2)
    line(   dataB.X, dataB.Y{i}, ...
        'LineStyle', fg.LineStyle{i}, ...
        'LineWidth', fg.LineWidth, ...
        'Color', fg.Color{i} );...
end

make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0 400];
ax(idx).YLim            = [0 1];
ax(idx).XTick           = [0 200 400];
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'time [a.u.]';
ax(idx).YLabel.String   = 'frequency of type A';


%% Panel C compare CA-CA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 3;
ax(idx)      =  axes(hFig);

dataC.Y   = dataA.PDT_end_mumax_r0;
dataC.X   = dataA.mumax;
xPos     = [1 1];

%plot data
fg = struct();
fg.LineStyle                    = {'-','-','-','-'};
fg.LineWidth                    = {1,1,1,1};
fg.Color                        = {colors.k, colors.k, colors.k, colors.k};
fg.Marker                       = {'none', 'none', 'none','none'};

fg.LineStyle2                    = {'none','none','none','none'};
fg.Marker2                       = {'o', 'o', 'o', 'o'};
fg.MarkerSize2                   = {25, 25, 25,25};
fg.MarkerFaceColor2              = {colors.green, colors.MagentaDark, colors.Blue, colors.red};
fg.MarkerFaceAlpha2              = {1, 1, 1, 1};
fg.MarkerEdgeColor2              = {'none','none','none','none'};

r0_idx=4;


for mm=r0_idx
    
    yData = dataC.Y(:,mm);
    xData = dataC.X;

    hold on
    plot(xData, yData, ...
        'LineStyle', fg.LineStyle{mm}, ...
        'LineWidth', fg.LineWidth{mm}, ...
        'Color', fg.Color{mm}, ...
        'Marker', fg.Marker{mm});
end

for mm=1:length(dataA.mumax)
    yData = dataC.Y(mm,r0_idx);
    xData = dataC.X(mm);
      
    hold on
    scatter(xData, yData, fg.MarkerSize2{mm}, ...
        'Marker', fg.Marker2{mm}, ...
        'MarkerFaceColor', fg.MarkerFaceColor2{mm}, ...
        'MarkerFaceAlpha', fg.MarkerFaceAlpha2{mm},...
        'MarkerEdgeColor', fg.MarkerEdgeColor2{mm});
    
end

make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).YLim            = [0 1];
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLim            = [0 7.3];
ax(idx).XTick           = [1 : 7];
ax(idx).XLabel.String   = {'\mu_B / \mu_A'}; ...'neighborhood size r_{A}=r_{B}'};
ax(idx).YLabel.String   = 'frequency of type A';

%% Panel D compare CA-CA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 4;
ax(idx)      =  axes(hFig);

dataD.Y   = dataB.PDT_end_mumax_r0;
dataD.X   = dataB.mumax;
xPos     = [1 1];

%plot data
fg = struct();
fg.LineStyle                    = {'-','-','-','-'};
fg.LineWidth                    = {1,1,1,1};
fg.Color                        = {colors.k, colors.k, colors.k, colors.k};
fg.Marker                       = {'none', 'none', 'none','none'};

fg.LineStyle2                    = {'none','none','none','none'};
fg.Marker2                       = {'v', 'v', 'v', 'v'};
fg.MarkerSize2                   = {30, 30, 30,30};
fg.MarkerFaceColor2              = {colors.green, colors.MagentaDark, colors.Blue, colors.red};
fg.MarkerFaceAlpha2              = {1, 1, 1, 1};
fg.MarkerEdgeColor2              = {'none','none','none','none'};

r0_idx=4;


for mm=r0_idx
    
    yData = dataD.Y(:,mm);
    xData = dataD.X;

    hold on
    plot(xData, yData, ...
        'LineStyle', fg.LineStyle{mm}, ...
        'LineWidth', fg.LineWidth{mm}, ...
        'Color', fg.Color{mm}, ...
        'Marker', fg.Marker{mm});
end

for mm=1:length(dataA.mumax)
    yData = dataD.Y(mm,r0_idx);
    xData = dataD.X(mm);
      
    hold on
    scatter(xData, yData, fg.MarkerSize2{mm}, ...
        'Marker', fg.Marker2{mm}, ...
        'MarkerFaceColor', fg.MarkerFaceColor2{mm}, ...
        'MarkerFaceAlpha', fg.MarkerFaceAlpha2{mm},...
        'MarkerEdgeColor', fg.MarkerEdgeColor2{mm});
    
end

make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).YLim            = [0 1];
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLim            = [0 7.3];
ax(idx).XTick           = [1 : 7];
ax(idx).XLabel.String   = {'\mu_B / \mu_A'}; ...'neighborhood size r_{A}=r_{B}'};
ax(idx).YLabel.String   = 'frequency of type A';


%% ADJUST FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fg.figureWidth  = totWidth;
fg.figureHeight = totHeight;
hFig.Position   = [1 1 totWidth totHeight];

hFig.PaperUnits = 'centimeters';
hFig.PaperSize  = [totWidth totHeight];

% hFig.Color = 'none';
% hFig.InvertHardcopy = 'off';

%% SAVE FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf, 'Renderer', 'painters');
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



