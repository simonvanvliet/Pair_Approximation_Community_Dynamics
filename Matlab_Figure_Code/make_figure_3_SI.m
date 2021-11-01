%% LOAD ANALYTICAL FUNCTIONS
analytical_functions
colors = make_figure_plot_colors();

%% FIGURE SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.filename               = 'Fig3_SI';
settings.saveFig                = 1;
settings.saveFolder             = '~/Dropbox/Pair_Approximation_Manuscript/MainText/PairApproximationLatex/MatlabFigures/';

settings.AxisFontColor          = [1 1 1]*0.15;
settings.AxisFontSize           = 7;
settings.LabelFontSize          = 7;

settings.alphaConfInt           = 0.3;
settings.legendLocation         = 'lowerright';
settings.legendLineWidth        = 0.1;


%% SET margins in cm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 column wide (20.5 picas / 3.42? / 8.7 cm) 
%1.5 columns wide (27 picas / 4.5? / 11.4 cm) 
%2 columns wide (42.125 picas / 7? / 17.8 cm)

%set size of figure
totHeight    = 4.5;
totWidth     = 17.8;
numCol       = 3;
numRow       = 1;

%set margins
hMargin     = 1.2;
vMargin     = 1.1;
lMargin     = 0.8;
rMargin     = 0.2;
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



%% Panel A Patch Size as function of range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataCA = load('simulation_CA_patchSize.mat');
dataCA = dataCA.data;

idx = 1;
ax(idx)      =  axes(hFig);

data = struct();

data.X{1} = dataCA.nb0;
data.X{2} = dataCA.nb0;
data.Y{1} = dataCA.patch0;
data.Y{2} = dataCA.patch1;

tickInt = 60;
maxRange = ceil(max(dataCA.nb0)/tickInt)*tickInt;
rangeTick = 0:tickInt:maxRange;


fg = struct();
fg.LineStyle                    = {'-','-','-'};
fg.LineWidth                    = {1.5,1.5,1.5}; 
fg.Color                        = {colors.DP,colors.DT}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Marker'  , 'o', ...
    'MarkerEdgeColor'  , 'none', ...
    'MarkerFaceColor'  , fg.Color{i}, ...
    'MarkerSize'  , 5, ...
    'Color', fg.Color{i})
end
       

make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0.5 maxRange+0.5];
ax(idx).YLim            = [0 10];
ax(idx).XTick           = rangeTick;    
ax(idx).YTick           = [0 5 10];
ax(idx).XLabel.String   = 'neighborhood size r';
ax(idx).YLabel.String   = 'patch size (grid points)';
%ax(idx).YLabel.String   = {'local/global fraction other','P(B|A,r_A)/P(B)'};


%% Panel C Patch size relative to interaction range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataCA = load('simulation_CA_patchSize.mat');
dataCA = dataCA.data;


idx = 2;
ax(idx)      =  axes(hFig);

data = struct();

data.X{1} = dataCA.nb0;
data.X{2} = dataCA.nb1;
data.Y{1} = dataCA.normalized_patch0;
data.Y{2} = dataCA.normalized_patch1;

tickInt = 60;
maxRange = ceil(max(dataCA.nb0)/tickInt)*tickInt;
rangeTick = 0:tickInt:maxRange;


fg = struct();
fg.LineStyle                    = {'-','-','-'};
fg.LineWidth                    = {1.5,1.5,1.5}; 
fg.Color                        = {colors.DP,colors.DT}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Marker'  , 'o', ...
    'MarkerEdgeColor'  , 'none', ...
    'MarkerFaceColor'  , fg.Color{i}, ...
    'MarkerSize'  , 5, ...
    'Color', fg.Color{i})
end
       
line([0 maxRange],[1 1],'LineStyle',':','Color','k','LineWidth',1)


make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0.5 maxRange+0.5];
ax(idx).YLim            = [0 6];
ax(idx).XTick           = rangeTick;    
ax(idx).YTick           = [0 3 6];
ax(idx).XLabel.String   = 'neighborhood size r';
ax(idx).YLabel.String   = 'patch size / interaction range';



%% Panel C Clustering simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataCA = load('simulation_CA_patchSize.mat');
dataCA = dataCA.data;


idx = 3;
ax(idx)      =  axes(hFig);

data = struct();

data.X{1} = dataCA.nb0;
data.X{2} = dataCA.nb1;
data.Y{1} = dataCA.clustering0;
data.Y{2} = dataCA.clustering1;

tickInt = 60;
maxRange = ceil(max(dataCA.nb0)/tickInt)*tickInt;
rangeTick = 0:tickInt:maxRange;


fg = struct();
fg.LineStyle                    = {'-','-','-'};
fg.LineWidth                    = {1.5,1.5,1.5}; 
fg.Color                        = {colors.DP,colors.DT}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Marker'  , 'o', ...
    'MarkerEdgeColor'  , 'none', ...
    'MarkerFaceColor'  , fg.Color{i}, ...
    'MarkerSize'  , 5, ...
    'Color', fg.Color{i})
end
       
line([0 maxRange],[1 1],'LineStyle',':','Color','k','LineWidth',1)


make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0.5 maxRange+0.5];
ax(idx).YLim            = [0.8 1.1];
ax(idx).XTick           = rangeTick;    
ax(idx).YTick           = [0.8 0.9 1 1.1];
ax(idx).XLabel.String   = 'neighborhood size r';
ax(idx).YLabel.String   = 'local/global frequency';



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



