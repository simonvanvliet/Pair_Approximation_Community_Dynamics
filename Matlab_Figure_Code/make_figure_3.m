%% LOAD ANALYTICAL FUNCTIONS
analytical_functions
colors = make_figure_plot_colors();

%% FIGURE SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.filename               = 'Fig3';
settings.saveFig                = 1;
settings.saveFolder             = '~/Dropbox/Pair_Approximation_Manuscript/MainText/PairApproximationLatex/MatlabFigures/';

settings.AxisFontColor          = [1 1 1]*0.15;
settings.AxisFontSize           = 7;
settings.LabelFontSize          = 7;

settings.alphaConfInt           = 0.3;
settings.legendLocation         = 'lowerright';
settings.legendLineWidth        = 0.1;

maxRange                        = 21;
rangeTick                       = [3:6:21];

%% SET margins in cm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1 column wide (20.5 picas / 3.42? / 8.7 cm) 
%1.5 columns wide (27 picas / 4.5? / 11.4 cm) 
%2 columns wide (42.125 picas / 7? / 17.8 cm)

%set size of figure
totHeight    = 7;
totWidth     = 8.7;
numCol       = 2;
numRow       = 2;

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



%% Panel A EQ freq as function of growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 1;
ax(idx)      =  axes(hFig);
hold on

relMu = logspace(-1,1,1000);
rVec = [3,10];

data = struct();
for rr = 1:length(rVec)
    
    pA_bound = pA(relMu, 1, rVec(rr), rVec(rr));
    in_bounds = pA_bound>0 & pA_bound<1;
    
    data.X{rr} = relMu(in_bounds);
    data.Y{rr} = pA_bound(in_bounds);
    
end

data.X{length(rVec)+1} = relMu;
data.Y{length(rVec)+1} = pAWM(relMu, 1);
   
fg = struct();
fg.LineStyle                    = {'-','-','-'};
fg.LineWidth                    = {1.5,1.5,1.5}; 
fg.Color                        = {colors.green, colors.blue, colors.grey}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i})
end

[leg, obj] = legend({'','',''});
       
make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)
make_figure_makeLegendNice(leg, obj, ax(idx), settings);


ax(idx).XScale            = 'log';

ax(idx).XLim            = [0.1 10];
ax(idx).YLim            = [0 1];
ax(idx).XTick           = [0.1 0.2 0.5 1 2 5 10];    
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'growth ratio             ';
ax(idx).YLabel.String   = 'frequency type A';


%% Panel B EQ freq as function of range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 2;
ax(idx)      =  axes(hFig);

relMu = [0.5 1 2];
rVec = 3:maxRange;
rA = 3; 

data = struct();
for mm = 1:length(relMu)
    data.X{mm} = rVec/rA;
    data.Y{mm} = pA(1, relMu(mm), rA, rVec);
    data.X{mm+length(relMu)} = rVec/rA;
    data.Y{mm+length(relMu)} = pAWM(1, relMu(mm))*ones(size(rVec));
    
    
    
    
end

fg = struct();
fg.LineStyle                    = {'-','-','-',':',':',':'};
fg.LineWidth                    = {1.5,1.5,1.5,1.5,1.5,1.5}; 
fg.Color                        = {colors.green, colors.grey, colors.red,...
                                    colors.green, colors.grey, colors.red}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i})
end

% [leg, obj] = legend(compose('\\mu_B/\\mu_A=%.1f',relMu));
% make_figure_makeLegendNice(leg, obj, ax(idx), settings)

make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)


ax(idx).XLim            = [3/rA maxRange/rA];
ax(idx).YLim            = [0 1];
ax(idx).XTick           = rangeTick/rA;    
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'neighborhood size ratio r_B/r_A';
ax(idx).YLabel.String   = 'frequency type A';


%% Panel C clustering as function of range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 3;
ax(idx)      =  axes(hFig);

clustering = @(r) (r - 2) ./ (r - 1);
rVec = 2:maxRange;

data = struct();

data.X{1} = rVec;
data.Y{1} = clustering(rVec);

fg = struct();
fg.LineStyle                    = {'-','-','-'};
fg.LineWidth                    = {1.5,1.5,1.5}; 
fg.Color                        = {colors.grey}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i})
end
       
line([0 maxRange],[1 1],'LineStyle',':','Color','k','LineWidth',1)


make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0 maxRange];
ax(idx).YLim            = [0 1.2];
ax(idx).XTick           = rangeTick;    
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'neighborhood size r';
ax(idx).YLabel.String   = 'local/global frequency';
%ax(idx).YLabel.String   = {'local/global fraction other','P(B|A,r_A)/P(B)'};



%% Panel D rel fitness as function of range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 4;
ax(idx)      =  axes(hFig);

rVec = 2:maxRange;



relMu = [0.5 1 2];

data = struct();
for mm = 1:length(relMu)
    data.X{mm} = rVec;
    data.Y{mm} = Wrel(1, relMu(mm), rVec, rVec);
end


data.X{4} = [0 max(rVec)];
data.Y{4} = [1 1];

fg = struct();
fg.LineStyle                    = {'-','-','--',':','-'};
fg.LineWidth                    = {1.5,1.5,1.5,1,1.5}; 
fg.Color                        = {colors.green, colors.grey, colors.red, colors.k}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i})
end
       
make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0 maxRange];
ax(idx).YLim            = [0 1.2];
ax(idx).XTick           = rangeTick;    
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'neighborhood size r';
ax(idx).YLabel.String   = {'community growth rate','realtive to well-mixed'};




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



