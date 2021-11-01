%% LOAD ANALYTICAL FUNCTIONS
analytical_functions
colors = make_figure_plot_colors();

%load PA predictions
parametersModel = load('parameters_DP_DT_community_from_literature.mat','parameters');
parametersModel = parametersModel.parameters;

%% FIGURE SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.filename               = 'Fig6';
settings.labelName              = 'A';

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
totHeight    = 6;
totWidth     = 17;
numCol       = 3;
numRow       = 2;

%set margins
hMargin     = 1.2;
vMargin     = 1;
lMargin     = 0.8;
rMargin     = 0.3;
tMargin     = 0.2;
bMargin     = 0.6;

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

%% Panel A Experimental data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 1;
ax(idx)      =  axes(hFig);

%load data
data = load('data_community_composition.mat');
XDataNames = 'timeBetweenFramesMinutes';
YDataNames ='dTfractionInChamberAtFrame';
                              
data.Y   = {data.(YDataNames), nanmean(data.(YDataNames))};
data.X   = {data.(XDataNames) * [-54:-54-1+size(data.Y{1},2)] / 60,...
            data.(XDataNames) * [-54:-54-1+size(data.Y{1},2)] / 60};

%find and store end point measurement
tEnd = 60; 
index = find(data.X{1} == tEnd);
endPointData = data.Y{1}(:,index);
    
fg = struct();
fg.LineStyle                    = {'-','-'};
fg.LineWidth                    = LineWidth; 
fg.Color                        = {[0.0 0.0 0.0 0.2], colors.red}; 
fg.Marker                       = {'none', 'none'}; 
fg.MarkerSize                   = {4, 4};
fg.MarkerFaceColor              = {'none', 'none'}; 
fg.MarkerEdgeColor              = {'k','k'}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i}, ...
    'Marker', fg.Marker{i}, ...
    'MarkerSize', fg.MarkerSize{i}, ...
    'MarkerFaceColor', fg.MarkerFaceColor{i}, ...
    'MarkerEdgeColor', fg.MarkerEdgeColor{i});
end
       
make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0 62.7];
ax(idx).YLim            = [0 1];
ax(idx).XTick           = [0 30 60];    
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'time [h]';
ax(idx).YLabel.String   = 'frequency \DeltaT';


%% Panel B PA time  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 2;
ax(idx)      =  axes(hFig);

%get data
maxT = 30;
data = dynamics_PA_numerical_solve(0.035:0.03:0.965, maxT);
                
endPointPA  = data.yData(end,:);
data.Y      = {data.yData; nanmean(data.yData, 2)};
data.X      = {data.xData, data.xData};
    
fg = struct();
fg.LineStyle                    = {'-','-'};
fg.LineWidth                    = LineWidth; 
fg.Color                        = {[0.0 0.0 0.0 0.4], colors.blue}; 
fg.Marker                       = {'none', 'none'}; 
fg.MarkerSize                   = {4, 4};
fg.MarkerFaceColor              = {'none', 'none'}; 
fg.MarkerEdgeColor              = {'k','k'}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i}, ...
    'Marker', fg.Marker{i}, ...
    'MarkerSize', fg.MarkerSize{i}, ...
    'MarkerFaceColor', fg.MarkerFaceColor{i}, ...
    'MarkerEdgeColor', fg.MarkerEdgeColor{i});
end
       
make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0 maxT];
ax(idx).YLim            = [0 1];
ax(idx).XTick           = [0 maxT/2 maxT];    
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'time [h]';
ax(idx).YLabel.String   = 'frequency \DeltaT';

%% Panel C CA data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 3;
ax(idx)      =  axes(hFig);

%load data
data = load('simulation_CA_dynamics_rev.mat');
data = data.data;

data.Y   = {data.DTinTime, nanmean(data.DTinTime, 2)};
data.X   = {((1:size(data.Y{1},1))-1)'/100,...
                ((1:size(data.Y{1},1))-1)'/100};
    
fg = struct();
fg.LineStyle                    = {'-','-'};
fg.LineWidth                    = LineWidth; 
fg.Color                        = {[0.0 0.0 0.0 0.4], colors.green}; 
fg.Marker                       = {'none', 'none'}; 
fg.MarkerSize                   = {4, 4};
fg.MarkerFaceColor              = {'none', 'none'}; 
fg.MarkerEdgeColor              = {'k','k'}; 

for i = 1:length(data.X)
    line(   data.X{i}, data.Y{i}, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i}, ...
    'Marker', fg.Marker{i}, ...
    'MarkerSize', fg.MarkerSize{i}, ...
    'MarkerFaceColor', fg.MarkerFaceColor{i}, ...
    'MarkerEdgeColor', fg.MarkerEdgeColor{i});
end
       
make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0 40];
ax(idx).YLim            = [0 1];
ax(idx).XTick           = [0 20 40];    
ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = 'time steps x1000';
ax(idx).YLabel.String   = 'frequency \DeltaT';



%% Panel D compare data and model  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 4;
ax(idx)      =  axes(hFig);

randSpread = 0.4;
medWidth = 0.4;

%get data     
%load data
dataCA = load('simulation_CA_steady_state_rev.mat');
dataCA = dataCA.data;

data.Y   = {endPointData, dataCA.PDT_end};
xPos     = [1 3];

%plot data
fg = struct();
fg.LineStyle                    = {'-','-','-'};
fg.LineWidth                    = {3, 3, 3}; 
fg.Color                        = {colors.red, colors.green, colors.blue}; 
fg.Marker                       = {'none', 'none', 'none'}; 

fg.LineStyle2                    = {'none','none','none'};
fg.Marker2                       = {'o', 'o', 'o'}; 
fg.MarkerSize2                   = {10, 10, 10};
fg.MarkerFaceColor2              = {colors.k, colors.k, colors.k}; 
fg.MarkerFaceAlpha2              = {0.4, 0.4, 0.4}; 
fg.MarkerEdgeColor2              = {'none','none','none'}; 

for i = 1:length(data.Y)
    yData = data.Y{i};
    xData = ones(size(yData))*xPos(i) + randSpread*(rand(size(yData))-0.5);
    
    hold on
    scatter(xData, yData, fg.MarkerSize2{i}, ...
        'Marker', fg.Marker2{i}, ...
        'MarkerFaceColor', fg.MarkerFaceColor2{i}, ...
        'MarkerFaceAlpha', fg.MarkerFaceAlpha2{i},...
        'MarkerEdgeColor', fg.MarkerEdgeColor2{i});

    xMed = [-medWidth medWidth] + xPos(i);
    yMed = [1 1] * mean(data.Y{i});
    
    line(xMed, yMed, ...
        'LineStyle', fg.LineStyle{i}, ...
        'LineWidth', fg.LineWidth{i}, ...
        'Color', fg.Color{i}, ...
        'Marker', fg.Marker{i});

end

%plot prediction of PA:
xPos = 2;
data.Y   = {parametersModel.fDT};

fg = struct();
fg.LineStyle                    = {'-','-','-'};
fg.LineWidth                    = {3, 0.5, 0.5}; 
fg.BarWidth                     = {0.4, 0.3, 0.3};
fg.Color                        = {colors.blue, colors.k, colors.k}; 
fg.Marker                       = {'none', 'none', 'none'}; 

for i = length(data.Y):-1:1
    hold on
    xMed = [-1 1]*fg.BarWidth{i} + xPos;
    yMed = [1 1]*data.Y{i};
    line(   xMed, yMed, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i}, ...
    'Marker', fg.Marker{i});
end
      
make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0.5 3.5];
ax(idx).YLim            = [0 1];
ax(idx).XTick           = [1 2 3];   
ax(idx).XTickLabel      = {'Data','PA','CA'};   

ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = '';
ax(idx).YLabel.String   = 'frequency \DeltaT';


%report results
dataSet = {endPointData, dataCA.PDT_end};
referenceSet = {parametersModel.fDT, parametersModel.fDT};
name = {'Data', 'CA'};

fprintf('Eq. Composition:\n')

for dd = 1:length(dataSet)

    ref = referenceSet{dd};
    data = dataSet{dd};
    meanData = mean(data);
    n = length(data);
    sem = std(data) / sqrt(n);
    [~, p] = ttest(data - ref);

    dataCI = meanData + sem * tinv([0.025 0.975],n-1);
    fprintf('%s=%.3g (CI: %.3g-%.3g), PA=%.3g, n=%i, p=%.2g\n',...
         name{dd}, meanData, dataCI(1), dataCI(2),...
         ref, n, p)
end


%% Panel E clustering  measured %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 5;
ax(idx)      =  axes(hFig);

%get measured data
load('data_cell_clustering.mat');

randSpread = 0.4;
medWidth = 0.4;
                          
data.Y      = {data.clusteringDT, data.clusteringDP};
xPos        = [1 3];

fg = struct();
fg.LineStyle                    = {'-','-','-','-'};
fg.LineWidth                    = {3, 3, 3, 3}; 
%fg.Color                        = {colors.DT, colors.DT, colors.DP, colors.DP}; 
fg.Color                        = {colors.red, colors.red, colors.green, colors.green}; 

fg.Marker                       = {'none', 'none', 'none', 'none'}; 

fg.LineStyle2                    = {'none','none','none', 'none'};
fg.Marker2                       = {'o', 'o', 'o', 'o'}; 
fg.MarkerSize2                   = {10, 10, 10, 10};
fg.MarkerFaceColor2              = {colors.k, colors.k, colors.k, colors.k}; 
fg.MarkerFaceAlpha2              = {0.4, 0.4, 0.4, 0.4}; 
fg.MarkerEdgeColor2              = {'none','none','none', 'none'}; 

for i = 1:length(data.Y)
    yData = data.Y{i};
    xData = ones(size(yData))*xPos(i) + randSpread*(rand(size(yData))-0.5);
    
    hold on
    scatter(xData, yData, fg.MarkerSize2{i}, ...
            'Marker', fg.Marker2{i}, ...
            'MarkerFaceColor', fg.MarkerFaceColor2{i}, ...
            'MarkerFaceAlpha', fg.MarkerFaceAlpha2{i},...
            'MarkerEdgeColor', fg.MarkerEdgeColor2{i});
    
    xMed = [-medWidth medWidth] + xPos(i);
    yMed = [1 1] * mean(data.Y{i});
    
    line(xMed, yMed, ...
        'LineStyle', fg.LineStyle{i}, ...
        'LineWidth', fg.LineWidth{i}, ...
        'Color', fg.Color{i}, ...
        'Marker', fg.Marker{i});

end
 
%plot prediction of PA:
xPos = [2 4];
data.Y     = {parametersModel.clusteringDT, parametersModel.clusteringDP};

fg = struct();
fg.LineStyle                    = {'-','-','-','-','-','-'};
fg.LineWidth                    = {3, 3, 0.5,3, 0.5, 0.5}; 
fg.BarWidth                     = {0.4, 0.4, 0.3,0.4, 0.3, 0.3};
fg.Color                        = {colors.blue, colors.blue, colors.k,colors.blue, colors.k, colors.k}; 
fg.Marker                       = {'none', 'none', 'none','none', 'none', 'none'}; 

for i = length(data.Y):-1:1
    hold on
    xMed = [-1 1]*fg.BarWidth{i} + xPos(i);
    yMed = [1 1]*data.Y{i};
    line(   xMed, yMed, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i}, ...
    'Marker', fg.Marker{i});
end


%add reference line
line([0 6.5],[1 1],'LineStyle',':','Color','k','LineWidth',1)

make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0.5 4.5];
ax(idx).YLim            = [0 1.1];
ax(idx).XTick           = [1 2 3 4 5 6];   
ax(idx).XTickLabel      = {'Data','PA','Data','PA'};   

ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = '';
ax(idx).YLabel.String   = 'local/global frequency';

%report results
dataSet = {data.clusteringDT, data.clusteringDP};
referenceSet = {parametersModel.clusteringDT, parametersModel.clusteringDP};
name = {'DT_Data', 'DP_Data','DT_CA', 'DP_CA'};

fprintf('Clustering:\n')

for dd = 1:length(dataSet)

    ref = referenceSet{dd};
    data = dataSet{dd};
    meanData = mean(data);
    n = length(data);
    sem = std(data) / sqrt(n);
    [~, p] = ttest(data - ref);

    dataCI = meanData + sem * tinv([0.025 0.975],n-1);
    fprintf('%s=%.3g (CI: %.3g-%.3g), PA=%.3g, n=%i, p=%.2g\n',...
         name{dd}, meanData, dataCI(1), dataCI(2),...
         ref, n, p)
end



%% Panel F fitness decrease measured %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 6;
ax(idx)      =  axes(hFig);

%get measured data
data = load('data_community_productivity.mat');


%get model prediction
randSpread = 0.4;
medWidth = 0.4;
              
xPos        = [1];     
data.Y      = {data.muReal ./ data.muRand_mean};

fg = struct();
fg.LineStyle                    = {'-','-',':'};
fg.LineWidth                    = {3, 3, 1}; 
fg.Color                        = {colors.red, colors.green, colors.k}; 
fg.Marker                       = {'none', 'none', 'none'}; 

fg.LineStyle2                    = {'none','none','none'};
fg.Marker2                       = {'o', 'o', 'o'}; 
fg.MarkerSize2                   = {10, 10, 0};
fg.MarkerFaceColor2              = {colors.k, colors.k, colors.k}; 
fg.MarkerFaceAlpha2              = {0.4, 0.4, 0.4}; 
fg.MarkerEdgeColor2              = {'none','none','none'}; 

for i = 1:length(data.Y)
    yData = data.Y{i};
    xData = ones(size(yData))*xPos(i) + randSpread*(rand(size(yData))-0.5);
    
    hold on
    scatter(xData, yData, fg.MarkerSize2{i}, ...
    'Marker', fg.Marker2{i}, ...
    'MarkerFaceColor', fg.MarkerFaceColor2{i}, ...
    'MarkerFaceAlpha', fg.MarkerFaceAlpha2{i},...
    'MarkerEdgeColor', fg.MarkerEdgeColor2{i});

    xMed = [-medWidth medWidth] + xPos(i);
    yMed = [1 1]*mean(data.Y{i});
    
    line(   xMed, yMed, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i}, ...
    'Marker', fg.Marker{i});

end

%plot prediction of PA:
xPos = 2;
data.Y     = {parametersModel.relFitness};

fg = struct();
fg.LineStyle                    = {'-','-','-'};
fg.LineWidth                    = {3, 0.5, 0.5}; 
fg.BarWidth                     = {0.4, 0.3, 0.3};
fg.Color                        = {colors.blue, colors.k, colors.k}; 
fg.Marker                       = {'none', 'none', 'none'}; 

for i = length(data.Y):-1:1
    hold on
    xMed = [-1 1]*fg.BarWidth{i} + xPos;
    yMed = [1 1]*data.Y{i};
    line(   xMed, yMed, ...
    'LineStyle', fg.LineStyle{i}, ...
    'LineWidth', fg.LineWidth{i}, ...
    'Color', fg.Color{i}, ...
    'Marker', fg.Marker{i});
end

%add reference line
line([0 4],[1 1],'LineStyle',':','Color','k','LineWidth',1)
       
make_figure_makeAxisNice(ax(idx), plotPos{idx}, settings)

ax(idx).XLim            = [0.5 2.5];
ax(idx).YLim            = [0 1.2];
ax(idx).XTick           = [1 2 3];   
ax(idx).XTickLabel      = {'Data','PA'};   

ax(idx).YTick           = [0 0.5 1];
ax(idx).XLabel.String   = '';
ax(idx).YLabel.String   = {'community productivity','realtive to well-mixed'};


%report results
dataSet = {data.muReal ./ data.muRand_mean};
referenceSet = {parametersModel.relFitness};
name = {'Data', 'CA'};

fprintf('Relative fitness:\n')

for dd = 1:length(dataSet)

    ref = referenceSet{dd};
    data = dataSet{dd};
    meanData = mean(data);
    n = length(data);
    sem = std(data) / sqrt(n);
    [~, p] = ttest(data - ref);

    dataCI = meanData + sem * tinv([0.025 0.975],n-1);
    fprintf('%s=%.3g (CI: %.3g-%.3g), PA=%.3g, n=%i, p=%.2g\n',...
         name{dd}, meanData, dataCI(1), dataCI(2),...
         ref, n, p)
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
set(gcf, 'Renderer', 'painters');
if settings.saveFig
    saveas( hFig, [settings.saveFolder settings.filename '_painters.pdf']);
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



