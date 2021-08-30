% Process 3D Model extract growth rates store to disk
% Input data generated with InteractionRange_2Dvs3D_RunModel
% Plot Output with InteractionRange_2Dvs3D_PlotModel
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver & Biozentrum Uni Basel, [2] Harvard University
%
% Initial development: 2021.02.03

parpool(18);

%settings
edgeToIgnore = 10;
radiusVec = 0:0.5:15;
pathname = './ClusteredGridScan_20210202/';

% handle metadata
load([pathname, 'metaData.mat'],'settings2D','settings3D','initFractions','clusterSizes');
save([pathname, 'metaDataGrowth.mat'],'settings2D','settings3D','initFractions','clusterSizes','edgeToIgnore','radiusVec');

numChambers = length(initFractions);
numCluster = length(clusterSizes);

parfor cc = 1:numChambers
    saveName = sprintf('%sGridDataChamber%03i.mat', pathname, cc);
    modelSolutions = load(saveName, 'modelSolutions');
    modelSolutions = modelSolutions.modelSolutions;

    dataScan = cell(numCluster,1);
    for rr = 1:numCluster
        data = struct();

        %process 3D
        if ~isempty(modelSolutions{rr, 1})
            muGrid = modelSolutions{rr, 1}.mu;
            typeGrid = modelSolutions{rr, 1}.gridCellType;
            %compute growth curve (mu vs other/total) for current chamber
            [mu1,fo1,mu2,fo2] = SteadyState_3D_Calc_LocalProp(muGrid, typeGrid, settings3D.gridSize, edgeToIgnore, radiusVec);

            data.mu1_3D = mu1;
            data.fo1_3D = fo1;
            data.mu2_3D = mu2;
            data.fo2_3D = fo2;
        else
            data.mu1_3D = [];
            data.fo1_3D = [];
            data.mu2_3D = [];
            data.fo2_3D = [];
        end



        %process 2D
        if ~isempty(modelSolutions{rr, 2})
            muGrid = modelSolutions{rr, 2}.mu;
            typeGrid = modelSolutions{rr, 2}.gridCellType;
            %compute growth curve (mu vs other/total) for current chamber
            [mu1,fo1,mu2,fo2] = SteadyState_3D_Calc_LocalProp(muGrid, typeGrid, settings2D.gridSize, edgeToIgnore, radiusVec);

            data.mu1_2D = mu1;
            data.fo1_2D = fo1;
            data.mu2_2D = mu2;
            data.fo2_2D = fo2;
          else
            data.mu1_2D = [];
            data.fo1_2D = [];
            data.mu2_2D = [];
            data.fo2_2D = [];
        end

        dataScan{rr} = data;
    end
    saveName = sprintf('%sGrowthDataChamber%03i.mat', pathname, cc);
    saveOutput(saveName, dataScan)
end

function saveOutput(saveName, dataScan)
    save(saveName, 'dataScan','-v7.3')
end
