% Plot output of 3D Model
% Input data generated with InteractionRange_2Dvs3D_ProcessModel
%
% Written by Simon van Vliet & Alma Dal Co
% Eawag & ETH Zurich
%
% Initial development: 12.6.2018
% Last update: 2021.01.30

% load default settings


savePath = './ClusteredGridScan_20210202/';
load([savePath, 'metaDataGrowth.mat'],'settings2D','settings3D','initFractions','clusterSizes','edgeToIgnore','radiusVec');

numChambers = length(initFractions);
numCluster = length(clusterSizes);
nRep = 1;
numChmbrPerRep = numChambers/nRep;

reprocess = 1;

if exist([savePath, 'interactionRange.mat'],'file')==2 && ~reprocess
    load([savePath, 'interactionRange.mat'])

else
    %init output
    intRange = zeros(numCluster, nRep, 4);
    intRangeAllData = cell(numCluster, nRep, 2);

    dataFiles = cell(numChambers,1);

    for cc = 1:numChambers
        saveName = sprintf('%sGrowthDataChamber%03i.mat', savePath, cc);
        dataFiles{cc} = matfile(saveName);
    end


    for rr = 1:numCluster
        
        for nn =1:nRep
            mu1_3D = [];
            mu2_3D = [];
            fo1_3D = [];
            fo2_3D = [];

            mu1_2D = [];
            mu2_2D = [];
            fo1_2D = [];
            fo2_2D = [];

            for cc = 1:numChmbrPerRep
                
                chmrId = (cc-1)*nRep + nn;
                
                fprintf('%i-%i\n',rr,chmrId);

                %load(saveName, 'dataScan');
                dataCell = dataFiles{chmrId}.dataScan(rr,1);
                data = dataCell{1,1};

                mu1_3D = cat(1,mu1_3D,data.mu1_3D);
                mu2_3D = cat(1,mu2_3D,data.mu2_3D);
                fo1_3D = cat(1,fo1_3D,data.fo1_3D);
                fo2_3D = cat(1,fo2_3D,data.fo2_3D);

                mu1_2D = cat(1,mu1_2D,data.mu1_2D);
                mu2_2D = cat(1,mu2_2D,data.mu2_2D);
                fo1_2D = cat(1,fo1_2D,data.fo1_2D);
                fo2_2D = cat(1,fo2_2D,data.fo2_2D);


            end

            IR_3D = SteadyState_3D_GetInteractionRange(mu1_3D, mu2_3D, fo1_3D, fo2_3D, radiusVec);
            IR_2D = SteadyState_3D_GetInteractionRange(mu1_2D, mu2_2D, fo1_2D, fo2_2D, radiusVec);

            intRange(rr, nn, :) = [IR_2D.interactionRange1 IR_2D.interactionRange2 ...
                IR_3D.interactionRange1 IR_3D.interactionRange2];

            intRangeAllData{rr,nn,1} = IR_2D;
            intRangeAllData{rr,nn,2} = IR_3D;
        end

    end

    save([savePath, 'interactionRange.mat'],'intRangeAllData','intRange')
end

%%
[rangeA, rangeB] = GR_analytical_fromSettings(settings3D);
cmap=lines(4);

figure(101)
clf
hold on
for i=1:4
    data = squeeze(intRange(:,:,i))'*settings3D.cellSpacing;
    if nRep>1
        [meanData,ciUp,ciDown,ciDev]=createConfInterval(data);
        p(i) = plot(clusterSizes, meanData,'-','LineWidth',2,'Color',cmap(i,:));
        x = [clusterSizes, fliplr(clusterSizes)];
        inBetween = [ciUp, fliplr(ciDown)];
        fill(x, inBetween,cmap(i,:),'faceAlpha',0.2,'EdgeColor','none');
    else
        p(i) = plot(clusterSizes, data,'-','LineWidth',2,'Color',cmap(i,:));
    end
end    
    
plot([0 max(clusterSizes)],settings3D.cellSpacing*rangeA*[1 1],'-k','LineWidth',2);
legend(p,{'2D, type A','2D, type B','3D, type A','3D, type B','Analytical pred.'})
xlabel('Cluster size [a.u.]')
ylabel('Range [um]')
ylim([0 15])
set(gca,'YTick',0:5:15,'XTick',0:20:80)
