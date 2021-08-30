%This code reads data files and estimates parameters from data
%It also outputs statistics from model prediction and from data
%Output used for Table S3

% Written by Simon van Vliet[4, 5] & Alma Dal Co[1,2,3]
% 1: Eawag / 2: ETH Zurich / 3: Harvard / 4: UBC Vancouver / 5: Biozentrum, University Basel
%
% Updated 2021.08.26

%SET Settings
umPerPixel = 0.041;

%init output
parameters = struct();

%% process cell geometry
%load data (if needed)
load('data_cell_geometry.mat')

%contains cell array regionProps, with cell length and width measurement for each chamber
%contains 4D matrix segImProcessed, dimension 4 = chambers
%dimension 3: 
%-1: segemnted image DP, 
%-2: segemnted image DT, 
%-3: segemnted image DP + DT, 
%-4: segemnted image of area occupied by cells (morphologicla close of layer 3) 

%init output
nExp = length(regionProps);
cellLengthVec = zeros(nExp,1);
cellWidthVec = zeros(nExp,1);
cellDensityVec = zeros(nExp,1);

%loop chambers
for pp=1:nExp
    %get number of cell
    numCel = length(regionProps{pp,1});
    %get area occupied by cells 
    areaIm = squeeze(segImProcessed(:,:,4,pp));
    areacurr = sum(areaIm(:)) * umPerPixel^2;
    
    %get cell poroperties
    lcurr = [regionProps{pp,1}.MajorAxisLength]*umPerPixel;
    wcurr = [regionProps{pp,1}.MinorAxisLength]*umPerPixel;
    
    %store average properties of chamber
    cellLengthVec(pp) = mean(lcurr);
    cellWidthVec(pp)  = mean(wcurr);
    cellDensityVec(pp) = numCel / areacurr;
end

%calc mean and sem for cell geometry properties
cell_l = mean(cellLengthVec);
se_cell_l = std(cellLengthVec) / sqrt(nExp);

cell_w = mean(cellWidthVec);
se_cell_w = std(cellWidthVec) / sqrt(nExp);

cell_dens = mean(cellDensityVec);
se_cell_dens = std(cellDensityVec) / sqrt(nExp);


%% estimate interaction range size 
%load measured interaction range
load('data_interaction_range.mat')
IRDTVec = IR_data(:,1);
IRDPVec = IR_data(:,2);

%get mean and sem over all replicates
numIRMeas = size(IR_data, 1);
IRDT = mean(IRDTVec);
se_IRDT = std(IRDTVec) / sqrt(numIRMeas);
IRDP = mean(IRDPVec);
se_IRDP = std(IRDPVec) / sqrt(numIRMeas);

%report parameters
fprintf('IR DT = %.2g +/- %.2g um\n', IRDT, se_IRDT)
fprintf('IR DP = %.2g +/- %.2g um\n', IRDP, se_IRDP)

%calc 95% CI
ci_IRP = IRDP + tinv([0.025,0.975], numIRMeas-1) * se_IRDP;
ci_IRT = IRDT + tinv([0.025,0.975], numIRMeas-1) * se_IRDT;

fprintf('IR DT = %.2g (%.2g, %.2g) um\n', IRDT, ci_IRT)
fprintf('IR DP = %.3g (%.3g, %.3g) um\n', IRDP, ci_IRP)

%% estimate maximum growth rate
%load data
load('data_growth_vs_arrangement.mat')

%SET maximum frequency upto which linear fit is done
fMax = 1;

%init output
muMax = [nan nan];
se_muMax = [nan nan];

df_mu = [nan nan];

%loop cell type, 1=DPro, 2=DTrp
for i = 1:2 
    %get data
    f = fractionVec{i};
    mu = growthVec{i};
    %clip high frequency of other
    fLow = f(f < fMax);
    muLow = mu(f < fMax);
    
    %do linear regression
    mdl = fitlm(fLow,muLow);
    c = mdl.Coefficients.Estimate;
    se_c = mdl.Coefficients.SE;
    
    ci = coefCI(mdl);
    
    ci_muMax(i,:) = ci(2,:);
    
    %assign muMax (prediction at f=1) and sem
    muMax(i) = sum(c);
    se_muMax(i) = sqrt(sum(se_c.^2));
    
    df_mu(i) = mdl.DFE;
    
end
%assign output
muMaxP = muMax(1);
muMaxT = muMax(2);

se_muMaxP = se_muMax(1);
se_muMaxT = se_muMax(2);

fprintf('mu_Max DT = %.2g +/- %.2g 1/h\n', muMaxT, se_muMaxT)
fprintf('mu_Max DP = %.2g +/- %.2g 1/h\n', muMaxP, se_muMaxP)

%calc 95% CI
ci_muMaxP = muMaxP + tinv([0.025,0.975], df_mu(1)) * se_muMaxP;
ci_muMaxT = muMaxT + tinv([0.025,0.975], df_mu(2)) * se_muMaxT;
fprintf('mu_Max DT = %.2g (%.2g, %.2g) 1/h\n', muMaxT, ci_muMaxT)
fprintf('mu_Max DP = %.2g (%.2g, %.2g) 1/h\n', muMaxP, ci_muMaxP)


growthRatio = muMaxP/muMaxT;
se_growthRatio = growthRatio * sqrt(se_muMaxP^2 / muMaxP^2 + se_muMaxT^2 / muMaxT^2);
ci_growthRatio = growthRatio + norminv([0.025,0.975]) * se_growthRatio;

fprintf('mu_Max DP / DT = %.2g (%.2g, %.2g) 1/h\n', growthRatio, ci_growthRatio)


    
%% estimate neighborhood size 
syms IR l w rho se_IR se_l se_w se_rho

%formula for number of neighbors from IR
%calc num neighbor, assume cell is rectangle with rounded caps 
%cell area = w*(l-w) + pi*(w/2)^2 
%total area = pi*(r+w/2)^2 + (2r+w)*(l-w)
%area NB = 2*r*(l-w) + pi*(r+w/2)^2 - pi*(w/2)^2 
%number neighbor is area NB * cells/area
numNB = (2*IR*(l - w) + pi*(IR + w/2)^2 - pi*(w/2)^2) * rho;

se_numNB = sqrt(diff(numNB,   l)^2 * se_l^2 + ...
                diff(numNB,   w)^2 * se_w^2 + ...
                diff(numNB, rho)^2 * se_rho^2 + ...
                diff(numNB,  IR)^2 * se_IR^2);

%set parameter substituitions
parSym = [IR l w rho...
    se_IR se_l se_w se_rho];
parValDT = [IRDT cell_l cell_w cell_dens...
    se_IRDT se_cell_l se_cell_w se_cell_dens];
parValDP = [IRDP cell_l cell_w cell_dens...
    se_IRDP se_cell_l se_cell_w se_cell_dens];

%calculate values and ranges for data
rDT = double(subs(numNB, parSym, parValDT));
se_rDT = double(subs(se_numNB, parSym, parValDT));
rDP = double(subs(numNB, parSym, parValDP));
se_rDP = double(subs(se_numNB, parSym, parValDP));

%report estimate
fprintf('#Neighbor DT = %.2g +/- %.2g\n', rDT, se_rDT)
fprintf('#Neighbor DP = %.3g +/- %.2g\n', rDP, se_rDP)


%calc 95% CI
ci_rDP = rDP + norminv([0.025,0.975]) * se_rDP;
ci_rDT = rDT + norminv([0.025,0.975]) * se_rDT;

fprintf('#Neighbor DT = %.2g (%.2g, %.2g)\n', rDT, ci_rDT)
fprintf('#Neighbor DP = %.2f (%.2f, %.2f)\n', rDP, ci_rDP)


%% estimate equilibrium frequency 
syms muA muB rA rB se_muA se_muB se_rA se_rB
fA = (muA/rA - muB/rB + muA*(rA - 2)/rA) / (muB*(rB-2)/rB + muA*(rA-2)/rA);
se_fA = sqrt(diff(fA, muA)^2 * se_muA^2 + ...
             diff(fA, muB)^2 * se_muB^2 + ...
             diff(fA,  rA)^2 * se_rA^2 + ...
             diff(fA,  rB)^2 * se_rB^2);

%set parameter substituitions
parSym = [muA muB rA rB...
    se_muA se_muB se_rA se_rB];
parVal = [muMaxT muMaxP rDT rDP...
    se_muMaxT se_muMaxP se_rDT se_rDP];

%calculate values and ranges for data
fDT = double(subs(fA, parSym, parVal));
se_fDT = double(subs(se_fA, parSym, parVal));

%report estimate
fprintf('frequency DT = %.2g +/- %.2g\n', fDT, se_fDT)

%calc 95% CI
ci_fDT = fDT + norminv([0.025,0.975]) * se_fDT;
fprintf('frequency DT = %.2g (%.2g, %.2g)\n', fDT, ci_fDT)


%% estimate equilibrium frequency WELL MIXED
syms muA muB se_muA se_muB
fA_WM = (muA) / (muB + muA);
se_fA_WM = sqrt(diff(fA_WM, muA)^2 * se_muA^2 + ...
             diff(fA_WM, muB)^2 * se_muB^2);

%set parameter substituitions
parSym = [muA muB ...
    se_muA se_muB];
parVal = [muMaxT muMaxP...
    se_muMaxT se_muMaxP];

%calculate values and ranges for data
fDT_WM = double(subs(fA_WM, parSym, parVal));
se_fDT_WM = double(subs(se_fA_WM, parSym, parVal));

%report estimate
fprintf('frequency DT = %.2g +/- %.2g\n', fDT_WM, se_fDT_WM)

%calc 95% CI
ci_fDT_WM = fDT_WM + norminv([0.025,0.975]) * se_fDT_WM;
fprintf('frequency DT Well-Mixed = %.2g (%.2g, %.2g)\n', fDT_WM, ci_fDT_WM)

%% EQ frequency data
data = load('data_community_composition.mat');
XDataNames = 'timeBetweenFramesMinutes';
YDataNames ='dTfractionInChamberAtFrame';
                              
frac_in_time   = data.dTfractionInChamberAtFrame;
time   = data.timeBetweenFramesMinutes * (-54:-54-1+size(frac_in_time,2)) / 60;

%find and store end point measurement
tEnd = 60; 
index = find(time == tEnd);
fDT_end = frac_in_time(:,index);

ci_f = @(data) mean(data) + tinv([0.025,0.975],length(data)) * std(data)/sqrt(length(data));

fDT_Data = mean(fDT_end);
ci_fDT_Data = ci_f(fDT_end);

%calc 95% CI
fprintf('relFitness data = %.2g (%.2g, %.2g)\n', fDT_Data, ci_fDT_Data)




%% estimate clustering  
syms r se_r
cluster = (r-2)/(r-1);
se_cluster = sqrt(diff(cluster, r)^2 * se_r^2);

parSym = [r se_r];
parValDT = [rDT se_rDT];
parValDP = [rDP se_rDP];

clusteringDT = double(subs(cluster, parSym, parValDT));
se_clusteringDT = double(subs(se_cluster, parSym, parValDT));
clusteringDP = double(subs(cluster, parSym, parValDP));
se_clusteringDP = double(subs(se_cluster, parSym, parValDP));

%report estimate
fprintf('clustering DT = %.2g +/- %.2g\n', clusteringDT, se_clusteringDT)
fprintf('clustering DP = %.2g +/- %.2g\n', clusteringDP, se_clusteringDP)

%calc 95% CI
ci_clusteringDT = clusteringDT + norminv([0.025,0.975]) * se_clusteringDT;
ci_clusteringDP = clusteringDP + norminv([0.025,0.975]) * se_clusteringDP;

fprintf('clustering DT = %.2g (%.2g, %.2g)\n', clusteringDT, ci_clusteringDT)
fprintf('clustering DP = %.2f (%.2f, %.2f)\n', clusteringDP, ci_clusteringDP)


%% report clustering data
load('data_cell_clustering.mat');

ci_f = @(data) mean(data) + tinv([0.025,0.975],length(data)) * std(data)/sqrt(length(data));

clusteringDT_Data = mean(data.clusteringDT);
clusteringDP_Data = mean(data.clusteringDP);


ci_clusteringDT_Data = ci_f(data.clusteringDT);
ci_clusteringDP_Data = ci_f(data.clusteringDP);

fprintf('clustering DT data = %.2g (%.2g, %.2g)\n', clusteringDT_Data, ci_clusteringDT_Data)
fprintf('clustering DP data = %.2f (%.2f, %.2f)\n', clusteringDP_Data, ci_clusteringDP_Data)



%% estimate relative fitness 
syms muA muB rA rB se_muA se_muB se_rA se_rB
%global and local frequency in space
fA = (muA/rA - muB/rB + muA*(rA - 2)/rA) / (muB*(rB-2)/rB + muA*(rA-2)/rA);
rBA = (rA - 2) / (rA - 1) * (1 - fA);
qAB = (rB - 2) / (rB - 1) * fA;
%fitness is freq. A * average payoff A + freq. B * average payoff B
W_Space = fA * muA * rBA + (1 - fA) * muB * qAB;

%frequency and fitness in wel-mixed system
fA_WM = muA / (muA + muB);
W_WM = fA_WM * muA * (1 - fA_WM) + (1 - fA_WM) * muB * fA_WM;

%relative fitness
W_rel = W_Space / W_WM;

se_W_rel = sqrt(diff(W_rel, muA)^2 * se_muA^2 + ...
                diff(W_rel, muB)^2 * se_muB^2 + ...
                diff(W_rel,  rA)^2 * se_rA^2 + ...
                diff(W_rel,  rB)^2 * se_rB^2);


%set parameter substituitions
parSym = [muA muB rA rB...
    se_muA se_muB se_rA se_rB];
parVal = [muMaxT muMaxP rDT rDP...
    se_muMaxT se_muMaxP se_rDT se_rDP];

%calculate values and ranges for data
relFitness = double(subs(W_rel, parSym, parVal));
se_relFitness = double(subs(se_W_rel, parSym, parVal));

%report estimate
fprintf('relFitness = %.2g +/- %.2g\n', relFitness, se_relFitness)

%calc 95% CI
ci_relFitness = relFitness + norminv([0.025,0.975]) * se_relFitness;
fprintf('relFitness = %.2g (%.2g, %.2g)\n', relFitness, ci_relFitness)

%process data
data = load('data_community_productivity.mat');   
relGrowth = data.muReal ./ data.muRand_mean;

ci_f = @(data) mean(data) + tinv([0.025,0.975],length(data)) * std(data)/sqrt(length(data));

relFitness_Data = mean(relGrowth);
ci_relFitness_Data = ci_f(relGrowth);

%calc 95% CI
fprintf('relFitness data = %.2g (%.2g, %.2g)\n', relFitness_Data, ci_relFitness_Data)




%% store paramters
varNames = {'cell_l','cell_w','cell_dens',...
    'IRDT', 'IRDP', 'rDT', 'rDP',...
    'muMaxT', 'muMaxP',...
    'fDT', 'clusteringDT','clusteringDP','relFitness'};
for pp = 1:length(varNames)
    parameters.(varNames{pp}) = eval(varNames{pp});
    seVarName = ['se_',varNames{pp}];
    parameters.(seVarName) = eval(seVarName);
end


save('parameters_DP_DT_community_from_data.mat','parameters')

