% Literature and measured paramters for DPro-DTrp mutualistic community  
%
% AA1=pro, AA2=trp
% Type0 = DTry,  Type1 = DPro
% type 0 produces pro, growth limited by try
% type 1 produces try, growth limited by pro
%
% Written by Simon van Vliet[4, 5] & Alma Dal Co[1,2,3]
% 1: Eawag / 2: ETH Zurich / 3: Harvard / 4: UBC Vancouver / 5: Biozentrum, University Basel
%
% Updated 2021.08.26


lastUpdate='2021.08.26';
fprintf('August 2021 settings\n...last update %s\n',lastUpdate);
load('parameters_DP_DT_community_from_data.mat')


settings=struct;

%non-scaled rates in sec
tD_s = 2790; %s (46.5min, mu=1.29)
ru_s = 7.0; % 1/s
rl_s = 3.1E-6; % 1/s
D0_s = 7.61E2; %um^2/s

%rescaled rates
ru = ru_s * tD_s;
rl = rl_s * tD_s;
D0 = D0_s * tD_s;

%% SET  model parameters
settings.Iconst1                = 20; %Concentration of produced amino-acid 
settings.Iconst2                = 20; %Concentration of produced amino-acid

settings.D0                     = D0; %Diffusion constant in empty space
settings.cellSpacing            = 1.5; %Average cell spacing, used to rescale D0

settings.ru                     = ru;
settings.rl                     = rl; 

settings.delta_u                = 11.8; %ru2/ru1 = delta_u
settings.delta_l                = 1/26.3; %rl2/rl1 = delta_u
settings.delta_D                = 0.75; %D2/D1 = delta_u

settings.tD_s                   = tD_s;
settings.mu_wt                  = 3600 / tD_s;

%% SET MEASURED PARAMETERS
settings.rho                    = 0.65; %Density
settings.IR_P                   = 12.1; %measured interaction range in um for DPro
settings.IR_T                   = 3.2; %measured interaction range in um for DTrp
settings.beta                  = 0.88; %slope between IR and GR

% fprintf('...rho=%#.2g ru=%#.2g rl=%#.2g d_D=%#.2g d_ru=%#.2g d_rl=%#.2g \n',...
%     settings.rho,settings.ru,settings.rl,settings.delta_D,settings.delta_u,settings.delta_l);


%% Caculate growth range and neighrbood size
growth_range = @(r0, IC, rl, beta) settings.beta * r0 .* log(1/2 * rl .*IC .* (1+sqrt(1 + 8 ./ (rl .* IC))) + 4);

Deff = (settings.D0 ...
    *(1-settings.rho)/(1+settings.rho/2)....
    *(1-settings.rho)/settings.rho);

rlP = settings.rl/sqrt(settings.delta_l);
rlT = settings.rl*sqrt(settings.delta_l);

ruP = settings.ru/sqrt(settings.delta_u);
ruT = settings.ru*sqrt(settings.delta_u);

DP = Deff/sqrt(settings.delta_D);
DT = Deff*sqrt(settings.delta_D);

r0P = sqrt(DP ./ ((ruP + rlP)));
r0T = sqrt(DT ./ ((ruT + rlT)));

IcT = settings.Iconst1;
IcP = settings.Iconst2;


settings.rangePsimple = growth_range(r0P, IcP, rlP, settings.beta);
settings.rangeTsimple = growth_range(r0T, IcT, rlT, settings.beta);

%Convert no numer of neighbors
l = parameters.cell_l;
w = parameters.cell_w;
rho = parameters.cell_dens;

NBf = @(IR) (2*IR*(l - w) + pi*(IR + w/2)^2 - pi*(w/2)^2) * rho;

rDP = NBf(settings.rangePsimple);
rDT = NBf(settings.rangeTsimple);


fprintf('Interaction range DT (predicted, measured) = %.2g, %.2g\n', settings.rangeTsimple, settings.IR_T)
fprintf('Interaction range DP (predicted, measured) = %.3g, %.3g\n', settings.rangePsimple, settings.IR_P)

fprintf('# neighbor DT (predicted, measured) = %.2g, %.2g\n', rDT, parameters.rDT)
fprintf('# neighbor DP (predicted, measured) = %.3g, %.3g\n', rDP, parameters.rDP)


%% Caculate max growth rate 
IofEps = @(rl, Eps) (rl * (Eps - 1) + sqrt(rl^2 * (Eps + 1)^2 + 4 * rl * Eps)) / (2 * (1 + rl));
muMaxf = @(rl, Ic) IofEps(rl, Ic) / (1 + IofEps(rl, Ic));
muMaxSf = @(rl, Ic) (rl * Ic / 2) * (sqrt(1 + 4 / (rl * Ic)) - 1);

settings.muMaxPfull = muMaxf(rlP, IcP)*settings.mu_wt ;
settings.muMaxTfull = muMaxf(rlT, IcT)*settings.mu_wt ;

muMaxP = muMaxSf(rlP, IcP)*settings.mu_wt ;
muMaxT = muMaxSf(rlT, IcT)*settings.mu_wt ;

fprintf('mu_Max DT (predicted, measured) = %.2g, %.2g  1/h\n', muMaxT, parameters.muMaxT)
fprintf('mu_Max DP (predicted, measured) = %.2g, %.2g  1/h\n', muMaxP, parameters.muMaxP)
fprintf('mu_Max DP/DT (predicted, measured) = %.2g, %.2g 1/h\n', muMaxP/muMaxT, parameters.muMaxP/parameters.muMaxT)



%% estimate equilibrium frequency 
syms muA muB rA rB
fA = (muA/rA - muB/rB + muA*(rA - 2)/rA) / (muB*(rB-2)/rB + muA*(rA-2)/rA);

%set parameter substituitions
parSym = [muA muB rA rB];
parVal = [muMaxT muMaxP rDT rDP];

%calculate values and ranges for data
fDT = double(subs(fA, parSym, parVal));

%report estimate
fprintf('frequency DT = %.2g\n', fDT)


%% Estimate EQ frequency WM model


fprintf('frequency DT well-mixed = %.2g\n', muMaxT/(muMaxT+muMaxP))


%% estimate clustering  
syms r 
cluster = (r-2)/(r-1);

parSym = [r];
parValDT = [rDT];
parValDP = [rDP];

clusteringDT = double(subs(cluster, parSym, parValDT));
clusteringDP = double(subs(cluster, parSym, parValDP));

%report estimate
fprintf('clustering DT = %.2g\n', clusteringDT)
fprintf('clustering DP = %.2g\n', clusteringDP)


%% estimate relative fitness 
syms muA muB rA rB 
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


%set parameter substituitions
parSym = [muA muB rA rB];
parVal = [muMaxT muMaxP rDT rDP];

%calculate values and ranges for data
relFitness = double(subs(W_rel, parSym, parVal));

%report estimate
fprintf('relFitness = %.2g\n', relFitness)


%% store paramters
parameters = struct();
varNames = {'rDT', 'rDP',...
    'muMaxT', 'muMaxP',...
    'fDT', 'clusteringDT','clusteringDP','relFitness'};
for pp = 1:length(varNames)
    parameters.(varNames{pp}) = eval(varNames{pp});
end


save('parameters_DP_DT_community_from_literature.mat','parameters')

