function data = dynamics_PA_numerical_solve(initFrequencyA, tMax)

% Simulates pair approximation timecourse for DPro-DTrp mutualistic community 

% Written by Simon van Vliet[4,5] & Alma Dal Co[1,2,3, 5: Biozentrum, University Basel]
% 1: Eawag / 2: ETH Zurich / 3: Harvard / 4: UBC Vancouver
%
% Initial development: 2020.01.31
% Updated 2021.08.26

%SET ODE solver settings
opts = odeset('RelTol',1e-8,'AbsTol',1e-4,'NonNegative', [1 2 3]);

%get parameters from model community
%parameters = ParameterizePairApproximation;
parametersModel = load('parameters_DP_DT_community_from_literature.mat','parameters');
parametersModel = parametersModel.parameters;

parameters = struct();
parameters.muMax_P = parametersModel.muMaxP;
parameters.muMax_T = parametersModel.muMaxT;
parameters.NB_P = round(parametersModel.rDP);
parameters.NB_T = round(parametersModel.rDT);

%get analytical functions
analytical_functions      

%set time axis
tSpan = [0 tMax];
numTPoints = 100;
tVec = linspace(0, tMax, numTPoints);

%init output
numRuns = length(initFrequencyA);
trajectories = nan(numTPoints, numRuns); 


for tt = 1:numRuns
    
    %[x,y,z] = [ P(DT); (DP|DT, NB_T); P(DT|DP, NB_P)]
    %get steady state prediction
    x_init = initFrequencyA(tt);
    %set (DP|DT, NB_T) and P(DT|DP, NB_P) to predicted values
    y_init = (1-x_init) * (parameters.NB_T - 2) / (parameters.NB_T - 1);
    z_init =  x_init * (parameters.NB_P - 2) / (parameters.NB_P - 1);
    initState = [x_init; y_init; z_init];
          
    %run ODE
    [tOut,stateOut] = ode113(@(t, y) ode_pair_approximation(t,y,parameters), tSpan, initState, opts);    
    
    %interpolate result
    pAOut = interp1(tOut, stateOut(:,1), tVec, 'pchip');                                      
    trajectories(:, tt) = pAOut;
    
end

    
data = struct();
data.yData = trajectories;
data.xData = tVec;

                            