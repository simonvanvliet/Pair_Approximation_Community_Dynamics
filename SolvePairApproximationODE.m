function data = SolvePairApproximationODE(initFrequencyA, tMax)
% Simulates pair approximation timecourse for DPro-DTrp mutualistic community 

% Written by Simon van Vliet[4] & Alma Dal Co[1,2,3]
% 1: Eawag / 2: ETH Zurich / 3: Harvard / 4: UBC Vancouver
%
% Initial development: 2020.01.31

%SET Initial condition
pA_init = 0.5;

%SET time range to solve
tMax = 100;
numTPoints = 100;

%SET Parameters of system
parameters = struct();
parameters.muMax_P = 1;
parameters.muMax_T = 2;
parameters.NB_P = 4;
parameters.NB_T = 4;

%SET ODE solver settings
opts = odeset('RelTol',1e-8,'AbsTol',1e-4,'NonNegative', [1 2 3]);

%set time axis
tSpan = [0 tMax];
tVec = linspace(0, tMax, numTPoints);

%[x,y,z] = [ P(DT); (DP|DT, NB_T); P(DT|DP, NB_P)]
%get steady state prediction
x_init = pA_init;
%set (DP|DT, NB_T) and P(DT|DP, NB_P) to predicted values
y_init = (1-x_init) * (parameters.NB_T - 2) / (parameters.NB_T - 1);
z_init =  x_init * (parameters.NB_P - 2) / (parameters.NB_P - 1);
initState = [x_init; y_init; z_init];

%run ODE
[tOut,stateOut] = ode113(@(t, y) ODEPairApproximation(t,y,parameters), tSpan, initState, opts);    

%interpolate result
pAOut = interp1(tOut, stateOut(:,1), tVec, 'pchip');   

%plot results
figure(101)
plot(tVec, pAOut)
ylim([0 1])
 

                            