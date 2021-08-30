function [output]  =  SteadyState_3D_RefineGrid(settings,input)
% Grid refinement routine for SOR solver
% Calculates steady state solution for 3D model of amino acid exchange
% Grid size is refined iterativaly until the solution converges
%
%% settings: model parameters
%% input: initial guess of solution
%% output: steady state solution,  
%% muAsFunctionOfGridSize: structure with growth profile for each grid size
%% finalGridScaling: scale factor for first grid at which solution is converged
%
% type 0 produces pro, growth limited by try
% type 1 produces try, growth limited by pro
% Boundary conditions: no flux on left, right, and top wall, [E] = 0 on bottom wall
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020

TalkToMe = 0;

%% first run at base grid
settings.gridScaling = settings.gridScalingBase;
gridCellTypeScaled = input.gridCellTypeScaled;
[output, gridE1, gridE2] = SteadyState_3D_SOR_Solver(settings, input, input.gridE1, input.gridE2, gridCellTypeScaled);

if ~isempty(output)
    %% Run solver, keep refining grid untill solution converges
    muConverged = 0;

    while ~muConverged % refine grid
        %store previous growth rate
        muOld = output.mu;

        %increase grid size
        settings.gridScaling = settings.gridScaleFactor*settings.gridScaling;

        %input previous solution as new initial guess
        input = output;

        %calculate output
        if settings.gridSize(3) == 1 %2D
            gridE1 = imresize(gridE1,settings.gridScaleFactor,'nearest');
            gridE2 = imresize(gridE2,settings.gridScaleFactor,'nearest');
            gridCellTypeScaled = imresize(gridCellTypeScaled,settings.gridScaleFactor,'nearest');
        else %3D
            gridE1 = imresize3(gridE1,settings.gridScaleFactor,'nearest');
            gridE2 = imresize3(gridE2,settings.gridScaleFactor,'nearest');
            gridCellTypeScaled = imresize3(gridCellTypeScaled,settings.gridScaleFactor,'nearest');
        end    

        %solve at finer grid
        [output, gridE1, gridE2] = SteadyState_3D_SOR_Solver(settings, input, gridE1, gridE2, gridCellTypeScaled);
        
        if ~isempty(output)
            %compare absolute error between old and new mu to muTolerance
            muNew = output.mu;
            muError = max(max(abs(muNew-muOld)));
            if muError>settings.muTolerance
                if TalkToMe; fprintf('SteadyState_3D_RefineGrid-> Increasing gridScaling to %i, max abs error mu = %.2g\n',2*settings.gridScaling,muError); end
                if settings.gridScaling >= settings.maxGridScaling
                    muConverged = 1;
                    fprintf('!!! SteadyState_3D_RefineGrid-> mu has not converged. ru = %.2g, rl = %.2g\n',settings.ru,settings.rl);
                end
            else
                muConverged = 1;
                if TalkToMe; fprintf('SteadyState_3D_RefineGrid-> mu has converged. GridScale used = %i, max abs error mu = %.2g\n',settings.gridScaling,muError);end
            
            end
        else
            muConverged = 1;
        end
    end
end
