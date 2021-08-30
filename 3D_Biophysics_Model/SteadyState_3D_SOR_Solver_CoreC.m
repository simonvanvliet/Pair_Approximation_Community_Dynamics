function [exGridE1,exGridE2,numIt,errorCode] = SteadyState_3D_SOR_Solver_CoreC...
    (exGridE1, exGridE2, extendedGridType, ...
    gridSize, dx, boundaryType,...
    ru1,ru2,rl1,rl2,D1,D2,Iconst1,Iconst2,alpha,...
    tolerance,omega,timeOut) %#codegen

% Core solver, called from all other codes that solve 3D steady state problem
% Calculates steady state solution for 3D model of amino acid exchange
%
% two cell types:
%   type 0 produces AA1, growth limited by AA2
%   type 1 produces AA2, growth limited by AA1
%   type 0 produces pro, growth limited by try: Delta try on left
%   type 1 produces try, growth limited by pro: Delta pro on right
% internal concentration of produced AA is set to Iconst
% AA uptake is linear with [E], leakage is passive diffusion
% uptake, leakage, and diffusion rates can be assymetric between AA
% effective diffusion constant  =  (1-p)/(1+p/2), where p is density
% time scaled with growth rate (1h doubling time)
% space scaled with cell size (settings.cellSpacing um)
% concentrations scaled with Km of growth for each AA
%
% Boundary conditions: no flux on left, right, and top wall, [E] = 0 on bottom wall
%
% Solved with iterative finite difference scheme using successive over-relaxation (SOR) method
%
% Written by Simon van Vliet[1] & Alma Dal Co[2]
% [1] UBC Vancouver, [2] Harvard University
%
% Initial development: 06.02.2020
% Last update: 06.02.2020

%% Run SOR
isNotConverged = 1;
numIt = 0;
errorCode = 0;

%Grid size [Width (East-West), Depth (North-South), Height (Top-Bottom)]
numGridPointW = gridSize(1);
numGridPointD = gridSize(2);
numGridPointH = gridSize(3);


%Grid size [Width (East-West), Depth (North-South), Height (Top-Bottom)]
    
%boundaryType: 
% 0 is closed boundary (zero-flux)
% 1 os open boundary (zero concentration)
% order is: North, East, South, West, Top, Bottom
%boundaryType = [0, 0, 1, 0, 0, 0];  


tic

while isNotConverged
    
    %update grids
    prevGridE1  =  exGridE1;
    prevGridE2  =  exGridE2;
    numIt  =  numIt+1;
    
    %Implement boundary conditions
    %North boundary
    if boundaryType(1) == 0 %closed boundary, zero flux
        exGridE1(1, :, :)  =  exGridE1(3, :, :); 
        exGridE2(1, :, :)  =  exGridE2(3, :, :); 
    elseif boundaryType(1) == 1 %open boundary, zero concentration
        exGridE1(1, :, :)  =  0; 
        exGridE2(1, :, :)  =  0; 
    end
    
    %East boundary
    if boundaryType(2) == 0 %closed boundary, zero flux
        exGridE1(:, end, :)  =  exGridE1(:, end-2, :); 
        exGridE2(:, end, :)  =  exGridE2(:, end-2, :); 
    elseif boundaryType(2) == 1 %open boundary, zero concentration
        exGridE1(:, end, :)  =  0; 
        exGridE2(:, end, :)  =  0; 
    end
    
    %South boundary
    if boundaryType(3) == 0 %closed boundary, zero flux
        exGridE1(end, :, :)  =  exGridE1(end-2, :, :); 
        exGridE2(end, :, :)  =  exGridE2(end-2, :, :); 
    elseif boundaryType(3) == 1 %open boundary, zero concentration
        exGridE1(end, :, :)  =  0; 
        exGridE2(end, :, :)  =  0; 
    end
    
    %West boundary
    if boundaryType(4) == 0 %closed boundary, zero flux
        exGridE1(:, 1, :)  =  exGridE1(:, 3, :); 
        exGridE2(:, 1, :)  =  exGridE2(:, 3, :); 
    elseif boundaryType(4) == 1 %open boundary, zero concentration
        exGridE1(:, 1, :)  =  0; 
        exGridE2(:, 1, :)  =  0; 
    end
    
    if numGridPointH > 2 %System is 3D
        %Top boundary
        if boundaryType(5) == 0 %closed boundary, zero flux
            exGridE1(:, :, end)  =  exGridE1(:, :, end-2); 
            exGridE2(:, :, end)  =  exGridE2(:, :, end-2); 
        elseif boundaryType(5) == 1 %open boundary, zero concentration
            exGridE1(:, :, end)  =  0; 
            exGridE2(:, :, end)  =  0; 
        end

        %Bottom boundary
        if boundaryType(6) == 0 %closed boundary, zero flux
            exGridE1(:, :, 1)  =  exGridE1(:, :, 3); 
            exGridE2(:, :, 1)  =  exGridE2(:, :, 3); 
        elseif boundaryType(6) == 1 %open boundary, zero concentration
            exGridE1(:, :, 1)  =  0; 
            exGridE2(:, :, 1)  =  0; 
        end
    end
    
    localConvergence = 1;
    localDivergence = 0;
        
    %loop over grid
    for xx = 2:numGridPointW+1 %loop over interior, note: grid is extended by 2 pixels for boundary conditions
        for yy = 2:numGridPointD+1
            for zz = 2:numGridPointH+1
                
                %approximate Nabda E
                if numGridPointH==1 %2D
                    zz = int32(1);
                    dDif1 = exGridE1(yy-1, xx) + exGridE1(yy+1, xx) +...
                        exGridE1(yy, xx-1) + exGridE1(yy, xx+1);
                    dDif2 = exGridE2(yy-1, xx) + exGridE2(yy+1, xx) +...
                        exGridE2(yy, xx-1) + exGridE2(yy, xx+1);
                    
                    NDfactor = 4;
                else %3D
                    dDif1 = exGridE1(yy-1, xx, zz) + exGridE1(yy+1, xx, zz) +...
                        exGridE1(yy, xx-1, zz) + exGridE1(yy, xx+1, zz) +...
                        exGridE1(yy, xx, zz-1) + exGridE1(yy, xx, zz+1);
                    dDif2 = exGridE2(yy-1, xx, zz) + exGridE2(yy+1, xx, zz) +...
                        exGridE2(yy, xx-1, zz) + exGridE2(yy, xx+1, zz) +...
                        exGridE2(yy, xx, zz-1) + exGridE2(yy, xx, zz+1);   
                    NDfactor = 6;
                end
                
                %make sure concentrations do not go negative            
                E1 = exGridE1(yy, xx, zz);
                if E1<0, E1 = 0; end

                E2 = exGridE2(yy, xx, zz);
                if E2<0, E2 = 0; end

                %Calculate internal concentrations 
                type = extendedGridType(yy, xx, zz);
                
                if type == 0
                    I1 = Iconst1;
                    I2 = (E2*(ru2+rl2)-rl2)/(2+2*rl2)...
                        + sqrt( (ru2+rl2)^2*E2^2 + 2*(rl2+2)*(ru2+rl2)*E2 + rl2^2 ) / (2+2*rl2);
                else
                    I1 = (E1*(ru1+rl1)-rl1)/(2+2*rl1)...
                        + sqrt( (ru1+rl1)^2*E1^2 + 2*(rl1+2)*(ru1+rl1)*E1 + rl1^2 ) / (2+2*rl1);
                    I2 = Iconst2;
                end
                
         
                sourceE1 = alpha * E1 * (ru1+rl1)/D1 - alpha * I1 * rl1/D1;
                sourceE2 = alpha * E2 * (ru2+rl2)/D2 - alpha * I2 * rl2/D2;

                %update grid point, SOR algorithm
                exGridE1(yy, xx, zz) = (1-omega)*exGridE1(yy, xx, zz) + omega*( dDif1 - dx^2*sourceE1)/NDfactor;
                exGridE2(yy, xx, zz) = (1-omega)*exGridE2(yy, xx, zz) + omega*( dDif2 - dx^2*sourceE2)/NDfactor;

                errorE1 = abs((prevGridE1(yy, xx, zz)-exGridE1(yy, xx, zz))/exGridE1(yy, xx, zz));
                errorE2 = abs((prevGridE2(yy, xx, zz)-exGridE2(yy, xx, zz))/exGridE2(yy, xx, zz));

                if errorE1>tolerance
                    localConvergence = 0;
                end
                if errorE2>tolerance
                    localConvergence = 0;
                end

                if errorE1>1e6 || errorE2>1e6
                    localDivergence = 1;
                end
            end
        end
    end
    
    if localConvergence == 1
        isNotConverged = 0; %solved
    end
    
    if localDivergence == 1
        isNotConverged = 0; %no convergence
        errorCode = 1;
    end
    
    curT = toc;
    
    if curT > timeOut
        %timeOut give error and stop run
        errorCode = 1;
        exGridE1 = [];
        exGridE2 = [];
        isNotConverged = 0; 
    end    
end


