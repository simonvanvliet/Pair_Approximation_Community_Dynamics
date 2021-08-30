# Code belonging to: Global dynamics of microbial communities emerge from local interaction rules
## Simon van Vliet, Christoph Hauert, Martin Ackermann, and Alma Dal Co
### simon.vanvliet@unibas.ch or alma.dalco@gmail.com

## Contents
### Main folder
* PairApprox_200308_General.nb 
  * Mathematica notebook containing derivation of analytical expressions
* explore_model_temporal_dynamics.m
  * Matlab code to numerically solve Pair Approximation dynamical equations

### Matlab_Figure_Code folder
Contains Matlab scripts (.m) and data-files (.mat) to reproduce all main-text figures
* make_figure_X.m 
  * run code to recreate figure X

### 3D_Biophysics_Model
Contains Matlab script to solve individual based biophysical model on 2D and 3D grid and compare the interaction range
Run code in following order to reproduce Fig S1
1. InteractionRange_2Dvs3D_RunModel.m
2. InteractionRange_2Dvs3D_ProcessModel.m
3. InteractionRange_2Dvs3D_PlotModel.m 

The remaining code ("SteadyState_3D_*.m") implements the 3D boundary-value problem that is used to find steady-state growth in 3D grid of cells