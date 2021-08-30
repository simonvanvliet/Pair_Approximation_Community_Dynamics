function dydt = ode_pair_approximation(t,yIn,par)
% Implemets pair approximation ODE 
% Note: time has linear rescaling
%
% Written by Simon van Vliet[4,5] & Alma Dal Co[1,2,3]
% 1: Eawag / 2: ETH Zurich / 3: Harvard / 4: UBC Vancouver / 5: Biozentrum, University Basel
%
% Initial development: 2020.01.31
% Updated 2021.08.26

%get variables
x = yIn(1); %P(DT)
y = yIn(2); %P(DP|DT, NB_T)
z = yIn(3); %P(DT|DP, NB_P)

%get constants
a = par.muMax_T;
b = par.muMax_P;
k = par.NB_T;
l = par.NB_P;


%calc change 

dx = k.^(-1).*l.^(-1).*x.*y.*(a.*l.*(1+((-1)+k).*y)+(-1).*b.*(k+k.*(( ...
 -1)+l).*z));

dy = k.^(-2).*y.*(a.*((-1)+x).^(-1).*(1+((-1)+k).*y).*((-2).*((-1)+x+ ...
 x.*y)+k.*((-1)+x+y+x.*y))+(-1).*b.*((-2)+k).*k.*l.^(-1).*((-1)+y) ...
  .*(1+((-1)+l).*z));
 
dz = k.^(-1).*l.^(-2).*((-1)+x).^(-1).*y.*(a.*((-2)+l).*l.*x.*(1+((-1) ...
  +k).*y).*((-1)+z)+(-1).*b.*k.*(1+((-1)+l).*z).*((-2).*((-1)+l).*z+ ...
  ((-2)+l).*x.*(1+z)));
 
 
dydt = [dx; dy; dz];

end