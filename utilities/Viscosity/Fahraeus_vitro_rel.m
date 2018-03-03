function visco_vitro=Fahraeus_vitro_rel (radius, hematocrit)
%% visco_vitro compute the apparent relative viscosity of blood [-] in vitro using the empirical formula found by Pries et al. 1992
%
%        radius=radius of the vessel [µm]
%        hematocrit=descharge hematocrit in the vessel [-]
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 23/09/2016
%   Contact: simone.digre@gmail.com  
%%
D=2*radius;

mu_ast= 220 * exp (-1.3*D) + 3.2 - 2.44 * exp (-0.06 * D^0.645);

C=(0.8 + exp(-0.075 * D)) * (-1 + 1 / (1+ 10*(D/10)^12)) + 1/(1 + 10*(D/10)^12);

visco_vitro= 1 + (mu_ast - 1) * ((1-hematocrit)^C - 1)/((1-0.45)^C - 1);
