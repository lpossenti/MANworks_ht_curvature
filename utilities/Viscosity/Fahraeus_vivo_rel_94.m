function visco_vivo=Fahraeus_vivo_rel_94 (radius, hematocrit)
%% visco_vivo compute the apparent relative viscosity of blood [-] in vivo using the empirical formula found by Pries et al. 1994
%
%        radius=radius of the vessel [µm]
%        hematocrit=descharge hematocrit in the vessel [-]
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 23/09/2016
%   Contact: simone.digre@gmail.com  
%%

D=2*radius; %D = diametro

mu_ast= 6 * exp (-0.085*D) + 3.2 - 2.44 * exp (-0.06 * D^0.645);

C=(0.8 + exp(-0.075 * D)) * (-1 + 1 / (1+ 10*(D/10)^12)) + 1/(1 + 10*(D/10)^12);

visco_vivo= (1 + (mu_ast - 1) * ((1-hematocrit)^C - 1)/((1-0.45)^C - 1)*(D/(D-1.1))^2)*(D/(D-1.1))^2;