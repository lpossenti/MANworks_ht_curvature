function visco_vivo=Fahraeus_vivo_rel (radius, hematocrit)
%% visco_vivo compute the apparent relative viscosity of blood [-] in vivo using the empirical formula found by Pries et al. 2005
%
%        radius=radius of the vessel [µm]
%        hematocrit=descharge hematocrit in the vessel [-]
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 23/09/2016
%   Contact: simone.digre@gmail.com  
%%

D=2*radius; %D = diametro
Doff=2.4;
Dcrit=10.5;
D_50=100;
Eamp=1.1;
Ewidth=0.03;
Epeak=0.6;
E_Ht=1.18;

    if(D <= Doff)
        Was=0;
    else
        Was= (D-Doff)/(D+D_50-2*Doff);
    end
    
    if(D<=Doff)
        Wpeak=0;
    elseif(D>Doff && D<=Dcrit)
        Wpeak=Eamp*(D-Doff)/(Dcrit-Doff);
    else
        Wpeak=Eamp*exp(-Ewidth*(D-Dcrit));
    end
    
Wph=Was+Wpeak*Epeak;

Weff=Was+Wpeak*(1+hematocrit*E_Ht);

Dph=D-2*Wph;

Deff=D-2*Weff;

mu_ast= 220 * exp (-1.3*Dph) + 3.2 - 2.44 * exp (-0.06 * Dph^0.645);

C=(0.8 + exp(-0.075 * Dph)) * (-1 + 1 / (1+ 10*(Dph/10)^12)) + 1/(1 + 10*(Dph/10)^12);

visco_vitro= 1 + (mu_ast - 1) * ((1-hematocrit)^C - 1)/((1-0.45)^C - 1);

visco_vivo=visco_vitro*(D/Deff)^4;