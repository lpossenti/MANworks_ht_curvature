clear all
close all
clc

R=[0.5:0.1:1000]; %radius in [µm]
T=37; %temperature in °C
visco_plasma=1.808/(1+0.0337*T+0.00022*T*T)*1.8; %plasma viscosity [cP] using the formula from Biomachines course

for i=1:size(R,2)
    %in vivo viscosity
    visco_vivo_15(i)=Fahraeus_vivo_rel_94 (R(i), 0.15);
    visco_vivo_30(i)=Fahraeus_vivo_rel_94 (R(i), 0.30);
    visco_vivo_35(i)=Fahraeus_vivo_rel_94 (R(i), 0.35); 
    visco_vivo_40(i)=Fahraeus_vivo_rel_94 (R(i), 0.40);
    visco_vivo_45(i)=Fahraeus_vivo_rel_94 (R(i), 0.45);
    visco_vivo_50(i)=Fahraeus_vivo_rel_94 (R(i), 0.50);
    visco_vivo_60(i)=Fahraeus_vivo_rel_94 (R(i), 0.60);  

end

D=R*2;
% lol=visco_vitro_45';
figure
semilogx(D,visco_vivo_15,D,visco_vivo_35,D,visco_vivo_45,D,visco_vivo_60);
Lhandle=legend('Ht=15','Ht=35','Ht=45','Ht=60');
axis([0 2000 0 16])
xhandle=xlabel('Diameter [µm]');
yhandle=ylabel('Relative Apparent Viscosity [-]');
set(xhandle,'Fontsize',12)
set(xhandle,'Fontname','Times New Roman')
set(yhandle,'Fontsize',12)
set(yhandle,'Fontname','Times New Roman')
set(Lhandle,'Fontsize',10)
set(Lhandle,'Fontname','Times New Roman')

Viscosity=Fahraeus_vivo_rel_94(0.003*5e-4*1e6, 0.45)*visco_plasma; %if you want to compute viscosity in cP
