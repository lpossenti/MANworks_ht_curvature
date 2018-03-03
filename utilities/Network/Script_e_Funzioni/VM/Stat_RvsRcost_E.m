clear all
clc
close all

radii=[];
N_Reti=100;
Rcost=8e-3;

for i=1:N_Reti

fileName = sprintf('./Geometrie_Test_Density/Esagono_4m/Radius_%d.txt',i-1);
rapp = importdata(fileName);
radii=[radii rapp.data];
end

filename_con=sprintf('./Geometrie_Test_Density/Esagono_4m/Esagono_Conn.txt');
filename_nod=sprintf('./Geometrie_Test_Density/Esagono_4m/Esagono_Nodes.txt');
connectivity=load(filename_con,'connectivity','-ascii');
nodes=load(filename_nod,'nodes','-ascii');

%find lengths of vessels
vessels_length=compute_vessels_length(connectivity, nodes);

%find R*L
for j=1:size(radii,2)
	for i=1:size(radii,1)
        R_L(i,j)=radii(i,j)*vessels_length(i,2);
	end
end
    
Sum_R_L=zeros(1,size(radii,2));
    
for j=1:size(radii,2)
	for i=1:size(radii,1)
        Sum_R_L(1,j)=R_L(i,j)+Sum_R_L(j);
	end
end

%find Rcost*L
Sum_Rcost_L=sum(vessels_length(:,2))*Rcost;

%find ratio R*L/(Rcost*L)
ratioRvsRcost=Sum_R_L./Sum_Rcost_L;

%percentage difference
ratioRvsRcost=(ratioRvsRcost-1)*100;

mean_ratio=mean(ratioRvsRcost);
sd_ratio=std(ratioRvsRcost);

%table
R = table(ratioRvsRcost);
fileName = sprintf('rapportiR_Esagono.txt');
writetable(R,fileName,'Delimiter','\t');

%figure + saving images
figure
histogram(ratioRvsRcost,'Normalization','probability') 
titleStr = sprintf('Variazione densità [%] - Esagono Singolo');
title(titleStr)
xlabel('Variazione della densità [%]')
savefig(titleStr)
figureName = sprintf('Rapporto_densità_Esagono.png');
saveas(gcf,figureName)
   
normalized_ratio=(ratioRvsRcost-mean_ratio)/sd_ratio;
norm_tot= kstest(normalized_ratio);
