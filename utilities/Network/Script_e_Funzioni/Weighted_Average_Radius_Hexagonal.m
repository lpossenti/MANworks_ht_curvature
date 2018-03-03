% script for evaluation of weighted average radius in hexagonal geometry after computation of radius of vessels.
clear all
clc
close all

value_of_radius='4';

radii=[];

for iter=1:100

    fileName = sprintf('./Geometrie_Test_Density/Esagono_%s/Radius_%d.txt',value_of_radius,iter-1);
    rapp = importdata(fileName);
    radii=[radii rapp.data];
    
end

filename_con=sprintf('./Geometrie_Test_Density/Esagono_%s/Esagono_Conn.txt',value_of_radius);
filename_nod=sprintf('./Geometrie_Test_Density/Esagono_%s/Esagono_Nodes.txt',value_of_radius);
connectivityTOT=load(filename_con,'connectivityTOT','-ascii');
nodesTOT=load(filename_nod,'nodesTOT','-ascii');

vessel_length=zeros(size(connectivityTOT,1),2);

vessel_length=compute_vessels_length(connectivityTOT, nodesTOT);

    for j=1:size(radii,2)
        for i=1:size(radii,1)

            R_L(i,j)=radii(i,j)*vessel_length(i,2); %product R*L

        end
    end
    
    
R_L_tot=sum(R_L,1);
       
Weighted_Average_Radius=R_L_tot/sum(vessel_length(:,2))*5E-4*1E6;

median_radius=median(Weighted_Average_Radius);
mean_radius=mean(Weighted_Average_Radius);
sd_radius=std(Weighted_Average_Radius);
        
normalized_radius=(Weighted_Average_Radius-mean_radius)/sd_radius;
norm_test= kstest(normalized_radius);
hist(Weighted_Average_Radius)

