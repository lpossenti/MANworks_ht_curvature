clear all
close all
clc
%% Script for the selection of the networks which has values of radius inside a prescribed range
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 10/07/2017
%   Contact: simone.digre@gmail.com  
%%

N_geo=99;
radius_string='425';
min_radius=4e-3;
max_radius=0.012;

minimo=zeros(1,N_geo+1);
massimo=zeros(1,N_geo+1);
    r=1;  
for geo=0:N_geo
    clear raggi
    clear R_L
    raggi=[];


    
    for iter=1:100

        fileName = sprintf('./Raggi_medi_Voronoi8/Reti_radius_%s/Rete_%d_Radius_%d.txt',radius_string,geo,iter-1);
        rapp = importdata(fileName);
        raggi=[raggi rapp.data];
        
    end
    
    filename_con=sprintf('./Raggi_medi_Voronoi8/Voronoi_Conn_%d.txt',geo);
    filename_nod=sprintf('./Raggi_medi_Voronoi8/Voronoi_Nodes_%d.txt',geo);
    connectivity=load(filename_con,'connectivityTOT','-ascii');
    nodes=load(filename_nod,'nodesTOT','-ascii');

    vessels_length=compute_vessels_length(connectivity, nodes);

    total_length=sum(vessels_length(:,2));

    m=0;
    M=0;
    for j=1:size(raggi,2)
        for i=1:size(raggi,1)

            R_L(i)=raggi(i,j)*vessels_length(i,2);
            ratio(i)=vessels_length(i,2)/raggi(i,j);

            if raggi(i,j)<min_radius
                m=1;
            end

            if raggi(i,j)>max_radius
                M=1;
            end

        end
    
        R_L_tot=sum(R_L);
        density=R_L_tot*2*pi/5E-4;

        Weighted_Average_Radius=R_L_tot/sum(vessels_length(:,2))*5E-4*1E6;

        if ~M && ~m
            table_row(r,1)= geo;
            table_row(r,2)= j-1;
            table_row(r,3)= total_length;
            table_row(r,4)= density;
            table_row(r,5)= min(ratio);
            table_row(r,6)= Weighted_Average_Radius;
            r=r+1;
        end

        m=0;
        M=0;
    end
    
end 

if r~=1
    R = table(table_row(:,1),table_row(:,2),table_row(:,3),table_row(:,4),table_row(:,5),table_row(:,6),'VariableNames', {'geometria' 'it' 'l_tot' 'dens' 'ratio' 'average_radius'});
    fileName = sprintf('Reti_Buone_%s.txt',radius_string);
    writetable(R,fileName,'Delimiter','\t');
end
