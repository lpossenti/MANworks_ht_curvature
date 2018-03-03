clear all
clc
close all
%% Script for creating a network made by N random Voronoi networks, chosen by previous selection criteria
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 10/07/2017
%   Contact: simone.digre@gmail.com  
%%
N_geo=99;
n_reti=18;
d=5E-4;

fileName = sprintf('RETI_FINALI.txt');
good_network_struct = importdata(fileName);
good_networks=good_network_struct.data;


starting_radius=[4 4.25 4.5];
n_good=zeros(N_geo+1,size(starting_radius,2));
nr_good=zeros(size(starting_radius,2),1);
buone=0;
n_r=0;
selezionate=[];

for raggi_in=1:size(starting_radius,2)
        for check=1:size(good_networks,1)
            if good_networks(check,1)==starting_radius(raggi_in)
            nr_good(raggi_in)=nr_good(raggi_in)+1;
            end
        end
end

a=1;
for r=1:size(starting_radius,2)
    for geo=0:N_geo
        for check=a:(nr_good(r)+a-1)
            if good_networks(check,2)==geo
            n_good(geo+1,r)=n_good(geo+1,r)+1;
            no_zero=1;
            end
        end
    end
    a=a+nr_good(r);
end

starting_radius_selection=randi([1,size(starting_radius,2)],[n_reti,1]); %choose n_reti values for starting_radius

for i=1:n_reti
    ok=0;

    while ~ok
        ok_first=0;
        geometrie(i)=(randi([1,N_geo+1],1)); %choose a random good network for the specified starting_radius
        if n_good(geometrie(i),starting_radius_selection(i))~=0
            ok_first=1;
        end
        
% % %  Found different networks
        j=1;
        trovato_uguale=0;
        if ok_first
            while j<i && ~trovato_uguale
                if geometrie(i)==geometrie(j)
                    trovato_uguale=1;
                end
                j=j+1;
            end
        end
        ok=~trovato_uguale && ok_first;
% % %  FINE TROVARE GEOMETRIE DIVERSE
    end
    
end

%%% SE VUOI AVERE GEOEMTRIE UGUALI, MA ALMENO RAGGI DIVERSI
% % for i=1:n_reti
% % position(i)=1;
% % trovato=0;
% % % different=0;
% % %     while~different
% %         while ~trovato
% %             if reti_buone(position(i),1)==raggi_iniziali(raggi_part(i))
% %                 trovato=1;
% %             else
% %                 position(i)=position(i)+1;
% %             end
% %         end
% %      trovato_geo=0;  
% %          while ~trovato_geo
% %             if reti_buone(position(i),2)==geometrie(i)-1
% %                 trovato_geo=1;
% %             else
% %                 position(i)=position(i)+1;
% %             end
% %          end
% % 
% %          iterazione(i)=randi([1,n_good(geometrie(i),raggi_part(i))],1);
% %          position(i)=position(i)+iterazione(i);
% % %          for j=1:i
% % %              if position(i)==position(j) && i~=j
% % %                  different=0;
% % %              else
% % %                  different=1;
% % %              end
% % %          end
% % %     end
% % end 

% % %  SE VUOI POTER AVERE ANCHE LA STESSA GEOEMTRIA CON LO STESSO RAGGIO,
% OPPURE SE HAI GIA' CONTROLLATO CHE LA GEOMETRIA E' DIVERSA
for i=1:n_reti
position(i)=1;
trovato=0;
        while ~trovato
            if good_networks(position(i),1)==starting_radius(starting_radius_selection(i))
                trovato=1;
            else
                position(i)=position(i)+1;
            end
        end
trovato_geo=0;
         while ~trovato_geo
            if good_networks(position(i),2)==geometrie(i)-1
                trovato_geo=1;
            else
                position(i)=position(i)+1;
            end
         end

         iterazione(i)=randi([1,n_good(geometrie(i),starting_radius_selection(i))],1);
         position(i)=position(i)+iterazione(i);
         

end 

Total_Lateral_Surface=0;

for i=1:n_reti
Total_Lateral_Surface=Total_Lateral_Surface+good_networks(position(i),end);
end

density=Total_Lateral_Surface*2*pi/d*1E-6/d;

DOMAIN=1;
pos=0;
delta=DOMAIN/(n_reti);
Nnodi_max=0;
connectivityTOT=zeros(1,3);
nodesTOT=zeros(1,6);
nodes3DTOT=zeros(1,6);
raggi=[];

for n=1:n_reti
    
    if starting_radius(starting_radius_selection(n))==4
        string='4';
    elseif starting_radius(starting_radius_selection(n))==4.25
        string='425';
    elseif starting_radius(starting_radius_selection(n))==4.5
        string='4m';
    end
        
        
filename_con=sprintf('./Raggi_medi_Voronoi8/Voronoi_Conn_%d.txt',geometrie(n)-1);
filename_nod=sprintf('./Raggi_medi_Voronoi8/Voronoi_Nodes_%d.txt',geometrie(n)-1);
filename_raggio=sprintf('./Raggi_medi_Voronoi8/Reti_radius_%s/Rete_%d_Radius_%d.txt',string,geometrie(n)-1,good_networks(position(n),3));
connectivity=load(filename_con,'connectivity','-ascii');
nodes=load(filename_nod,'nodes','-ascii');
radius_data = importdata(filename_raggio);
raggi=[raggi; radius_data.data];

if max(raggi) > 1.2E-2 || min(raggi)<4E-3
    disp ('ERROR')
end

    if n==1
        pos=pos+delta/2;
    else
        pos=pos+delta;
    end

    for i=1:size(connectivity,1)
        connectivity(i,2)=connectivity(i,2)+Nnodi_max;
        connectivity(i,3)=connectivity(i,3)+Nnodi_max;
    end

    nodes3D=nodes;

    for i=1:size(nodes,1)
        nodes(i,1)=nodes(i,1)+Nnodi_max;
        nodes(i,4)=pos;
    end
    
    %Transformation from 2D to 3D in a range (-delta/2,+delta/2)
    nodes3D=nodes;
    r = randi([ceil((pos-1/2*delta)*1000) ceil((pos+1/2*delta)*1000)],size(nodes3D,1),1)/1000;
    nodes3D(:,4)=r;


    connectivityTOT=[connectivityTOT; connectivity];
    nodesTOT=[nodesTOT; nodes];
    nodes3DTOT=[nodes3DTOT; nodes3D];

    Nnodi_max=size(nodes,1)+Nnodi_max;
end

connectivityTOT=connectivityTOT(2:size(connectivityTOT,1),:);
nodesTOT=nodesTOT(2:size(nodesTOT,1),:);
nodes3DTOT=nodes3DTOT(2:size(nodes3DTOT,1),:);

if max(nodes3DTOT(:,2)) > 1 || max(nodes3DTOT(:,3)) > 1 ||max(nodes3DTOT(:,4)) > 1 || min(nodes3DTOT(:,2))<0 ||min(nodes3DTOT(:,3))<0||min(nodes3DTOT(:,4))<0
    disp ('ERROR in length, one or more vertices are outside the domain')
end

for i=1:size(connectivityTOT,1)
    connectivityTOT(i,1)=i;
end

vessels_3Dlength=compute_vessels_length(connectivityTOT, nodes3DTOT);
total_length=sum(vessels_3Dlength(:,2));
R_L3D=vessels_3Dlength(:,2).*raggi(:,1);
Density3D=sum(R_L3D)*2*pi/d;
vessels_2Dlength=compute_vessels_length(connectivityTOT, nodesTOT);
total_length2d=sum(vessels_2Dlength(:,2));
R_L2D=vessels_2Dlength(:,2).*raggi(:,1);
Density2D=sum(R_L2D)*2*pi/d;
 
filename_con='Voronoi_Conn_2D_18.txt';
filename_nod='Voronoi_Nodes_2D_18.txt';
filename_nod3D='Voronoi_Nodes_3D_18.txt';
save(filename_con,'connectivityTOT','-ascii');
save(filename_nod,'nodesTOT','-ascii');
save(filename_nod3D,'nodes3DTOT','-ascii');

figure
printNetwork(nodesTOT,connectivityTOT,20,'Voronoi2D_18.pts',0,1);
figure
printNetwork(nodes3DTOT,connectivityTOT,20,'Voronoi3D_18.pts',0,1);

filenameR='Radius_Voronoi_18_dens.pts';
fileID=fopen(filenameR,'w');
fprintf(fileID,'BEGIN_LIST\n');
for i=1:size(raggi,1)
    fprintf(fileID, '%d\n',raggi(i,1));
end
fprintf(fileID,'END_LIST\n');
fclose(fileID);


