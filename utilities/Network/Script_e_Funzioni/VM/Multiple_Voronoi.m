clear all
clc
close all

%% Script for the creation of a network made by the superposition of a defined number of Voronoi networks
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 10/07/2017
%   Contact: simone.digre@gmail.com  
%%
N=2;
DOMAIN=1;
position=0;
delta=DOMAIN/(N);
connectivityTOT=zeros(1,3);
nodesTOT=zeros(1,6);
nodes3DTOT=zeros(1,6);
Nnodi_max=0;
n_points=8;
DIR_IN=32;
DIR_OUT=15;
triD_dim=0;
plot=1;

for n_reti=1:N;
    
    clear nodes
    clear connectivity
    
    if n_reti==1
position=position+delta/2;
    else
position=position+delta;
    end

     [connectivity, nodes]=Create_Voronoi_Network(n_points,DIR_IN,DIR_OUT,triD_dim,plot);

     for i=1:size(connectivity,1)
        connectivity(i,2)=connectivity(i,2)+Nnodi_max;
        connectivity(i,3)=connectivity(i,3)+Nnodi_max;
     end
     
    for i=1:size(nodes,1)
    nodes(i,1)=nodes(i,1)+Nnodi_max;
    nodes(i,4)=position;
    end

    nodes3D=nodes;
    r = randi([ceil((position-1/4*delta)*1000) ceil((position+1/4*delta)*1000)],size(nodes3D,1),1)/1000;
    nodes3D(:,4)=r;

    connectivityTOT=[connectivityTOT; connectivity];
    nodesTOT=[nodesTOT; nodes];
    nodes3DTOT=[nodes3DTOT; nodes3D];

    Nnodi_max=size(nodes,1)+Nnodi_max;

end

connectivityTOT=connectivityTOT(2:size(connectivityTOT,1),:);
nodesTOT=nodesTOT(2:size(nodesTOT,1),:);
nodes3DTOT=nodes3DTOT(2:size(nodes3DTOT,1),:);

for i=1:size(connectivityTOT,1)
    connectivityTOT(i,1)=i;
end

filename_con=sprintf('Voronoi_Conn_2D_%d.txt',N);
filename_nod=sprintf('Voronoi_Nodes_2D_%d.txt',N);
filename_nod3D=sprintf('Voronoi_Nodes_3D_%d.txt',N);
filename_network3D=sprintf('Voronoi3D_%d.pts',N);
filename_network2D=sprintf('Voronoi2D_%d.pts',N);
figure
printNetwork(nodesTOT,connectivityTOT,20,filename_network2D,0,1);
figure
printNetwork(nodes3DTOT,connectivityTOT,20,filename_network3D,0,1);
save(filename_con,'connectivityTOT','-ascii');
save(filename_nod,'nodesTOT','-ascii');
save(filename_nod3D,'nodes3DTOT','-ascii');
