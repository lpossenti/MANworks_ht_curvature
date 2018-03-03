clear all
clc
close all

%% Script for the creation of a defined number of Voronoi networks
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 10/07/2017
%   Contact: simone.digre@gmail.com  
%%
NRETI=10;
n_point=8;
triD_dim=0;
DIR_IN=32;
DIR_OUT=15;
DIR_IN_R=5E-6;
DIR_OUT_R=5E-6;
DIR_IN_R_4=4E-6; 
DIR_OUT_R_4=4E-6;
DIR_IN_R_4m=4.5E-6; %condizione Dirichlet in ingresso
DIR_OUT_R_4m=4.5E-6;
plot=1;

for n_reti=1:NRETI
    
    [connectivity, nodes]=Create_Voronoi_Network(n_point,DIR_IN,DIR_OUT,triD_dim,plot);
    nodes_R=Change_BC_Voronoi(nodes,DIR_IN_R,DIR_OUT_R,DIR_IN,DIR_OUT);
    nodes_R_4=Change_BC_Voronoi(nodes,DIR_IN_R_4,DIR_OUT_R_4,DIR_IN,DIR_OUT);
    nodes_R_4m=Change_BC_Voronoi(nodes,DIR_IN_R_4m,DIR_OUT_R_4m,DIR_IN,DIR_OUT);
    
    filename_con=sprintf('Voronoi_Conn_%d.txt',n_reti-1);
    filename_nodes=sprintf('Voronoi_Nodes_%d.txt',n_reti-1);
    filename_rete = sprintf('Voronoi_%d.pts',n_reti-1);
    filename_rete_R = sprintf('Voronoi_%d_R.pts',n_reti-1);
    printNetwork(nodes,connectivity,4,filename_rete,0,1);
    printNetwork(nodes_R,connectivity,4,filename_rete_R,0,0);
    filename_rete_R4 = sprintf('Voronoi_%d_R4.pts',n_reti-1);
    printNetwork(nodes_R_4,connectivity,4,filename_rete_R4,0,0);
    filename_rete_R4m = sprintf('Voronoi_%d_R4m.pts',n_reti-1);
    printNetwork(nodes_R_4m,connectivity,4,filename_rete_R4m,0,0);
    save(filename_con,'connectivity','-ascii');
    save(filename_nodes,'nodes','-ascii');

end