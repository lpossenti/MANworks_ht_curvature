clear all
close all
clc

%% Script for the creation of a Voronoi Network (also with different BC and 3D )
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 10/07/2017
%   Contact: simone.digre@gmail.com  
%%
N_points=8;
DIR_IN=32;
DIR_OUT=15;
triD_dim=0;

[connectivity, nodes]=Create_Voronoi_Network(N_points,DIR_IN,DIR_OUT,triD_dim,1);
% 
% filename_con='rete_Conn.txt';
% filename_nod='rete_Nodes.txt';
% figure
% printNetwork(nodes,connectivity,20,'rete.pts',0,1);
% save(filename_con,'connectivity','-ascii');
% save(filename_nod,'nodes','-ascii');

DIR_IN_new=4E-6;
DIR_OUT_new=5E-6;
[new_nodes]=Change_BC_Voronoi(nodes,DIR_IN_new,DIR_OUT_new,DIR_IN,DIR_OUT);

% filename_con='rete_Conn.txt';
% filename_nod='rete_Nodes.txt';
% figure
% printNetwork(new_nodes,connectivity,20,'rete2.pts',0,1);
% save(filename_con,'connectivity','-ascii');
% save(filename_nod,'nodes','-ascii');

delta=0.1;
[nodes3D]=Change_From2DTo3D_Delta(nodes,connectivity,delta,1);

% filename_con='rete3D_Conn.txt';
% filename_nod='rete3D_Nodes.txt';
% figure
% printNetwork(nodes3D,connectivity,20,'rete3D.pts',0,1);
% save(filename_con,'connectivity','-ascii');
% save(filename_nod,'nodes','-ascii');