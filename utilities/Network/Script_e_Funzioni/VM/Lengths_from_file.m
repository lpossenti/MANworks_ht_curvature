clc
close all
clear all

filename_con='.Raggi_medi_Voronoi8\Voronoi_Conn_85.txt';
filename_nod='.Raggi_medi_Voronoi8\Voronoi_Nodes_85.txt';
connectivity=load(filename_con,'connectivity','-ascii');
nodes=load(filename_nod,'nodes','-ascii');
% nodes(:,4)=0.5;

lung=calcola_lunghezze(connectivity, nodes);
total_length=sum(lung(:,2));
max_length=max(lung(:,2));
min_length=min(lung(:,2));
mean_length=mean(lung(:,2));
sd_length=std(lung(:,2));
median_length = median(lung(:,2));
first_qua_length = quantile((lung(:,2)),0.25);
third_qua_length = quantile((lung(:,2)),0.75);
