% clear all
% clc
% close all

filename_con='.\Geometrie_Simulazioni\Geometrie_7000m-1\Esagoni_46_dens_Conn.txt';
filename_nod='.\Geometrie_Simulazioni\Geometrie_7000m-1\Esagoni_46_dens_Nodes.txt';
connectivity=load(filename_con,'connectivity','-ascii');
nodes=load(filename_nod,'nodes','-ascii');

%PLOT THE NETWORK
figure
hold on
    for segment=1:size(connectivity,1)
        startingNode=find_coord_nodes(nodes,connectivity(segment,2));
        endingNode=find_coord_nodes(nodes,connectivity(segment,3));
         plot3([startingNode(1) endingNode(1)],[startingNode(2) endingNode(2)],[startingNode(3) endingNode(3)],'k');
    end
axis([0 1 0 1 0 1])
xlabel('x');
ylabel('y');
zlabel('z');
hold off

%CREATE FILE .pts
figure
printNetwork(nodes,connectivity,20,'rete.pts',0,1);
