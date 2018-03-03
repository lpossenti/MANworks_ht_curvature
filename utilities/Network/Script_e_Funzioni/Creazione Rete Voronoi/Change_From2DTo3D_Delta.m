function [nodes3D]=Change_From2DTo3D_Delta(nodes2D,connectivity,delta,plot)
%% Change_From2DTo3D_Delta make a transformation from a 2D geometry to a 3D one setting the third coordinate randomly in a range (-delta,+ delta)
%
%        nodes2D = nodes of the 2D geometry
%        connectivity = connectivity of geometry (mandatory if plot=1)
%        delta = displacement of third coordinate
%        plot = flag to plot the new 3D geometry
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 08/05/2017
%   Contact: simone.digre@gmail.com  
%%
nodes3D=nodes2D;

for i=1:size(nodes2D,1)
    position=nodes2D(i,4);
    nodes3D(i,4)= randi([ceil((position-delta)*1000) ceil((position+delta)*1000)],1,1)/1000;
end

if plot
figure
hold on
    for segment=1:size(connectivity,1)
        startingNode=find_coord_nodes(nodes3D,connectivity(segment,2));
        endingNode=find_coord_nodes(nodes3D,connectivity(segment,3));
         plot3([startingNode(1) endingNode(1)],[startingNode(2) endingNode(2)],[startingNode(3) endingNode(3)],'k');
    end
axis([0 1 0 1 0 1])
xlabel('x');
ylabel('y');
zlabel('z');
hold off
end