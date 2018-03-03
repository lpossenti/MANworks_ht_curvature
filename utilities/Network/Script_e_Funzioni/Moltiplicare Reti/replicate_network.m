function [nodes, connectivity]=replicate_network(Single_Nodes,Single_Connectivity, N_points_Y1, N_points_Y2, N_points_Z, SizeY, SizeZ, Inversion, Plot)
%% compute_starting_point_position finds the coordinates of the points from which the pattern networks will start
% 
% 
%               Single_Nodes = Nodes of the pattern (id coordinates...)
%               Single_Connectivity = Connectivity of the pattern
%               N_points_Y1 = Number of points on the odd index planes
%               N_points_Y2 = Number of points on the even index planes
%               N_points_Z = Half number of planes (it can be integer e.g. 5
%               or an integer + 0.5 e.g. 5.5)
%               SizeY = size of domain along y direction
%               SizeZ = size of domain along z direction
%               Inversion = flag to invert the directions of the planes
%               Plot = flag to plot the network
%    
% You will obtain
% (N_points_Y1+N_points_Y2)*mod(N_points_Z,1)+N_points_Y1*(N_points_Z-mod(N_points_Z,1))
% 
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 07/07/2017
%   Contact: simone.digre@gmail.com  

starting_points=compute_starting_point_position(N_points_Y1, N_points_Y2, N_points_Z, SizeY, SizeZ);

nodes=multiply_network(starting_points, Single_Nodes);

if Inversion
 nodes=perpendicular_planes(nodes,N_points_Z,N_points_Y1, N_points_Y2, Single_Nodes);
end

connectivity=mult_connec(Single_Nodes,Single_Connectivity,starting_points);

if N_points_Y1==0 || N_points_Y2==0
    for i=1:size(nodes,1)
        nodes(i,4)=0.5; 
    end
end


if Plot
figure
hold on
    for segment=1:size(connectivity,1)
        startingNode=find_coord_nodes(nodes,connectivity(segment,2));
        endingNode=find_coord_nodes(nodes,connectivity(segment,3));
         plot3([startingNode(1) endingNode(1)],[startingNode(2) endingNode(2)],[startingNode(3) endingNode(3)],'k');
    end
axis([0 1 0 1 0 1])
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
hold off
end

%% OTHER FIGURES
% figure
% plot(starting_points(:,2),starting_points(:,3),'k*')
% axis equal
% axis([0 1 0 1])
% xhandle=xlabel('z');
% yhandle=ylabel('y');
% set(xhandle,'Fontsize',20)
% set(xhandle,'Fontname','Times New Roman')
% set(yhandle,'Fontsize',20)
% set(yhandle,'Fontname','Times New Roman')
% set(gca,'XTick',[0 1],'Fontname','Times New Roman','Fontsize',20)
% set(gca,'YTick',[0 1],'Fontname','Times New Roman','Fontsize',20)
% hold off

% grid on
% figure
% % plot3(zeros(5,1),starting_points(1:5,2),starting_points(1:5,3),'k*',starting_points(6:10,2),zeros(5,1),starting_points(6:10,3),'k*')
% hold on
% plot3([0 0], [0 1], [0 0] ,'k')
% plot3([0 1], [1 1], [0 0] ,'k')
% plot3([1 1], [0 1], [0 0] ,'k')
% plot3([0 1], [0 0], [0 0] ,'k')
% plot3([0 0], [0 0], [0 1] ,'k')
% plot3([0 1], [1 1], [1 1] ,'k')
% plot3([0 0], [0 1], [1 1] ,'k')
% plot3([1 1], [0 1], [1 1] ,'k')
% plot3([0 1], [0 0], [1 1] ,'k')
% plot3([1 1], [1 1], [0 1] ,'k')
% plot3([0 0], [0 1], [1 1] ,'k')
% plot3([1 1], [0 0], [0 1] ,'k')
% plot3([0 0], [1 1], [1 0] ,'k')
% axis([0 1 0 1 0 1])
% axis equal
% xhandle=xlabel('x');
% yhandle=ylabel('y');
% zhandle=zlabel('z');
% set(xhandle,'Fontsize',20)
% set(xhandle,'Fontname','Times New Roman')
% set(yhandle,'Fontsize',20)
% set(yhandle,'Fontname','Times New Roman')
% set(zhandle,'Fontsize',20)
% set(zhandle,'Fontname','Times New Roman')
% set(gca,'XTick',[0 1],'Fontname','Times New Roman','Fontsize',20)
% set(gca,'YTick',[0 1],'Fontname','Times New Roman','Fontsize',20)
% set(gca,'ZTick',[0 1],'Fontname','Times New Roman','Fontsize',20)
% hold off
