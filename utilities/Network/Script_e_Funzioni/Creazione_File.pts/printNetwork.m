
function printNetwork(nodes,connettivity,N,filename,plot_dots,plot_segm)
%network(nodes,connettivity,N,filename,plot) print the network in a
%Gambit-like format.
%
%   nodes is the matrix of the nodes: [ID x y z BCtype BCvalue]
%   BCtype= 0 -> DIR =1->MIX 2-> INT
%   connettivity is [#segment startingNode endingNode]
%   N number of discretization per each segment
%   filename is the string with the file name
%   plot is a flag 1=plot the network
%
%   Author: Luca Possenti 
%   Politecnico di Milano, 30/09/2016
%   Contact: luca.possenti@polimi.it   
%% opening
fileID=fopen(filename,'w');
fprintf(fileID,'BEGIN_LIST\n');

if nargin <4
    plot_dots = 0;
    plot_segm = 0;
end
if nargin <5
    plot_segm = 0;
end
        if plot_segm
        figure
        end
%% network computation
for segment=1:size(connettivity,1)
    % retireve nodes
    startingNode=findNode(nodes,connettivity(segment,2));
    endingNode=findNode(nodes,connettivity(segment,3));
    %get points
    segmentPoints=getPoints(startingNode(2:4),endingNode(2:4),N,connettivity(segment,1));
    if plot_dots
        plotSegment(segmentPoints);
    end
    switch startingNode(5) 
        case 0 
            startBC=sprintf('DIR %d',startingNode(end));
        case 1
            startBC=sprintf('MIX %d',startingNode(end));
        case 2
            startBC=sprintf('INT');
        otherwise
            startBC=sprintf('INT');
            disp('Error: BC label in network.m');  
    end
    switch endingNode(5)
        case 0
            endBC=sprintf('DIR %d',endingNode(end));
        case 1
            endBC=sprintf('MIX %d',endingNode(end));
        case 2
            endBC=sprintf('INT');
        otherwise
            endBC=sprintf('INT');
            disp('Error: BC label in network.m');      
    end
    printSegment(segmentPoints,filename,0,startBC,endBC,fileID);
end

%% closing
fprintf(fileID,'END_LIST\n');
fclose(fileID);

%% plotNetwork
hold on
if plot_segm
    for segment=1:size(connettivity,1)
        startingNode=findNode(nodes,connettivity(segment,2));
        endingNode=findNode(nodes,connettivity(segment,3));
        plot3([startingNode(2) endingNode(2)],[startingNode(3) endingNode(3)],[startingNode(4) endingNode(4)],'r-');
    end
end
xlabel('x');
ylabel('y');
zlabel('z');
hold off

end

%%auxialiary function
function [node]=findNode(nodes,ID)
    
    for i=1:size(nodes,1)
        if ID == nodes(i,1)
            node = nodes(i,:);
            return;
        end
    end
    
end