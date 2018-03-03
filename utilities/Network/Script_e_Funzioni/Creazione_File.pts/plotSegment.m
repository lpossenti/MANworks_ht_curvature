function plotSegment(pointMatrix)
%plotSegment(pointMatrix) plot the segment with getPoints output
%
%   to generate pointMatrix see getPoints.m
%
%   Author: Luca Possenti 
%   Politecnico di Milano, 30/09/2016
%   Contact: luca.possenti@polimi.it   

%%ordering
orderedPoints=zeros(size(pointMatrix,1),4);
for i=1:size(pointMatrix,1)

    switch i
        case 1
            orderedPoints(i,:)=pointMatrix(i,:);
        case 2
            orderedPoints(end,:)=pointMatrix(i,:);
        otherwise 
            orderedPoints(i-1,:)=pointMatrix(i,:);
    end
            
end

%% Plot
figure
plot3(orderedPoints(:,2),orderedPoints(:,3),orderedPoints(:,3),'ro')
xlabel('x');
ylabel('y');
zlabel('z');

end
