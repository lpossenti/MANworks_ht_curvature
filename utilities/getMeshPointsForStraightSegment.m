function [points,l] = getMeshPointsForStraightSegment(x0,xEnd,y0,yEnd, Npoints, ID)
% Author: Luca Possenti <luca.possenti@polimi.it>
% function to calculate uniform mesh of a curve vessel knowing initial and 
% ending point, curvature radius and # of points
% Initial point = (x0,y0)
% Ending point = (xEnd, yEnd)
% N points is the number of nodes desired
% ID is the id of the segment (to be printed as first number of each row)

% element projection
dx = (xEnd - x0)/(Npoints-1);
dy = (yEnd - y0)/(Npoints-1);

% parametrization
points = zeros(Npoints,4);
for n=1:Npoints
        points(n,1) = ID; % x
        points(n,2) = x0 + (n-1)*dx; % x
        points(n,3) = y0 + (n-1)*dy; % y
        points(n,4) = 0.5; %z in the middle of the domain
end

%lenght
l = ((xEnd - x0)^2 + (yEnd -y0)^2)^0.5;

end