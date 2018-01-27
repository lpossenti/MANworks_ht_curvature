function [points,l] = getMeshPointsForSegment(x0,xEnd,y0,yEnd,curvatureRadius, curvatureType, Npoints, ID, R)
% Author: Luca Possenti <luca.possenti@polimi.it>
% function to calculate uniform mesh of a curve vessel knowing initial and 
% ending point, curvature radius and # of points
% Initial point = (x0,y0)
% Ending point = (xEnd, yEnd)
% Curvature radius = 1/curvature
% Curvature type 1 -> curved vessel above the segment different from 1
% (usually -1) for curved vessel below the segment
% N points is the number of nodes desired
% ID is the id of the segment (to be printed as first number of each row)
% Radius of the vessel (R) is need to check kR parameter


%input checking
if ~curvatureRadius
    disp('Invalid null radius')
    return
end
curvature = 1/curvatureRadius;
disp('curvature ')
disp(curvature)
disp('kR ')
disp(curvature * R)

%cord lenght
segment = ((xEnd - x0)^2 + ( yEnd - y0)^2)^0.5;
% center angle for curved vessel
beta = 2 * asin(segment/2/curvatureRadius);

% mean point of XY
Mx = (x0 + xEnd)/2;
My = (y0 + yEnd)/2;

% MC lenght
MC = curvatureRadius * cos(beta/2);

% Obtaing versor perpendicular to cord
if(curvatureType==1)
    Vy = (xEnd - x0)/(yEnd - y0);
    Vx = -1;
else
    Vy = -(xEnd - x0)/(yEnd - y0);
    Vx = 1;
end
%Vx = 1;
mod = (Vx^2 + Vy^2)^0.5;
Vx = Vx/mod;
Vy = Vy/mod;

% center coordinate
Cx = Mx + Vx*MC;
Cy = My + Vy*MC;

% teta0
if(curvatureType==1)
    teta0 = acos((xEnd - Cx)/curvatureRadius);
else
    teta0 = -acos((x0 - Cx)/curvatureRadius);
end


%parametrization
points = zeros(Npoints,4);
if(curvatureType==1)
    for n=1:Npoints
        teta = teta0 + beta - beta/(Npoints - 1) * (n - 1);
        points(n,1) = ID; % x
        points(n,2) = Cx + curvatureRadius*cos(teta); % x
        points(n,3) = Cy + curvatureRadius*sin(teta); % y
        points(n,4) = 0.5; %z in the middle of the domain
    end
else
    for n=1:Npoints
        teta = beta/(Npoints - 1) * (n - 1) + teta0;
        points(n,1) = ID; % x
        points(n,2) = Cx + curvatureRadius*cos(teta); % x
        points(n,3) = Cy + curvatureRadius*sin(teta); % y
        points(n,4) = 0.5; %z in the middle of the domain
    end
end

% lenght
l = beta * curvatureRadius;
disp('Lenght: ')
disp(l)

%plotting
plot([x0 xEnd],[y0 yEnd],'--')
hold on
plot(Mx,My,'o')
plot([Mx Mx+Vx],[My My+Vy],'->')
plot(Cx,Cy,'o')
plot(points(:,2),points(:,3),'o')
hold off
axis square
axis equal
end