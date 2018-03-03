function [valore]=nella_retta(x1, x2, y1, y2,xy,coordinata)
%% nella_retta find the x (or y) value of a point whose only y (or x) coordinate is known 
%
%        (x1,y1), (x2,y2) = points belonging to the line of equation (x-x1)/(x2-x1)=(y-y1)/(y2-y1)
%                      xy = flag to define which one coordinate is unknown
%              coordinata = known coordinate of the point
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 07/10/2016
%   Contact: simone.digre@gmail.com  
%%

%se xy=0 cerco la coordinata x e coordinata è la coordinata y del punto che
%mi interessa, se xy=1 cerco la coordinata y e coordinata è la coordinata x
%del punto che mi interessa.
%Questa funzione prende in ingresso x1,x2, y1,y2 e costruisce la retta di
%equazione (x-x1)/(x2-x1)=(y-y1)/(y2-y1)

if xy==0;
    
    valore=(coordinata-y1)/(y2-y1)*(x2-x1)+x1;
elseif xy==1;
    
    valore=(coordinata-x1)/(x2-x1)*(y2-y1)+y1;
else
    error('error xy must be equal to either 0 or 1')
    
end