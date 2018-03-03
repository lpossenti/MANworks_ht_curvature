function [pointMatrix]=getPoints(ptA,ptB,N,ID)
%getPoints(ptA,ptB,N,ID) get the points of the segment for Gambit-like file
%for FEM
%
%   ptA and ptB are start and end point of the segment respectively. Please
%   note that the problem is 3D.
%   N is the number of point desired (start and end included). please note
%   that N>=3.
%   ID is the segment ID
%   pointMatrix is [IDi xi yi zi]. row 1: start || row 2: end || other rows:
%   points
%
%   Author: Luca Possenti 
%   Politecnico di Milano, 30/09/2016
%   Contact: luca.possenti@polimi.it             

%% input check
switch nargin
    case 1
        disp('Warning: ptB automatically set to (0,0,0).');
        ptB=[0 0 0];
        disp('Warning: N automatically set to 3.');
        N=3;
        disp('Warning: ID automatically set to 1.');
        ID=1;
    case 2
        disp('Warning: N automatically set to 3.');
        N=3;
        disp('Warning: ID automatically set to 1.');
        ID=1;
    case 3
        disp('Warning: ID automatically set to 1.');
        ID=1;
end
if (length(ptA)~=length(ptB)) || length(ptA)~=3 
    disp('Error: ptA or ptB dimension are not correct. Please check in the command window.');
    pointMatrix=zeros(2);
    return;
end
if N<3
    disp('Error: N must be higher or equal to 3.');
    pointMatrix=zeros(2);
    return;
end
if ID<=0
    disp('Warning: negative ID set. ID automatically set to 1.');
    ID=1;
end

%% vector calcuation
n=ptB-ptA;

%%Matrix computation
for i=1:N+1
    switch i
        case 1
            pointMatrix(1,:)=[ID ptA(1) ptA(2) ptA(3)];
        case 2
            pointMatrix(2,:)=[ID ptB(1) ptB(2) ptB(3)];
        otherwise 
            pointMatrix(i,:)=[ID ptA(1)+(i-2)/N*n(1) ptA(2)+(i-2)/N*n(2) ptA(3)+(i-2)/N*n(3)];
    end
end

end

