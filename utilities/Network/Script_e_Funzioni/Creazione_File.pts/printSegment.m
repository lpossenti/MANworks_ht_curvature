function printSegment(pointMatix,filename,open,BC1,BC2,fileID_opened)
%printSegment(pointMatrix,filename,open) print the segment with getPoints
%output in a Gambit-like format
%
%   to generate pointMatrix see getPoints.m
%   filename is the string with the file name
%   open is a flag 0 = the file is already opened, 1=open the file
%   BC1 and BC2 are stings with the BC labels and values
%
%   Author: Luca Possenti 
%   Politecnico di Milano, 30/09/2016
%   Contact: luca.possenti@polimi.it   

%% input check
switch nargin
    case 1
        disp('Warning: filename automatically set to "segment.pts".');
        filename='segment.pts';
        disp('Warning: open flag automatically set to 1. The file will be opend');
        open=1;
        disp('Warning: BC1 automatically set to "DIR 1.0"');
        BC1='DIR 1.0';
        disp('Warning: BC2 automatically set to "DIR 0.0"');
        BC2='DIR 0.0';
    case 2
        disp('Warning: open flag automatically set to 1. The file will be opend');
        open=1;
        disp('Warning: BC1 automatically set to "DIR 1.0"');
        BC1='DIR 1.0';
        disp('Warning: BC2 automatically set to "DIR 0.0"');
        BC2='DIR 0.0';
    case 3
        disp('Warning: BC1 automatically set to "DIR 1.0"');
        BC1='DIR 1.0';
        disp('Warning: BC2 automatically set to "DIR 0.0"');
        BC2='DIR 0.0';
    case 4
        disp('Warning: BC2 automatically set to "DIR 0.0"');
        BC2='DIR 0.0';
end
if ~open && nargin<6
    disp('Error: please insert the fileId or open a new file.');
    return
end
%% opening
if open
    fileID = fopen(filename,'w');
else 
    fileID =fileID_opened;
end
if fileID == -1
    disp('Warning: filename automatically set to "segment.pts".');
        filename='segment.pts';
        fileID = fopen(filename,'w');
end

%% printing
fprintf(fileID,'BEGIN_ARC\n');
fprintf(fileID,'BC %s\n',BC1);
fprintf(fileID,'BC %s\n',BC2);
for i=1:size(pointMatix,1)
    switch i
        case 1
            fprintf(fileID,'%d %d %d %d start\n',pointMatix(i,1),pointMatix(i,2),pointMatix(i,3),pointMatix(i,4));
        case 2
            fprintf(fileID,'%d %d %d %d end\n',pointMatix(i,1),pointMatix(i,2),pointMatix(i,3),pointMatix(i,4));
        otherwise
            fprintf(fileID,'%d %d %d %d point\n',pointMatix(i,1),pointMatix(i,2),pointMatix(i,3),pointMatix(i,4));    
    end
end

fprintf(fileID,'END_ARC\n');
%% closing
if open
    fclose(fileID);
end
