%script for repetion of a pattern inside the domain

clear all
close all
clc
%CHOOSE THE NETWORK TO REPEAT (the y and z coordinates of the center
%point of the pattern must be equal to 1, while the other must be relative
%to the point (x,1,1)

%SINGLEBRANCH
singola_rete=[  1	0	1   1   0	32;
                2   1   1   1   0   15];
singola_connettivita=[1 1 2];

%ESAGONO DIMO
% % % singola_rete=[  1   0   1   1   0   32;
% % %                 2   0.3 1   1   2   0;
% % %                 3   0.35 0.95 1   2   0;
% % %                 4   0.65 0.95 1   2   0;
% % %                 5   0.35 1.05 1   2   0;
% % %                 6   0.65 1.05 1   2   0;
% % %                 7   0.7 1   1   2   0;
% % %                 8   1   1   1   0   15];
% % % singola_connettivita=[1 1   2;
% % %                       2 2   3;
% % %                       3 2   5;
% % %                       4 3   4;
% % %                       5 5   6;
% % %                       6 4   7;
% % %                       7 6   7;
% % %                       8 7   8];

%ESAGONO DIMO_ SENS
% singola_rete=[  1   0   0.975   0.99   0   32;
%                 2   0.3 0.975   0.99   2   0;
%                 3   0.35 1.025 0.99   2   0;
%                 4   0.65 1.025 0.99   2   0;
%                 5   0.35 0.925 0.99   2   0;
%                 6   0.65 0.925 0.99   2   0;
%                 7   0.7 0.975   0.99   2   0;
%                 8   1   0.975   0.99   0   15;
%                 9   0   1.025   1.01   0   32;
%                 10   0.3 1.025   1.01   2   0;
%                 11   0.35 1.075 1.01   2   0;
%                 12   0.65 1.075 1.01   2   0;
%                 13   0.35 0.975 1.01   2   0;
%                 14   0.65 0.975 1.01   2   0;
%                 15   0.7 1.025   1.01   2   0;
%                 16   1   1.025   1.01   0   15];
%             
% singola_connettivita=[1 1   2;
%               2 2   3;
%               3 3   4;
%               4 4   7;
%               5 2   5;
%               6 5  6;
%               7 6   7;
%               8 7   8;
%               9    9   10;
%               10    10  11;
%               11    11  12;
%               12    12  15;
%               13    10  13;
%               14    13  14;
%               15    14  15;
%               16    15  16;];
                  
% % % % % % % % % MESH
% % % % %  singola_rete=[1	0       1     1	0	32;
% % % % %         2   0.35    1     1	2	0;
% % % % %         3   0.45    0.9     1 2   0;
% % % % %         4   0.5     0.95    1 2   0;
% % % % %         5   0.55    0.9     1 2   0;
% % % % %         6   0.5     0.85    1 2   0;
% % % % %         7   0.65    1     1 2   0;
% % % % %         8   0.45    1.1     1 2   0;
% % % % %         9   0.50    1.15    1 2   0;
% % % % %         10  0.55    1.10    1 2   0;
% % % % %         11  0.50    1.05    1 2   0;
% % % % %         12  1.00    1     1 0   15]
% % % % % 
% % % % % singola_connettivita=[1 1   2;
% % % % %               2 2   8;
% % % % %               3 2   3;
% % % % %               4 3   4;
% % % % %               5 3   6;
% % % % %               6 4   5;
% % % % %               7 6   5;
% % % % %               8 5   7;
% % % % %               9 8   9;
% % % % %               10    8   11;
% % % % %               11    9   10;
% % % % %               12    11  10;
% % % % %               13    10  7;
% % % % %               14    7   12];
% %           
% %singola_rete=[1	0       1     1	0	32;
% %         2   0.3    1     1	2	0;
% %         3   0.4    1.1     1 2   0;
% %         4   0.5     1.1    1 2   0;
% %         5   0.6    1.1     1 2   0;
% %         6   0.4     0.9    1 2   0;
% %         7   0.5    0.9     1 2   0;
% %         8   0.6    0.9     1 2   0;
% %         9   0.7    1    1 2   0;
% %         10  1    1    1 0  15];
% % 
% % singola_connettivita=[1 1   2;
% %               2 2   3;
% %               3 3   4;
% %               4 4   5;
% %               5 5   9;
% %               6 2   6;
% %               7 6   7;
% %               8 7   8;
% %               9 8   9;
% %               10    9   10;];

N_points_Y1=3;
N_points_Y2=5;
N_points_Z=2;
DimensioneY=1;
DimensioneZ=1;

[nodes, connectivity]=replicate_network(singola_rete,singola_connettivita, N_points_Y1, N_points_Y2, N_points_Z, DimensioneY, DimensioneZ, 1, 1);

% IF YOU WANT FILE .pts
% figure
% printNetwork(nodes,connectivity,20,'rete.pts',0,1);
% filename_con='rete_conn.txt';
% filename_nod='rete_nod.txt';
% save(filename_con,'connectivity','-ascii');
% save(filename_nod,'nodes','-ascii');
