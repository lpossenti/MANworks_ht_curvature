function [new_nodes]=Change_BC_Voronoi(nodes,DIR_IN_new,DIR_OUT_new,DIR_IN_old,DIR_OUT_old)
%% Change_BC_Voronoi changes the BC of a generic geometry
%
%        nodes = nodes of the geometry
%        DIR_..._old = old BC
%        DIR_..._new = new BC
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 10/07/2017
%   Contact: simone.digre@gmail.com  
%%
new_nodes=nodes;
for j=1:size(nodes,1)
      if nodes(j,5)==0
            if nodes(j,6)== DIR_IN_old
                new_nodes(j,6)=DIR_IN_new;
            elseif nodes(j,6)== DIR_OUT_old
                new_nodes(j,6)=DIR_OUT_new;
            else
                disp('error in changing BC')
            end
      end
end