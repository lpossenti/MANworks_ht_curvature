function coordinates=find_coord_nodes(nodes,id_node)
%% find_coord_nodes find the coordinates of a node with known id
%
%        nodes=coordinates of nodes of network vessels
%        id_node=id of the node
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 07/07/2017
%   Contact: simone.digre@gmail.com  
%%
j=1;
found=0;
   while found==0 && j~=size(nodes,1);
   
        if nodes(j,1)==id_node
            found=1;
        else
            j=j+1;
        end %end if
       
   end %end while
   
   coordinates=nodes(j,2:4);