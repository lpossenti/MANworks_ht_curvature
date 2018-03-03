function multi=multiply_network(start_points, single_network)
%% multiply_network takes a network and repeat it in the domain starting from pre-established points
%
%        start_points=starting points of the repeated networks
%        single_network=pattern of the future network
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 09/02/2017
%   Contact: simone.digre@gmail.com  
%%
multi=zeros( size(start_points,1)*size(single_network,1),6);
j=0;
for i=1:size(start_points,1)
            multi((1+(i-1)*size(single_network,1)):i*size(single_network,1),:)=single_network;
            for j=1:size(single_network,1)
            multi((1+(i-1)*size(single_network,1))+j-1,3)=multi((1+(i-1)*size(single_network,1))+j-1,3)+(start_points(i,2)-1);
            multi((1+(i-1)*size(single_network,1))+j-1,4)=multi((1+(i-1)*size(single_network,1))+j-1,4)+(start_points(i,3)-1);
            end
end

for i=1:size(multi,1)
    multi(i,1)=i;
end
