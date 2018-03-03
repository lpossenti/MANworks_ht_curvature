function global_connect=mult_connec(single_nodes,single_connectivity,starting_points)

global_connect=zeros(size(single_connectivity,1)*size(starting_points,1),size(single_connectivity,2));

           for i=1:size(starting_points,1)
            global_connect((1+(i-1)*size(single_connectivity,1)):i*size(single_connectivity,1),:)=single_connectivity;
           for j=1:size(single_connectivity,1)
           global_connect((1+(i-1)*size(single_connectivity,1))+j-1,2)=global_connect((1+(i-1)*size(single_connectivity,1))+j-1,2)+(i-1)*size(single_nodes,1);
global_connect((1+(i-1)*size(single_connectivity,1))+j-1,3)=global_connect((1+(i-1)*size(single_connectivity,1))+j-1,3)+(i-1)*size(single_nodes,1);
          end
end

for i=1:size(global_connect,1)
    global_connect(i,1)=i;
end