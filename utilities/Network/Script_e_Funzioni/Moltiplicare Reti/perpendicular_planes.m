function nodes=perpendicular_planes(nodes,N_puntiz,N_puntiy1, N_puntiy2,single_nodes)


for i=1:fix(N_puntiz)
    if ~(mod(i,2))
    temp=nodes((i-1)*(N_puntiy1 + N_puntiy2)*size(single_nodes,1)+1:i*(N_puntiy1 + N_puntiy2)*size(single_nodes,1),2);
    nodes((i-1)*(N_puntiy1 + N_puntiy2)*size(single_nodes,1)+1:i*(N_puntiy1 + N_puntiy2)*size(single_nodes,1),2)=nodes((i-1)*(N_puntiy1 + N_puntiy2)*size(single_nodes,1)+1:i*(N_puntiy1 + N_puntiy2)*size(single_nodes,1),3);
    nodes((i-1)*(N_puntiy1 + N_puntiy2)*size(single_nodes,1)+1:i*(N_puntiy1 + N_puntiy2)*size(single_nodes,1),3)=temp;
    end
end