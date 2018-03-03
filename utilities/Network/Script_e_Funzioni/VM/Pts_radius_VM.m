clear all
clc
close all

Valore_BC_Radius=4E-6;
Max_N_Voronoi=100;

for i=0:(Max_N_Voronoi-1)
filename_con=sprintf('Voronoi_Conn_%d.txt',i);
filename_nod=sprintf('Voronoi_Nodes_%d.txt',i);

connectivity=load(filename_con,'connectivity','-ascii');
nodes=load(filename_nod,'nodes','-ascii');

for j=1:size(nodes,1)
      if nodes(j,5)==0
          nodes(j,6)=Valore_BC;
      end
end

figure
filename_file=sprintf('Voronoi_%d_R4.pts',i);
printNetwork(nodes,connectivity,20,filename_file,0,1);

end