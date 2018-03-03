clear all
clc
close all
%% Script for the creation of a Voronoi Network STARTING FROM NORMAL DISTRIBUITED POINTS
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 10/07/2017
%   Contact: simone.digre@gmail.com  
%%
numero_punti=250;
x_notNorm = randn(numero_punti,1); %scelgo delle coordinate x casuali
y_notNorm = randn(length(x_notNorm), 1); %scelgo delle coordinate x casuali

x_pos = x_notNorm - (min(x_notNorm)* (min(x_notNorm)<0));
y_pos = y_notNorm - (min(y_notNorm)* (min(y_notNorm)<0));
x = x_pos./max(x_pos);
y = y_pos./max(y_pos);


%Get Voronoi
[Vx,Vy] = voronoi(x,y); %creo il diagramma di voronoi 2D
plot(x,y,'k.',Vx,Vy,'k-');
save('MM_Coordinates.txt','Vx','Vy','-ascii')
axis equal;
figure
plot(x,y,'k.',Vx,Vy,'k-');
hold on
plot([0,0],[0,1],'r',[0,1],[1,1],'r',[1,1],[1,0],'r',[1,0],[0,0],'r')
axis equal;
% axis([0,1,0,1]);
[vx,vy]=taglia_segmenti(Vx,Vy); %circoscrivo il diagramma di Voronoi al dominio [0,1] x [0,1]
figure
plot(x,y,'r+',vx,vy,'b-');
plot(x,y,'k.',vx,vy,'k-');
axis equal;
%trovo tutti i vertici del diagramma evitando ripetizioni e ponendo le
%coordinate in una matrice Nx3 con N numero dei vertici e le colonne sono
%ID x y
full_zero=ones(size(vx,2),3)*5;
s=1;
ok=0;
for i=1:size(vx,2)
    for j=1:size(vx,2)
        if full_zero(j,2)==vx(1,i) && full_zero(j,3)==vy(1,i)
            ok=1;
        end
    end
        if ok==0
            full_zero(s,2)=vx(1,i);
            full_zero(s,3)=vy(1,i);
            full_zero(s,1)=s;
            s=s+1;
        end
        
        ok=0;
end
for i=1:size(vx,2)
    for j=1:size(vx,2)
        if full_zero(j,2)==vx(2,i) && full_zero(j,3)==vy(2,i)
            ok=1;
        end
    end
        if ok==0
            full_zero(s,2)=vx(2,i);
            full_zero(s,3)=vy(2,i);
            full_zero(s,1)=s;
            s=s+1;
        end
        
        ok=0;
end

full=full_zero(1:s-1,1:3);
connectivity=zeros(3,size(full,2));
s=1;
%Creo la matrice di connettività prendendo per ogni segmento gli id dei
%vertici
for i=1:size(vx,2)
    for j=1:size(full,1)
        if vx(1,i)==full(j,2) && vy(1,i)==full(j,3)
            connectivity(s,2)=full(j,1);
                    for k=1:size(full,1)
                      if vx(2,i)==full(k,2) && vy(2,i)==full(k,3)
                       connectivity(s,3)=full(k,1);
                      connectivity(s,1)=s;
                         s=s+1;
                      end
                end
    end
end
end    
% per la funzione Network i nodi devono essere posti con ID x y z BClabel
% BCvalue

nodes=zeros(size(full,1),6);
nodes(:,1:3)=full;
nodes(:,4)=0.5; %pongo z=costante=0.5
DIR_IN=30; %condizione Dirichlet in ingresso
DIR_OUT=10; %condizione Dirichlet in uscita
%i nodi che si trovano in x=0 o y=0 sono considerati vertici arteriosi,
%quindi DIR_IN, mentre i vertici con x=1 o y=1 sono considerati vertici
%venosi quindi DIR_OUT; gli altri sono nodi interni.
a=1;
for i=1:size(nodes,1)
    if nodes(i,2)==0 || nodes(i,3)==0
        nodes(i,5)=0;
        nodes(i,6)=DIR_IN;
        array_in(a)=i;
        a=a+1;
    elseif nodes(i,2)==1 || nodes(i,3)==1
        nodes(i,5)=0;
        nodes(i,6)=DIR_OUT;
    else
        nodes(i,5)=2;
        nodes(i,6)=0;
    end
end

s=1;
for a=1:size(array_in,2)
 [riga, colonna]=find(connectivity(:,2:3)==array_in(a));
tmp=connectivity(s,:);
connectivity(s,:)=connectivity(riga,:);
connectivity(riga,:)=tmp;
if (colonna==2)
    tmp=connectivity(s,2);
    connectivity(s,2)=connectivity(s,3);
    connectivity(s,3)=tmp;
end
s=s+1;
end

for i=1:size(connectivity,1)
    connectivity(i,1)=i;
end

filename_con='Voronoi_Conn.txt';
filename_nod='Voronoi_Nodes.txt';
figure
printNetwork(nodes,connectivity,20,'Voronoi2D.pts',0,1);
save(filename_con,'connectivity','-ascii');
save(filename_nod,'nodes','-ascii');



% set(gcf, 'Position', get(0,'Screensize')); % Enlarge figure to full
