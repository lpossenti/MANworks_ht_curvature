function [new_vx, new_vy]=taglia_segmenti(vx,vy)
%% taglia_segmenti rescale the segment whose vertices are outside the domain [0.1] x [0.1] x [0.1]
%
%        [vx,vy]=outputs of voronoi.m
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 07/10/2016
%   Contact: simone.digre@gmail.com  
%%
%la funzione taglia i segmenti del diagramma di Voronoi che partono all'interno del dominio e poi escono dal
%dominio [0.1] x [0.1] x [0.1]
s=1;

%tutti i segmenti che non iniziano e non finiscono all'interno del dominio
%vengono eliminati

for i=1:size(vx,2)
    if ~((vx(1,i)<0 && vx(2,i)<0) || (vx(1,i)>1 && vx(2,i)>1) || (vy(1,i)<0 && vy(2,i)<0) || (vy(1,i)>1 && vy(2,i)>1))
        Vx(:,s)=vx(:,i);
        Vy(:,s)=vy(:,i);
        s=s+1;
    end
end

vx=Vx;
vy=Vy;
%
%il dominio 2D è un quadrato di lati xA, xB, yA, yB che corrispondono a xA:
%x=0, xB: x=1, yA: y=0, yB: y=1;
%
%           _____yB____
%          |           |
%       xA |           | xB
%          |           |
%          |___________| 
%               yA
xa=1;
xb=1;

for i=1:size(vx,2)
    if vx(1,i)<0 || vx(2,i)<0
        nodi_xA(xa)=i; %salvo le posizioni dei punti che si trovano a sinistra di xA (x<0)
        xa=xa+1;
    end
    if vx(1,i)>1 || vx(2,i)>1 %salvo le posizioni dei punti che si trovano a destra di xB (x>1)
        nodi_xB(xb)=i;
        xb=xb+1;
    end
end

for j=1:xa-1
new_y=nella_retta(vx(1,nodi_xA(j)), vx(2,nodi_xA(j)), vy(1,nodi_xA(j)), vy(2,nodi_xA(j)),1,0); %trovo le coordinate del punto(0,y) del punto appartente alla retta che collega un punto interno al dominio con uno esterno e che passa per x=0;
if vx(1,nodi_xA(j))<0
    position=1;
else
    position=2;
end
vx(position,nodi_xA(j))=0;
vy(position,nodi_xA(j))=new_y;
end

for j=1:xb-1
new_y=nella_retta(vx(1,nodi_xB(j)), vx(2,nodi_xB(j)), vy(1,nodi_xB(j)), vy(2,nodi_xB(j)),1,1); %trovo le coordinate del punto(1,y) del punto appartente alla retta che collega un punto interno al dominio con uno esterno e che passa per x=1;
if vx(1,nodi_xB(j))>1
    position=1;
else
    position=2;
end
vx(position,nodi_xB(j))=1;
vy(position,nodi_xB(j))=new_y;
end

ya=1;
yb=1;

for i=1:size(vy,2)
    if vy(1,i)<0 || vy(2,i)<0 %salvo le posizioni dei punti che si trovano sotto a yA (y<0)
        nodi_yA(ya)=i;
        ya=ya+1;
    end
    if vy(1,i)>1 || vy(2,i)>1 %salvo le posizioni dei punti che si trovano sopra a yB (y>1)
        nodi_yB(yb)=i;
        yb=yb+1;
    end
end

for j=1:ya-1
new_x=nella_retta(vx(1,nodi_yA(j)), vx(2,nodi_yA(j)), vy(1,nodi_yA(j)), vy(2,nodi_yA(j)),0,0); %trovo le coordinate del punto(x,0) del punto appartente alla retta che collega un punto interno al dominio con uno esterno e che passa per y=0;
if vy(1,nodi_yA(j))<0
    position=1;
else
    position=2;
end
vx(position,nodi_yA(j))=new_x;
vy(position,nodi_yA(j))=0;
end

for j=1:yb-1
new_x=nella_retta(vx(1,nodi_yB(j)), vx(2,nodi_yB(j)), vy(1,nodi_yB(j)), vy(2,nodi_yB(j)),0,1); %trovo le coordinate del punto(x,1) del punto appartente alla retta che collega un punto interno al dominio con uno esterno e che passa per y=1;

if vy(1,nodi_yB(j))>1
    position=1;
else
    position=2;
end
vx(position,nodi_yB(j))=new_x;
vy(position,nodi_yB(j))=1;
end

new_vx=vx;
new_vy=vy;
    