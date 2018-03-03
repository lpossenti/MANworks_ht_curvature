clear all
close all
clc

%parameter for plot
pt_plot_min=-5;
pt_plot_max=10;
Pt=[pt_plot_min:0.01:pt_plot_max];

%parameter for building lymphatic drainage
Q_work=8.55E-7; % lymphatic flow in physiological condition
Qmax=1.7E-5; % maximum lymphatic flow
A=Qmax*1.01; %coefficient A=101%*Qmax
Qmean_= (Qmax+Q_work) / 2; % mean lymphatic flow
p_work=-1; % interstitial fluid pressure in physiological condition
p_max=5; % maximum interstitial fluid pressure in physiological condition
p_mean=(p_max+p_work) / 2; % mean interstitial fluid pressure


%parameter for the iterative method
res=1000;
Dold=0.5;
Cold=0.719;
Bold=2E-7;

%iterative process
iter=0;

while res > 0.0001 && iter<10000

B=(A-Q_work)*(1+exp((p_work+Dold)/Cold)); 
C=(p_max+Dold)/log((A-B-Qmax)/(Qmax-A));
D=C*log((Qmean_+B-A)/(A-Qmean_))-p_mean;

res= max([abs((C-Cold)/C),abs((D-Dold)/D),abs((B-Bold)/B)]);
iter=iter+1;

Bold=B;
Cold=C;
Dold=D;

end

%building sigmoid lymphatic function
QLF=A-B./(1+exp((Pt+D)/C));

%building Chamney curve
Q_derivative=(Qmax-Q_work)/(p_max-p_work);
QLF_retta=Q_derivative.*([-1:0.01:5]+1)+Q_work;
dim_max=size([5:0.01:pt_plot_max],2)-1;
dim_min=size([pt_plot_min:0.01:-1],2)-1;
MIN=ones(1,dim_min)*Q_work;
MAX=ones(1,dim_max)*Qmax;
QL_Cha=[MIN,QLF_retta,MAX];

%PLOT COMPARISON WITH Chamney et al 1999
figure
h =plot(Pt,QLF,'r',Pt,QL_Cha,'b');
Lhandle=legend('Sigmoide','[Chamney et al 1999] modificata','Location','NorthWest');
xhandle=xlabel('Pt [mmHg]');
yhandle=ylabel('Portata del Sistema Linfatico [1/s]');
set(xhandle,'Fontsize',12)
set(xhandle,'Fontname','Times New Roman')
set(yhandle,'Fontsize',12)
set(yhandle,'Fontname','Times New Roman')
set(Lhandle,'Fontsize',10)
set(Lhandle,'Fontname','Times New Roman')
set(h, 'Linewidth', 2);

%PLOT only lymphatic system
figure
h =plot(Pt,QLF,'b');
xhandle=xlabel('Pt [mmHg]');
yhandle=ylabel('Portata del Sistema Linfatico [1/s]');
set(xhandle,'Fontsize',12)
set(xhandle,'Fontname','Times New Roman')
set(yhandle,'Fontsize',12)
set(yhandle,'Fontname','Times New Roman')
set(h, 'Linewidth', 2);

%Obtain Lymphedema
A_50=A/2;
B_50=B/2;
QLF_50=A_50-B_50./(1+exp((Pt+D)/C));
A_75=A/4;
B_75=B/4;
QLF_75=A_75-B_75./(1+exp((Pt+D)/C));
QLF_100=zeros(1,size(QLF,2));

%PLOT Lymphedema
figure
h =plot(Pt,QLF,'r',Pt,QLF_50,'b',Pt,QLF_75,'g',Pt,QLF_100,'y');
Lhandle=legend('Sigmoide Fisiologica','Sigmoide Patologica (-50%)','Sigmoide Patologica (-75%)','Sigmoide Patologica (-100%)','Location','NorthWest');
xhandle=xlabel('Pt [mmHg]');
yhandle=ylabel('Portata Linfatico [1/s]');
set(xhandle,'Fontsize',12)
set(xhandle,'Fontname','Times New Roman')
set(yhandle,'Fontsize',12)
set(yhandle,'Fontname','Times New Roman')
set(Lhandle,'Fontsize',10)
set(Lhandle,'Fontname','Times New Roman')
set(h, 'Linewidth', 2);