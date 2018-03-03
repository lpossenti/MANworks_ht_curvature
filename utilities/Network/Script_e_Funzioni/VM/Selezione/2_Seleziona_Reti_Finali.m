clear all
close all
clc
%% Script for the selection of the networks which has values of average radius in a range of variability
%
%   Author: Simone Di Gregorio
%   Politecnico di Milano, 10/07/2017
%   Contact: simone.digre@gmail.com  
%%
radius_4 = '4';
radius_425 = '425';
radius_4m = '4m';
start_radius=[4 4.25 4.5];

percentage=5; % percentage of variation from the wanted average radius
wanted_average_radius=4; % µm
raggio_medio_Sup=wanted_average_radius*(100+percentage)/100;
raggio_medio_Inf=wanted_average_radius*(100-percentage)/100;

k=1;

for i=1:size(start_radius,2)
    if i==1
    fileName = sprintf('Reti_Buone_%s.txt',radius_4);
    elseif i==2
    fileName = sprintf('Reti_Buone_%s.txt',radius_425);
    elseif i==3
    fileName = sprintf('Reti_Buone_%s.txt',radius_4m);
    end
    
    good_networks_struct = importdata(fileName);
    good_networks=[good_networks_struct.data(:,1:3) good_networks_struct.data(:,end)];

    for j=1:size(good_networks,1)
            if good_networks(j,4)< raggio_medio_Sup && good_networks(j,4)>raggio_medio_Inf
                selection(k,2:5)=good_networks(j,:);
                selection(k,1)=start_radius(i);
                selection(k,6)=selection(k,4).*selection(k,5); %R*L
                k=k+1;
            end
    end

end

% area_media=mean(selection(:,end));
% area_sd=std(selection(:,end));
% area_mediana=median(selection(:,end));
% 
% normalize_area=(selection(:,end)-mean(selection(:,end)))/std(selection(:,end));
% norm_area=kstest(normalize_area);
% histogram(selection(:,end));
% 
% raggio_medio=mean(selection(:,end-1));
% raggio_sd=std(selection(:,end-1));
% raggio_mediana=median(selection(:,end-1));
% 
% normalize_raggio=(selection(:,end-1)-mean(selection(:,end-1)))/std(selection(:,end-1));
% [norm_raggio,p]=kstest(normalize_raggio);
% figure
% histogram(selection(:,end-1));
% incertezza_raggio=(max(selection(:,end-1))-min(selection(:,end-1)))/2;
% 
% area_laterale=selection(:,4).*selection(:,5);
% figure
% histogram(area_laterale);
% figure
% histogram(normalize_area);
% incertezza_area=(max(area_laterale)-min(area_laterale))/2;
% delta_perc=incertezza_area/mean(area_laterale);

 T = table(selection(:,1),selection(:,2),selection(:,3),selection(:,4),selection(:,5),selection(:,6),'VariableNames', {'starting_radius' 'N_geo' 'N_iter_radius' 'total_length' 'average_radius' 'lateral_surface'});
 fileName = sprintf('RETI_FINALI.txt');
 writetable(T,fileName,'Delimiter','\t');