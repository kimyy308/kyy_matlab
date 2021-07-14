
clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
result_path3 = '/home/ygkim/ipcc/data/tos/new';
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%--------STEP1. Load result data-----%
%processing result file
load(strcat(result_path,'/processing_result_heatcontent_NWP.mat'));
%longitude,latitude meshgride data                 : lon2,lat2
%basin parameter(NA:1 SA:2 NP:3 SP:4 no:5 SO:6 AO:7: new_basin,ocean_num
%observation heat content(time,lat,lon)            : obs_heat2
%model heat content(model,time,lat,lon)            : model_heat2
%model name information(i-3)th                     : root_path2(i).name
%RMSE,bias (model,lat,lon)                         : RMSE,bias
%total area mean of RMSE,bias                      :total_RMSE,total_bias
%total area mean of RMSE's std                     :total_std
%bias (model,time,lat,lon),mean bias               :bias2,total_bias2
%each ocean mean of RMSE,bias                      :ocean_RMSE,ocean_bias
%each ocean mean of RMSE's std                     :ocean_std
%total trend of obs,mod                   :total_trend_obs,total_trend_mod
%ocean trend of obs,mod                   :ocean_trend_obs,ocean_trend_mod
%each grid trend of obs,mod               :each_trend_obs,each_trend_mod
%std of RMSE model performance(model)              :se,se_ocean
%ensemble result                                   :en



%find model name
for i = 4:length(root_path2)
    model_name(i-3) = {root_path2(i).name};
end
model_name(end+1) = {'Ensemble'};



%Talyer diagram example 

%Need variables : Std, RMSE, Corr
ens_total_mod_heat_new = nanmean(total_mod_heat_new);
data = [total_obs_heat_new;total_mod_heat_new;ens_total_mod_heat_new];
%data = [total_obs_tos_new;total_mod_tos_new];

%models , time

%----1.calculate std of data ----%
std_data = zeros(length(data),1);
for i = 1:length(std_data)
    mod_data = data(i,:);
    std_data(i) = sqrt(sum((mod_data-mean(mod_data)).^2)/numel(mod_data));
end


%----2.calculate RMSE of data ----%
RMSE_data = zeros(length(std_data),1);
for i = 1:length(std_data)
    obs_data = data(1,:);
    mod_data = data(i,:);
    RMSE_data(i) = sqrt(sum( ((mod_data-mean(mod_data))-(obs_data-mean(obs_data))).^2)/numel(obs_data));
end

%----3.calculate corr of data ----%
cor_data = zeros(length(std_data),1);
for i = 1:length(std_data)
    obs_data = data(1,:);
    mod_data = data(i,:);
    cor_data(i) = sum( (mod_data-mean(mod_data)).*(obs_data-mean(obs_data)))/numel(obs_data)/(std_data(1)*std_data(i));
end

check = RMSE_data - sqrt( std_data.^2 + std_data(1)^2 - 2.*std_data*std_data(1).*cor_data);

%Make Taylor diagram
std_data = std_data/10^8;
RMSE_data = RMSE_data/10^8;


figure1 = figure;hold on;
[hp, ht, axl] = taylordiag(std_data,RMSE_data,cor_data,'titleRMS',10^8,'titleSTD',10^8,...
                'tickRMS',[0.5:0.5:2],'tickRMSangle',100,'limSTD',2.2);
markerstyle = '+o*.xsd^v><ph';            
for i =1:length(ht)
    set(ht(i),'fontsize',0.1);
    set(hp(i),'linestyle','none');
    if i ==1
       set(ht(i),'fontsize',8,'string','observation','color',[0 0 0]);
       set(hp(i),'markersize',12,'color',[0 0 0]);
    elseif i==50
        set(ht(i),'fontsize',8,'string','Ens','color',[0 0 0]);
        set(hp(i),'markersize',12,'color',[0 0 0]);
    elseif i <= 10
        set(ht(i),'string',num2str(i-1),'color',[1 0 0]*i/10);
        set(hp(i),'marker',markerstyle(i),'markersize',4,'color',[1 0 0],'markerfacecolor',[1 0 0]);
    elseif i<= 20
        set(ht(i),'string',num2str(i-1),'color',[1 1 0]*(i-10)/10);
        set(hp(i),'marker',markerstyle(i-10),'markersize',4,'color',[0.9 0.8 0],'markerfacecolor',[0.8 0.5 0]);
    elseif i <=30
        set(ht(i),'string',num2str(i-1),'color',[0.5 1 0.5]*(i-20)/10);
        set(hp(i),'marker',markerstyle(i-20),'markersize',4,'color',[1 0.5 0],'markerfacecolor',[1 0.5 0]);
    elseif i <=40
        set(ht(i),'string',num2str(i-1),'color',[0.5 1 1]*(i-30)/10);
        set(hp(i),'marker',markerstyle(i-30),'markersize',4,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0]);
    else
        set(ht(i),'string',num2str(i-1),'color',[0 0 1]*(i-40)/10);
        set(hp(i),'marker',markerstyle(i-40),'markersize',4,'color',[0 0 1],'markerfacecolor',[0 0 1]);
    end
end

tt = 'Taylor Diagram for UOHC';
title(tt);
model_name2 = {'a','b','c'};
model_name2 = [model_name2,model_name];
legend(model_name2);

ht = axl(2).handle;
for i = 1:length(ht)
    set(ht(i),'fontsize',10,'fontweight','normal');
end
set(axl(1).handle,'fontweight','normal');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;
