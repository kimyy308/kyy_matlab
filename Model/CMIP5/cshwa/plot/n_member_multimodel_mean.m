clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
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
%ensemble result                                   :ens_*

%------STEP1.1. Determine best number ------%
% ens_num = 5;
% 
% %Load data - each component
% load(strcat('/home/ygkim/ipcc/script/heatcontent/plot/rmse_bias_new/','result_model_performance.mat'));
% each_box_data = box_data(1:48,:); %individual model box data
% each_sort_num = sort_num;
% each_sort_num(sort_num==49) = []; %remove ensemble results

%Load data 
load(strcat('./','result_model_performance_add_SSH.mat'));

%----Plot Ensemble number vs Model Performance ----%
figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);
hold on;

plot(ens_I,'linewidth',3);hold on;
%plot(box_data,'linewidth',1);
%plot(I+se*1.96);plot(I-se*1.96);
plot(zeros(41,1)+min(ens_I),':');
set(gca,'FontSize',16,'xlim',[1 41]);
% Create ylabel
xlabel('Ensemble Model Number','fontsize',25);
tt = strcat('Ensemble Model Performance');
title(tt,'fontsize',30,'fontweight','bold');
legend('Performance Index','RMSE OHC','trend OHC','RMSE SST','trend SST','Location','southeast');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;



% 
% 
% %----Plot Model Performance vs each model performance
% figure1 = figure;
% set(figure1,'OuterPosition',[0,0,1000,400]);
% set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);
% hold on;
% 
% plot(ens_I,'linewidth',3);hold on;
% plot(each_box_data(each_sort_num,:),'linewidth',1);
% %plot(I+se*1.96);plot(I-se*1.96);
% plot(zeros(48,1)+min(ens_I),':');
% set(gca,'FontSize',16,'xlim',[1 48]);
% % Create ylabel
% xlabel('Ensemble Model Number','fontsize',25);
% tt = strcat('Ensemble Model Performance(each component)');
% title(tt,'fontsize',30,'fontweight','bold');
% legend('Performance Index','RMSE OHC','trend OHC','RMSE SST','trend SST','Location','southeast');
% saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
% close all;

% 
% %------STEP2. Make Table ------
% % row - each region // colum : obs_mean Mean(std) bias ratio (best models)
% 
% data = zeros(8,3);
% %total region
%  %obs_mean
% data(1,1) = nanmean(total_obs_heat_new(:));
%  %Mean(std)
% data(1,2) = mean(ens_total_RMSE_new2(ens_num,:));
% %data(1,3) = std(total_RMSE_new(ens_num,:));
%  %bias
% data(1,3) = mean(ens_total_bias_new2(ens_num,:));
% 
% 
% for j = 1:ocean_num
%     idx = find(new_basin == j);
%      %obs_mean
%     obs = ocean_obs_heat(j,:);
%     data(j+1,1) = nanmean(obs(:));
%      %Mean(std)
%     data(j+1,2) = mean(ens_ocean_RMSE2(ens_num,j));
% %    data(j+1,3) = std(ocean_RMSE(:,j));
%      %bias
%     data(j+1,3) = mean(ens_ocean_bias2(ens_num,j));
% end
% RMSE_data = data;
% 
% %Make trend table
% 
% data = zeros(8,2);
% %total region
%  %obs_mean
% data(1,1) = total_trend_obs_new;
%  %Mean(std)
% data(1,2) = mean(ens_total_trend_mod_new2(ens_num,:));
% %data(1,3) = std(total_trend_mod_new);
%  
% 
% for j = 1:ocean_num
%     idx = find(new_basin == j);
%      %obs_mean
%     obs = ocean_trend_obs(j);
%     model = ens_ocean_trend_mod2(ens_num,j);
%     data(j+1,1) = obs;
%      %Mean(std)
%     data(j+1,2) = mean(model);
%     %data(j+1,3) = std(model);
% end
% trend_data = data;
