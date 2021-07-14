
clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%--------STEP1.1 Load result data(SSH)-----%
%processing result file
load(strcat('/home/ygkim/ipcc/data/zos/new','/','processing_result_v2_zos_NWP_csh3.mat'));
RMSE_SSH = [ens_ocean_RMSE2(:,2)];
trend_SSH = [ens_ocean_trend_mod2(:,2)];
trend_SSH_obs = ocean_trend_obs(2,1);
trend_SSH = abs(trend_SSH-trend_SSH_obs);

%------STEP2. Normalization Process(Take mean)------%
 %load mean data
 load(strcat('./','mean_cshwa_24_3.mat'));

RMSE_SSH_nor = RMSE_SSH/mean_RMSE_SSH;
trend_SSH_nor = trend_SSH/mean_trend_SSH;

box_data = [RMSE_SSH_nor,trend_SSH_nor];

%------STEP3. make performance index ------%
I = zeros(1,24);%kyy
I = mean(box_data,2);

ens_I = I;

plot(I,'linew',1);
xlabel('Ensemble Model Number','fontsize',25);
ylabel('Performance Index','fontsize',25);
set(gca,'fontsize',25);
set(gca,'linew',1);



%------STEP4. save data ------%
% save('result_model_performance_ensemble_3_zos.mat','ens_I','box_data');


