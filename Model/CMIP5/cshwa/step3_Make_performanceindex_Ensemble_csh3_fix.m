
clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%--------STEP1.1 Load result data(OHC)-----%
%processing result file
load(strcat(result_path,'/processing_result_v2_heatcontent_NWP_csh3.mat'));

RMSE_OHC = [ens_ocean_RMSE2(:,2)];
trend_OHC = [ens_ocean_trend_mod2(:,2)];
trend_OHC_obs = ocean_trend_obs(2,1);
trend_OHC = abs(trend_OHC-trend_OHC_obs);



%--------STEP1.2 Load result data(SST)-----%
%processing result file
load(strcat('/home/ygkim/ipcc/data/tos/new','/','processing_result_v2_tos_NWP_csh3.mat')); 

RMSE_SST = [ens_ocean_RMSE2(:,2)];
trend_SST = [ens_ocean_trend_mod2(:,2)];
trend_SST_obs = ocean_trend_obs(2,1);
trend_SST = abs(trend_SST-trend_SST_obs);



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

RMSE_OHC_nor = RMSE_OHC/mean_RMSE_OHC;
trend_OHC_nor = trend_OHC/mean_trend_OHC;
RMSE_SST_nor = RMSE_SST/mean_RMSE_SST;
trend_SST_nor = trend_SST/mean_trend_SST;
RMSE_SSH_nor = RMSE_SSH/mean_RMSE_SSH;
trend_SSH_nor = trend_SSH/mean_trend_SSH;

box_data = [RMSE_OHC_nor,trend_OHC_nor,RMSE_SST_nor,trend_SST_nor,RMSE_SSH_nor,trend_SSH_nor];


%------STEP3. make performance index ------%
I = zeros(1,24);%kyy
I = mean(box_data,2);

ens_I = I;

plot(I,'linew',1);set(gca,'linew',1);
xlabel('Ensemble Model Number','fontsize',25);
ylabel('Performance Index','fontsize',25);
set(gca,'fontsize',25);

%------STEP4. save data ------%
save('result_model_performance_ensemble.mat','ens_I','box_data');
