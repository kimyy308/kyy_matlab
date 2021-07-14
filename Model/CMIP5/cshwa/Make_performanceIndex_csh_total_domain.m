
clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%--------STEP1.1 Load result data(OHC)-----%
%processing result file
load(strcat(result_path,'/processing_result_v2_heatcontent_NWP_csh2.mat'));

RMSE_OHC = [ens_ocean_RMSE2(:,1)];
trend_OHC = [ens_ocean_trend_mod2(:,1)];
trend_OHC_obs = ocean_trend_obs(2,1);
trend_OHC = abs(trend_OHC-trend_OHC_obs);



%--------STEP1.2 Load result data(SST)-----%
%processing result file
load(strcat('/home/ygkim/ipcc/data/tos/new','/','processing_result_v2_tos_NWP_csh2.mat')); 

RMSE_SST = [ens_ocean_RMSE2(:,1)];
trend_SST = [ens_ocean_trend_mod2(:,1)];
trend_SST_obs = ocean_trend_obs(2,1);
trend_SST = abs(trend_SST-trend_SST_obs);



%--------STEP1.1 Load result data(SSH)-----%
%processing result file
load(strcat('/home/ygkim/ipcc/data/zos/new','/','processing_result_v2_zos_NWP_csh2.mat'));
RMSE_SSH = [ens_ocean_RMSE2(:,1)];
trend_SSH = [ens_ocean_trend_mod2(:,1)];
trend_SSH_obs = ocean_trend_obs(2,1);
trend_SSH = abs(trend_SSH-trend_SSH_obs);



%------STEP2. Normalization Process(Take mean)------%
 %load mean data
%  load(strcat('./','result_model_performance_add_SSH.mat'));
% mean_RMSE_OHC = mean(RMSE_OHC);
% mean_trend_OHC = mean(trend_OHC);
% mean_RMSE_SST = mean(RMSE_SST);
% mean_trend_SST = mean(trend_SST);
% mean_RMSE_SSH = mean(RMSE_SSH);
% mean_trend_SSH = mean(trend_SSH);
load(strcat('./','mean_cshwa_24_total_domain.mat'));

RMSE_OHC_nor = RMSE_OHC/mean_RMSE_OHC;
trend_OHC_nor = trend_OHC/mean_trend_OHC;
RMSE_SST_nor = RMSE_SST/mean_RMSE_SST;
trend_SST_nor = trend_SST/mean_trend_SST;
RMSE_SSH_nor = RMSE_SSH/mean_RMSE_SSH;
trend_SSH_nor = trend_SSH/mean_trend_SSH;

box_data = [RMSE_OHC_nor,trend_OHC_nor,RMSE_SST_nor,trend_SST_nor,RMSE_SSH_nor,trend_SSH_nor];


%------STEP3. make performance index ------%
% I = zeros(1,41);%kyy
% I = mean(box_data,2);
% cshwa_num = [7,9,12,16,18,20,21,22,25,26,29,34,35,36,38,39,41,42,43,45,47,48,49,50]; %cshwa
% cshwa_num =cshwa_num - 3;
I = zeros(1,24);
I = mean(box_data,2);

%---STEP4.2. Bar plot for performance index---%

%find model name

%model_name(end+1) = {'Ensemble'};
name24 = model_name(sort_num(csh_index));

% [sort_idx, sort_num] = sort(I) ;
I = I;

% select24 = [1, 3, 5, 6, 7, 8, 9, 11, 16, 17, 19, 20 ,21, 26, 27, 29, 31, 32, 34, 35, 36, 37, 39, 41];
% select_24 = sort(select24)
% sort_num_24 = sort_num(select_24);

% I_24 = I(sort_num_24);
% [sort_idx_csh, sort_num_csh] = sort(I_24) 
 
figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',name24(sort_num),...
    'XTick',[1:24],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');

bar(1:8,I(sort_num(1:8)),'BaseValue',0,'FaceColor','blue');
bar(9:16,I(sort_num(9:16)),'BaseValue',0,'FaceColor','red');
bar(17:24,I(sort_num(17:24)),'BaseValue',0,'FaceColor','green');


set(gca,'FontSize',16,'xlim',[0 25]);
% Create ylabel
xlabel('Model Name','fontsize',25);
ylabel('Performance Index','fontsize',25);
tt = strcat('Model Performance Result(only SSHA)');
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% multi_model mean %%%%
% ens_I = I(sort_num);

figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);
hold on;

plot(I(sort_num(csh_index)),'linewidth',3);hold on;
%plot(box_data,'linewidth',1);
%plot(I+se*1.96);plot(I-se*1.96);
plot(zeros(24,1)+min(I(sort_num(csh_index))),':');
set(gca,'FontSize',16,'xlim',[1 24]);
% Create ylabel
xlabel('Ensemble Model Number','fontsize',25);
tt = strcat('Ensemble Model Performance');
title(tt,'fontsize',30,'fontweight','bold');
%legend('Performance Index','RMSE OHC','trend OHC','RMSE SST','trend SST','Location','southeast');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;

ens_cal = 0
for i = 1:length(ens_I)
    if i == 1
        ens_cal(1) = ens_I(1)
    else
    ens_cal(i) = (sum(ens_I(1:i)))/i
    end
end

plot(ens_cal)

ens_cal = 0
for i = 1:length(ens_I)
    if i == 1
        ens_cal(1) = ens_I(1);
    else
    ens_cal(i) = (mean(ens_I(1:i)));
    end
end

plot(ens_cal)



%------STEP4. save data ------%
save('result_model_performance_ensemble_csh.mat','ens_I','box_data');
