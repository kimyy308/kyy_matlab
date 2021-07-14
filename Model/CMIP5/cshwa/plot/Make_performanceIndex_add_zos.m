
clc;clear all;close all;
%consider NWP
load('result_model_performance_add_SSH.mat');
sort_num_NWP = sort_num;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 


%--------STEP1.1 Load result data(OHC)-----%
%processing result file
load(strcat(result_path2,'/processing_result_v4_area_ratio.mat'));

RMSE_OHC = [total_RMSE];
trend_OHC = [total_trend_mod];
trend_OHC_obs = total_trend_obs;
trend_OHC = abs(trend_OHC-trend_OHC_obs);


%-------STEP1.2 Load result data(SST) ------%

 %Load SST data
load(strcat('/home/ygkim/ipcc/data/tos/new','/','processing_result_v4_area_ratio.mat')); 
RMSE_SST = [total_RMSE];
trend_SST = [total_trend_mod];	%Pacific Ocean
trend_SST_obs = total_trend_obs;
trend_SST = abs(trend_SST-trend_SST_obs);

%-------STEP1.3 Load result data(zos) ------%

 %Load SST data
load(strcat('/home/ygkim/ipcc/data/zos/new','/','processing_result_v4_area_ratio.mat')); 
RMSE_SSH = [total_RMSE];
trend_SSH = [total_trend_mod];
trend_SSH_obs = total_trend_obs;
trend_SSH = abs(trend_SSH-trend_SSH_obs);

%------STEP1.4 Select models(included all component)
%find zos file
zos_idx = [];
zos_num = 0;
for i = 4:length(root_path2);
    %check zos folder
    file_zos =dir(strcat(root_path,'/',root_path2(i).name,'/','*zos_*'));
    [m,n] = size(file_zos);
    if m~=0
        zos_num = zos_num +1;
        zos_idx(end+1) = i;
    end
end

RMSE_OHC = RMSE_OHC(zos_idx-3);
trend_OHC = trend_OHC(zos_idx-3);
RMSE_SST = RMSE_SST(zos_idx-3);
trend_SST = trend_SST(zos_idx-3);


%------STEP2. Normalization Process(Take mean)------%
mean_RMSE_OHC = mean(RMSE_OHC);
mean_trend_OHC = mean(trend_OHC);
mean_RMSE_SST = mean(RMSE_SST);
mean_trend_SST = mean(trend_SST);
mean_RMSE_SSH = mean(RMSE_SSH);
mean_trend_SSH = mean(trend_SSH);

RMSE_OHC_nor = RMSE_OHC/mean_RMSE_OHC;
trend_OHC_nor = trend_OHC/mean_trend_OHC;
RMSE_SST_nor = RMSE_SST/mean_RMSE_SST;
trend_SST_nor = trend_SST/mean_trend_SST;
RMSE_SSH_nor = RMSE_SSH/mean_RMSE_SSH;
trend_SSH_nor = trend_SSH/mean_trend_SSH;

box_data = [RMSE_OHC_nor,trend_OHC_nor,RMSE_SST_nor,trend_SST_nor,RMSE_SSH_nor,trend_SSH_nor];


%------STEP3. make performance index ------%
I = zeros(1,length(zos_idx));
I = mean(box_data,2);


%------STEP4. Plot data ------%

%---STEP4.1. box plot for each data ---%
 %set variables
label = {'RMSE of OHC','Trend of OHC','RMSE of SST','Trend of SST','RMSE of SSH','Trend of SSH'};
data = box_data;

figure;hold on;
boxplot(data,'Notch','on','Labels',label);
set(gca,'FontSize',20);
% Create titel,label
xlabel('Variables','fontsize',25);
ylabel('Normalized Value','fontsize',25);
tt ='Normalized results(Global)';
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;

%---STEP4.2. Bar plot for performance index---%

%find model name
for i = 4:length(root_path2)
    model_name(i-3) = {root_path2(i).name};
end
%model_name(end+1) = {'Ensemble'};
model_name = model_name(zos_idx-3);

[sort_idx, sort_num] = sort(I) ;
I = I;



figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',model_name(sort_num),...
    'XTick',[1:41],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');

%bar plot(blue to red)
%bar(1,I(sort_num(1)),'BaseValue',1,'FaceColor',rgb('dimgray'));hold on;
% bar(1:11,I(sort_num(1:11)),'BaseValue',0,'FaceColor',rgb('indigo'));
% bar(12:21,I(sort_num(12:21)),'BaseValue',0,'FaceColor',rgb('royalblue'));
% bar(22:31,I(sort_num(22:31)),'BaseValue',0,'FaceColor',rgb('ForestGreen'));
% bar(32:41,I(sort_num(32:41)),'BaseValue',0,'FaceColor',rgb('darkorange'));
%bar(42:48,I(sort_num(42:48)),'BaseValue',1,'FaceColor',rgb('crimson'));
bar(1:11,I(sort_num(1:11)),'BaseValue',0,'FaceColor','blue');
bar(12:21,I(sort_num(12:21)),'BaseValue',0,'FaceColor','red');
bar(22:31,I(sort_num(22:31)),'BaseValue',0,'FaceColor','green');
bar(32:41,I(sort_num(32:41)),'BaseValue',0,'FaceColor','cyan');

set(gca,'FontSize',16,'xlim',[0 42]);
% Create ylabel
xlabel('Model Name','fontsize',25);
ylabel('Performance Index','fontsize',25);
tt = strcat('Model Performance Result(Global)');
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;

%save('result_model_performance_add_SSH.mat','I','sort_idx','sort_num',...
 %   'mean_RMSE_OHC','mean_trend_OHC','mean_RMSE_SST','mean_trend_SST','mean_RMSE_SSH','mean_trend_SSH');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cshwa 24 model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold off;
select24 = [1, 2, 3, 7, 8, 9, 13, 14, 15, 23, 25, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40, 41, 11];
select_24 = sort(select24)
sort_num_24 = sort_num(select_24);

% I_24 = I(sort_num_24);
% [sort_idx_csh, sort_num_csh] = sort(I_24) 

figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',model_name(sort_num_24),...
    'XTick',[1:24],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');

bar(1:8,I(sort_num_24(1:8)),'BaseValue',0,'FaceColor','blue');
bar(9:16,I(sort_num_24(9:16)),'BaseValue',0,'FaceColor','red');
bar(17:24,I(sort_num_24(17:24)),'BaseValue',0,'FaceColor','green');


set(gca,'FontSize',16,'xlim',[0 25]);
% Create ylabel
xlabel('Model Name','fontsize',25);
ylabel('Performance Index','fontsize',25);
tt = strcat('Model Performance Result(Global)');
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;

