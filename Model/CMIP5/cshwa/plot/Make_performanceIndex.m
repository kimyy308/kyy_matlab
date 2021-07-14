
clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%--------STEP1.1 Load result data(OHC)-----%
%processing result file
load(strcat(result_path,'/processing_result_heatcontent_NWP.mat'));
idx_ocean = 2;

RMSE_OHC = [ocean_RMSE(:,idx_ocean)];
trend_OHC = [ocean_trend_mod(:,idx_ocean)];
trend_OHC_obs = ocean_trend_obs(idx_ocean);
trend_OHC = abs(trend_OHC-trend_OHC_obs);


%-------STEP1.2 Load result data(SST) ------%

%RMSE-RMSE
 %Load SST data
load(strcat('/home/ygkim/ipcc/data/tos/new','/','processing_result_tos_NWP.mat')); 
RMSE_SST = [ocean_RMSE(:,idx_ocean)];
trend_SST = [ocean_trend_mod(:,idx_ocean)];
trend_SST_obs = ocean_trend_obs(idx_ocean);
trend_SST = abs(trend_SST-trend_SST_obs);

%------STEP2. Normalization Process(Take mean)------%
mean_RMSE_OHC = mean(RMSE_OHC);
mean_trend_OHC = mean(trend_OHC);
mean_RMSE_SST = mean(RMSE_SST);
mean_trend_SST = mean(trend_SST);

RMSE_OHC_nor = RMSE_OHC/mean_RMSE_OHC;
trend_OHC_nor = trend_OHC/mean_trend_OHC;
RMSE_SST_nor = RMSE_SST/mean_RMSE_SST;
trend_SST_nor = trend_SST/mean_trend_SST;

box_data = [RMSE_OHC_nor,trend_OHC_nor,RMSE_SST_nor,trend_SST_nor];


%------STEP3. make performance index ------%
I = zeros(1,48);
I = mean(box_data,2);


%------STEP4. Plot data ------%

%---STEP4.1. box plot for each data ---%
 %set variables
label = {'RMSE of OHC','Trend of OHC','RMSE of SST','Trend of SST'};
data = box_data;

figure;hold on;
boxplot(data,'Notch','on','Labels',label);
set(gca,'FontSize',20);
% Create titel,label
xlabel('Variables','fontsize',25);
ylabel('Normalized Value','fontsize',25);
tt ='Normalized results';
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;

%---STEP4.2. Bar plot for performance index---%

%find model name
for i = 4:length(root_path2)
    model_name(i-3) = {root_path2(i).name};
end
%model_name(end+1) = {'Ensemble'};

[sort_idx, sort_num] = sort(I) ;
I = I;

figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',model_name(sort_num),...
    'XTick',[1:48],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');

%bar plot(blue to red)
%bar(1,I(sort_num(1)),'BaseValue',1,'FaceColor',rgb('dimgray'));hold on;
bar(1:11,I(sort_num(1:11)),'BaseValue',1,'FaceColor',rgb('indigo'));
bar(12:21,I(sort_num(12:21)),'BaseValue',1,'FaceColor',rgb('royalblue'));
bar(22:31,I(sort_num(22:31)),'BaseValue',1,'FaceColor',rgb('ForestGreen'));
bar(32:41,I(sort_num(32:41)),'BaseValue',1,'FaceColor',rgb('darkorange'));
bar(42:48,I(sort_num(42:48)),'BaseValue',1,'FaceColor',rgb('crimson'));

set(gca,'FontSize',16,'xlim',[0 50]);
% Create ylabel
xlabel('Model Name','fontsize',25);
ylabel('Performance Index','fontsize',25);
tt = strcat('Model Performance Result(Boundary Area)');
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;

save('result_model_performance.mat','I','sort_idx','sort_num',...
    'mean_RMSE_OHC','mean_trend_OHC','mean_RMSE_SST','mean_trend_SST');
