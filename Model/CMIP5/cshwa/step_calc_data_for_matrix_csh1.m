clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

load('result_model_performance_add_SSH.mat','sort_num') %for getting sort_num
load(strcat(result_path,'/','processing_result_v2_heatcontent_NWP_csh2.mat'),'csh_index');

%--------STEP1.1 Load result data(OHC)-----%
%processing result file
load(strcat(result_path,'/processing_result_heatcontent_NWP.mat'));
idx_ocean = 2;

RMSE_OHC = [ocean_RMSE(:,idx_ocean)];
trend_OHC = [ocean_trend_mod(:,idx_ocean)];
trend_OHC_obs = ocean_trend_obs(idx_ocean);
trend_OHC = abs(trend_OHC-trend_OHC_obs);


%-------STEP1.2 Load result data(SST) ------%

 %Load SST data
load(strcat('/home/ygkim/ipcc/data/tos/new','/','processing_result_tos_NWP.mat')); 
RMSE_SST = [ocean_RMSE(:,idx_ocean)];
trend_SST = [ocean_trend_mod(:,idx_ocean)];
trend_SST_obs = ocean_trend_obs(idx_ocean);
trend_SST = abs(trend_SST-trend_SST_obs);

%-------STEP1.3 Load result data(zos) ------%

 %Load SST data
load(strcat('/home/ygkim/ipcc/data/zos/new','/','processing_result_zos_NWP.mat')); 
RMSE_SSH = [ocean_RMSE(:,idx_ocean)];
trend_SSH = [ocean_trend_mod(:,idx_ocean)];
trend_SSH_obs = ocean_trend_obs(idx_ocean);
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

%find model name
for i = 4:length(root_path2)
    model_name(i-3) = {root_path2(i).name};
end

model_name = model_name(zos_idx-3);


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



%% Plot data

%plot data
I = zeros(1,41);
I = mean(box_data,2);

% [sort_idx, sort_num] = sort(I) ;

%---STEP4.2. Bar plot for performance index---%

%find model name

%model_name(end+1) = {'Ensemble'};
name14 = model_name(sort_num);

figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',name14,...
    'XTick',[1:41],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');

% bar(1:8,I(sort_num(csh_index(1:8))),'BaseValue',0,'FaceColor','blue');
% bar(9:16,I(sort_num(csh_index(9:16))),'BaseValue',0,'FaceColor','red');
% bar(17:24,I(sort_num(csh_index(17:24))),'BaseValue',0,'FaceColor','green');
% bar(1:8,I(csh_index(1:8)),'BaseValue',0,'FaceColor','blue');
% bar(9:16,I(sort_num(9:16)),'BaseValue',0,'FaceColor','red');
% bar(17:24,I(sort_num(17:24)),'BaseValue',0,'FaceColor','green');
bar(1:10,I(sort_num(1:10)),'BaseValue',0,'FaceColor','blue');
bar(11:20,I(sort_num(11:20)),'BaseValue',0,'FaceColor','red');
bar(21:30,I(sort_num(21:30)),'BaseValue',0,'FaceColor','green');
bar(31:41,I(sort_num(31:41)),'BaseValue',0,'FaceColor','cyan');
set(gca,'FontSize',16,'xlim',[0 42]);

%%% 24 model
name14 = model_name(sort_num);
figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',name14(csh_index),...
    'XTick',[1:24],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');

% bar(1:8,I(sort_num(csh_index(1:8))),'BaseValue',0,'FaceColor','blue');
% bar(9:16,I(sort_num(csh_index(9:16))),'BaseValue',0,'FaceColor','red');
% bar(17:24,I(sort_num(csh_index(17:24))),'BaseValue',0,'FaceColor','green');
% bar(1:8,I(csh_index(1:8)),'BaseValue',0,'FaceColor','blue');
% bar(9:16,I(sort_num(9:16)),'BaseValue',0,'FaceColor','red');
% bar(17:24,I(sort_num(17:24)),'BaseValue',0,'FaceColor','green');
bar(1:8,I(sort_num(csh_index(1:8))),'BaseValue',0,'FaceColor','blue');
bar(9:16,I(sort_num(csh_index(9:16))),'BaseValue',0,'FaceColor','red');
bar(17:24,I(sort_num(csh_index(17:24))),'BaseValue',0,'FaceColor','green');
% bar(31:41,I(sort_num(31:41)),'BaseValue',0,'FaceColor','cyan');
set(gca,'FontSize',16,'xlim',[0 25]);



%---STEP4.1. box plot for each data ---%
 %set variables
 
%label_y = {'RMSE of OHC','bias of OHC','Trend of OHC'};
label_y = {'RMSE of OHC','Trend of OHC','RMSE of SST','Trend of SST','RMSE of SSHA','Trend of SSHA'};
label_x = model_name_zos;
data = box_data;
[m,n] = size(data);
x = [1:m];
y = [1:n];

figure1 = figure;hold on;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

imagesc(x,y,data');
caxis([0 2]);axis('tight');
colorbar('EastOutSide');
set(gca,'XTick',x,'XTickLabel',label_x,'XTicklabelrotation',90,'YTick',y,'YTickLabel',label_y);
set(gca,'xgrid','on','ygrid','on','gridlinestyle','-','xcolor','k','ycolor','k');

tt = strcat('Matrix results');
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;