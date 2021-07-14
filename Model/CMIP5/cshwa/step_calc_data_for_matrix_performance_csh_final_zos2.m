clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

load('result_model_performance_add_SSH.mat','sort_num') %for getting sort_num
load(strcat(result_path,'/','processing_result_v2_heatcontent_NWP_csh2.mat'),'csh_index');

%-------STEP1.3 Load result data(zos) ------%

 %Load SST data
idx_ocean = 2;
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


%------STEP2. Normalization Process(Take mean)------%
mean_RMSE_SSH = mean(RMSE_SSH);
mean_trend_SSH = mean(trend_SSH);

RMSE_SSH_nor = RMSE_SSH/mean_RMSE_SSH;
trend_SSH_nor = trend_SSH/mean_trend_SSH;

box_data = [RMSE_SSH_nor,trend_SSH_nor];


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

bar(1:10,I(sort_num(1:10)),'BaseValue',0,'FaceColor','blue');
bar(11:20,I(sort_num(11:20)),'BaseValue',0,'FaceColor','red');
bar(21:30,I(sort_num(21:30)),'BaseValue',0,'FaceColor','green');
bar(31:41,I(sort_num(31:41)),'BaseValue',0,'FaceColor','cyan');
set(gca,'FontSize',16,'xlim',[0 42]);

%%% 24 model
[sort_idx1, sort_num1] = sort(I) ; 
name14 = model_name(sort_num1);
%match the 24 model name
csh_name = {'CCSM4', 'CESM1-CAM5', 'CSIRO-Mk3-6-0', 'GFDL-ESM2G', 'GFDL-ESM2M','HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC5', 'MIROC-ESM-CHEM','NorESM1-M', 'NorESM1-ME', 'bcc-csm1-1', 'bcc-csm1-1-m', 'CESM1-WACCM', 'CNRM-CM5','CanESM2', 'EC-EARTH', 'FGOALS-g2', 'GISS-E2-R', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3'};
for i = 1:24
 srt_index(i) = find(strcmp(csh_name(i),model_name(sort_num1)) == 1);
end
csh_index = sort(srt_index);

figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',name14(csh_index),...
    'XTick',[1:24],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');

bar(1:8,I(sort_num1(csh_index(1:8))),'BaseValue',0,'FaceColor','blue');
bar(9:16,I(sort_num1(csh_index(9:16))),'BaseValue',0,'FaceColor','red');
bar(17:24,I(sort_num1(csh_index(17:24))),'BaseValue',0,'FaceColor','green');
% bar(31:41,I(sort_num(31:41)),'BaseValue',0,'FaceColor','cyan');
xlabel('Model Name','fontsize',25);
ylabel('Performance Index','fontsize',25);
set(gca,'FontSize',16,'xlim',[0 25]);


%%%%% 14model
[sort_idx1, sort_num1] = sort(I) ; 
name14 = model_name(sort_num1);
%match the 24 model name
csh_name = {'CCSM4', 'CESM1-CAM5', 'CSIRO-Mk3-6-0', 'GFDL-ESM2G', 'GFDL-ESM2M','HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC5', 'MIROC-ESM-CHEM','NorESM1-M', 'NorESM1-ME', 'bcc-csm1-1', 'bcc-csm1-1-m', 'CESM1-WACCM', 'CNRM-CM5','CanESM2', 'EC-EARTH', 'FGOALS-g2', 'GISS-E2-R', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3'};
for i = 1:24
 srt_index(i) = find(strcmp(csh_name(i),model_name(sort_num1)) == 1);
end
csh_index = sort(srt_index);

in_index = [5,6,7,10,11,12,13,14,15,16,18,20,22,23];
in_index1 = [5,6,7,10,11];
in_index2 = [12,13,14,15,16];
in_index3 = [18,20,22,23];

figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);
% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',name14(csh_index(in_index)),...
    'XTick',[1:14],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');

bar(1:5,I(sort_num1(csh_index(in_index1))),'BaseValue',0,'FaceColor','blue');
bar(6:10,I(sort_num1(csh_index(in_index2))),'BaseValue',0,'FaceColor','red');
bar(11:14,I(sort_num1(csh_index(in_index3))),'BaseValue',0,'FaceColor','green');
set(gca,'FontSize',16,'xlim',[0 15]);
xlabel('Model Name','fontsize',25);
ylabel('Performance Index','fontsize',25);
set(gca,'linew',1);


%%%% re-plot
[sort_idx1, sort_num1] = sort(I) ; 
name14 = model_name(sort_num1);
%match the 24 model name
csh_name = {'CCSM4', 'CESM1-CAM5', 'CSIRO-Mk3-6-0', 'GFDL-ESM2G', 'GFDL-ESM2M','HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC5', 'MIROC-ESM-CHEM','NorESM1-M', 'NorESM1-ME', 'bcc-csm1-1', 'bcc-csm1-1-m', 'CESM1-WACCM', 'CNRM-CM5','CanESM2', 'EC-EARTH', 'FGOALS-g2', 'GISS-E2-R', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3'};
for i = 1:24
 srt_index(i) = find(strcmp(csh_name(i),model_name(sort_num1)) == 1);
end
csh_index = sort(srt_index);

in_index = [5,6,7,10,11,12,13,14,15,16,18,20,22,23];
in_index1 = [5,6,7,10,11];
in_index2 = [12,13,14,15,16];
in_index3 = [18,20,22,23];

subplot(2,10,[1 9]); hold on;
box(gca,'on');
set(gca,'XTick',[1:14]);
set(gca,'XTickLabel',name14(csh_index(in_index)));
set(gca,'XTickLabelRotation',45);
bar(1:5,I(sort_num1(csh_index(in_index1))),'BaseValue',0,'FaceColor','blue');
bar(6:10,I(sort_num1(csh_index(in_index2))),'BaseValue',0,'FaceColor','red');
bar(11:14,I(sort_num1(csh_index(in_index3))),'BaseValue',0,'FaceColor','green');
set(gca,'FontSize',10,'xlim',[0 15]);
xlabel('Model Name','fontsize',10);
ylabel('Performance Index','fontsize',10);
set(gca,'linew',1);


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