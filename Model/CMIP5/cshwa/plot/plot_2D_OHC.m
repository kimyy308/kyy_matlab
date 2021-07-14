clc;clear all;close all;

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%--------STEP1. Load result data-----%
%processing result file
load(strcat(result_path2,'/processing_result_v4_new.mat'));
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

%% ------ Set the map lon,lat ------
lon_idx = [115.5,164.5]; %lon : 115 ~164
lat_idx = [15.5,52.5];   %lat : 15 ~ 52





%% ------ Plot Map of each grid RMSE, bias,Trend

%RMSE data
data = ens_RMSE;

figure1 = figure;
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 2 1]*6);
m_proj('mercator','lon',lon_idx,'lat',lat_idx);
m_gshhs_l('color','k');
%m_gshhs_l('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
m_pcolor(lon2,lat2,data);shading flat;
set(gca,'FontSize',14);
tt = strcat('OHC RMSE of ensemble model(each grid)');
title(tt,'fontsize',20,'fontweight','bold');  
colorbar('EastOutSide');colormap jet;
axis tight;
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;


%bias data
data = ens_bias;

figure1 = figure;
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 2 1]*6);
m_proj('mercator','lon',lon_idx,'lat',lat_idx);
m_gshhs_l('color','k');
%m_gshhs_l('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
m_pcolor(lon2,lat2,data);shading flat;
set(gca,'FontSize',14);
tt = strcat('OHC Bias of ensemble model(each grid)');
title(tt,'fontsize',20,'fontweight','bold');  
colorbar('EastOutSide');colormap jet;
axis tight;
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;

%trend data
data = ens_trend;
data(data==0) = NaN;

figure1 = figure;
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 2 1]*6);
m_proj('mercator','lon',lon_idx,'lat',lat_idx);
m_gshhs_l('color','k');
%m_gshhs_l('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
m_pcolor(lon2,lat2,data);shading flat;
set(gca,'FontSize',14);
tt = strcat('OHC Trend of ensemble model(each grid)');
title(tt,'fontsize',20,'fontweight','bold');  
colorbar('EastOutSide');colormap jet;
axis tight;caxis([-50 50]*10^6);
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;


%trend data
data = trend_obs;
data(data==0) = NaN;

 %smoothing data
xx = [0.5:4:360];
yy = [-89.5:4:90];
[xx,yy]=meshgrid(xx,yy);
data2 = interp2(lon2,lat2,data,xx,yy);

figure1 = figure;
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 2 1]*6);
m_proj('mercator','lon',lon_idx,'lat',lat_idx);
m_gshhs_l('color','k');
%m_gshhs_l('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
m_pcolor(lon2,lat2,data);shading flat;
%m_pcolor(xx,yy,data2);shading flat;
set(gca,'FontSize',14);
tt = strcat('OHC Trend of observation(each grid)');
title(tt,'fontsize',20,'fontweight','bold');  
colorbar('EastOutSide');colormap jet;
axis tight;caxis([-50 50]*10^6);
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;


%trend data
data = ens_trend - trend_obs;
data(data==0) = NaN;

 %smoothing data
xx = [0.5:4:360];
yy = [-89.5:4:90];
[xx,yy]=meshgrid(xx,yy);
data2 = interp2(lon2,lat2,data,xx,yy);

figure1 = figure;
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 2 1]*6);
m_proj('mercator','lon',lon_idx,'lat',lat_idx);
m_gshhs_l('color','k');
%m_gshhs_l('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;
m_pcolor(lon2,lat2,data);shading flat;
%m_pcolor(xx,yy,data2);shading flat;
set(gca,'FontSize',14);
tt = strcat('OHC - Trend difference of observation(each grid)');
title(tt,'fontsize',20,'fontweight','bold');  
colorbar('EastOutSide');colormap jet;
axis tight;%caxis([-0.1 0.1]*0.03);
saveas(gca,strcat('./figure/',tt,'.fig'),'fig');
close all;

