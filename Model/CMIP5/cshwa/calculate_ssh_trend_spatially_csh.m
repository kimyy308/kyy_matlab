% THis code is used for calculating RMSE of NWP
clc;clear all;close all;
%time range : 1976.01 ~ 2005.12
%area : global

% result_path = '/home/ygkim/temp_work'; %result path -cshwa
result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/zos/new';
result_path3 = '/home/ygkim/ipcc/data/zos';
root_path ='/data/ygkim/ipcc/'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%--------STEP1. Load processing result data(Original, all variables)-----%
%--------------Additional mask process ----------------------------------%

% load(strcat(result_path2,'/','processing_result_v4_area_ratio.mat'));
load(strcat(result_path2,'/','recon_ssha_1976_2005_yearly.mat'));
%longitude,latitude meshgride data                 : lon2,lat2
%observation heat content(time,lat,lon)            : obs_heat2
%model heat content(model,time,lat,lon)            : model_heat2
%model name information(i-3)th                     : root_path2(i).name
%RMSE,bias (model,lat,lon)                         : RMSE,bia
%bias (model,time,lat,lon)                         : bias2


%------STEP2. Calculate RMSE and Trend of NWP

 %---2.1. Determine NWP domain
%  lon_idx = [115.5,164.5]; %lon : 115 ~164
%  lat_idx = [15.5,52.5];   %lat : 15 ~ 52
  lon_idx = [115,164]; %lon : 115 ~164
 lat_idx = [15,52];   %lat : 15 ~ 52
 
 
% nwp_idx{1} = find(lon2>=lon_idx(1) & lon2<=lon_idx(2) & lat2>=lat_idx(1) & lat2<=lat_idx(2) );
nwp_idx{1} = find(lon2==lon_idx(1) & lon2==lon_idx(2) & lat2==lat_idx(1) & lat2==lat_idx(2) );

% aa = find(lon2>=lon_idx(1) & lon2<=lon_idx(2) & lat2>=lat_idx(1)-1 & lat2<= lat_idx(1)+1); %east bc
% bb = find(lon2>=lon_idx(2)-1 & lon2<=lon_idx(2)+1 & lat2>=lat_idx(1) & lat2<= lat_idx(2)); %south bc
% cc = find(lon2>=lon_idx(1)-1 & lon2<=lon_idx(1)+1 & lat2>=lat_idx(1) & lat2<=lat_idx(2) ); %west bc
aa = find(lon2==lon_idx(1) & lon2==lon_idx(2) & lat2==lat_idx(1)-1 & lat2== lat_idx(1)+1); %east bc
bb = find(lon2==lon_idx(2)-1 & lon2==lon_idx(2)+1 & lat2==lat_idx(1) & lat2== lat_idx(2)); %south bc
cc = find(lon2==lon_idx(1)-1 & lon2==lon_idx(1)+1 & lat2==lat_idx(1) & lat2==lat_idx(2) ); %west bc
nwp_idx{2} = sort(unique([aa',bb', cc']));

ocean_num = 2;

 %---2.3. calculate RMSE and trend results

% spatial trend of obs
 x = [1:30]';
 ocean_trend_obs = NaN(180,360);
 ocean_trend_obs2 = NaN(180,360);
 for i = 1:180
    for j= 1:360
    data = squeeze(obs_zos3(:,i,j));
    a = polyfit(x,data,1);
    ocean_trend_obs(i,j) = a(1);
    ocean_trend_obs2(i,j) = a(2);
end
end

 
 % spatial trend of model
 x = [1:30]; %time index
 ocean_trend_model = NaN(41,180,360);
 ocean_trend_model2 = NaN(41,180,360);
for model = 1:41
 for i = 1:180
    for j= 1:360
    data = squeeze(model_zos3(model,:,i,j));
    a = polyfit(x,data,1);
    ocean_trend_model(model,i,j) = a(1);
    ocean_trend_model2(model,i,j) = a(2);
end
end
end

save('compare_result_zos.mat');

