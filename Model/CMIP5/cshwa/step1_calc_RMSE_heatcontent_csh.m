% THis code is used for calculating RMSE of NWP
clc;clear all;close all;
%time range : 1976.01 ~ 2005.12
%area : global

result_path = '/home/ygkim/temp_work'; %result path -cshwa
% result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new';
result_path3 = '/home/ygkim/ipcc/data/zos';
root_path ='/data/ygkim/ipcc/'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%--------STEP1. Load processing result data(Original, all variables)-----%
%--------------Additional mask process ----------------------------------%

load(strcat(result_path2,'/','processing_result_v4_area_ratio.mat'));
%longitude,latitude meshgride data                 : lon2,lat2
%observation heat content(time,lat,lon)            : obs_heat2
%model heat content(model,time,lat,lon)            : model_heat2
%model name information(i-3)th                     : root_path2(i).name
%RMSE,bias (model,lat,lon)                         : RMSE,bias
%bias (model,time,lat,lon)                         : bias2


%------STEP2. Calculate RMSE and Trend of NWP

 %---2.1. Determine NWP domain
 lon_idx = [115.5,164.5]; %lon : 115 ~164
 lat_idx = [15.5,52.5];   %lat : 15 ~ 52
 
nwp_idx{1} = find(lon2>=lon_idx(1) & lon2<=lon_idx(2) & lat2>=lat_idx(1) & lat2<=lat_idx(2) );

aa = find(lon2>=lon_idx(1) & lon2<=lon_idx(2) & lat2>=lat_idx(1)-1 & lat2<= lat_idx(1)+1); %east bc
bb = find(lon2>=lon_idx(2)-1 & lon2<=lon_idx(2)+1 & lat2>=lat_idx(1) & lat2<= lat_idx(2)); %south bc
cc = find(lon2>=lon_idx(1)-1 & lon2<=lon_idx(1)+1 & lat2>=lat_idx(1) & lat2<=lat_idx(2) ); %west bc
nwp_idx{2} = sort(unique([aa',bb', cc']));

ocean_num = 2;
cshwa_num = [7,9,12,16,18,20,21,22,25,26,29,34,35,36,38,39,41,42,43,45,47,48,49,50]; %cshwa
cshwa_num = cshwa_num - 3;
 %---2.2. calculate NWP heat time series
 %obs
ocean_obs_heat = zeros(ocean_num,30);
ocean_area = zeros(ocean_num,1);
for i = 1:ocean_num
    idx = nwp_idx{i};
    ocean_area(i) = len*len*nansum(ratio(idx));
    for j = 1:30
        ocean_obs_heat(i,j) = nansum(obs_heat3(j,idx))/ocean_area(i);
    end
end
%mod
ocean_mod_heat = zeros(24,ocean_num,nn);
for k = 1:24
    for i = 1:ocean_num
        for j = 1:30
            idx = nwp_idx{i};
            ocean_mod_heat(k,i,j) = nansum(model_heat3(cshwa_num(k),j,idx))/ocean_area(i); 
        end
    end
end

 %---2.3. calculate NWP RMSE and trend results

%calculate RMSE and bias
ocean_RMSE = zeros(24,ocean_num);
ocean_bias = zeros(24,ocean_num);
for i = 1:24
    for j = 1:ocean_num
        ocean_RMSE(i,j) = sqrt(squeeze(nanmean((squeeze(ocean_mod_heat(i,j,:))-ocean_obs_heat(j,:)').^2)));
        ocean_bias(i,j)    = squeeze(nanmean(squeeze(ocean_mod_heat(i,j,:))-ocean_obs_heat(j,:)')) ;
    end
end

x = [1:30];
%obs_trend
 ocean_trend_obs = zeros(ocean_num,1);
 for i = 1:ocean_num
    data = squeeze(ocean_obs_heat(i,:));
    [a,b] = polyfit(x,data,1);
    ocean_trend_obs(i) = a(1);
 end
 

 %mod trend
 x = [1:30]';
ocean_trend_mod = zeros(24,ocean_num);
for i = 1:24
     for j = 1:ocean_num
         data = squeeze(ocean_mod_heat(i,j,:));
         [a,b] = polyfit(x,data,1);
         ocean_trend_mod(i,j) = a(1);
     end
end


clear nc;
save(strcat(result_path,'/','processing_result_heatcontent_NWP_csh.mat'));

