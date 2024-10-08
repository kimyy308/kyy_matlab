clc;clear all;close all;
%time range : 1971.01 ~ 2005.12
%area : global

result_path = '/home/ygkim/ipcc/data/heatcontent'; %result path
result_path2 = '/home/ygkim/ipcc/data/heatcontent/rmse_bias_new'; %result path
root_path ='/data/ygkim/ipcc'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 


%--------STEP1. Load result data-----%
%processing result file
load(strcat(result_path,'/processing_result_heatcontent_NWP.mat'));
%longitude,latitude meshgride data                 : lon2,lat2
%basin parameter(NA:1 SA:2 NP:3 SP:4 no:5 SO:6 AO:7: new_basin,ocean_num
%observation heat content(time,lat,lon)            : obs_heat2
%model heat content(model,time,lat,lon)            : model_heat2
%model name information(i-3)th                     : root_path2(i).name
%RMSE,bias (model,lat,lon)                         : RMSE,bias
%total area mean of RMSE,bias             

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

%select models (SSH exists)
% ocean_mod_heat = ocean_mod_heat(zos_idx-3,:,:);
% 

%------STEP2. ensemble mean process ------%
ens_ocean_mod_heat2 = zeros(24,2,30);
ens_ocean_RMSE2 = zeros(24,2);
ens_ocean_bias2 = zeros(24,2);
ens_ocean_trend_mod2 = zeros(24,2);
% ens_ocean_mod_heat2 = zeros(length(sort_num),2,30);
% ens_ocean_RMSE2 = zeros(length(sort_num),2);
% ens_ocean_bias2 = zeros(length(sort_num),2);         :total_RMSE,total_bias
%total area mean of RMSE's std                     :total_std
%bias (model,time,lat,lon),mean bias               :bias2,total_bias2
%each ocean mean of RMSE,bias                      :ocean_RMSE,ocean_bias
%each ocean mean of RMSE's std                     :ocean_std
%total trend of obs,mod                   :total_trend_obs,total_trend_mod
%ocean trend of obs,mod                   :ocean_trend_obs,ocean_trend_mod
%each grid trend of obs,mod               :each_trend_obs,each_trend_mod
%std of RMSE model performance(model)              :se,se_ocean

%---Load Model performance result ---%
%load(strcat('/home/ygkim/ipcc/script/NWP/','result_model_performance_add_SSH.mat'));
% cshwa_num = [7,9,12,16,18,20,21,22,25,26,29,34,35,36,38,39,41,42,43,45,47,48,49,50]; %cshwa
cshwa_num = [2,3,5,8,11,12,16,17,18,19,23,24,26,27,28,29,30,31,34,37,38,39,40,45]; % cshwa from model_name(sort_num)
% I : model performance
%[sort_idx, sort_num] = sort(total_RMSE_new) ;   %set sort parameter
% ens_ocean_trend_mod2 = zeros(length(sort_num),2);
[sort_idx, sort_num] = sort(total_RMSE_new) ;   %set sort parameter
ens_ocean_trend_mod2 = zeros(length(cshwa_num),2);

%%%%%%%%% cshwa 24 model index pick up
k=0
% sort_cal = sort_num +3;
% cshwa_num = cshwa_num -3 ;
for j = 1:length(cshwa_num)
    for i = 1:length(sort_num)
        if sort_num(i) == cshwa_num(j);
            k = k+1;
            csh_index(k) = i;
        end
    end
end
for i = 4:length(root_path2)
    model_name(i-3) = {root_path2(i).name};
end

csh_name = {'CCSM4', 'CESM1-CAM5', 'CSIRO-Mk3-6-0', 'GFDL-ESM2G', 'GFDL-ESM2M','HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'MIROC-ESM', 'MIROC5', 'MIROC-ESM-CHEM','NorESM1-M', 'NorESM1-ME', 'bcc-csm1-1', 'bcc-csm1-1-m', 'CESM1-WACCM', 'CNRM-CM5','CanESM2', 'EC-EARTH', 'FGOALS-g2', 'GISS-E2-R', 'MPI-ESM-LR', 'MPI-ESM-MR', 'MRI-CGCM3'};
for i = 1:24
 srt_index(i) = find(strcmp(csh_name(i),model_name(sort_num)) == 1);
end
sort(srt_index);
%model_name(end+1) = {'Ensemble'};
% model_name = model_name(zos_idx-3);


for k = 1:length(cshwa_num) %kyy
%     [m,n] = size(file_zos);
%     if m~=0
%         zos_num = zos_num +1;
%         zos_idx(end+1) = i;
%     end
% sort_num : sorted model number
% sort_idx : sorted model result

%------STEP1.4 Select models(included all component)
%find zos file
    %ensemble mean process
    nn = k; %ensemble number
    fin_idx = sort_num(cshwa_num(1:nn));

    ens_ocean_mod_heat = squeeze(mean(ocean_mod_heat(fin_idx,:,:),1));
    %----3.2. calc RMSE----%
     %calculate RMSE and bias

    for i = 1:2
        ens_ocean_RMSE(i) = sqrt(squeeze(nanmean((squeeze(ens_ocean_mod_heat(i,:))-ocean_obs_heat(i,:)).^2)));
        ens_ocean_bias(i)    = squeeze(nanmean(squeeze(ens_ocean_mod_heat(i,:))-ocean_obs_heat(i,:))) ;
    end


    %----3.3. calc Trend----%


     %mod trend
     x = [1:30];
    ens_ocean_trend_mod = zeros(1,2);
         for j = 1:2
             data = squeeze(ens_ocean_mod_heat(j,:));
             [a,b] = polyfit(x,data,1);
             ens_ocean_trend_mod(j) = a(1);
         end
    
    %save each ens results
    ens_ocean_mod_heat2(k,:) = ens_ocean_mod_heat(:);
    ens_ocean_RMSE2(k,:) = ens_ocean_RMSE(:);
    ens_ocean_bias2(k,:) = ens_ocean_bias(:);
    ens_ocean_trend_mod2(k,:) = ens_ocean_trend_mod(:);
end


%save data
save(strcat(result_path,'/','processing_result_v2_heatcontent_NWP_csh.mat'));
