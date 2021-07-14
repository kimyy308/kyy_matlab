clc;clear all;close all;
%time range : 1971.01 ~ 2005.12
%area : global

%Setting variable and path
result_path = '/home/ygkim/ipcc/data/zos'; %result path
result_path2= '/home/ygkim/ipcc/data/zos/new'; %save data in this path
root_path ='/data/ygkim/ipcc/'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 
time = 35; %35 year


%---- load mask file
load(strcat('basin_mask.mat'),'basin_mask');
idx_nan = zeros(30,180,360);
for i =1:30
    idx_nan(i,:) = basin_mask(:);
end
idx_nan = find(isnan(idx_nan) ==1);

%--------STEP1. Load obs climatology data-----%
%zos_obs file
load(strcat(result_path,'/','recon_ssha_1976_2005_yearly.mat'));

%make anomaly data
recon_mean = nanmean(recon_ssh);
obs_zos = zeros(size(recon_ssh));
for i = 1: 30
    obs_zos(i,:) = recon_ssh(i,:) - recon_mean(1,:);
end

%mask process
obs_zos(idx_nan) = NaN;

%plot_obs
for i = 1:30
mean_obs(i) = nanmean(nanmean(obs_zos(i,:,:)));
end

figure
plot(mean_obs)
title(strcat('NOAA recon. SSHA'))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')

load(strcat(result_path2,'/','processing_result.mat'));
nanmax(nanmax(nanmax(nanmax(model_zos))))
nanmin(nanmin(nanmin(nanmin(model_zos))))

x = [1:30]';
ocean_trend_obs = NaN(180,360);
ocean_trend_obs2 = NaN(180,360);
 for i = 1:180
    for j= 1:360
    data = squeeze(obs_zos(:,i,j));
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
    data = squeeze(model_zos(model,:,i,j));
    a = polyfit(x,data,1);
    ocean_trend_model(model,i,j) = a(1);
    ocean_trend_model2(model,i,j) = a(2);
end
end
end

save('compare_result_zos_fix.mat');

%% summation for ensemble
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

csh_name = {'CanESM2','NorESM1-M','MPI-ESM-MR','CNRM-CM5','bcc-csm1-1-m','FGOALS-g2','GFDL-ESM2G','MRI-CGCM3','CSIRO-Mk3-6-0','IPSL-CM5A-MR','MIROC5','MIROC-ESM-CHEM','MPI-ESM-LR','IPSL-CM5A-LR'};
for i = 1:14
 srt_index(i) = find(strcmp(csh_name(i),model_name) == 1);
end
% csh_index = sort(srt_index);

for i = 1:14
if i==1
    temp(i,:,:,:) = squeeze(model_zos(srt_index(1),:,:,:));
else
tem = squeeze(temp(i-1,:,:,:)) + squeeze(model_zos(srt_index(i),:,:,:)); %ensemble mean for total performance order
temp(i,:,:,:) = tem ./ i;
i
end
end
ens_model_zos = temp;

x = [1:30]; %time index
ens_trend_model = NaN(14,180,360);
ens_trend_model2 = NaN(14,180,360);
for model = 1:14
 for i = 1:180
    for j= 1:360
    data = squeeze(ens_model_zos(model,:,i,j));   
    a = polyfit(x,data,1);
    ens_trend_model(model,i,j) = a(1);
    ens_trend_model2(model,i,j) = a(2);
end
end
end
save('compare_result_zos_fix_1.mat');

nanmax(nanmax(nanmax(nanmax(ens_model_zos))))
nanmin(nanmin(nanmin(nanmin(ens_model_zos))))

nanmax(nanmax(nanmax(ens_trend_model)))
nanmin(nanmin(nanmin(ens_trend_model)))

nanmax(nanmax(nanmax(ocean_trend_model)))
nanmin(nanmin(nanmin(ocean_trend_model)))

nanmax(nanmax(nanmax(ocean_trend_obs)))
nanmin(nanmin(nanmin(ocean_trend_obs)))

for i = 1:14
figure
pcolor(lon_ref2,lat_ref2,squeeze(ocean_trend_model(srt_index(i),:,:)))
shading flat
colorbar
caxis([-0.004 0.004]);
title(model_name(srt_index(i)))
set(gca,'fontsize',20,'fontweight','bold');
xlim([115.5 164.5])
ylim([15.5 52.5])
end

figure
pcolor(lon_obs,lat_obs,ocean_trend_obs(:,:))
shading flat
colorbar
caxis([-0.004 0.004]);
title('NOAA-REC.SSHA')
set(gca,'fontsize',20,'fontweight','bold');
xlim([115.5 164.5])
ylim([15.5 52.5 ])

for i = 1:14
figure
pcolor(lon_ref2,lat_ref2,squeeze(ens_trend_model(i,:,:)))
shading flat
colorbar
caxis([-0.004 0.004]);
title(strcat('Ensemble - ',num2str(i),' model'))
set(gca,'fontsize',20,'fontweight','bold');
xlim([115.5 164.5])
ylim([15.5 52.5])
saveas(gca,strcat('Ensemble_',num2str(i),'_model.png'))
end

%%%%%%%%%%%%%%%%%% estimate zos only order 

csh_name1 = {'MPI-ESM-LR','NorESM1-M','IPSL-CM5A-LR','CanESM2','FGOALS-g2','CSIRO-Mk3-6-0','IPSL-CM5A-MR','MRI-CGCM3','MIROC5','bcc-csm1-1-m','MPI-ESM-MR','MIROC-ESM-CHEM','GFDL-ESM2G','CNRM-CM5'};
for i = 1:14
 srt_index1(i) = find(strcmp(csh_name1(i),model_name) == 1);
end

for i = 1:14
if i==1
    temp(i,:,:,:) = squeeze(model_zos(srt_index1(1),:,:,:));
else
tem = squeeze(temp(i-1,:,:,:)) + squeeze(model_zos(srt_index1(i),:,:,:)); %ensemble mean for total performance order
temp(i,:,:,:) = tem ./ i;
i
end
end
ens_model_zos1 = temp;

x = [1:30]; %time index
ens_trend_model1 = NaN(14,180,360);
ens_trend_model12 = NaN(14,180,360);
for model = 1:14
 for i = 1:180
    for j= 1:360
    data = squeeze(ens_model_zos1(model,:,i,j));   
    a = polyfit(x,data,1);
    ens_trend_model1(model,i,j) = a(1);
    ens_trend_model12(model,i,j) = a(2);
end
end
end
save('compare_result_zos_fix_2.mat');

%% 4 non-boussinesq model
load 'compare_result_zos_fix_2.mat'
% csh_name1 = {'NorESM1-M','NorESM1-ME','GISS-E2-R','GISS-E2-R-CC'};
csh_name_3 = {'NorESM1-M','GISS-E2-R','NorESM1-ME'}; %% GISS-E2-R-CC get only RCP 4.5, 8.5 
for i = 1:3
 srt_index_3(i) = find(strcmp(csh_name_3(i),model_name) == 1);
end

clearvars temp tem
for i = 1:3
if i==1
    temp(i,:,:,:) = squeeze(model_zos(srt_index_3(1),:,:,:));
else
tem = squeeze(temp(i-1,:,:,:)) + squeeze(model_zos(srt_index_3(i),:,:,:)); %ensemble mean for total performance order
temp(i,:,:,:) = tem ./ i;
i
end
end
ens_model_zos3 = temp;

nanmax(nanmax(nanmax()))
nanmin(nanmin(nanmin())

x = [1:30]; %time index
ens_trend_model3 = NaN(3,180,360);
ens_trend_model32 = NaN(3,180,360);
for model = 1:3
 for i = 1:180
    for j= 1:360
    data = squeeze(ens_model_zos3(model,:,i,j));   
    a = polyfit(x,data,1);
    ens_trend_model3(model,i,j) = a(1);
    ens_trend_model32(model,i,j) = a(2);
end
end
end

for i = 1:3
figure
pcolor(lon_ref2,lat_ref2,squeeze(ens_trend_model3(i,:,:)))
shading flat
colorbar
caxis([-0.004 0.004]);
title(strcat('Ensemble v3- ',num2str(i),' model'))
set(gca,'fontsize',20,'fontweight','bold');
xlim([115.5 164.5])
ylim([15.5 52.5])
saveas(gca,strcat('Ensemble_v3_',num2str(i),'_model.png'))
end


%%% not ensemble
x = [1:30]; %time index
ens_trend_model3 = NaN(3,180,360);
ens_trend_model32 = NaN(3,180,360);
for model = 1:3
 for i = 1:180
    for j= 1:360
    data = squeeze(model_zos(srt_index_3(model),:,i,j));   
    a = polyfit(x,data,1);
    ens_trend_model31(model,i,j) = a(1);
    ens_trend_model312(model,i,j) = a(2);
end
end
end

for i = 1:3
figure
pcolor(lon_ref2,lat_ref2,squeeze(ens_trend_model31(i,:,:)))
shading flat
colorbar
caxis([-0.004 0.004]);
title(model_name(srt_index_3(i)));
set(gca,'fontsize',20,'fontweight','bold');
xlim([115.5 164.5])
ylim([15.5 52.5])
saveas(gca,strcat('mpd_v3_',num2str(i),'_model.png'))
end

% ssha time serise (for 3 model) 
for i = 1:3
    temp= squeeze(nanmean(nanmean(model_zos(srt_index_3(i),:,:,:),4),3));
    mod_zos_time_3(i,:) = temp;
end

for i = 1:3
figure
plot(mod_zos_time_3(i,:))
title(model_name(srt_index_3(i)))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.04 0.04])
saveas(gca,strcat('time_mod_v3_',num2str(i),'.png'))
end

%%%%%%%%%%%%%%%%%%
% ssha time serise 
for i = 1:14
    temp= squeeze(nanmean(nanmean(model_zos(srt_index(i),:,:,:),4),3));
    mod_zos_time(i,:) = temp;
end

for i = 1:14
figure
plot(mod_zos_time(i,:))
title(model_name(srt_index(i)))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.04 0.04])
saveas(gca,strcat('time_mod_',num2str(i),'.png'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
csh_name24 = {'EC-EARTH','CanESM2', 'CESM1-WACCM', 'NorESM1-M','GISS-E2-R','CESM1-CAM5',...
    'MPI-ESM-MR','CNRM-CM5','bcc-csm1-1-m','FGOALS-g2','CCSM4','GFDL-ESM2G','MRI-CGCM3',...
    'CSIRO-Mk3-6-0','NorESM1-ME','MIROC-ESM','IPSL-CM5A-MR','HadGEM2-ES','MIROC5','GFDL-ESM2M','MIROC-ESM-CHEM',...
    'MPI-ESM-LR','bcc-csm1-1','IPSL-CM5A-LR'};

for i = 1:24
 srt_index_24(i) = find(strcmp(csh_name24(i),model_name) == 1);
end

clearvars temp tem
for i = 1:24
if i==1
    temp(i,:,:,:) = squeeze(model_zos(srt_index_24(1),:,:,:));
else
tem = squeeze(temp(i-1,:,:,:)) + squeeze(model_zos(srt_index_24(i),:,:,:)); %ensemble mean for total performance order
temp(i,:,:,:) = tem ./ i;
i
end
end
ens_model_zos24 = temp;

nanmax(nanmax(nanmax()))
nanmin(nanmin(nanmin())

x = [1:30]; %time index
ens_trend_model_24 = NaN(24,180,360);
ens_trend_model_242 = NaN(24,180,360);
for model = 1:24
 for i = 1:180
    for j= 1:360
    data = squeeze(ens_model_zos24(model,:,i,j));   
    a = polyfit(x,data,1);
    ens_trend_model_24(model,i,j) = a(1);
    ens_trend_model_242(model,i,j) = a(2);
end
end
end

for i = 1:24
figure
pcolor(lon_ref2,lat_ref2,squeeze(ens_trend_model_24(i,:,:)))
shading flat
colorbar
caxis([-0.004 0.004]);
title(strcat('Ensemble v3- ',num2str(i),' model'))
set(gca,'fontsize',20,'fontweight','bold');
xlim([115.5 164.5])
ylim([15.5 52.5])
saveas(gca,strcat('Ensemble_v24_',num2str(i),'_model.png'))
end

% ssha time serise 
for i = 1:24
    temp24= squeeze(nanmean(nanmean(model_zos(srt_index_24(i),:,:,:),4),3));
    mod_zos_time24(i,:) = temp24;
end

for i = 1:24
figure
plot(mod_zos_time24(i,:))
title(model_name(srt_index_24(i)))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.04 0.04])
saveas(gca,strcat('time_mod_v24_',num2str(i),'.png'))
end

for i = 1:41
    temp41= squeeze(nanmean(nanmean(model_zos(i,:,:,:),4),3));
    mod_zos_time41(i,:) = temp41;
end
for i = 1:41
figure
plot(mod_zos_time41(i,:))
title(model_name(i))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.04 0.04])
saveas(gca,strcat('time_mod_v41_',num2str(i),'.png'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%NWP
load(strcat(result_path2,'/','processing_result_v4_area_ratio.mat'),'ratio','len');

%---2.1. Determine NWP domain
lon_idx = [115.5,164.5]; %lon : 115 ~164
lat_idx = [15.5,52.5];   %lat : 15 ~ 52
 
nwp_idx{1} = find(lon_ref2>=lon_idx(1) & lon_ref2<=lon_idx(2) & lat_ref2>=lat_idx(1) & lat_ref2<=lat_idx(2) );

aa = find(lon_ref2>=lon_idx(1) & lon_ref2<=lon_idx(2) & lat_ref2>=lat_idx(1)-1 & lat_ref2<= lat_idx(1)+1); %east bc
bb = find(lon_ref2>=lon_idx(2)-1 & lon_ref2<=lon_idx(2)+1 & lat_ref2>=lat_idx(1) & lat_ref2<= lat_idx(2)); %south bc
cc = find(lon_ref2>=lon_idx(1)-1 & lon_ref2<=lon_idx(1)+1 & lat_ref2>=lat_idx(1) & lat_ref2<=lat_idx(2) ); %west bc
nwp_idx{2} = sort(unique([aa',bb', cc']));

nwp_idx{2}

%---2.2. calculate NWP zos time series
%obs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% consider real length on globe
% % len = 114.64*1000 ; %length of 1 degree %%
% % ratio = cos(lat2*pi/180); %area ratio  %%%
[mm,nn,qq,pp] = size(model_zos);
marine_num = 2;
marine_obs_zos = zeros(marine_num,30);
marine_area = zeros(marine_num,1);
for i = 1:marine_num
    idx = nwp_idx{i};
    marine_area(i) = len*len*nansum(ratio(idx));
    for j = 1:30
        marine_obs_zos(i,j) = nansum(obs_zos(j,idx))/marine_area(i);
    end
end
%mod
marine_mod_zos = zeros(mm,marine_num,nn);
for k = 1:mm
    for i = 1:marine_num
        for j = 1:30
            idx = nwp_idx{i};
            marine_mod_zos(k,i,j) = nansum(model_zos(k,j,idx))/marine_area(i); 
        end
    end
end

nanmax(nanmax(nanmax(marine_mod_zos)))
nanmax(nanmax(nanmax(marine_obs_zos)))
nanmin(nanmin(nanmin(marine_mod_zos)))
nanmin(nanmin(nanmin(marine_obs_zos)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just consider it on grid  

[mm,nn,qq,pp] = size(model_zos);
marine_num = 2;
marine_obs_zos = zeros(marine_num,30);
marine_area = zeros(marine_num,1);
for i = 1:marine_num
    idx = nwp_idx{i};
    marine_area(i) = len*len*nansum(ratio(idx));
    for j = 1:30
        marine_obs_zos(i,j) = nansum(obs_zos(j,idx))/length(idx);
    end
end
%mod
marine_mod_zos = zeros(mm,marine_num,nn);
for k = 1:mm
    for i = 1:marine_num
        for j = 1:30
            idx = nwp_idx{i};
            marine_mod_zos(k,i,j) = nansum(model_zos(k,j,idx))/length(idx); 
        end
    end
end

nanmax(nanmax(nanmax(marine_mod_zos)))
nanmax(nanmax(nanmax(marine_obs_zos)))
nanmin(nanmin(nanmin(marine_mod_zos)))
nanmin(nanmin(nanmin(marine_obs_zos)))


figure
plot(squeeze(marine_obs_zos(1,:))) %% 1 is NWP, 2 is BC
title(strcat('NOAA recon. SSHA'))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
saveas(gca,strcat('time_obs_NWP.png'))
close all;

figure
plot(squeeze(marine_obs_zos(2,:))) %% 1 is NWP, 2 is BC
title(strcat('NOAA recon. SSHA'))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
saveas(gca,strcat('time_obs_BC.png'))
close all;

for i = 1:41
figure
plot(squeeze(marine_mod_zos(i,2,:))) %% 1 is NWP, 2 is BC
title(model_name(i))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
saveas(gca,strcat('time_mod_BC_',num2str(i),'.png'))
end
close all;

for i = 1:41
figure
plot(squeeze(marine_mod_zos(i,1,:))) %% 1 is NWP, 2 is BC
title(model_name(i))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
saveas(gca,strcat('time_mod_NWP_',num2str(i),'.png'))
end
close all;


%%%%%% 24 same with upper plot

for i = 1:24
figure
plot(squeeze(marine_mod_zos(srt_index_24(i),2,:))) %% 1 is NWP, 2 is BC
title(model_name(srt_index_24(i)))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
saveas(gca,strcat('time_mod_BC_',num2str(i),'.png'))
end
close all;

for i = 1:24
figure
plot(squeeze(marine_mod_zos(srt_index_24(i),1,:))) %% 1 is NWP, 2 is BC
title(model_name(srt_index_24(i)))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
saveas(gca,strcat('time_mod_NWP_',num2str(i),'.png'))
end
close all;

% 3 model
for i = 1:3
figure
plot(squeeze(marine_mod_zos(srt_index_3(i),2,:))) %% 1 is NWP, 2 is BC
title(model_name(srt_index_3(i)))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
saveas(gca,strcat('3_time_mod_BC_',num2str(i),'.png'))
end
close all;

for i = 1:3
figure
plot(squeeze(marine_mod_zos(srt_index_3(i),1,:))) %% 1 is NWP, 2 is BC
title(model_name(srt_index_3(i)))
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
saveas(gca,strcat('3_time_mod_NWP_',num2str(i),'.png'))
end
close all;

%NWP-24
plot(squeeze(marine_obs_zos(1,:)),char(GetRandomLineStyleForPlot))
hold on;
for i = 1:24
plot(squeeze(marine_mod_zos(srt_index_24(i),1,:)),char(GetRandomLineStyleForPlot))
hold on;
end
title('NWP SSHA Compare')
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05])
set(gca,'ytick',[-0.05:0.01:0.05])
name_fig(1,2:25)=model_name(srt_index_24);
name_fig(1,1)= cellstr('NOAA-OBS');
legend(name_fig,'Fontsize',5)

%BC-24
plot(squeeze(marine_obs_zos(2,:)),char(GetRandomLineStyleForPlot))
hold on;
for i = 1:24
plot(squeeze(marine_mod_zos(srt_index_24(i),2,:)),char(GetRandomLineStyleForPlot))
hold on;
end
title('BC SSHA Compare')
set(gca,'fontsize',20,'fontweight','bold');
xlabel('year')
ylabel('meter')
ylim([-0.05 0.05]) 
set(gca,'ytick',[-0.05:0.01:0.05])
name_fig(1,2:25)=model_name(srt_index_24);
name_fig(1,1)= cellstr('NOAA-OBS');
legend(name_fig,'Fontsize',5)

