close all; clear all;  clc;
% %  Updated 16-Oct-2021 by Yong-Yub Kim,  historical + ssp585

% % % configuration of RCM
% RCM_info.name={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
% RCM_info.name={'test2107', 'test2108', 'test2109'};
% RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};

% RCM_info.name={  'test2107'};
% RCM_info.abbs = {  'RCM-CNRM'};

RCM_info.model = 'nwp_1_20';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';
RCM_info.region = {'AKP4'};
RCM_info.years = 2015:2050;
RCM_info.months = 1:12;
RCM_grid.dl = 1/20;

% % % % % configuration of GCM
GCM_info.name={'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
GCM_info.abbs = {'GCM-CNE', 'GCM-ECV', 'GCM-ACC', 'GCM-CMH', 'GCM-CMC'};
% GCM_info.name={'CNRM-ESM2-1'};
% GCM_info.abbs = {  'GCM-CNRM'};
GCM_info.model = GCM_info.name;
GCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'CMIP6', filesep, 'NWP', filesep];
GCM_info.saveroot = GCM_info.dataroot;
GCM_info.phase = RCM_info.phase;
GCM_info.region = RCM_info.region;
GCM_info.years = RCM_info.years;
GCM_info.months = RCM_info.months;
GCM_grid.dl = 1/2;

% % % configuration of CMEMS
CMEMS_info.filedir = 'D:\Data\Observation\CMEMS\';
CMEMS_grid.dl = 1/4;

% % % configuration of flags (make file)
flags.make_std_mat = 1;
flags.make_std_nc = 1;

tmp.fs=filesep;
tmp.variable = 'zeta';
regionind=1;
RCM_info.regionname = RCM_info.region{regionind};

RCM_info.savedir = [RCM_info.dataroot, RCM_info.name{end}, tmp.fs, 'run', tmp.fs, ...
            tmp.variable, tmp.fs];
RCM_info.matname_mean_data_trend = [RCM_info.savedir,'ENSg','_',RCM_info.regionname, '_RCM_ssh_mean_data_trend_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
GCM_info.matname_mean_data_trend = [RCM_info.savedir,'ENSg','_',RCM_info.regionname, '_GCM_ssh_mean_data_trend_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];


RCM_info.tidal_st= [ 126.5919, 37.4528; ...  % Incheon (Yellow Sea)
                126.1322, 36.6436; ...  % Anheung
                126.4858, 36.4065; ...  % Boryung
                126.5431, 35.9756; ...  % Gunsan
                126.3035, 35.6177; ...  % Wuido
                126.3756, 34.7297; ...  % Mokpo
                125.4356, 34.6842; ...  % Heuksando
                126.3003, 33.9619; ...  % Chujado (South Sea)
                126.7597, 34.2756; ...  % Wando
                127.8056, 34.7472; ...  % Yeosu
                128.4547, 34.8278; ...  % Tongyeong
                128.8008, 35.0242; ...  % Gadeokdo
                129.0453, 35.0564; ...  % Busan
                126.5431, 33.5675; ...  % Jeju
                126.5617, 33.2400; ...  % Seoguipo
                127.3089, 33.9983; ...  % Geomundo
                129.3872, 35.4719; ...  % Ulsan (East Sea)
                129.3839, 36.0372; ...  % Pohang
                129.1164, 37.6003; ...  % Mukho
                128.5942, 38.2072; ...  % Sokcho
                130.9136, 37.5614];  % Ulleungdo

RCM_info.tidal_name = { 'West-01-Incheon', 'West-02-Anheung', 'West-03-Boryung', 'West-04-Gunsan', 'West-05-Wuido', 'West-06-Mokpo', 'West-07-Heuksando', ...
               'South-08-Chujado', 'South-09-Wando', 'South-10-Yeosu', 'South-11-Tongyeong', 'South-12-Gadeokdo',  ...
               'South-13-Busan', 'South-14-Jeju', 'South-15-Seoguipo', 'South-16-Geomundo' ...
               'East-17-Ulsan', 'East-18-Pohang', 'East-19-Mukho', 'East-20-Sokcho', 'East-21-Ulleungdo'};
RCM_info.tidal_name_korean = { '서해안-01-인천', '서해안-02-안흥', '서해안-03-보령', '서해안-04-군산', '서해안-05-위도', '서해안-06-목포', '서해안-07-흑산도', ...
                '남해안-08-추자도', '남해안-09-완도', '남해안-10-여수', '남해안-11-통영', '남해안-12-가덕도', ...
                '남해안-13-부산', '남해안-14-제주', '남해안-15-서귀포', '남해안-16-거문도', ...
                '동해안-17-울산', '동해안-18-포항', '동해안-19-묵호', '동해안-20-속초', '동해안-21-울릉도'};
load(RCM_info.matname_mean_data_trend);
load(GCM_info.matname_mean_data_trend);


RCM_info_ssp585=RCM_info;
RCM_mean_data_ssp585=RCM_mean_data;
RCM_mean_trend_ssp585=RCM_mean_trend;
GCM_info_ssp585=GCM_info;
GCM_mean_data_ssp585=GCM_mean_data;
GCM_mean_trend_ssp585=GCM_mean_trend;


% RCM_info.name={'test2102', 'test2103', 'test2104'};
RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};

RCM_info.years = 1985:2014;
RCM_info.savedir = [RCM_info.dataroot, RCM_info.name{end}, tmp.fs, 'run', tmp.fs, ...
            tmp.variable, tmp.fs];
RCM_info.matname_mean_data_trend = [RCM_info.savedir,'ENSg','_',RCM_info.regionname, '_RCM_ssh_mean_data_trend_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
GCM_info.matname_mean_data_trend = [RCM_info.savedir,'ENSg','_',RCM_info.regionname, '_GCM_ssh_mean_data_trend_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
load(RCM_info.matname_mean_data_trend);
load(GCM_info.matname_mean_data_trend);

RCM_info_historical=RCM_info;
RCM_mean_data_historical=RCM_mean_data;
RCM_mean_trend_historical=RCM_mean_trend;
GCM_info_historical=GCM_info;
GCM_mean_data_historical=GCM_mean_data;
GCM_mean_trend_historical=GCM_mean_trend;

RCM_mean_data_all.ens.(RCM_info.regionname)(1:length(RCM_info_historical.years)) = ...
    RCM_mean_data_historical.ens.(RCM_info.regionname);
RCM_mean_data_all.ens.(RCM_info.regionname)(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
    RCM_mean_data_ssp585.ens.(RCM_info.regionname);
GCM_mean_data_all.ens.(RCM_info.regionname)(1:length(RCM_info_historical.years)) = ...
    GCM_mean_data_historical.ens.(RCM_info.regionname);
GCM_mean_data_all.ens.(RCM_info.regionname)(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
    GCM_mean_data_ssp585.ens.(RCM_info.regionname);

RCM_mean_data_all.ens.([RCM_info.regionname, '_95_14'])= mean(RCM_mean_data_all.ens.(RCM_info.regionname)(11:30));
GCM_mean_data_all.ens.([RCM_info.regionname, '_95_14'])= mean(GCM_mean_data_all.ens.(RCM_info.regionname)(11:30));
RCM_mean_data_all.ens.([RCM_info.regionname, '_slr_50'])= RCM_mean_data_all.ens.(RCM_info.regionname)(end) - RCM_mean_data_all.ens.([RCM_info.regionname, '_95_14']);
GCM_mean_data_all.ens.([RCM_info.regionname, '_slr_50'])= GCM_mean_data_all.ens.(RCM_info.regionname)(end) - GCM_mean_data_all.ens.([RCM_info.regionname, '_95_14']);

RCM_info_all.subregions = {'ES_KHOA', 'YS_KHOA', 'SS_KHOA'};
for i=1:length(RCM_info_all.subregions)
    tmp.subregion=RCM_info_all.subregions{i};
    RCM_mean_data_all.ens.(tmp.subregion)(1:length(RCM_info_historical.years)) = ...
        RCM_mean_data_historical.ens.(tmp.subregion);
    RCM_mean_data_all.ens.(tmp.subregion)(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
        RCM_mean_data_ssp585.ens.(tmp.subregion);
    GCM_mean_data_all.ens.(tmp.subregion)(1:length(RCM_info_historical.years)) = ...
        GCM_mean_data_historical.ens.(tmp.subregion);
    GCM_mean_data_all.ens.(tmp.subregion)(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
        GCM_mean_data_ssp585.ens.(tmp.subregion);
    RCM_mean_data_all.ens.([tmp.subregion, '_95_14'])= mean(RCM_mean_data_all.ens.(tmp.subregion)(11:30));
    GCM_mean_data_all.ens.([tmp.subregion, '_95_14'])= mean(GCM_mean_data_all.ens.(tmp.subregion)(11:30));
    RCM_mean_data_all.ens.([tmp.subregion, '_slr_50'])= RCM_mean_data_all.ens.(tmp.subregion)(end) - RCM_mean_data_all.ens.([tmp.subregion, '_95_14']);
    GCM_mean_data_all.ens.([tmp.subregion, '_slr_50'])= GCM_mean_data_all.ens.(tmp.subregion)(end) - GCM_mean_data_all.ens.([tmp.subregion, '_95_14']);
end


RCM_mean_data_all.ens.([RCM_info.regionname, '_rand_lower_limit'])(1:length(RCM_info_historical.years)) = ...
     RCM_mean_data_historical.ens.([RCM_info.regionname, '_rand_lower_limit']);
RCM_mean_data_all.ens.([RCM_info.regionname, '_rand_upper_limit'])(1:length(RCM_info_historical.years)) = ...
     RCM_mean_data_historical.ens.([RCM_info.regionname, '_rand_upper_limit']);
RCM_mean_data_all.ens.([RCM_info.regionname, '_rand_lower_limit'])(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
     RCM_mean_data_ssp585.ens.([RCM_info.regionname, '_rand_lower_limit']);
RCM_mean_data_all.ens.([RCM_info.regionname, '_rand_upper_limit'])(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
     RCM_mean_data_ssp585.ens.([RCM_info.regionname, '_rand_upper_limit']);
 
GCM_mean_data_all.ens.([RCM_info.regionname, '_rand_lower_limit'])(1:length(RCM_info_historical.years)) = ...
     GCM_mean_data_historical.ens.([RCM_info.regionname, '_rand_lower_limit']);
GCM_mean_data_all.ens.([RCM_info.regionname, '_rand_upper_limit'])(1:length(RCM_info_historical.years)) = ...
     GCM_mean_data_historical.ens.([RCM_info.regionname, '_rand_upper_limit']);
GCM_mean_data_all.ens.([RCM_info.regionname, '_rand_lower_limit'])(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
     GCM_mean_data_ssp585.ens.([RCM_info.regionname, '_rand_lower_limit']);
GCM_mean_data_all.ens.([RCM_info.regionname, '_rand_upper_limit'])(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
     GCM_mean_data_ssp585.ens.([RCM_info.regionname, '_rand_upper_limit']); 
 
 for i=1:length(RCM_info_all.subregions)
    tmp.subregion=RCM_info_all.subregions{i};
    RCM_mean_data_all.ens.([tmp.subregion, '_rand_lower_limit'])(1:length(RCM_info_historical.years)) = ...
         RCM_mean_data_historical.ens.([tmp.subregion, '_rand_lower_limit']);
    RCM_mean_data_all.ens.([tmp.subregion, '_rand_upper_limit'])(1:length(RCM_info_historical.years)) = ...
         RCM_mean_data_historical.ens.([tmp.subregion, '_rand_upper_limit']);
    RCM_mean_data_all.ens.([tmp.subregion, '_rand_lower_limit'])(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
         RCM_mean_data_ssp585.ens.([tmp.subregion, '_rand_lower_limit']);
    RCM_mean_data_all.ens.([tmp.subregion, '_rand_upper_limit'])(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
         RCM_mean_data_ssp585.ens.([tmp.subregion, '_rand_upper_limit']);
    GCM_mean_data_all.ens.([tmp.subregion, '_rand_lower_limit'])(1:length(RCM_info_historical.years)) = ...
         GCM_mean_data_historical.ens.([tmp.subregion, '_rand_lower_limit']);
    GCM_mean_data_all.ens.([tmp.subregion, '_rand_upper_limit'])(1:length(RCM_info_historical.years)) = ...
         GCM_mean_data_historical.ens.([tmp.subregion, '_rand_upper_limit']);
    GCM_mean_data_all.ens.([tmp.subregion, '_rand_lower_limit'])(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
         GCM_mean_data_ssp585.ens.([tmp.subregion, '_rand_lower_limit']);
    GCM_mean_data_all.ens.([tmp.subregion, '_rand_upper_limit'])(length(RCM_info_historical.years)+1:length(RCM_info_historical.years)+length(RCM_info_ssp585.years)) = ...
         GCM_mean_data_ssp585.ens.([tmp.subregion, '_rand_upper_limit']); 
 end

 
RCM_info_all.years=min(RCM_info_historical.years) : max(RCM_info_ssp585.years);

RCM_mean_data_all.ens.([RCM_info.regionname, '_slr_upper_50'])= ...
    RCM_mean_data_all.ens.([RCM_info.regionname, '_rand_upper_limit'])(end) - ...
    RCM_mean_data_all.ens.([RCM_info.regionname, '_95_14']);
RCM_mean_data_all.ens.([RCM_info.regionname, '_slr_lower_50'])= ...
    RCM_mean_data_all.ens.([RCM_info.regionname, '_rand_lower_limit'])(end) - ...
    RCM_mean_data_all.ens.([RCM_info.regionname, '_95_14']);
GCM_mean_data_all.ens.([RCM_info.regionname, '_slr_upper_50'])= ...
    GCM_mean_data_all.ens.([RCM_info.regionname, '_rand_upper_limit'])(end) - ...
    GCM_mean_data_all.ens.([RCM_info.regionname, '_95_14']);
GCM_mean_data_all.ens.([RCM_info.regionname, '_slr_lower_50'])= ...
    GCM_mean_data_all.ens.([RCM_info.regionname, '_rand_lower_limit'])(end) - ...
    GCM_mean_data_all.ens.([RCM_info.regionname, '_95_14']);

for i=1:length(RCM_info_all.subregions)
    tmp.subregion=RCM_info_all.subregions{i};
     RCM_mean_data_all.ens.([tmp.subregion, '_slr_upper_50'])= ...
        RCM_mean_data_all.ens.([tmp.subregion, '_rand_upper_limit'])(end) - ...
        RCM_mean_data_all.ens.([tmp.subregion, '_95_14']);
    RCM_mean_data_all.ens.([tmp.subregion, '_slr_lower_50'])= ...
        RCM_mean_data_all.ens.([tmp.subregion, '_rand_lower_limit'])(end) - ...
        RCM_mean_data_all.ens.([tmp.subregion, '_95_14']);
    GCM_mean_data_all.ens.([tmp.subregion, '_slr_upper_50'])= ...
        GCM_mean_data_all.ens.([tmp.subregion, '_rand_upper_limit'])(end) - ...
        GCM_mean_data_all.ens.([tmp.subregion, '_95_14']);
    GCM_mean_data_all.ens.([tmp.subregion, '_slr_lower_50'])= ...
        GCM_mean_data_all.ens.([tmp.subregion, '_rand_lower_limit'])(end) - ...
        GCM_mean_data_all.ens.([tmp.subregion, '_95_14']);
end


val_transparent = 0.2;
figure_ts.RCM_cmap = [0,0,1]; % blue  [1 0 0] -> red
figure_ts.RCM_cmap_range = rgb2hsv(figure_ts.RCM_cmap);
figure_ts.RCM_cmap_range(:,2) =  val_transparent;
figure_ts.RCM_cmap_range = hsv2rgb(figure_ts.RCM_cmap_range);

figure_ts.GCM_cmap = [1,0, 0]; % blue  [1 0 0] -> red
figure_ts.GCM_cmap_range = rgb2hsv(figure_ts.GCM_cmap);
figure_ts.GCM_cmap_range(:,2) =  val_transparent;
figure_ts.GCM_cmap_range = hsv2rgb(figure_ts.GCM_cmap_range);



% % %  RCM time series with confident range
figure_ts.fig_RCM_range=fill([RCM_info_ssp585.years, flip(RCM_info_ssp585.years)], ...
    [RCM_mean_data_ssp585.ens.AKP4_rand_lower_limit, flip(RCM_mean_data_ssp585.ens.AKP4_rand_upper_limit)].*100.0, figure_ts.RCM_cmap_range);
figure_ts.fig_RCM_range.FaceColor = figure_ts.RCM_cmap_range;
figure_ts.fig_RCM_range.EdgeColor = 'none';
hold on
figure_ts.fig_RCM=plot(RCM_info_all.years, RCM_mean_data_all.ens.AKP4.*100.0, 'color', figure_ts.RCM_cmap);
figure_ts.fig_RCM.LineWidth= 2;
hold off
xlabel('Year')
ylabel('Sea-level (cm)')
set(gca,'fontsize',15)
ylim([40 100])
xlim([min(RCM_info_all.years) max(RCM_info_all.years)])
tifname= ['D:\MEPL\project\SSH\6th_year\figure\nwp_1_20\all\ts\RCM_ssh_ts_hist_', 'ssp585', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.tif'];
saveas(gcf,tifname,'tif'); RemoveWhiteSpace([], 'file', tifname);

% % %  GCM time series with confident range
figure_ts.fig_GCM_range=fill([RCM_info_ssp585.years, flip(RCM_info_ssp585.years)], ...
    [GCM_mean_data_ssp585.ens.AKP4_rand_lower_limit, flip(GCM_mean_data_ssp585.ens.AKP4_rand_upper_limit)].*100.0, figure_ts.RCM_cmap_range);
figure_ts.fig_GCM_range.FaceColor = figure_ts.GCM_cmap_range;
figure_ts.fig_GCM_range.EdgeColor = 'none';
hold on
figure_ts.fig_GCM=plot(RCM_info_all.years, GCM_mean_data_all.ens.AKP4.*100.0, 'color', figure_ts.GCM_cmap);
figure_ts.fig_GCM.LineWidth= 2;
hold off
xlabel('Year')
ylabel('Sea-level (cm)')
set(gca,'fontsize',15)
ylim([40 100])
xlim([min(RCM_info_all.years) max(RCM_info_all.years)])
tifname= ['D:\MEPL\project\SSH\6th_year\figure\nwp_1_20\all\ts\GCM_ssh_ts_hist_', 'ssp585', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.tif'];
saveas(gcf,tifname,'tif'); RemoveWhiteSpace([], 'file', tifname);
close all;
