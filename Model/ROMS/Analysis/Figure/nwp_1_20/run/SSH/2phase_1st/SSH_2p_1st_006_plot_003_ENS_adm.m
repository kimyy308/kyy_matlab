close all; clear all;  clc;
% %  Updated 14-Oct-2021 by Yong-Yub Kim, structure
% %  Updated 06-Dec-2021 by Yong-Yub Kim, structure

% % SSP only

 % %     set dropbox path
tmp.fs=filesep;
tmp.dropboxpath = 'C:\Users\User\Dropbox';
addpath(genpath([tmp.dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'function']));

[tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
    tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
    'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));

% % % configuration of RCM
% RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};

% RCM_info.name={'test2107', 'test2108', 'test2109'};
RCM_info.abbs = {'RCM-CNE', 'RCM-ECV', 'RCM-ACC', 'RCM-CNH', 'RCM-CMC'};
% RCM_info.name={  'test2107'};
% RCM_info.abbs = {  'RCM-CNRM'};

RCM_info.model = 'nwp_1_20';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';
% RCM_info.region = {'NWP', 'AKP4'};
RCM_info.region = {'AKP4'};
% RCM_info.years = 1985:2014;
% RCM_info.years = 1989:2014;
RCM_info.years = 2015:2050;
RCM_info.months = 1:12;
RCM_grid.dl = 1/20;

% % % % % configuration of GCM
GCM_info.name={'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
GCM_info.abbs = {'GCM-CNE', 'GCM-ECV', 'GCM-ACC', 'GCM-CNH', 'GCM-CMC'};
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
RCM_info.matname_mean_data_trend = [RCM_info.savedir,'ENS_adm','_',RCM_info.regionname, '_RCM_ssh_mean_data_trend_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
GCM_info.matname_mean_data_trend = [RCM_info.savedir,'ENS_adm','_',RCM_info.regionname, '_GCM_ssh_mean_data_trend_', ...
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

[RCM_info.scenario, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(RCM_info.name{end});

% % %  AKP4 (domain)
% RCM_fig.tot_trend_rand=[RCM_mean_trend.ens.([RCM_info.regionname,'_rand'])'];
% GCM_fig.tot_trend_rand=[GCM_mean_trend.ens.([RCM_info.regionname,'_rand'])'];

% % adm_div_all
RCM_fig.tot_trend_rand=[RCM_mean_trend.ens.(['adm_div_all','_rand'])'];
GCM_fig.tot_trend_rand=[GCM_mean_trend.ens.(['adm_div_all','_rand'])'];


RCM_info.subregions = {'adm_div_YS', 'adm_div_SS', 'adm_div_ES'};
RCM_info.subregions_korean = {'황해', '남해', '동해'};
for i=1:length(RCM_info.subregions)
    tmp.subregion=RCM_info.subregions{i};
    RCM_fig.tot_trend_rand = [RCM_fig.tot_trend_rand;  RCM_mean_trend.ens.([tmp.subregion,'_rand'])' ];
    GCM_fig.tot_trend_rand = [GCM_fig.tot_trend_rand;  GCM_mean_trend.ens.([tmp.subregion,'_rand'])' ];
end

for ind_sta=1:size(RCM_info.tidal_st,1)
    RCM_fig.tot_trend_rand = [RCM_fig.tot_trend_rand; squeeze(RCM_mean_trend.ens.tidal_station_rand(:,ind_sta)) ];
    GCM_fig.tot_trend_rand = [GCM_fig.tot_trend_rand; squeeze(GCM_mean_trend.ens.tidal_station_rand(:,ind_sta)) ];
end


RCM_fig.trend_rand_name = ...
    repmat({'한반도 주변'},length(RCM_mean_trend.ens.([RCM_info.regionname,'_rand'])),1);
RCM_fig.tot_trend_rand_name = RCM_fig.trend_rand_name;
for i=1:length(RCM_info.subregions)
    tmp.subregion=RCM_info.subregions_korean{i};
    RCM_fig.trend_rand_name = ...
        repmat({tmp.subregion},length(RCM_mean_trend.ens.([RCM_info.regionname,'_rand'])),1);
    RCM_fig.tot_trend_rand_name = [RCM_fig.tot_trend_rand_name; RCM_fig.trend_rand_name];
end

for ind_sta=1:size(RCM_info.tidal_st,1)
    tmp.tidal_station_name = RCM_info.tidal_name_korean{ind_sta};
    RCM_fig.trend_rand_name = ...
        repmat({tmp.tidal_station_name},length(RCM_mean_trend.ens.([RCM_info.regionname,'_rand'])),1);
    RCM_fig.tot_trend_rand_name = [RCM_fig.tot_trend_rand_name; RCM_fig.trend_rand_name];
end

q3=norminv(0.75);
q95=norminv(0.975);
w95=(q95-q3)/(2*q3);

boxplot(RCM_fig.tot_trend_rand, RCM_fig.tot_trend_rand_name, 'BoxStyle', 'outline', 'Colors', 'krrrbbbbbbbbbbbbbbbbbbbbb', 'ColorGroup', RCM_fig.tot_trend_rand_name, ...
    'Whisker', w95, 'symbol','')
% plotlim=[-4 10];
plotlim=[0 12];
ylim(plotlim)

set(gca, 'fontsize', 20)
ylabel('Trend (mm/yr)')
xlabel('Region (station)');
set(gcf,'PaperPosition', [0 0 30 15]) 
hold off
tifname= ['D:\MEPL\project\SSH\6th_year\figure\nwp_1_20\all\station\RCM_station_trend_ENS_adm_hist_',...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.tif'];
saveas(gcf,tifname,'tif'); RemoveWhiteSpace([], 'file', tifname);
close all;

boxplot(GCM_fig.tot_trend_rand, RCM_fig.tot_trend_rand_name, 'BoxStyle', 'outline', 'Colors', 'krrrbbbbbbbbbbbbbbbbbbbbb', 'ColorGroup', RCM_fig.tot_trend_rand_name, ...
    'Whisker', w95, 'symbol','')
ylim(plotlim)

set(gca, 'fontsize', 20)
ylabel('Trend (mm/yr)')
xlabel('Region (station)');
set(gcf,'PaperPosition', [0 0 30 15]) 
hold off
tifname= ['D:\MEPL\project\SSH\6th_year\figure\nwp_1_20\all\station\GCM_station_trend_ENS_adm_hist_',...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.tif'];
saveas(gcf,tifname,'tif'); RemoveWhiteSpace([], 'file', tifname);
close all;


% % % h=boxplot(RCM_fig.tot_trend_rand, RCM_fig.tot_trend_rand_name, 'BoxStyle', 'outline', 'Colors', 'krrrbbbbbbbbbbbbbbbbbbbbb', 'ColorGroup', RCM_fig.tot_trend_rand_name, ...
% % %     'Whisker', w95, 'symbol','')
% % % set(h,{'linew'},{2})
% % % set(gca, 'fontsize', 20)
% % % ylabel('Trend (mm/yr)')
% % % xlabel('Region (station)');
% % % set(gcf,'PaperPosition', [0 0 100 15]) 
% % % tifname= ['Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\all\station\RCM_station_trend_4_1_',...RCM_fig
% % %     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.tif'];
% % % ylim([4 9])
% % % saveas(gcf,tifname,'tif'); RemoveWhiteSpace([], 'file', tifname);
% close all;
% 
% quantile(RCM_mean_trend.ens.(['AKP4_rand']), 0.975)
% quantile(RCM_mean_trend.ens.(['AKP4_rand']), 0.025)
% quantile(GCM_mean_trend.ens.(['AKP4_rand']), 0.975)
% quantile(GCM_mean_trend.ens.(['AKP4_rand']), 0.025)

