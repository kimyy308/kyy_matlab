close all; clear all;  clc;
% %  Updated 05-Jul-2021 by Yong-Yub Kim, structure

% % % configuration of RCM
% RCM_info.name={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
% RCM_info.name={'test2107', 'test2108', 'test2109', 'test2110', 'test2111'};
RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};

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
RCM_info.years = 2015:2050;
RCM_info.months = 1:12;
RCM_grid.dl = 1/20;

% % % % % configuration of GCM
% GCM_info.name={'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
% GCM_info.abbs = {'GCM-CNRM', 'GCM-EC-Veg', 'GCM-ACC', 'GCM-CNRM-HR', 'GCM-CMCC'};
GCM_info.name={'CNRM-ESM2-1'};
GCM_info.abbs = {  'GCM-CNRM'};
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


% %  working
for testnameind=1:length(RCM_info.name)
    for regionind=1:length(RCM_info.region)
        clearvars '*' -except regionind testnameind flags flg_lev param ...
            RCM_info RCM_grid GCM_info GCM_grid CMEMS_info CMEMS_grid ...
            RCM_mean_data RCM_mean_trend GCM_mean_data GCM_mean_trend
        
        tmp.variable = 'zeta';
        tmp.variable_GCM = 'zos';
        tmp.variable_CMEMS = 'sla';
        tmp.fs = filesep; % file separator win = '\', linux = '/'

        % %     set dropbox path
        tmp.dropboxpath = 'C:\Users\User\Dropbox';
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));

% %     set temporary variables (testname, regionname, filesep, ...)
        RCM_info.testname = RCM_info.name{testnameind};
        RCM_info.regionname = RCM_info.region{regionind};
        RCM_info.abb = RCM_info.abbs{testnameind};
        [RCM_info.scenario, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(RCM_info.testname);
        
        GCM_info.testname = GCM_info.name{testnameind};
        GCM_info.regionname = RCM_info.regionname;
        GCM_info.abb = GCM_info.abbs{testnameind};
        GCM_info.scenario = RCM_info.scenario;
                
        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
            tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
            'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));
        
        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(RCM_info.regionname);

        tmp.param_script = [tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', ...
            tmp.fs, 'Model', tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', ...
            tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, 'fig_param', tmp.fs, ...
            'fig_param2_kyy_', RCM_info.regionname, '.m'];
        RCM_info.filedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs, 'run', tmp.fs, ...
            tmp.variable, tmp.fs];
        RCM_info.savedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs, 'run', tmp.fs, ...
            tmp.variable, tmp.fs];
        GCM_info.filedir = [GCM_info.dataroot, tmp.variable_GCM, tmp.fs, GCM_info.scenario, tmp.fs, ...
            'Omon', tmp.fs, GCM_info.testname, tmp.fs];
        
        RCM_info.matname_mean_data_trend = [RCM_info.savedir,'all','_',RCM_info.regionname, '_RCM_ssh_mean_data_trend_', ...
            num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
        GCM_info.matname_mean_data_trend = [RCM_info.savedir,'all','_',RCM_info.regionname, '_GCM_ssh_mean_data_trend_', ...
            num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
        
% % %         time set
        for folding=1:1
            tmp.tind=1;
            for yearij = 1:length(RCM_info.years)
                for month=1:length(RCM_info.months) 
                    tmp.year = RCM_info.years(yearij);
%                     tmp.month = RCM_info.months(monthij);
                    RCM_time.ftime(tmp.tind) = datenum(tmp.year,month,15) - datenum(1900,12,31);
                    tmp.tind=tmp.tind+1;
                end
            end
            for month=1:length(RCM_info.months)  
                tmp.year = RCM_info.years(yearij);
                RCM_time.climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
            end

            for i =1:length(RCM_info.years) 
                tmp.year = RCM_info.years(yearij);
                for month=1:length(RCM_info.months)  
                    RCM_time.xData((12*(i-1))+month) = datenum([num2str(tmp.year),'-',num2str(month,'%02i'),'-01',]); 
                end
            end

            RCM_time.trendtime=RCM_info.years(1):1/length(RCM_info.months) : RCM_info.years(end)+1-1/length(RCM_info.months) ;
            RCM_time.trendtime_yearly=RCM_info.years(1) : RCM_info.years(end);
        end     
        
        RCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_ssh_', ...
            num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
        RCM_info.matname_interped = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_cmems_interped_ssh_', ...
            num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
        GCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_ssh_', ...
            num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
        GCM_info.matname_interped = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_cmems_interped_ssh_', ...
            num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
        CMEMS_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_CMEMS_ssh_', ...
            num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
        if (exist(RCM_info.matname , 'file') == 2 || flags.fig_switch(1)~=2)   
            load(RCM_info.matname); 
            load(GCM_info.matname);
            load(GCM_info.matname_interped);
            load(CMEMS_info.matname);
        end
        
% % %        reference set : yearly_mean - mean(RCM(95-14)) + mean(CMEMS(95-14));
        RCM_data.yearly_mean = RCM_data.yearly_mean ...
            - Func_0017_SSH_correction_for_CMIP6_RMSE(RCM_info.testname) + Func_0017_SSH_correction_for_CMIP6_RMSE('CMEMS');
        GCM_data.yearly_mean = GCM_data.yearly_mean ...
            - Func_0017_SSH_correction_for_CMIP6_RMSE(GCM_info.testname) + Func_0017_SSH_correction_for_CMIP6_RMSE('CMEMS');
        GCM_data_interped.yearly_mean = GCM_data_interped.yearly_mean ...
            - Func_0017_SSH_correction_for_CMIP6_RMSE(GCM_info.testname) + Func_0017_SSH_correction_for_CMIP6_RMSE('CMEMS');
        
% % %         get mean sea-level (domain, (AKP4))
        RCM_time.tlen=length(RCM_time.trendtime_yearly);
        for i=1:RCM_time.tlen
% % %             ex) RCM_mean_data.test2107.AKP4(i) = mean(RCM_data.yearly_mean(:)) 
            [RCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)(i), tmp.error_status] = ...
                Func_0011_get_area_weighted_mean(RCM_data.yearly_mean(:,:,i), RCM_grid.lon_rho, RCM_grid.lat_rho); 
            [GCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)(i), tmp.error_status] = ...
                Func_0011_get_area_weighted_mean(GCM_data.yearly_mean(:,:,i), GCM_grid.lon, GCM_grid.lat);
            [GCM_mean_data_interped.([RCM_info.testname]).(RCM_info.regionname)(i), tmp.error_status] = ...
                Func_0011_get_area_weighted_mean(GCM_data_interped.yearly_mean(:,:,i), CMEMS_grid.lon2, CMEMS_grid.lat2);
        end
        if testnameind==1
            RCM_mean_data.ens.(RCM_info.regionname)=zeros(1,length(RCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)));
            RCM_mean_data.ens.(RCM_info.regionname)= ...
                RCM_mean_data.ens.(RCM_info.regionname) + RCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)/length(RCM_info.name);
            GCM_mean_data.ens.(RCM_info.regionname)=zeros(1,length(GCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)));            
            GCM_mean_data.ens.(RCM_info.regionname)= ...
                GCM_mean_data.ens.(RCM_info.regionname) + GCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)/length(RCM_info.name);
        else
            RCM_mean_data.ens.(RCM_info.regionname)= ...
                RCM_mean_data.ens.(RCM_info.regionname) + RCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)/length(RCM_info.name);
            GCM_mean_data.ens.(RCM_info.regionname)= ...
                GCM_mean_data.ens.(RCM_info.regionname) + GCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)/length(RCM_info.name);
        end
        
% % %         get trend (AKP4)
        [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
             squeeze(RCM_mean_data.([RCM_info.testname]).(RCM_info.regionname))'.*1000.0, 'poly1'); %% m/y -> mm/y
        RCM_mean_trend.([RCM_info.testname]).(RCM_info.regionname) = tmp_fitobject.p1;
        [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
             squeeze(GCM_mean_data.([RCM_info.testname]).(RCM_info.regionname))'.*1000.0, 'poly1'); %% m/y -> mm/y
        GCM_mean_trend.([RCM_info.testname]).(RCM_info.regionname) = tmp_fitobject.p1;


% % %         get mean sea-level (ES_KHOA, YS_KHOA, SS_KHOA)
        RCM_info.subregions = {'ES_KHOA', 'YS_KHOA', 'SS_KHOA'};
        for i=1:length(RCM_info.subregions)
            tmp.subregion=RCM_info.subregions{i};
            [RCM_grid.(['polygon_', tmp.subregion]), RCM_grid.(['domain_', tmp.subregion]), tmp.error_status] = ...
                Func_0007_get_polygon_data_from_regionname(tmp.subregion);
            RCM_grid.(['mask_', tmp.subregion]) = ...
                double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho, ...
                RCM_grid.(['polygon_', tmp.subregion])(:,1), ...
                RCM_grid.(['polygon_', tmp.subregion])(:,2)));
            RCM_grid.(['mask_', tmp.subregion])(RCM_grid.(['mask_', tmp.subregion])==0)=NaN;
            GCM_grid.(['mask_', tmp.subregion]) = ...
                double(inpolygon(GCM_grid.lon,GCM_grid.lat, ...
                RCM_grid.(['polygon_', tmp.subregion])(:,1), ...
                RCM_grid.(['polygon_', tmp.subregion])(:,2)));
            GCM_grid.(['mask_', tmp.subregion])(GCM_grid.(['mask_', tmp.subregion])==0)=NaN;
            CMEMS_grid.(['mask_', tmp.subregion]) = ...
                double(inpolygon(CMEMS_grid.lon2,CMEMS_grid.lat2, ...
                RCM_grid.(['polygon_', tmp.subregion])(:,1), ...
                RCM_grid.(['polygon_', tmp.subregion])(:,2)));
            CMEMS_grid.(['mask_', tmp.subregion])(CMEMS_grid.(['mask_', tmp.subregion])==0)=NaN;
            for j=1:RCM_time.tlen
                [RCM_mean_data.([RCM_info.testname]).(tmp.subregion)(j), tmp.error_status] = ...
                    Func_0011_get_area_weighted_mean(RCM_data.yearly_mean(:,:,j) .* RCM_grid.(['mask_', tmp.subregion]), ...
                    RCM_grid.lon_rho, RCM_grid.lat_rho);
                [GCM_mean_data.([RCM_info.testname]).(tmp.subregion)(j), tmp.error_status] = ...
                    Func_0011_get_area_weighted_mean(GCM_data.yearly_mean(:,:,j) .* GCM_grid.(['mask_', tmp.subregion]), ...
                    GCM_grid.lon, GCM_grid.lat);
                [GCM_mean_data_interped.([RCM_info.testname]).(tmp.subregion)(j), tmp.error_status] = ...
                    Func_0011_get_area_weighted_mean(GCM_data_interped.yearly_mean(:,:,j) .* CMEMS_grid.(['mask_', tmp.subregion]), ...
                    CMEMS_grid.lon2, CMEMS_grid.lat2);
            end
            
            if testnameind==1
                RCM_mean_data.ens.(tmp.subregion)=zeros(1,length(RCM_mean_data.([RCM_info.testname]).(tmp.subregion)));
                RCM_mean_data.ens.(tmp.subregion)= ...
                    RCM_mean_data.ens.(tmp.subregion) + RCM_mean_data.([RCM_info.testname]).(tmp.subregion)/length(RCM_info.name);
                GCM_mean_data.ens.(tmp.subregion)=zeros(1,length(GCM_mean_data.([RCM_info.testname]).(tmp.subregion)));            
                GCM_mean_data.ens.(tmp.subregion)= ...
                    GCM_mean_data.ens.(tmp.subregion) + GCM_mean_data.([RCM_info.testname]).(tmp.subregion)/length(RCM_info.name);
            else
                RCM_mean_data.ens.(tmp.subregion)= ...
                    RCM_mean_data.ens.(tmp.subregion) + RCM_mean_data.([RCM_info.testname]).(tmp.subregion)/length(RCM_info.name);
                GCM_mean_data.ens.(tmp.subregion)= ...
                    GCM_mean_data.ens.(tmp.subregion) + GCM_mean_data.([RCM_info.testname]).(tmp.subregion)/length(RCM_info.name);
            end
            
            [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
                 squeeze(RCM_mean_data.([RCM_info.testname]).(tmp.subregion))'.*1000.0, 'poly1'); %% m/y -> mm/y
            RCM_mean_trend.([RCM_info.testname]).(tmp.subregion) = tmp_fitobject.p1;
            [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
                 squeeze(GCM_mean_data.([RCM_info.testname]).(tmp.subregion))'.*1000.0, 'poly1'); %% m/y -> mm/y
            GCM_mean_trend.([RCM_info.testname]).(tmp.subregion) = tmp_fitobject.p1;
        end
% % %         ex) East Sea
%         [RCM_grid.polygon_ES_KHOA, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('ES_KHOA');
%         RCM_grid.mask_model_ES_KHOA = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,RCM_grid.polygon_ES_KHOA(:,1),RCM_grid.polygon_ES_KHOA(:,2)));
%         for i=1:RCM_time.tlen
%             [RCM_mean_data.([RCM_info.testname]).ES_KHOA(i), tmp.error_status] = ...
%                 Func_0011_get_area_weighted_mean(RCM_data.yearly_mean(:,:,i), RCM_grid.lon_rho, RCM_grid.lat_rho);
%         end
        
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

        RCM_grid.(['lon_rho_', RCM_info.regionname])=RCM_grid.lon_rho.*RCM_grid.mask_ocean;
        RCM_grid.(['lat_rho_', RCM_info.regionname])=RCM_grid.lat_rho.*RCM_grid.mask_ocean;
        GCM_grid.(['lon_', RCM_info.regionname])=GCM_grid.lon.*GCM_grid.mask_ocean;
        GCM_grid.(['lat_', RCM_info.regionname])=GCM_grid.lat.*GCM_grid.mask_ocean;
        GCM_grid.(['lon_interped_', RCM_info.regionname])=CMEMS_grid.lon2.*GCM_grid.mask_ocean_interped;
        GCM_grid.(['lat_interped_', RCM_info.regionname])=CMEMS_grid.lat2.*GCM_grid.mask_ocean_interped;

% % %             get mean sea-level (tidal stations)
        for ind_sta=1:size(RCM_info.tidal_st,1)
% % %                 RCM
            if (testnameind==1)
                for i=1:size(RCM_grid.(['lon_rho_', RCM_info.regionname]),1)
                    for j=1:size(RCM_grid.(['lon_rho_', RCM_info.regionname]),2)
                        dist(i,j)=m_lldist([RCM_grid.(['lon_rho_', RCM_info.regionname])(i,j), ...
                            RCM_info.tidal_st(ind_sta,1)], [RCM_grid.(['lat_rho_', RCM_info.regionname])(i,j), RCM_info.tidal_st(ind_sta,2)]);
                    end
                end
                RCM_info.tidal_ind_sta(ind_sta)=find(dist(:)==min(dist(:)));
                RCM_mean_data.ind_sta_model_lon(ind_sta)= ...
                    mod(RCM_info.tidal_ind_sta(ind_sta),size(RCM_grid.(['lon_rho_', RCM_info.regionname]),1));
                RCM_mean_data.ind_sta_model_lat(ind_sta)= ...
                    floor(RCM_info.tidal_ind_sta(ind_sta)/size(RCM_grid.(['lon_rho_', RCM_info.regionname]),1))+1;
            end
            RCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:) = ...
                RCM_data.yearly_mean(RCM_mean_data.ind_sta_model_lon(ind_sta), RCM_mean_data.ind_sta_model_lat(ind_sta), :);
            RCM_mean_data.([RCM_info.testname]).tidal_station_lon(ind_sta,:) = ...
                RCM_grid.lon_rho(RCM_mean_data.ind_sta_model_lon(ind_sta), RCM_mean_data.ind_sta_model_lat(ind_sta));
            RCM_mean_data.([RCM_info.testname]).tidal_station_lat(ind_sta,:) = ...
                RCM_grid.lat_rho(RCM_mean_data.ind_sta_model_lon(ind_sta), RCM_mean_data.ind_sta_model_lat(ind_sta));
            clear dist
            
            [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
                 squeeze(RCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:))'.*1000.0, 'poly1'); %% m/y -> mm/y
            RCM_mean_trend.([RCM_info.testname]).tidal_station(ind_sta) = tmp_fitobject.p1;
            
% % %                 GCM
            for i=1:size(GCM_grid.(['lon_', RCM_info.regionname]),1)
                for j=1:size(GCM_grid.(['lon_', RCM_info.regionname]),2)
                    dist(i,j)=m_lldist([GCM_grid.(['lon_', RCM_info.regionname])(i,j), ...
                        RCM_info.tidal_st(ind_sta,1)], [GCM_grid.(['lat_', RCM_info.regionname])(i,j), RCM_info.tidal_st(ind_sta,2)]);
                end
            end
            GCM_info.tidal_ind_sta(ind_sta)=find(dist(:)==min(dist(:)));
            GCM_mean_data.([RCM_info.testname]).ind_sta_model_lon(ind_sta)= ...
                mod(GCM_info.tidal_ind_sta(ind_sta),size(GCM_grid.(['lon_', RCM_info.regionname]),1));
            GCM_mean_data.([RCM_info.testname]).ind_sta_model_lat(ind_sta)= ...
                floor(GCM_info.tidal_ind_sta(ind_sta)/size(GCM_grid.(['lat_', RCM_info.regionname]),1))+1;
            GCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:) = ...
                GCM_data.yearly_mean(GCM_mean_data.([RCM_info.testname]).ind_sta_model_lon(ind_sta), GCM_mean_data.([RCM_info.testname]).ind_sta_model_lat(ind_sta), :);
            GCM_mean_data.([RCM_info.testname]).tidal_station_lon(ind_sta,:) = ...
                GCM_grid.lon(GCM_mean_data.([RCM_info.testname]).ind_sta_model_lon(ind_sta), GCM_mean_data.([RCM_info.testname]).ind_sta_model_lat(ind_sta));
            GCM_mean_data.([RCM_info.testname]).tidal_station_lat(ind_sta,:) = ...
                GCM_grid.lat(GCM_mean_data.([RCM_info.testname]).ind_sta_model_lon(ind_sta), GCM_mean_data.([RCM_info.testname]).ind_sta_model_lat(ind_sta));
            clear dist
            [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
                 squeeze(GCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:))'.*1000.0, 'poly1'); %% m/y -> mm/y
            GCM_mean_trend.([RCM_info.testname]).tidal_station(ind_sta) = tmp_fitobject.p1;
            
            if testnameind==1
                RCM_mean_data.ens.tidal_station(ind_sta,:)=zeros(1,length(RCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:)));
                RCM_mean_data.ens.tidal_station(ind_sta,:)= ...
                    RCM_mean_data.ens.tidal_station(ind_sta,:) + RCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:)/length(RCM_info.name);
                GCM_mean_data.ens.tidal_station(ind_sta,:)=zeros(1,length(GCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:)));            
                GCM_mean_data.ens.tidal_station(ind_sta,:)= ...
                    GCM_mean_data.ens.tidal_station(ind_sta,:) + GCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:)/length(RCM_info.name);
            else
                RCM_mean_data.ens.tidal_station(ind_sta,:)= ...
                    RCM_mean_data.ens.tidal_station(ind_sta,:) + RCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:)/length(RCM_info.name);
                GCM_mean_data.ens.tidal_station(ind_sta,:)= ...
                    GCM_mean_data.ens.tidal_station(ind_sta,:) + GCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,:)/length(RCM_info.name);
            end
        end
    end
end


% % %         get ens trend (domain)
[tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
     squeeze(RCM_mean_data.ens.(RCM_info.regionname))'.*1000.0, 'poly1'); %% m/y -> mm/y
RCM_mean_trend.ens.(RCM_info.regionname) = tmp_fitobject.p1;
[tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
     squeeze(GCM_mean_data.ens.(RCM_info.regionname))'.*1000.0, 'poly1'); %% m/y -> mm/y
GCM_mean_trend.ens.(RCM_info.regionname) = tmp_fitobject.p1;
% % %         get ens trend (subregion)
for i=1:length(RCM_info.subregions)
    tmp.subregion=RCM_info.subregions{i};
    [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
         squeeze(RCM_mean_data.ens.(tmp.subregion))'.*1000.0, 'poly1'); %% m/y -> mm/y
    RCM_mean_trend.ens.(tmp.subregion) = tmp_fitobject.p1;
    [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
         squeeze(GCM_mean_data.ens.(tmp.subregion))'.*1000.0, 'poly1'); %% m/y -> mm/y
    GCM_mean_trend.ens.(tmp.subregion) = tmp_fitobject.p1;
end
% % %         get ens trend (tidal station)
for ind_sta=1:size(RCM_info.tidal_st,1)
    [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
         squeeze(RCM_mean_data.ens.tidal_station(ind_sta,:))'.*1000.0, 'poly1'); %% m/y -> mm/y
    RCM_mean_trend.ens.tidal_station(ind_sta) = tmp_fitobject.p1;
    [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
         squeeze(GCM_mean_data.ens.tidal_station(ind_sta,:))'.*1000.0, 'poly1'); %% m/y -> mm/y
    GCM_mean_trend.ens.tidal_station(ind_sta) = tmp_fitobject.p1;
end
a=1

% % % get std of ens members (domain)
for i=1:RCM_time.tlen
    for testnameind=1:length(RCM_info.name)
        RCM_info.testname = RCM_info.name{testnameind};
        if testnameind==1
            tmp.RCM_ens_members=RCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)(i);
            tmp.GCM_ens_members=GCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)(i);
        else
            tmp.RCM_ens_members=[tmp.RCM_ens_members, RCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)(i)];
            tmp.GCM_ens_members=[tmp.GCM_ens_members, GCM_mean_data.([RCM_info.testname]).(RCM_info.regionname)(i)];
        end
        RCM_mean_data.ens.([RCM_info.regionname, '_std'])(i) = std(tmp.RCM_ens_members);
        GCM_mean_data.ens.([RCM_info.regionname, '_std'])(i) = std(tmp.GCM_ens_members);
    end
end

% % % get std of ens members (subregion)
for i=1:length(RCM_info.subregions)
    tmp.subregion=RCM_info.subregions{i};
    for j=1:RCM_time.tlen
        for testnameind=1:length(RCM_info.name)
            RCM_info.testname = RCM_info.name{testnameind};
            if testnameind==1
                tmp.RCM_ens_members=RCM_mean_data.([RCM_info.testname]).(tmp.subregion)(j);
                tmp.GCM_ens_members=GCM_mean_data.([RCM_info.testname]).(tmp.subregion)(j);
            else
                tmp.RCM_ens_members=[tmp.RCM_ens_members, RCM_mean_data.([RCM_info.testname]).(tmp.subregion)(j)];
                tmp.GCM_ens_members=[tmp.GCM_ens_members, GCM_mean_data.([RCM_info.testname]).(tmp.subregion)(j)];
            end
            RCM_mean_data.ens.([tmp.subregion, '_std'])(j) = std(tmp.RCM_ens_members);
            GCM_mean_data.ens.([tmp.subregion, '_std'])(j) = std(tmp.GCM_ens_members);
        end
    end
end


% % % get std of ens members (tidal station)
for ind_sta=1:size(RCM_info.tidal_st,1)
    for i=1:RCM_time.tlen
        for testnameind=1:length(RCM_info.name)
            RCM_info.testname = RCM_info.name{testnameind};
            if testnameind==1
                tmp.RCM_ens_members=RCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,i);
                tmp.GCM_ens_members=GCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,i);
            else
                tmp.RCM_ens_members=[tmp.RCM_ens_members, RCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,i)];
                tmp.GCM_ens_members=[tmp.GCM_ens_members, GCM_mean_data.([RCM_info.testname]).tidal_station(ind_sta,i)];
            end
            RCM_mean_data.ens.(['tidal_station_std'])(ind_sta,i) = std(tmp.RCM_ens_members);
            GCM_mean_data.ens.(['tidal_station_std'])(ind_sta,i) = std(tmp.GCM_ens_members);
        end
    end
end

% RCM_mean_data.ens=rmfield(RCM_mean_data.ens, 'AKP4_rand')
% for k=1:1000
%     for i=1:RCM_time.tlen
%         for j=1:1500 %% if > 1500, mean error of std = -0.0043, max err = 0.3687, min err = -0.3824 (1,000 tries)
%             RCM_mean_data.ens.([RCM_info.regionname, '_rand'])(j,i) = ...
%                 normrnd(RCM_mean_data.ens.(RCM_info.regionname)(i), RCM_mean_data.ens.([(RCM_info.regionname),'_std'])(i));
%         end
%         RCM_mean_data.ens.([RCM_info.regionname, '_rand_std'])(i)= std(RCM_mean_data.ens.([RCM_info.regionname, '_rand'])(:,i));
%     end
%     err(k)=mean(RCM_mean_data.ens.AKP4_rand_std - RCM_mean_data.ens.AKP4_std).*1000;
% end


RCM_mean_data.ens_lap = 1500;

% % % get population from avg and std; get trend from population (domain)
for i=1:RCM_time.tlen
    for j=1:RCM_mean_data.ens_lap
        RCM_mean_data.ens.([RCM_info.regionname, '_rand'])(j,i) = ...
            normrnd(RCM_mean_data.ens.(RCM_info.regionname)(i), RCM_mean_data.ens.([(RCM_info.regionname),'_std'])(i));
        GCM_mean_data.ens.([RCM_info.regionname, '_rand'])(j,i) = ...
            normrnd(GCM_mean_data.ens.(RCM_info.regionname)(i), GCM_mean_data.ens.([(RCM_info.regionname),'_std'])(i));
    end
    RCM_mean_data.ens.([RCM_info.regionname, '_rand_std'])(i)= std(RCM_mean_data.ens.([RCM_info.regionname, '_rand'])(:,i));
    GCM_mean_data.ens.([RCM_info.regionname, '_rand_std'])(i)= std(GCM_mean_data.ens.([RCM_info.regionname, '_rand'])(:,i));
end
for j=1:RCM_mean_data.ens_lap
    [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
         squeeze(RCM_mean_data.ens.([RCM_info.regionname, '_rand'])(j,:))'.*1000.0, 'poly1'); %% m/y -> mm/y
    RCM_mean_trend.ens.([RCM_info.regionname, '_rand'])(j) = tmp_fitobject.p1;
    [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
         squeeze(GCM_mean_data.ens.([RCM_info.regionname, '_rand'])(j,:))'.*1000.0, 'poly1'); %% m/y -> mm/y
    GCM_mean_trend.ens.([RCM_info.regionname, '_rand'])(j) = tmp_fitobject.p1;
end

% % quantile(RCM_mean_trend.ens.([RCM_info.regionname, '_rand']), 0.83)

% % % get population from avg and std; get trend from population (subregion)
for k=1:length(RCM_info.subregions)
    tmp.subregion=RCM_info.subregions{k};
    for i=1:RCM_time.tlen
        for j=1:RCM_mean_data.ens_lap
            RCM_mean_data.ens.([tmp.subregion, '_rand'])(j,i) = ...
                normrnd(RCM_mean_data.ens.(tmp.subregion)(i), RCM_mean_data.ens.([(tmp.subregion),'_std'])(i));
            GCM_mean_data.ens.([tmp.subregion, '_rand'])(j,i) = ...
                normrnd(GCM_mean_data.ens.(tmp.subregion)(i), GCM_mean_data.ens.([(tmp.subregion),'_std'])(i));
        end
        RCM_mean_data.ens.([tmp.subregion, '_rand_std'])(i)= std(RCM_mean_data.ens.([tmp.subregion, '_rand'])(:,i));
        GCM_mean_data.ens.([tmp.subregion, '_rand_std'])(i)= std(GCM_mean_data.ens.([tmp.subregion, '_rand'])(:,i));
    end
    for j=1:RCM_mean_data.ens_lap
        [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
             squeeze(RCM_mean_data.ens.([tmp.subregion, '_rand'])(j,:))'.*1000.0, 'poly1'); %% m/y -> mm/y
        RCM_mean_trend.ens.([tmp.subregion, '_rand'])(j) = tmp_fitobject.p1;
        [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
             squeeze(GCM_mean_data.ens.([tmp.subregion, '_rand'])(j,:))'.*1000.0, 'poly1'); %% m/y -> mm/y
        GCM_mean_trend.ens.([tmp.subregion, '_rand'])(j) = tmp_fitobject.p1;
    end
end

% % % get population from avg and std; get trend from population (tidal stations)
for ind_sta=1:size(RCM_info.tidal_st,1)
    for i=1:RCM_time.tlen
        for j=1:RCM_mean_data.ens_lap
            RCM_mean_data.ens.(['tidal_station_rand'])(j,ind_sta,i) = ...
                normrnd(RCM_mean_data.ens.tidal_station(ind_sta,i), ...
                RCM_mean_data.ens.('tidal_station_std')(ind_sta,i));
            GCM_mean_data.ens.(['tidal_station_rand'])(j,ind_sta,i) = ...
                normrnd(GCM_mean_data.ens.tidal_station(ind_sta,i), ...
                GCM_mean_data.ens.('tidal_station_std')(ind_sta,i));
        end
        RCM_mean_data.ens.(['tidal_station_rand_std'])(ind_sta,i)= ...
            std(RCM_mean_data.ens.('tidal_station_rand')(:,ind_sta,i));
        GCM_mean_data.ens.(['tidal_station_rand_std'])(ind_sta,i)= ...
            std(GCM_mean_data.ens.('tidal_station_rand')(:,ind_sta,i));
    end
    for j=1:RCM_mean_data.ens_lap
        [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
             squeeze(RCM_mean_data.ens.(['tidal_station_rand'])(j,ind_sta,:)).*1000.0, 'poly1'); %% m/y -> mm/y
        RCM_mean_trend.ens.(['tidal_station_rand'])(j, ind_sta) = tmp_fitobject.p1;
        [tmp_fitobject, tmp_gof] = fit(RCM_time.trendtime_yearly', ...
             squeeze(GCM_mean_data.ens.(['tidal_station_rand'])(j,ind_sta,:)).*1000.0, 'poly1'); %% m/y -> mm/y
        GCM_mean_trend.ens.(['tidal_station_rand'])(j, ind_sta) = tmp_fitobject.p1;
    end
end

save(RCM_info.matname_mean_data_trend, 'RCM_mean_data', 'RCM_mean_trend', '-v7.3');
save(GCM_info.matname_mean_data_trend, 'GCM_mean_data', 'GCM_mean_trend', '-v7.3');



% % plot(RCM_mean_data.test2107.(RCM_info.regionname));
% % hold on
% % for i=2:length(RCM_info.name)
% %     RCM_info.testname = RCM_info.name{i}
% %     plot(RCM_mean_data.([RCM_info.testname]).(RCM_info.regionname));
% % end
% % hold off
% % 
% % plot(GCM_mean_data.test2107.(RCM_info.regionname));
% % hold on
% % for i=2:length(RCM_info.name)
% %     RCM_info.testname = RCM_info.name{i}
% %     plot(GCM_mean_data.([RCM_info.testname]).(RCM_info.regionname));
% % end
% % hold off







% f = fieldnames(RCM_mean_data);
% v = getfield(s,f{1},'value');

% RCM_mean_trend.ens.(RCM_info.regionname)=mean(


% %  checking trends of test2111
% % % % RCM_info.testname = 'test2110';
% RCM_info.savedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs, 'run', tmp.fs, ...
%             tmp.variable, tmp.fs];
% % RCM_info.matname_trends = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_ssh_trend_', ...
% %     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
% RCM_info.matname_trends = [RCM_info.savedir,RCM_info.testname,'_','NWP', '_RCM_ssh_trend_', ...
%     num2str(min(2015),'%04i'),'_',num2str(max(2050),'%04i'),'.mat'];
% load(RCM_info.matname_trends)
% pcolor(RCM_data_trend.yearly_trend'); hold on; contour(RCM_data_trend.yearly_trend', 'k', 'ShowText','on'); shading flat; colorbar; caxis([0 10]); hold off
% 
% % % % GCM_info.testname = 'CNRM-CM6-1-HR';
% % GCM_info.matname_trends = [RCM_info.savedir,GCM_info.testname,'_',GCM_info.regionname, '_GCM_ssh_trend_', ...
% %     num2str(min(GCM_info.years),'%04i'),'_',num2str(max(GCM_info.years),'%04i'),'.mat'];
% GCM_info.matname_trends = [RCM_info.savedir,GCM_info.testname,'_','NWP', '_GCM_ssh_trend_', ...
%     num2str(min(2015),'%04i'),'_',num2str(max(2050),'%04i'),'.mat'];
% load(GCM_info.matname_trends)
% pcolor(GCM_data_trend.yearly_trend'); hold on; contour(GCM_data_trend.yearly_trend', 'k', 'ShowText','on'); shading flat; colorbar; caxis([0 10]); hold off
% 


% clear comb_b
% for i=1:5
%     comb_b(i)=nchoosek(5,i)
% end
% sum(comb_b)