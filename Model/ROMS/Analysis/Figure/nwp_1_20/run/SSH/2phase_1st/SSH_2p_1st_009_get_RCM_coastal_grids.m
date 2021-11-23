close all; clear all;  clc;
% %  Updated 16-Nov-2021 by Yong-Yub Kim

% % % configuration of RCM
% RCM_info.name={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
% RCM_info.name={'test2107', 'test2108', 'test2109', 'test2110', 'test2111'};
RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
RCM_info.abbs = {'RCM-CNRM', 'RCM-EC-Veg', 'RCM-ACC', 'RCM-CNRM-HR', 'RCM-CMCC'};
% RCM_info.name={  'test2107'};
% RCM_info.abbs = {  'RCM-CNRM'};

RCM_info.model = 'nwp_1_20';
% RCM_info.dataroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';
% RCM_info.region = {'NWP', 'AKP4'};
RCM_info.region = {'AKP4'};
RCM_info.years = 1985:2014;
% RCM_info.years = 2015:2050;
RCM_info.months = 1:12;
RCM_info.coastal_distance = 30; %(km)
RCM_grid.dl = 1/20;
RCM_grid.adm_areas={'SK_coastal', 'SK_EEZ', 'adm_div_all', 'adm_div_YS', 'adm_div_SS', 'adm_div_ES' ...
                        'adm_div_GGD', 'adm_div_CCND', 'adm_div_JBD', 'adm_div_JND', 'adm_div_JJD' ...
                        'adm_div_GSND', 'adm_div_GSBD', 'adm_div_GSBD_coastal', 'adm_div_ULD' ...
                        'adm_div_ULD_only', 'adm_div_DD_only', 'adm_div_GWD'};



% % % configuration of flags (make file)
flags.make_std_mat = 1;
flags.make_std_nc = 1;

% % % configuration of figure levels
param.fig_lev_shad = [-2 2];
param.fig_lev_shad_trend = [0 4];
param.fig_lev_shad_bias = [-4 4];
param.fig_lev_con = 0:5:35;
param.fig_lev_rms = [0 4];

[tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

% %  working
for testnameind=1:length(RCM_info.name)
    for regionind=1:length(RCM_info.region)
        clearvars '*' -except regionind testnameind flags flg_lev param ...
            RCM_info RCM_grid GCM_info GCM_grid CMEMS_info CMEMS_grid
        
        tmp.variable = 'zeta';
        tmp.variable_GCM = 'zos';
        tmp.variable_CMEMS = 'sla';
        tmp.fs = filesep; % file separator win = '\', linux = '/'

        % %     set dropbox path
        if (strcmp(computer,'PCWIN64'))
            tmp.dropboxpath = 'C:\Users\User\Dropbox';
        else
            tmp.dropboxpath = '/home/kimyy/Dropbox';
        end
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));

% %     set temporary variables (testname, regionname, filesep, ...)
        tmp.testname = RCM_info.name{testnameind};
        tmp.regionname = RCM_info.region{regionind};
        tmp.abb = RCM_info.abbs{testnameind};
        [tmp.scenario, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(tmp.testname);
                
        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
            tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
            'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));
        
        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);
        
        for admi=1:length(RCM_grid.adm_areas)
            [RCM_grid.([RCM_grid.adm_areas{admi}, '_polygon']), RCM_grid.(['domain_',RCM_grid.adm_areas{admi}]), tmp.error_status] = ...
                Func_0007_get_polygon_data_from_regionname(RCM_grid.adm_areas{admi});
        end
%         [RCM_grid.SK_coastal_polygon, RCM_grid.domain_SK_coastal, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('SK_coastal');
%         [RCM_grid.SK_EEZ_polygon, RCM_grid.domain_SK_EEZ, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('SK_EEZ');
%         [RCM_grid.adm_div_all_polygon, RCM_grid.domain_adm_div_all, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('adm_div_all');
%         [RCM_grid.adm_div_YS_polygon, RCM_grid.domain_adm_div_YS, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('adm_div_YS');
%         [RCM_grid.adm_div_SS_polygon, RCM_grid.domain_adm_div_SS, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('adm_div_SS');
%         [RCM_grid.adm_div_ES_polygon, RCM_grid.domain_adm_div_SS, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('adm_div_ES');

        tmp.param_script = [tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', ...
            tmp.fs, 'Model', tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', ...
            tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, 'fig_param', tmp.fs, ...
            'fig_param2_kyy_', tmp.regionname, '.m'];
        RCM_info.filedir = [RCM_info.dataroot, tmp.testname, tmp.fs, 'run', tmp.fs, ...
            tmp.variable, tmp.fs];
        RCM_info.savedir = [RCM_info.dataroot, tmp.testname, tmp.fs, 'run', tmp.fs, ...
            tmp.variable, tmp.fs];

        
% % %         flag configuration (process)
        for folding=1:1
            flags.fig_name{1}='get coastal data from RCM';

            for flagi=1:length(flags.fig_name)
                flags.fig_switch(flagi)=0;
            end
            flags.fig_switch(1)=2;
            flags.fig_switch(2)=1;
        end
        
        if flags.fig_switch(1) > 0
            flags.fig_tmp = flags.fig_switch(1);
            SSH_2p_1st_009_sub_001_get_coastal_data;
        end
        
        if flags.fig_switch(2) > 0
            flags.fig_tmp = flags.fig_switch(2);
            SSH_2p_1st_009_sub_002_cal_slr_adm_region;
        end
        

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
        

        

    end
end
