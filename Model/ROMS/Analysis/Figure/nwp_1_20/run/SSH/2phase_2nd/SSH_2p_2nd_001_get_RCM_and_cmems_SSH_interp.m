close all; clear all;  clc;
% %  Updated 05-Jul-2021 by Yong-Yub Kim, structure
% %  Updated 23-Jun-2022 by Yong-Yub Kim, minor correction

% % % configuration of RCM
% RCM_info.name={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
% RCM_info.name={'test2107', 'test2108', 'test2109', 'test2110', 'test2111'};
% RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
RCM_info.abbs = {'RCM-CNE', 'RCM-ECV', 'RCM-ACC', 'RCM-CNH', 'RCM-CMC'};
% RCM_info.name={  'test2129'};
% RCM_info.abbs = { 'RCM-ACC'};

RCM_info.model = 'nwp_1_20';
% RCM_info.dataroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';
% RCM_info.region = {'EKWC2'};
% RCM_info.region = {'AKP4'};
RCM_info.region = {'NWP'};

% RCM_info.years = 1985:2014;
% RCM_info.years = 1989:2014;
RCM_info.years = 2015:2100;
% RCM_info.months = [12, 1, 2];
RCM_info.season = 'all';
RCM_grid.dl = 1/20;

% % % configuration of GCM
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
GCM_grid.dl = 1/2;

% % % configuration of CMEMS
CMEMS_info.filedir = 'D:\Data\Observation\CMEMS\';
CMEMS_grid.dl = 1/4;

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

% RCM_info2=RCM_info;
% GCM_info2=GCM_info;

% %  working
for testnameind=1:length(RCM_info.name)
    for regionind=1:length(RCM_info.region)
        clearvars '*' -except regionind testnameind flags flg_lev param ...
            RCM_info RCM_grid GCM_info GCM_grid CMEMS_info CMEMS_grid
        [RCM_info.months,error_status] = Func_0019_get_month_from_season(RCM_info.season);
        GCM_info.months = RCM_info.months;
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
        
%         RCM_info=RCM_info2;
%         GCM_info=GCM_info2;
        
%% set temporary variables (testname, regionname, filesep, dir, ...)
        RCM_info.testname = RCM_info.name{testnameind};
        RCM_info.regionname = RCM_info.region{regionind};
        [RCM_info.abb, tmp.error_status] = Func_0018_abbname(RCM_info.testname);
        [RCM_info.scenario, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(RCM_info.testname);
        
        GCM_info.testname = Func_0004_get_GCMname_from_RCM(RCM_info.testname);
        GCM_info.regionname = RCM_info.regionname;
        [GCM_info.abb, tmp.error_status] = Func_0018_abbname(GCM_info.testname);
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
        
%% flag configuration (process)
        for folding=1:1
            flags.fig_name{1}='get model (RCM, GCM) data and cmems data';
            flags.fig_name{2}='calculation of RCM sea level trend (raw grid)';
            flags.fig_name{3}='calculation of GCM sea level trend (raw grid)';
            flags.fig_name{4}='interped sea level trend analysis';
            flags.fig_name{5}='interped sea level correlation analysis';
            flags.fig_name{6}='low pass filtered interped sea level trend analysis';
            flags.fig_name{7}='low pass filtered interped sea level correlation analysis';
            flags.fig_name{8}='detrended sea level correlation analysis';
            flags.fig_name{9}='moving averaged interped sea level trend analysis';
            flags.fig_name{10}='moving averaged interped sea level correlation analysis';

            for flagi=1:length(flags.fig_name)
                flags.fig_switch(flagi)=0;
            end
            flags.fig_switch(1)=1;
            flags.fig_switch(2)=0;
            flags.fig_switch(3)=0;
            flags.fig_switch(4)=0;
            flags.fig_switch(5)=0;
        end
        
%% time set
        SSH_2p_2nd_001_sub_000_time_set;
        
%% get meridional velocity  
        if flags.fig_switch(1) > 0
            flags.fig_tmp = flags.fig_switch(1);
            SSH_2p_2nd_001_sub_001_get_SSH;
        end

        if flags.fig_switch(2) > 0
            flags.fig_tmp = flags.fig_switch(2);
            SSH_2p_2nd_001_sub_002_cal_SSH_trends;
        end
        
        if flags.fig_switch(3) > 0
            flags.fig_tmp = flags.fig_switch(3);
            SSH_2p_2nd_001_sub_003_cal_GCM_SSH_trends;
        end
        
        if flags.fig_switch(4) > 0
            flags.fig_tmp = flags.fig_switch(4);
            SSH_2p_2nd_001_sub_004_cal_CMEMS_SSH_trends;
        end
    end
end


% clear comb_b
% for i=1:5
%     comb_b(i)=nchoosek(5,i)
% end
% sum(comb_b)