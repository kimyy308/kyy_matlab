close all; clear all;  clc;
% %  Updated 03-Feb-2022 by Yong-Yub Kim, structure

%% configuration of RCM
% RCM_info.testnames={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
% RCM_info.testnames={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
RCM_info.testnames={'test2117', 'test2118', 'test2119'};

RCM_info.model = 'nwp_1_20';
RCM_info.figroot = ['D:', filesep, 'Research', filesep, 'Ph_D_course', filesep, ...
    '2022_winter_EKWC_inter_model_difference', filesep, 'figure', filesep];
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'cut_ES_daily', filesep];
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'cut_ES_daily', filesep, 'analysis', filesep];
RCM_info.atmroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.transroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'transport_barot', filesep];
RCM_info.phase = 'run';
RCM_info.region = {'EKWC2'};
RCM_info.years = 1985:1987;
% RCM_info.years = 2015:2050;
% RCM_info.months = [12, 1, 2];
RCM_info.season ='winter';
[RCM_info.months, tmp.error_status]=Func_0019_get_month_from_season(RCM_info.season);
RCM_grid.filename= 'D:\Data\Model\ROMS\nwp_1_20\input\test2117\roms_grid_nwp_1_20_test2117.nc';
RCM_grid.dl = 1/20;
RCM_grid.varnames = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v', 'lon_psi', 'lat_psi', 'pm', 'pn', 'f'};

% % % configuration of GCM
% GCM_info.name={'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
% GCM_info.name={'CNRM-ESM2-1', 'EC-Earth3-Veg'};
% GCM_info.name={'CNRM-ESM2-1'};
% GCM_info.model = GCM_info.name;
GCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'CMIP6', filesep, 'NWP', filesep];
GCM_info.saveroot = GCM_info.dataroot;
GCM_info.phase = RCM_info.phase;
GCM_info.region = RCM_info.region;
GCM_info.years = RCM_info.years;
GCM_info.season = RCM_info.season;
% GCM_info.months = RCM_info.months;
GCM_grid.dl = 1/2;

%% configuration of CMEMS
CMEMS_info.filedir = 'D:\Data\Observation\CMEMS\';
CMEMS_grid.dl = 1/4;

%% configuration of flags (make file)
flags.make_std_mat = 1;
flags.make_std_nc = 1;

%% configuration of figure levels
param.fig_lev_shad = [-2 2];
param.fig_lev_shad_trend = [0 4];
param.fig_lev_shad_bias = [-4 4];
param.fig_lev_con = 0:5:35;
param.fig_lev_rms = [0 4];

%%  working
for testnameind=1:length(RCM_info.testnames)
    for regionind=1:length(RCM_info.region)
        clearvars '*' -except regionind testnameind flags flg_lev param ...
            RCM_info RCM_grid GCM_info GCM_grid CMEMS_info CMEMS_grid
        
        tmp.variable = 'v';
        tmp.variable_GCM = 'vo';
%         tmp.variable = 'u';
%         tmp.variable_GCM = 'uo';
%         tmp.variable = 'temp';
%         tmp.variable_GCM = 'thetao';
%         tmp.variable = 'salt';
%         tmp.variable_GCM = 'so';
%         
        tmp.variable_CMEMS = 'sla';

        tmp.fs = filesep; % file separator win = '\', linux = '/'

        % %     set dropbox path
        if (strcmp(computer,'PCWIN64'))
            tmp.dropboxpath = 'C:\Users\User\Dropbox';
        else
            tmp.dropboxpath = '/home/kimyy/Dropbox';
        end
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%         RCM_info=RCM_info2;
%         GCM_info=GCM_info2;
        
%% set temporary variables (testname, regionname, filesep, dir, ...)
        RCM_info.testname = RCM_info.testnames{testnameind};
        RCM_info.regionname = RCM_info.region{regionind};
        [RCM_info.abb, tmp.error_status] = Func_0018_abbname(RCM_info.testname);
        [RCM_info.scenario, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(RCM_info.testname);
        
%         GCM_info.testname = GCM_info.name{testnameind};
        [GCM_info.testname] =  Func_0004_get_GCMname_from_RCM(RCM_info.testname);
        GCM_info.regionname = RCM_info.regionname;
        [GCM_info.abb, tmp.error_status] = Func_0018_abbname(GCM_info.testname);
        GCM_info.scenario = RCM_info.scenario;
                
        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
            tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
            'run', tmp.fs, 'SSH', tmp.fs, '2phase_2nd', tmp.fs, 'subroutine']));
        
        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(RCM_info.regionname);

        tmp.param_script = [tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', ...
            tmp.fs, 'Model', tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', ...
            tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, 'fig_param', tmp.fs, ...
            'fig_param2_kyy_', RCM_info.regionname, '.m'];
        RCM_info.filedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs];
        RCM_info.atmfiledir = [RCM_info.atmroot, RCM_info.testname, tmp.fs, 'run', tmp.fs];
        RCM_info.transfiledir = [RCM_info.transroot, RCM_info.testname, tmp.fs];
        RCM_info.savedir = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            tmp.variable, tmp.fs];
        RCM_info.windsavedir = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'wind', tmp.fs];
        GCM_info.filedir = [GCM_info.dataroot, tmp.variable_GCM, tmp.fs, GCM_info.scenario, tmp.fs, ...
            'Omon', tmp.fs, GCM_info.testname, tmp.fs];
        if (exist(RCM_info.savedir , 'dir') ~= 7)
            mkdir(RCM_info.savedir);
        end
        if (exist(RCM_info.windsavedir , 'dir') ~= 7)
            mkdir(RCM_info.windsavedir);
        end
%% flag configuration (process)
        for folding=1:1
            flags.fig_name{1}='get model (RCM, GCM) data and cmems data';
            flags.fig_name{2}='get wind data';
            flags.fig_name{3}='get geostrophic & ageostrophic current';
            flags.fig_name{4}='correlation between wind and the EKWC';
            flags.fig_name{5}='get rho';

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
            SSH_2p_2nd_003_sub_001_get_v;
        end 

        if flags.fig_switch(2) > 0
            flags.fig_tmp = flags.fig_switch(2);
            SSH_2p_2nd_003_sub_002_get_wind;
        end
        
        if flags.fig_switch(3) > 0
            flags.fig_tmp = flags.fig_switch(3);
            SSH_2p_2nd_003_sub_003_cal_geo_vel;
        end
        
        if flags.fig_switch(4) > 0
            flags.fig_tmp = flags.fig_switch(4);
            SSH_2p_2nd_003_sub_004_cal_EKWC;
        end    
        
        if flags.fig_switch(5) > 0
            flags.fig_tmp = flags.fig_switch(5);
            SSH_2p_2nd_003_sub_005_cal_rho;
        end        

    end
end


% clear comb_b
% for i=1:5
%     comb_b(i)=nchoosek(5,i)
% end
% sum(comb_b)