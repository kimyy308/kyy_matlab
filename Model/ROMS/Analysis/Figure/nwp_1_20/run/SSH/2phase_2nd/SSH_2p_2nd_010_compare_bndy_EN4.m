close all; clc; clear all;



    error_status=1;


% % % configuration of RCM
 
% RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
RCM_info.name={'test2117'};

RCM_info.model = 'nwp_1_20';

% RCM_info.dataroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';  % run or spinup

% RCM_info.vars = {'SSS'};  % 'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'shflux', 'u', 'v'

% RCM_info.vert_vars = {'salt'};
% RCM_info.vert_section = [129, 131, 37, 37, -200, 0]; % lonmin, lonmax, latmin, latmax, depthmin, depthmax

% RCM_info.vars = {'wstrcurl'};
RCM_info.bndy_directions={'north', 'east', 'south', 'west'};

% RCM_info.vars = {'Uwind', 'Vwind'};

% RCM_info.vars = {'SSH'};
% RCM_info.years = 1985:2014;  
RCM_info.years = 1993:2014;  
% RCM_info.years = 1995:2014;  
%         RCM_info.years = 2030:2030;  
% RCM_info.years = years_group(years_groupi);  

%         RCM_info.months =1:12;  

% seasons_group={'all'};

% seasons_group={ 'all', 'February', 'August', 'July'};
seasons_group={ 'all'};

RCM_grid.dl = 1/20;
RCM_grid.gridname = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v', 'lon_psi', 'lat_psi', 'pm', 'pn', 'f', 'h', ...
    'mask_rho'};


for seasons_groupi=1:length(seasons_group)
    for testnameind2=1:length(RCM_info.name)
        close all;
        clearvars '*' -except RCM_info RCM_grid testnameind2 regionind2 years_groupi years_group seasons_group seasons_groupi season
        filesep=filesep;  

        % % % 
        % %     set dropbox path
        if (strcmp(computer,'PCWIN64'))
            tmp.dropboxpath = 'C:\Users\User\Dropbox';
        else
            tmp.dropboxpath = '/home/kimyy/Dropbox';
        end
        addpath(genpath([tmp.dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'function']));
        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Model' ...
            filesep, 'ROMS', filesep, 'Analysis', filesep, 'Figure', filesep, 'nwp_1_20', filesep ...
            'run', filesep, 'SSH', filesep, '2phase_2nd', filesep, 'subroutine']));
        RCM_grid.stddepth=[0, -10,-20,-30, -40, -50, ...
        -60, -70, -80, -90, -100, ...
        -110, -120, -130, -140, -150, ...
        -160, -170, -180, -190, -200, ...
        -220, -240, -260, -280, -300, ...
        -350, -400,-500, ...
        -700,-1000,-1250,-1500,...
        -1750,-2000,-2250,-2500,-3000,-3500, -4000, -4500, -5000];
        
        RCM_info.season=seasons_group{seasons_groupi};
        [RCM_info.months,error_status] = Func_0019_get_month_from_season(RCM_info.season);

%         tmp.shadlev = [0 35];
%         tmp.rms_shadlev = [0 4];
%         tmp.trendlev = [-10 10];  %% trend lev
%         tmp.abstrendlev =[4 7];
%         tmp.reltrendlev =[-5 5];
%         tmp.conlev  = 0:5:35;
%         tmp.meanplotlev =[-0.3 0.3];
%         tmp.trendplotlev = [3 7];
%         tmp.sshlev =[-0.7 1.3];
%         tmp.sshdifflev = [40 70];

        % for snu_desktop
        tmp.testname=RCM_info.name{testnameind2};   % % need to change
%             tmp.abb=RCM_info.abbs{testnameind2};    % % need to change
        [tmp.abb, tmp.error_status] = Func_0018_abbname(tmp.testname);
        [tmp.scenname, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(tmp.testname);
        [tmp.testname_GCM] = Func_0004_get_GCMname_from_RCM(tmp.testname);
        [tmp.testname_bndy, tmp.error_status] = Func_0022_get_RCMbndyname_from_GCM(tmp.testname_GCM, tmp.scenname);
%         RCM_info.years = [1985:2014]; % % put year which you want to plot [year year ...]
%           RCM_info.months = [1:12]; % % put month which you want to plot [month month ...]

        flags.fig_name{2}='GCM_bndy_get_data_stddepth';

        for flagi=1:7
            fig_flags{flagi,2}=0;
        end
        flags.fig_switch(1)=0;  %1 or 2
        flags.fig_switch(2)=0;  %1 or 2
        flags.fig_switch(3)=0;  %1 or 2
        flags.fig_switch(4)=0;  %1 or 2
        flags.fig_switch(5)=2;  %1 or 2

        tmp.variable ='salt';
        tmp.variable_GCM = 'so';
        
%         tmp.variable ='temp';
%         tmp.variable_GCM = 'thetao';
        
%         run('nwp_polygon_point.m');
%         tmp.regionname=RCM_info.region{regionind2};
        tmp.regionname='NWP';

        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

        [cmaps.bwrmap, tmp.error_status] = Func_0009_get_colormaps('bwr', tmp.dropboxpath);
        cmaps.wrmap = cmaps.bwrmap(51:100,:);
        [cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);
        [cmaps.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr', tmp.dropboxpath);

        dirs.figrawdir =strcat('D:\MEPL\project\SSH\7th_year(2022)\figure\nwp_1_20\'); % % where figure files will be saved
        dirs.figrawdir_GLORYS =strcat('D:\MEPL\project\SSH\7th_year(2022)\figure\GLORYS\'); % % where figure files will be saved
        tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param2_kyy_', tmp.regionname, '.m'];
%         tmp.param_script =['/home/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/fig_param/fig_param2_kyy_', tmp.regionname, '.m'];
        dirs.filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\backup_surf\', tmp.testname, '\run\'); % % where data files are          
        dirs.matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\mean\');
%         dirs.griddir = strcat('D:\Data\Model\ROMS\nwp_1_20\backup_surf\'); % % where grid data files are     
        dirs.griddir = strcat('D:\Data\Model\ROMS\nwp_1_20\input\test2117\'); % % where grid data files are     
        dirs.vert_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\cut_ES_stddepth\', tmp.testname, '\'); % % where data files are          
        dirs.bndyfiledir = strcat('D:\Data\Model\ROMS\nwp_1_20\input\OBC\', tmp.scenname, '\', tmp.testname_GCM, '\'); % % where data files are          
%         dirs.glorysdir = ['/data2/ykang/data/MYOCEAN/'];
        dirs.glorysdir = ['D:\Data\Reanalysis\MyOcean'];
        RCM_grid.filename = [dirs.griddir, filesep, 'roms_grid_nwp_1_20_test2117.nc'];
        for gridi=1:length(RCM_grid.gridname)
            RCM_grid.(['filename_', RCM_grid.gridname{gridi}])=[dirs.griddir, 'NWP_pck_ocean_', RCM_grid.gridname{gridi}, '_NWP.nc'];
        end

        run(tmp.param_script);

        if (exist(strcat(dirs.matdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.matdir));
        end 

        if flags.fig_switch(1) > 0
            flags.fig_tmp = flags.fig_switch(1);
            SSH_2p_2nd_010_sub_001_get_bndy;
        end
        
        if flags.fig_switch(2) > 0
            flags.fig_tmp = flags.fig_switch(2);
            SSH_2p_2nd_010_sub_002_get_bndy_GLORYS;
        end
        
        if flags.fig_switch(3) > 0
            flags.fig_tmp = flags.fig_switch(3);
            SSH_2p_2nd_010_sub_003_get_RMS_bndy_GLORYS;
        end
        
        if flags.fig_switch(4) > 0
            flags.fig_tmp = flags.fig_switch(4);
            SSH_2p_2nd_010_sub_004_plot_bndy_RCM_GLORYS;
        end
        
        if flags.fig_switch(5) > 0
            flags.fig_tmp = flags.fig_switch(5);
            SSH_2p_2nd_010_sub_005_plot_bndy_GLORYS;
        end
        
    end
end