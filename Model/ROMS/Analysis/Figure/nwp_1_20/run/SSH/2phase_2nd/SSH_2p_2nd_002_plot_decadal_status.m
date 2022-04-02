% %  Updated 09-Oct-2021 by Yong-Yub Kim, structure
% %  Updated 10-Nov-2021 by Yong-Yub Kim, added wind plot module
% %  Updated 07-Dec-2021 by Yong-Yub Kim, added wind plot module

close all; clear all;  clc;
warning off;

% if(isempty(gcp('nocreate')))
%     parpool(4);
% end

% % % configuration of RCM
%         RCM_info.name={ 'test2117', 'test2120'};
%         RCM_info.name={ 'test2118', 'test2121'};   
%         RCM_info.name={ 'test2119'};      
        RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121', 'v04'};
% RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
%         RCM_info.name={'v04', 'v05'};
%         RCM_info.name={'v04'};

% RCM_info.abbs = {'RCM-CNE', 'RCM-ECV', 'RCM-ACC', 'RCM-CNH', 'RCM-CMC'};
% RCM_info.abbs ={'RCM-SOD', 'RCM-GRS'};
RCM_info.model = 'nwp_1_20';
% RCM_info.dataroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';  % run or spinup
% RCM_info.region = {'NES2'}; % NWP, AKP4, ES_KHOA, YS, ...
RCM_info.region = {'ES_KHOA'}; % NWP, AKP4, ES_KHOA, YS, ...
% RCM_info.region = {'EKWC2'}; % NWP, AKP4, ES_KHOA, YS, ...

% RCM_info.vars = {'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'shflux', 'u', 'v'};  % 'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'shflux', 'u', 'v'
% RCM_info.vars = {'SSS'};  % 'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'shflux', 'u', 'v'

% RCM_info.vert_vars = {'salt'};
RCM_info.vert_section = [129, 131, 37, 37, -200, 0]; % lonmin, lonmax, latmin, latmax, depthmin, depthmax

% RCM_info.vars = {'u', 'v'};
% RCM_info.vars = { 'Uwind', 'Vwind'};
% RCM_info.vars = {'wstrcurl', 'SST', 'SSS', 'SSH'};
% RCM_info.vars = {'SSS'};
RCM_info.vars = {'wstrcurl'};

% RCM_info.vars = {'Uwind', 'Vwind'};

% RCM_info.vars = {'SSH'};
RCM_info.years = 1985:2014;  
% RCM_info.years = 1993:2014;  
% RCM_info.years = 1995:2014;  
%         RCM_info.years = 2030:2030;  
% RCM_info.years = years_group(years_groupi);  

%         RCM_info.months =1:12;  

% seasons_group={'all'};
% seasons_group={'all', 'spring', 'summer', 'fall', 'winter'};
% seasons_group={'winter'};
% seasons_group={'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'};
% seasons_group={ 'February', 'August'};
seasons_group={ 'all', 'February', 'August', 'July'};
% seasons_group={ 'all'};


RCM_grid.dl = 1/20;
RCM_grid.gridname = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v', 'lon_psi', 'lat_psi', 'pm', 'pn', 'f'};

% RCM_info.testnames = {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};

% all_region2 ={'NWP', 'YS', 'AKP2'}

% all_region2 ={'NWP'};

% RCM_info.vars = {'SST', 'SSH', 'SSS'};
% RCM_info.vars = {'SSH'};

% all_region2 ={'NWP'}


        
for seasons_groupi=1:length(seasons_group)
    for testnameind2=1:length(RCM_info.name)
        for regionind2=1:length(RCM_info.region)
            close all;
            clearvars '*' -except RCM_info RCM_grid testnameind2 regionind2 years_groupi years_group seasons_group seasons_groupi season
            tmp.fs=filesep;  
            
            % % % 
            % %     set dropbox path
            if (strcmp(computer,'PCWIN64'))
                tmp.dropboxpath = 'C:\Users\User\Dropbox';
            else
                tmp.dropboxpath = '/home/kimyy/Dropbox';
            end
            addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
            addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
                tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
                'run', tmp.fs, 'SSH', tmp.fs, '2phase_2nd', tmp.fs, 'subroutine']));
            
            RCM_info.season=seasons_group{seasons_groupi};
            [RCM_info.months,error_status] = Func_0019_get_month_from_season(RCM_info.season);
            
            tmp.shadlev = [0 35];
            tmp.rms_shadlev = [0 4];
            tmp.trendlev = [-10 10];  %% trend lev
            tmp.abstrendlev =[4 7];
            tmp.reltrendlev =[-5 5];
            tmp.conlev  = 0:5:35;
            tmp.meanplotlev =[-0.3 0.3];
            tmp.trendplotlev = [3 7];
            tmp.sshlev =[-0.7 1.3];
            tmp.sshdifflev = [40 70];

            % for snu_desktop
            tmp.testname=RCM_info.name{testnameind2};   % % need to change
%             tmp.abb=RCM_info.abbs{testnameind2};    % % need to change
            [tmp.abb, tmp.error_status] = Func_0018_abbname(tmp.testname);

    %         RCM_info.years = [1985:2014]; % % put year which you want to plot [year year ...]
    %           RCM_info.months = [1:12]; % % put month which you want to plot [month month ...]

            flags.fig_name{1}='earlier decadal current plot';
            flags.fig_name{2}='later decadal current plot';
            flags.fig_name{3}='earlier decadal SST, SSS plot';
            flags.fig_name{4}='later decadal SST, SSS plot';
            flags.fig_name{5}='earlier decadal YSBCW plot';
            flags.fig_name{6}='later decadal YSBCW plot';
            flags.fig_name{7}='decadal wind plot';
            flags.fig_name{8}='decadal ekman transport plot';
            flags.fig_name{9}='vertical variables plot';

            for flagi=1:7
                fig_flags{flagi,2}=0;
            end
            flags.fig_switch(1)=0;  %1 or 2
            flags.fig_switch(2)=0; 
            flags.fig_switch(3)=0;
            flags.fig_switch(4)=0;
            flags.fig_switch(5)=0;
            flags.fig_switch(6)=0;
            flags.fig_switch(7)=0;
            flags.fig_switch(8)=2;
            flags.fig_switch(9)=0;
            
            tmp.variable ='zeta';
    %         run('nwp_polygon_point.m');
            tmp.regionname=RCM_info.region{regionind2};
            [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

            [cmaps.bwrmap, tmp.error_status] = Func_0009_get_colormaps('bwr', tmp.dropboxpath);
            cmaps.wrmap = cmaps.bwrmap(51:100,:);
            [cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);
            [cmaps.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr', tmp.dropboxpath);
            
            dirs.figrawdir =strcat('D:\MEPL\project\SSH\7th_year(2022)\figure\nwp_1_20\'); % % where figure files will be saved
            tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param2_kyy_', tmp.regionname, '.m'];
            dirs.filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\backup_surf\', tmp.testname, '\run\'); % % where data files are          
            dirs.matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\mean\');
            dirs.griddir = strcat('D:\Data\Model\ROMS\nwp_1_20\backup_surf\'); % % where grid data files are            
            dirs.vert_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\cut_ES_stddepth\', tmp.testname, '\'); % % where data files are          


            for gridi=1:length(RCM_grid.gridname)
                RCM_grid.(['filename_', RCM_grid.gridname{gridi}])=[dirs.griddir, 'NWP_pck_ocean_', RCM_grid.gridname{gridi}, '_NWP.nc'];
            end

            run(tmp.param_script);

            if (exist(strcat(dirs.matdir) , 'dir') ~= 7)
                mkdir(strcat(dirs.matdir));
            end 

            if flags.fig_switch(1) > 0
                flags.fig_tmp = flags.fig_switch(1);
                SSH_2p_2nd_002_sub_001_current_plot;
            end

            if flags.fig_switch(3) > 0
                flags.fig_tmp = flags.fig_switch(3);
                SSH_2p_2nd_002_sub_003_surface_variable_plot;
            end

            if flags.fig_switch(7) > 0
                flags.fig_tmp = flags.fig_switch(7);
                SSH_2p_2nd_002_sub_007_wind_plot;
            end
            
            if flags.fig_switch(8) > 0
                flags.fig_tmp = flags.fig_switch(8);
                SSH_2p_2nd_002_sub_008_ekman_transport_plot;
            end

        end
    end
end