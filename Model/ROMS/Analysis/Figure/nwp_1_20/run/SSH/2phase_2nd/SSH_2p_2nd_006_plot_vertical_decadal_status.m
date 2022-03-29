% %  Updated 22-Feb-2022 by Yong-Yub Kim, make

close all; clear all;  clc;
warning off;

% if(isempty(gcp('nocreate')))
%     parpool(4);
% end

% % % configuration of RCM
        RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
% RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
%         RCM_info.name={'v04', 'v05'};

% RCM_info.abbs = {'RCM-CNE', 'RCM-ECV', 'RCM-ACC', 'RCM-CNH', 'RCM-CMC'};
% RCM_info.abbs ={'RCM-SOD', 'RCM-GRS'};
RCM_info.model = 'nwp_1_20';
% RCM_info.dataroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'cut_ES_stddepth', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'cut_ES_stddepth', filesep];
RCM_info.phase = 'run';  % run or spinup
RCM_info.region = 'ES_KHOA'; % NWP, AKP4, ES_KHOA, YS, ...
% RCM_info.vert_vars = {'vert_temp', 'vert_salt', 'vert_u', 'vert_v'};
RCM_info.vert_vars = {'vert_u','vert_v', 'vert_temp', 'vert_salt', 'vert_temp_nogrd', 'vert_salt_nogrd',};

% lonmin, lonmax, latmin, latmax, depthmin, depthmax
RCM_info.vert_sections = [129, 131, 35, 35, -200, 0;
                          129, 131, 36, 36, -200, 0;
                          129, 131, 37, 37, -200, 0;
                          128, 131, 38, 38, -200, 0;
                          127, 131, 39, 39, -200, 0;
                          127, 131, 40, 40, -200, 0;
                          127, 131, 41, 41, -200, 0;
                          138, 141, 43, 43, -300, 0;
                          138, 142, 47, 47, -300, 0
                          130, 130, 40, 41, -300, 0;
                          131, 131, 41, 43, -300, 0;
                          132, 132, 41, 44, -300, 0;
                          133, 133, 41, 44, -300, 0;
                          134, 134, 41, 44, -300, 0;
                          135, 135, 41, 44, -300, 0;
                          136, 136, 42, 45, -300, 0;
                          137, 137, 43, 46, -300, 0;
                          138, 138, 43, 46, -300, 0; 
                          141.9, 141.9, 45, 47, -100, 0];
RCM_info.years = 1985:2014;  

% seasons_group={'all'};
% seasons_group={'all', 'spring', 'summer', 'fall', 'winter'};
% seasons_group={'winter'};
% seasons_group={'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'};
seasons_group={'February', 'August'};


RCM_grid.dl = 1/20;
% RCM_grid.gridname = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v', 'lon_psi', 'lat_psi', 'pm', 'pn', 'f'};
% RCM_info.testnames = {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};

% all_region2 ={'NWP', 'YS', 'AKP2'}

% all_region2 ={'NWP'};

% RCM_info.vars = {'SST', 'SSH', 'SSS'};
% RCM_info.vars = {'SSH'};

% all_region2 ={'NWP'}
for sectionind2 = 1:size(RCM_info.vert_sections, 1)
    for seasons_groupi=1:length(seasons_group)
        for testnameind2=1:length(RCM_info.name)
                close all;
                clearvars '*' -except RCM_info RCM_grid testnameind2 sectionind2 ...
                    years_groupi years_group seasons_group seasons_groupi season
                tmp.fs=filesep;
                
                RCM_info.vert_section = RCM_info.vert_sections(sectionind2,:);
                RCM_info.season=seasons_group{seasons_groupi};
                [RCM_info.months,error_status] = Func_0019_get_month_from_season(RCM_info.season);

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

                tmp.shadlev = [0 35];
                tmp.abstrendlev =[4 7];
                tmp.reltrendlev =[-5 5];
                tmp.conlev  = 0:5:35;
                

                % for snu_desktop
                tmp.testname=RCM_info.name{testnameind2};   % % need to change
    %             tmp.abb=RCM_info.abbs{testnameind2};    % % need to change
                [tmp.abb, tmp.error_status] = Func_0018_abbname(tmp.testname);

        %         RCM_info.years = [1985:2014]; % % put year which you want to plot [year year ...]
        %           RCM_info.months = [1:12]; % % put month which you want to plot [month month ...]

                flags.fig_name{1}='vertical variables plot';

                for flagi=1:7
                    fig_flags{flagi,2}=0;
                end
                flags.fig_switch(1)=1;  %1 or 2

                tmp.variable ='zeta';
        %         run('nwp_polygon_point.m');
                tmp.regionname=RCM_info.region;
                [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

                [cmaps.bwrmap, tmp.error_status] = Func_0009_get_colormaps('bwr', tmp.dropboxpath);
                cmaps.wrmap = cmaps.bwrmap(51:100,:);
                [cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);
                [cmaps.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr', tmp.dropboxpath);

                dirs.figrawdir =strcat('D:\MEPL\project\SSH\7th_year(2022)\figure\nwp_1_20\'); % % where figure files will be saved
                tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param2_kyy_', tmp.regionname, '.m'];
                dirs.filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\cut_ES_stddepth\', tmp.testname, '\run\'); % % where data files are          
                dirs.matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\vert_mean\');
                dirs.griddir = strcat('D:\Data\Model\ROMS\nwp_1_20\cut_ES_stddepth\'); % % where grid data files are            
                dirs.vert_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\cut_ES_stddepth\', tmp.testname, '\'); % % where data files are          
                
%                 for gridi=1:length(RCM_grid.gridname)
%                     RCM_grid.(['filename_', RCM_grid.gridname{gridi}])=[dirs.griddir, 'NWP_pck_ocean_', RCM_grid.gridname{gridi}, '_NWP.nc'];
%                 end

                run(tmp.param_script);

                if (exist(strcat(dirs.matdir) , 'dir') ~= 7)
                    mkdir(strcat(dirs.matdir));
                end 

                if flags.fig_switch(1) > 0
                    flags.fig_tmp = flags.fig_switch(1);
                    SSH_2p_2nd_006_sub_001_vertical_variable_plot;
                end
        end %test
    end % season
end %section