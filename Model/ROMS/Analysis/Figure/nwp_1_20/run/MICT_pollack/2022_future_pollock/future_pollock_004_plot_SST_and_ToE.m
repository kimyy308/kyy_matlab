% %  Updated 21-Apr-2021 by Yong-Yub Kim, 

close all; clear all;  clc;
warning off;

% RCM_info.name={ 'test06'};
% RCM_info.model = 'nwp_1_10';

% RCM_info.name={ 'test2127'};
% RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121', 'ens2202'};
RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131', 'ens2201'};
% RCM_info.name={'test2128', 'test2129', 'test2131'};
% RCM_info.name={'test2131'};
% RCM_info.name={'ens2201'};

% 
RCM_info.model = 'nwp_1_20';

RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, RCM_info.model, filesep, 'backup_surf', filesep];
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, RCM_info.model, filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';  % run or spinup
% RCM_info.region = {'EKB2'}; % NWP, AKP4, ES_KHOA, YS, ...
RCM_info.region = {'pollock_egg3'}; % NWP, AKP4, ES_KHOA, YS, ...

RCM_info.vars = {'SST'};
% RCM_info.vars = {'v'};
% RCM_info.vars = {'SSS', 'SSH', 'u', 'v', 'Uwind', 'Vwind', 'shflux', 'SST', 'swrad', 'lwrad', 'sensible', 'latent'};
% RCM_info.vars = {'SSS', 'SSH', 'Uwind', 'Vwind'};
% RCM_info.vars = {'swrad', 'shflux'};
% RCM_info.vars = {'swrad', 'shflux', 'lwrad', 'sensible', 'latent'};
% RCM_info.vars = {'wstrcurl'};

% RCM_info.years = 1983:2021;  
% RCM_info.years = [2015:2050, 2081:2100];  
% RCM_info.years = [1983:1987];  
% RCM_info.years = [1988:1992];  
% RCM_info.years = [1995:2014];  
% RCM_info.years = [1993:2021];  
RCM_info.years = [2081:2100];  
% RCM_info.years = [2100];  
% RCM_info.years = [2015:2100];  

RCM_info.years_his=[1995:2014];


% seasons_group={'February', 'January', 'JF-'};
% seasons_group={'JF-'};
seasons_group={'all'};

RCM_grid.dl = 1/20;
RCM_grid.gridname = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v', 'lon_psi', 'lat_psi', 'pm', 'pn', 'f'};

for seasons_groupi=1:length(seasons_group)
    for testnameind2=1:length(RCM_info.name)
        for regionind2=1:length(RCM_info.region)
            
            close all;
            clearvars '*' -except RCM_info RCM_grid testnameind2 regionind2 years_groupi years_group seasons_group seasons_groupi season
            tmp.fs=filesep;  
            
            %%     set dropbox path
            addpath(genpath('C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\MICT_pollack\2022_future_pollock\subroutine\'))
            tmp.dropboxpath = 'C:\Users\User\Dropbox';
            addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
            addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
                tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
                'run', tmp.fs, 'SSH', tmp.fs, '2phase_2nd', tmp.fs, 'subroutine']));
            %% get colormaps
            [cmaps.byrmap3, tmp.error_status] = Func_0009_get_colormaps('byr3', tmp.dropboxpath);
            [cmaps.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr2', tmp.dropboxpath);        
            [cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);    
            
            RCM_info.season=seasons_group{seasons_groupi};
            [RCM_info.months,error_status] = Func_0019_get_month_from_season(RCM_info.season);

            tmp.testname=RCM_info.name{testnameind2};   % % need to change
            [tmp.abb, tmp.error_status] = Func_0018_abbname(tmp.testname);
            tmp.regionname=RCM_info.region{regionind2};
            [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

            tmp.variable ='temp';
            
%             dirs.figrawdir =strcat('D:\MEPL\project\SSH\7th_year(2022)\figure\nwp_1_20\'); % % where figure files will be saved
            dirs.figrawdir =strcat('D:\Research\Ph_D_course\2022_pollock_future\figure',filesep, RCM_info.model, filesep); % % where figure files will be saved            
            tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\', RCM_info.model, '\run\fig_param\fig_param2_kyy_', tmp.regionname, '.m'];
            dirs.filedir = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\backup_surf\', tmp.testname, '\run\'); % % where data files are          
            dirs.matdir = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\', tmp.testname, '\run\mean\');
            dirs.griddir = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\backup_surf\'); % % where grid data files are            
            
            dirs.matdir_reana = strcat('D:\Data\Model\ROMS\', 'nwp_1_10', '\', 'test06', '\run\mean\');
            
            for gridi=1:length(RCM_grid.gridname)
                RCM_grid.(['filename_', RCM_grid.gridname{gridi}])=[dirs.griddir, 'NWP_pck_ocean_', RCM_grid.gridname{gridi}, '_NWP.nc'];
            end
            
            flags.fig_switch(1)=2; % get data
            flags.fig_switch(2)=1; % temporal mean pcolor
            flags.fig_switch(3)=0; % yearly spatial mean time series
            flags.fig_switch(4)=0; % temporal mean vec
            flags.fig_switch(5)=0; % get RMSE and bias
            flags.fig_switch(6)=2; % temporal mean pcolor diff
            flags.fig_switch(7)=0; % temporal mean vec diff
            flags.fig_switch(8)=0; % temporal mean windvec diff

%             1/10 SST read
%             comb
%             make yearly
%             get trend
%             detrend 
%             get std
            

%         fig_flag=fig_flags{2,2};
%         while (fig_flag)
%              for checkti=1:length(checktime)
%                 temp_checktime=checktime(checkti);
%                 jpgname=strcat(outfile, '_', testname,'_',regionname, '_loc_prob_', num2str(temp_checktime, '%02i'),'days', ...
%                     num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
%                     num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
%                 future_pollock_003_subroutine_002;
%              end
%             fig_flag=0;
%         end
            
         fig_flag=flags.fig_switch(1);
         while (fig_flag)
            future_pollock_004_subroutine_001
            fig_flag=0;
         end
         
         fig_flag=flags.fig_switch(2);
         while (fig_flag)
            future_pollock_004_subroutine_002
            fig_flag=0;
         end
         
         fig_flag=flags.fig_switch(3);
         while (fig_flag)
            future_pollock_004_subroutine_003
            fig_flag=0;
         end
         
         fig_flag=flags.fig_switch(4);
         while (fig_flag)
            future_pollock_004_subroutine_004
            fig_flag=0;
         end
         
         fig_flag=flags.fig_switch(5);
         while (fig_flag)
            future_pollock_004_subroutine_005
            fig_flag=0;
         end
         
        fig_flag=flags.fig_switch(6);
        if fig_flag>0
            future_pollock_004_subroutine_006
            fig_flag=0;
        end
         
        fig_flag=flags.fig_switch(7);
        if fig_flag>0
            future_pollock_004_subroutine_007
            fig_flag=0;
        end
        
        fig_flag=flags.fig_switch(8);
        if fig_flag>0
            future_pollock_004_subroutine_008
            fig_flag=0;
        end
        
        end
    end
end