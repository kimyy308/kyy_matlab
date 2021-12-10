
% %  Updated 07-Dec-2021 by Yong-Yub Kim, 

close all; clear all;  clc;
warning off;

% years_group=[1990, 2000, 2010];
% years_group=[2015, 2050];
years_group=[2015, 2020, 2030, 2040, 2050];
seasons_group={'all', 'spring', 'summer', 'fall', 'winter'};

for years_groupi=1:length(years_group)
    for seasons_groupi=1:length(seasons_group)


        % % % % % % % % Ensemble of 5 members (GCM initial)
        % RCM_info.name = {'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
        RCM_info.name = {'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};

        % scenname='historical';
%         scenname='ssp585';

        RCM_info.region ={'NWP', 'AKP4'};
        % RCM_info.region ={'YS'};

        % RCM_info.vars = {'SST', 'SSH', 'SSS'};
        % RCM_info.vars = {'SSH'};
        % RCM_info.vars = {'SST', 'SSS'};
        RCM_info.vars = {'SST', 'SSH', 'SSS'};
        % RCM_info.vars = {'BT'};

        % RCM_info.region ={'NWP'}
        for regionind2=1:length(RCM_info.region)
            close all;
        %     clearvars '*' -except regionind2 testnameind2 RCM_info.region RCM_info.name RCM_info.vars
            clearvars '*' -except RCM_info RCM_grid testnameind2 regionind2 years_groupi years_group seasons_group seasons_groupi
            tmp.fs=filesep;

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
                'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));

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
            RCM_info.years = years_group(years_groupi);  

            season=seasons_group{seasons_groupi};
            switch(season)
            case 'all'
                RCM_info.months =1:12;  
            case 'spring'
                RCM_info.months =[3,4,5];  
            case 'summer'
                RCM_info.months =[6,7,8];  
            case 'fall'
                RCM_info.months =[9,10,11];  
            case 'winter'
                RCM_info.months =[12,1,2];  
            end         

            RCM_grid.dl = 1/20;
            RCM_grid.gridname = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v'};
        
            flags.fig_name{1}='earlier decadal current plot';
            flags.fig_name{2}='later decadal current plot';
            flags.fig_name{3}='earlier decadal SST, SSS plot';
            flags.fig_name{4}='later decadal SST, SSS plot';
            flags.fig_name{5}='earlier decadal YSBCW plot';
            flags.fig_name{6}='later decadal YSBCW plot';
            flags.fig_name{7}='earlier decadal wind plot';

            for flagi=1:7
                fig_flags{flagi,2}=0;
            end
            flags.fig_switch(1)=2;  %1 or 2
            flags.fig_switch(2)=0; 
            flags.fig_switch(3)=2;
            flags.fig_switch(4)=0;
            flags.fig_switch(5)=0;
            flags.fig_switch(6)=0;
            flags.fig_switch(7)=0;

            tmp.regionname=RCM_info.region{regionind2};
            [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

            [cmaps.bwrmap, tmp.error_status] = Func_0009_get_colormaps('bwr', tmp.dropboxpath);
            cmaps.wrmap = cmaps.bwrmap(51:100,:);
            [cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);

            dirs.figrawdir =strcat('D:\MEPL\project\SSH\6th_year\figure\nwp_1_20\'); % % where figure files will be saved           
            tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param2_kyy_', tmp.regionname, '.m'];
            dirs.enssavedir =strcat('D:\Data\Model\ROMS\nwp_1_20\2phase_1st\RCM_ENSg\mean\');
            dirs.griddir = strcat('D:\Data\Model\ROMS\nwp_1_20\backup_surf\'); % % where grid data files are            
            tmp.abb="RCM-ENS";
            
            
            for gridi=1:length(RCM_grid.gridname)
                RCM_grid.(['filename_', RCM_grid.gridname{gridi}])=[dirs.griddir, 'NWP_pck_ocean_', RCM_grid.gridname{gridi}, '_NWP.nc'];
            end
            
        %     tmp.testname=testname;
%             RCM_info.years=inputyear1;
%             RCM_info.months=inputmonth;

            run(tmp.param_script);
            
            if (exist(strcat(dirs.enssavedir) , 'dir') ~= 7)
                mkdir(strcat(dirs.enssavedir));
            end 

            if flags.fig_switch(1) > 0
                flags.fig_tmp = flags.fig_switch(1);
                SSH_2p_1st_013_sub_001_current_plot;
            end
            
            if flags.fig_switch(3) > 0
                flags.fig_tmp = flags.fig_switch(3);
                SSH_2p_1st_013_sub_003_surface_variable_plot;
            end
        end

    end
end

