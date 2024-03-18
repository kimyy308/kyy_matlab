close all; clear all;  clc;
% %  get cmemsstructed SSH. compare model and reSSH. save. 

% RCM_info.name = {'prob_ens2203'}; 
% RCM_info.name = {'ens2203'}; 

% RCM_info.name = {'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
RCM_info.name = {'test2127', 'test2128', 'test2129', 'test2130'};
% RCM_info.name = {'test2127', 'test2130'};
% RCM_info.name = {'t17s_t27c', 't18s_t28c', 't19s_t29c', 't20s_t30c', 't21s_t31c'};

% RCM_info.name = {'test2128'};

% RCM_info.all_regions ={'ES_KHOA'};
% RCM_info.all_regions ={'pollock_egg', 'pollock_egg2'};
% RCM_info.all_regions ={'ES'};
% RCM_info.all_regions ={'ES_KHOA','YS_KHOA', 'SS_KHOA'};
RCM_info.all_regions ={'pollock_egg3'};

for testnameind=1:length(RCM_info.name)
    for regionind=1:length(RCM_info.all_regions)
        clearvars '*' -except regionind testnameind RCM_info
        % % % 
%         tmp.dropboxpath='C:/users/user/Dropbox/';
%         addpath(genpath('C:/Users/User/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/MICT_pollack/2022_future_pollock/subroutine/'))
%         addpath(genpath('C:/Users/User/Dropbox/source/matlab/function/'))

        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
        tmp.fs=filesep;
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model', tmp.fs, ...
            'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, ...
            'MICT_pollack', tmp.fs, '2022_future_pollock', tmp.fs, 'subroutine', tmp.fs]))
        [tmp.dropboxpath, error_status] = Func_0008_set_dropbox_path(computer);
        
        dl=1/10;
        
        [cmaps.byrmap3, tmp.error_status] = Func_0009_get_colormaps('byr3', tmp.dropboxpath);
        [cmaps.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr2', tmp.dropboxpath);   
        [cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);        
        [cmaps.yrmap3, tmp.error_status] = Func_0009_get_colormaps('yr3', tmp.dropboxpath);        
        [cmaps.bymap3, tmp.error_status] = Func_0009_get_colormaps('by3', tmp.dropboxpath);        
        [cmaps.bwr_10, tmp.error_status] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        [cmaps.bw_10, tmp.error_status] = Func_0009_get_colormaps('bw_10', tmp.dropboxpath);

        % for snu_desktop
        tmp.testname_ssp=RCM_info.name{testnameind};
        [tmp.testname_his, tmp.error_status] = Func_0023_RCM_CMIP6_testname_his(tmp.testname_ssp);
        disp(['testname_ssp is ', tmp.testname_ssp])
        disp(['testname_his is ', tmp.testname_his])
        
        RCM_info.years_ssp = [2081:2100]; % % put year which you want to plot [year year ...]
%         RCM_info.years = [1995:2014];
%         RCM_info.years = [2000:2009]; % % put year which you want to plot [year year ...]
%         RCM_info.years = [2010:2019]; % % put year which you want to plot [year year ...]
%         RCM_info.years_ssp = [0]; % % put year which you want to plot [year year ...]

        RCM_info.years_his = [1995:2014];
%                 RCM_info.years_his = [2081:2100];
%         RCM_info.years_his = [0];

        RCM_info.days = [1:1000];
%         refyear = [1983];
%         refyear =[1995:2014];
%         refyear =[2081:2100];
        
%         allyear=[1995:2014];
%           allyear=[2081:2100];
%         allyear = [refyear, RCM_info.years];

%         RCM_info.months = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        RCM_info.season='JF-';
        [RCM_info.months,tmp.error_status] = Func_0019_get_month_from_season(RCM_info.season);
        %         RCM_info.months = [1,2,3,12]; % % put month which you want to plot [month month ...]
        RCM_info.checktime=[15,30];
        tmp.varname ='zeta';
        tmp.variable='SST';
        tmp.regionname=RCM_info.all_regions{regionind};
%         run('nwp_polygon_point.m');
        
% % %         switch region
        [tmp.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

%         tmp.param_script =['C:/users/user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/fig_param/fig_param2_kyy_', tmp.regionname, '.m'];
        tmp.param_script =[tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model', tmp.fs, ...
            'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, ...
            'fig_param', tmp.fs, 'fig_param2_kyy_', tmp.regionname, '.m'];

        dirs.pollock_fdir='/Users/kimyy/Desktop/backup/Research/Ph_D_course/2022_pollock_future/';
%         dirs.pollock_fdir='D:/Research/Ph_D_course/2022_pollock_future/';
        dirs.figrawdir =strcat(dirs.pollock_fdir, 'figure', tmp.fs, tmp.testname_ssp, tmp.fs); % % where figure files will be saved

        dirs.nwproot='/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/';
%         dirs.nwproot='D:/Data/Model/ROMS/nwp_1_20/';
        dirs.filedir_ssp = strcat(dirs.nwproot, tmp.testname_ssp, tmp.fs, 'pollock', tmp.fs); % % where data files are
        dirs.savedir_ssp = strcat(dirs.nwproot, tmp.testname_ssp, tmp.fs, 'pollock', tmp.fs);
        
        dirs.filedir_his = strcat(dirs.nwproot, tmp.testname_his, tmp.fs, 'pollock', tmp.fs); % % where data files are
        dirs.savedir_his = strcat(dirs.nwproot, tmp.testname_his, tmp.fs, 'pollock', tmp.fs);
        
%         dirs.figrawdir_allmean = strcat(dirs.pollock_fdir, 'figure', tmp.fs, 'allmean', tmp.fs); % % where figure files will be saved
        dirs.figrawdir_allmean = strcat(dirs.pollock_fdir, 'figure', tmp.fs, 'allmean4', tmp.fs); % % where figure files will be saved
        
%         inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
%         mkdir(savedir);
        tmp.LTRANS_testname='pollock1';
        
        % % %         flag configuration

%             for flagi=1:100
%                 flags.fig_switch(flagi)=0;
%             end

            flags.fig_switch(1)=0;  % 'probability plot for pollock spwaning ground, compared to total period';
            flags.fig_switch(2)=0;  %'# of individual plot at each grid (elapsed ?? days)';
            flags.fig_switch(3)=0;  %'spawning ground time series (regime)';
            flags.fig_switch(4)=0;  %'plot ensemble spawning probability at certain time (from std of RCMs)';
            flags.fig_switch(5)=0;  %'plot SSPR mean of all test';
            flags.fig_switch(6)=2;  %'plot # of individual plot of all test';
            flags.fig_switch(7)=0;  %'plot intermodel std of individual plot of all test';
            flags.fig_switch(8)=0;  %'plot intermodel std of SSPR of all test';
            flags.fig_switch(9)=0;  %'calculation of longitudinal/latitudinal movement of particles of all test';

%             flags.fig_switch(1)=2;  % 'probability plot for pollock spwaning ground';
%             flags.fig_switch(2)=1;  % 'probability plot for pollock location (elapsed ?? days)';
% % %             flags.fig_switch(3)=1;  % 'temp diff plot between 80s and later';
%             flags.fig_switch(4)=1;  % 'vec diff plot between 80s and later';
%             flags.fig_switch(5)=1;  % 'egg probability (elapsed ?? days) diff plot between 80s and later';
% % %             flags.fig_switch(6)=1;  % 'spawning ground diff plot between 80s and later';
%             flags.fig_switch(7)=1;  % 'wind diff plot between 80s and later';
%             flags.fig_switch(8)=0;  % 'tair diff plot between 80s and later';
%             flags.fig_switch(9)=0;  % 'shflux diff plot between 80s and later';
%             flags.fig_switch(10)=0;  % 'relative vorticity plot';
%             flags.fig_switch(11)=0;  % 'wind stress curl plot';
%             flags.fig_switch(12)=0;  % 'vorticity diff plot between 80s and later';
%             flags.fig_switch(13)=0;  % 'wind stress curl diff plot between 80s and later';
%             flags.fig_switch(14)=1;  % 'vec plot';
% %             flags.fig_switch(15)=1;  % 'probability plot for pollock spwaning ground (hovmoller)';
%             flags.fig_switch(16)=1;  % 'wind plot';
%             flags.fig_switch(17)=0;  %'test';
%             flags.fig_switch(18)=0;  %'Southern East Korean Bay nwv & numegg(15d) time series';
%             flags.fig_switch(19)=0;  %'Around Vladivostok nev & numegg(15d) time series'
%             flags.fig_switch(20)=0;  %'Southern East Korean Bay nwv & easterly time series';
%             flags.fig_switch(21)=0;  %'Around Vladivostok nev & southerly time series'
%             flags.fig_switch(22)=0;  %'Southern East Korean Bay2 nwv & numegg(15d) time series';
%             flags.fig_switch(23)=0;  %'Southern East Korean Bay2 nwv & easterly time series';
%             flags.fig_switch(24)=0;  %'Southern East Korean Bay2 nwv & East Korean Bay curl time series';
%             flags.fig_switch(25)=0;  %'Southern East Korean Bay2 nwv & East Korean Bay curl(all) time series'
%             flags.fig_switch(26)=0;  %'Around Vladivostok2 nev & numegg(15d) time series';
%             flags.fig_switch(27)=0;  %'Around Vladivostok2 nev & southerly time series';
%             flags.fig_switch(28)=0;  %'Around Vladivostok southerly & AOI time series';
%             flags.fig_switch(29)=0;  %'Around Southern East Korean Bay2 easterly & AOI time series';
%             flags.fig_switch(30)=0;  %'pair diff plot between 80s and later';
%             flags.fig_switch(31)=2;  %'SST plot';
%             flags.fig_switch(32)=1;  %'OISST plot';
%             flags.fig_switch(33)=0;  %'CMEMS geostrophic vec plot';
%             flags.fig_switch(34)=0;  %'Southern East Korean Bay2 spawning ground time series (regime)';
%             flags.fig_switch(35)=0;  %'Around Vladivostok spawning ground time series (regime)';
%             flags.fig_switch(36)=0;  %'Southern East Korean Bay2 numegg(??d) time series (regime)';
%             flags.fig_switch(37)=0;  %'Around Vladivostok numegg(15d) time series (regime)';
%             flags.fig_switch(38)=0;  %'Around Southern East Korean Bay2 numegg(15d) & AOI time series';
%             flags.fig_switch(39)=0;  %'Northern East Sea spawning ground time series (regime)';
%             flags.fig_switch(40)=0;  %'Northern East Sea numegg(15d) time series (regime)';
%             flags.fig_switch(41)=0;  %'Russian EEZ spawning ground time series (regime)';
%             flags.fig_switch(42)=0;  %'Russian EEZ numegg(15d) time series (regime)';
%             flags.fig_switch(43)=0;  %'each EEZ spawning ground time series (regime)';
%             flags.fig_switch(44)=0;  %'each EEZ numegg(??d) time series (regime)';
%             flags.fig_switch(45)=0;  %'each EEZ numegg(30d) time series (regime)';
%             flags.fig_switch(46)=0;  %'EEZ plot';
%             flags.fig_switch(47)=0;  %'lon lat distance time series (not completed)';
%             flags.fig_switch(48)=0;  % 'egg probability (elapsed 30 days) diff plot between 80s and later';
% % %             flags.fig_switch(49)=1;  %'Southern Korea EEZ spawning ground time series (regime)';
% % %             flags.fig_switch(50)=1;  %'Southern Korea EEZ numegg(??d) time series (regime)';
%             flags.fig_switch(51)=0;  %'Southern Korea EEZ numegg(30d) time series (regime)';
%             flags.fig_switch(52)=0;  %'probability plot for pollock location (elapsed 30 days)';
%             flags.fig_switch(53)=1;  %'particle moved distance-lat time series (regime)';
% % %             flags.fig_switch(54)=1;  %'latitudinal particle density (horizontal bar)';
%             flags.fig_switch(55)=1;  %'latitudinal particle density (elapsed ??days, horizontal bar)';
%             flags.fig_switch(56)=1;  %'SST RMSE plot';
%             flags.fig_switch(57)=1;  %'latitudinal particle density difference (horizontal bar)';
%             flags.fig_switch(58)=1;  %'latitudinal particle density difference (elapsed ??days, horizontal bar)';
%             flags.fig_switch(59)=1;  %'spawning ground, all period, dot plot';
%             flags.fig_switch(60)=1;  %'AO Index plot (regime)';
%             flags.fig_switch(61)=1;  %'latitudinal particle density difference (elapsed ??days, horizontal bar)';
%             flags.fig_switch(62)=1;  %'vec + isotherm(2,5) plot';
%             flags.fig_switch(63)=1;  %'particle latitude boxplot (elapsed ?? days)';
%             flags.fig_switch(64)=1;  %'EAWMI plot (regime)';
% % %             flags.fig_switch(65)=2;  %'# of individual plot at each grid (elapsed ?? days)';
%             flags.fig_switch(66)=2;  %'# of individual (elapsed ?? days) diff plot between 80s and later';
% % %             flags.fig_switch(67)=1;  %'probability plot for pollock spwaning ground (~41N)';
%             flags.fig_switch(68)=0;  %'particle moved distance-lon time series (regime)';
%             flags.fig_switch(69)=1;  %'coastal area (<500m) numegg(??d) time series (regime)';
%             flags.fig_switch(70)=1;  %'southern coastal area (<500m, 39N) numegg(??d) time series (regime)';
%         
%         for flagi=1:100
%             flags.fig_switch(flagi)=0;
%         end
%         flags.fig_switch(63)=2; 
%         flags.fig_switch(43)=2; 

        


        dirs.figdir=[dirs.figrawdir,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/spawn/'];
        if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.figdir));
        end 
        
        dirs.spdiff_figdir=[dirs.figrawdir,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/spawn_diff/'];
        if (exist(strcat(dirs.spdiff_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.spdiff_figdir));
        end 
        
        dirs.climfigdir=[dirs.figrawdir,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/clim/'];
        if (exist(strcat(dirs.climfigdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.climfigdir));
        end 
        
        dirs.clim_atfigdir=[dirs.figrawdir,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/clim_at/'];
        if (exist(strcat(dirs.clim_atfigdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.clim_atfigdir));
        end 
        
        dirs.clim_diff_figdir=[dirs.figrawdir,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/clim_diff/'];
        if (exist(strcat(dirs.clim_diff_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.clim_diff_figdir));
        end 
       
        dirs.clim_atdiff_figdir=[dirs.figrawdir,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/clim_atdiff/'];
        if (exist(strcat(dirs.clim_atdiff_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.clim_atdiff_figdir));
        end 
        
        dirs.all_figdir=[dirs.figrawdir,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/all/'];
        if (exist(strcat(dirs.all_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.all_figdir));
        end 
        
        dirs.regime_figdir=[dirs.figrawdir,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/regime_shift/'];
        if (exist(strcat(dirs.regime_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.regime_figdir));
        end         
        
        dirs.figdir_allmean=[dirs.figrawdir_allmean,tmp.LTRANS_testname, tmp.fs, tmp.regionname, '/spawn/'];
        if (exist(strcat(dirs.figdir_allmean) , 'dir') ~= 7)
            mkdir(strcat(dirs.figdir_allmean));
        end 
        
        
        
% % %         run
        for flagi=1:9
            if flags.fig_switch(flagi)>=1
                run(['future_pollock_003_subroutine_', num2str(flagi, '%03i')]);
            end
        end
        
        
        
% % % % %       'probability plot for pollock spwaning ground'
% %         fig_flag=flags.fig_switch(1);
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_spawning_ground_probability_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             future_pollock_003_subroutine_001;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'probability plot for pollock location (elapsed ?? days)'
% %         fig_flag=flags.fig_switch(2);
% %         while (fig_flag)
% %              for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_loc_prob_', num2str(tmp.checktime, '%02i'),'days', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 future_pollock_003_subroutine_002;
% %              end
% %             fig_flag=0;
% %         end
% % 
% % % % %       'temp diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(3);
% %         while (fig_flag)
% %             jpgname=strcat(clim_diff_outfile, '_', tmp.testname,'_',tmp.regionname, '_temp_diff_fut_his_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 RCM_info.season, '.tif'); %% ~_year_month.jpg
% %             future_pollock_003_subroutine_003;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'vec diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(4);
% %         fig_name=flags.fig_switch(4,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_diff_outfile, '_', tmp.testname,'_',tmp.regionname, '_vec_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_004;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'egg probability (elapsed 15 days) diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(5);
% %         fig_name=flags.fig_switch(5,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(spdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_egg_prob_', num2str(tmp.checktime, '%02i'), 'd_diff_80s_', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_005;
% %             end
% %             fig_flag=0;
% %         end
% % 
% % % % %       'spwaning ground diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(6);
% %         fig_name=flags.fig_switch(6,1};
% %         while (fig_flag)
% %             jpgname=strcat(spdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_spawn_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             future_pollock_003_subroutine_006;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'wind diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(7);
% %         fig_name=flags.fig_switch(7,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_atdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_wind_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_007;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'tair diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(8);
% %         fig_name=flags.fig_switch(8,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_atdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_tair_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_008;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'shflux diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(9);
% %         fig_name=flags.fig_switch(9,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_atdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_shflux_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_009;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'vort plot'
% %         fig_flag=flags.fig_switch(10);
% %         fig_name=flags.fig_switch(10,1};
% %         while (fig_flag)
% %             jpgname=strcat(climoutfile, '_', tmp.testname,'_',tmp.regionname, '_vort_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_010;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'wind stress curl plot'
% %         fig_flag=flags.fig_switch(11);
% %         fig_name=flags.fig_switch(11,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_at_outfile, '_', tmp.testname,'_',tmp.regionname, '_wind_stress_curl_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             future_pollock_003_subroutine_011;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'vort diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(12);
% %         fig_name=flags.fig_switch(12,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_diff_outfile, '_', tmp.testname,'_',tmp.regionname, '_vort_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_012;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'wstr curl diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(13);
% %         fig_name=flags.fig_switch(13,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_atdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_wstr_curl_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_013;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'vec plot'
% %         fig_flag=flags.fig_switch(14);
% %         fig_name=flags.fig_switch(14,1};
% %         while (fig_flag)
% %             jpgname=strcat(climoutfile, '_', tmp.testname,'_',tmp.regionname, '_vec_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_014;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'probability plot for pollock spwaning ground (hovmoller)'
% %         fig_flag=flags.fig_switch(15);
% %         fig_name=flags.fig_switch(15,1};
% %         while (fig_flag)
% %             jpgname=strcat(all_outfile, '_', tmp.testname,'_',tmp.regionname, '_sp_ground_prob_hovm_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_015;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'wind plot'
% %         fig_flag=flags.fig_switch(16);
% %         fig_name=flags.fig_switch(16,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_at_outfile, '_', tmp.testname,'_',tmp.regionname, '_wind_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_016;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'test'
% %         fig_flag=flags.fig_switch(17);
% %         fig_name=flags.fig_switch(17,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_at_outfile, '_', tmp.testname,'_',tmp.regionname, '_wsc_ts_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_017;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'Southern East Korean Bay nwv & numegg(15d) time series'
% %         fig_flag=flags.fig_switch(18);
% %         fig_name=flags.fig_switch(18,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB_ts_egg_density_15d_nwv_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_018;
% %             fig_flag=0;
% %         end       
% %   
% % % % %       'Around Vladivostok nev & numegg(15d) time series'
% %         fig_flag=flags.fig_switch(19);
% %         fig_name=flags.fig_switch(19,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SVLA_ts_egg_density_15d_nev_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_019;
% %             fig_flag=0;
% %         end   
% % 
% % % % %       'Southern East Korean Bay nwv & easterly time series'
% %         fig_flag=flags.fig_switch(20);
% %         fig_name=flags.fig_switch(20,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB_ts_easterly_nwv_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_020;
% %             fig_flag=0;
% %         end    
% % 
% % % % %       'Around Vladivostok nev & southerly time series'
% %         fig_flag=flags.fig_switch(21);
% %         fig_name=flags.fig_switch(21,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SVLA_ts_southerly_nev_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_021;
% %             fig_flag=0;
% %         end  
% % 
% % % % %       'Southern East Korean Bay2 nwv & numegg(15d) time series'
% %         fig_flag=flags.fig_switch(22);
% %         fig_name=flags.fig_switch(22,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_ts_egg_density_15d_nwv_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_022;
% %             fig_flag=0;
% %         end    
% % 
% % % % %       'Southern East Korean Bay2 nwv & easterly time series'
% %         fig_flag=flags.fig_switch(23);
% %         fig_name=flags.fig_switch(23,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_ts_easterly_nwv_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_023;
% %             fig_flag=0;
% %         end   
% % 
% % % % %       'Southern East Korean Bay2 nwv & East Korean Bay curl time series'
% %         fig_flag=flags.fig_switch(24);
% %         fig_name=flags.fig_switch(24,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_ts_wstr_curl_nwv_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_024;
% %             fig_flag=0;
% %         end   
% % 
% % % % %       'Southern East Korean Bay2 nwv & East Korean Bay curl(all) time series'
% %         fig_flag=flags.fig_switch(25);
% %         fig_name=flags.fig_switch(25,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_ts_allwstr_curl_nwv_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_025;
% %             fig_flag=0;
% %         end 
% %         
% % % % %       'Around Vladivostok2 nev & numegg(15d) time series'
% %         fig_flag=flags.fig_switch(26);
% %         fig_name=flags.fig_switch(26,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SVLA2_ts_egg_density_15d_nev_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_026;
% %             fig_flag=0;
% %         end   
% %         
% % % % %       'Around Vladivostok2 nev & southerly time series'
% %         fig_flag=flags.fig_switch(27);
% %         fig_name=flags.fig_switch(27,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SVLA2_ts_southerly_nev_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_027;
% %             fig_flag=0;
% %         end  
% %  
% % % % %       'Around Vladivostok southerly & AOI time series'
% %         fig_flag=flags.fig_switch(28);
% %         fig_name=flags.fig_switch(28,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SVLA_ts_southerly_AOI_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_028;
% %             fig_flag=0;
% %         end  
% % 
% % % % %       'Around Southern East Korean Bay2 easterly & AOI time series'
% %         fig_flag=flags.fig_switch(29);
% %         fig_name=flags.fig_switch(29,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_ts_easterly_AOI_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_029;
% %             fig_flag=0;
% %         end  
% % 
% % % % %       'pair diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(30);
% %         fig_name=flags.fig_switch(30,1};
% %         while (fig_flag)
% %             jpgname=strcat(clim_atdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_pair_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_030;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'SST plot'
% %         fig_flag=flags.fig_switch(31);
% %         fig_name=flags.fig_switch(31,1};
% %         while (fig_flag)
% %             jpgname=strcat(climoutfile, '_', tmp.testname,'_',tmp.regionname, '_SST_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_031;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'OISST plot, from 003 code'
% %         fig_flag=flags.fig_switch(32);
% %         fig_name=flags.fig_switch(32,1};
% %         while (fig_flag)
% %             jpgname=strcat(climoutfile, '_', tmp.testname,'_',tmp.regionname, '_OISST_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_032;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'CMEMS geostrophic vec plot, from CMEMS data (DAMO)'
% %         fig_flag=flags.fig_switch(33);
% %         fig_name=flags.fig_switch(33,1};
% %         while (fig_flag)
% %             jpgname=strcat(climoutfile, '_', tmp.testname,'_',tmp.regionname, '_cmems_vec_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_033;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'Southern East Korean Bay2 spawning ground time series (regime)';
% %         fig_flag=flags.fig_switch(34);
% %         fig_name=flags.fig_switch(34,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_regime_ts_sp_ground_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_034;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'Around Vladivostok spawning ground time series (regime)';
% %         fig_flag=flags.fig_switch(35);
% %         fig_name=flags.fig_switch(35,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_SVLA_regime_ts_sp_ground_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_035;
% %             fig_flag=0;
% %         end    
% % 
% % % % %       'Southern East Korean Bay2 numegg(15d) time series (regime)';
% %         fig_flag=flags.fig_switch(36);
% %         fig_name=flags.fig_switch(36,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_regime_ts_egg_', num2str(tmp.checktime, '%02i'),'d_', ...
% %                     num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_036;
% %             end
% %             fig_flag=0;
% %         end    
% % 
% % % % %       'Around Vladivostok numegg(15d) time series (regime)';
% %         fig_flag=flags.fig_switch(37);
% %         fig_name=flags.fig_switch(37,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_SVLA_regime_ts_egg_15d_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_037;
% %             fig_flag=0;
% %         end   
% % 
% % % % %       'Around Southern East Korean Bay2 numegg(15d) & AOI time series'
% %         fig_flag=flags.fig_switch(38);
% %         fig_name=flags.fig_switch(38,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_ts_egg_15d_AOI_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_038;
% %             fig_flag=0;
% %         end  
% % 
% % % % %       'Northern East Sea spawning ground time series (regime)';
% %         fig_flag=flags.fig_switch(39);
% %         fig_name=flags.fig_switch(39,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_NES_regime_ts_sp_ground_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_039;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'Northern East Sea numegg(15d) time series (regime)';
% %         fig_flag=flags.fig_switch(40);
% %         fig_name=flags.fig_switch(40,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_NES_regime_ts_egg_15d_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_040;
% %             fig_flag=0;
% %         end    
% % 
% % % % %       'Russian EEZ spawning ground time series (regime)';
% %         fig_flag=flags.fig_switch(41);
% %         fig_name=flags.fig_switch(41,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_REEZ_regime_ts_sp_ground_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_041;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'Russian EEZ numegg(15d) time series (regime)';
% %         fig_flag=flags.fig_switch(42);
% %         fig_name=flags.fig_switch(42,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_REEZ_regime_ts_egg_15d_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_042;
% %             fig_flag=0;
% %         end    
% % 
% % % % %       'each EEZ spawning ground time series (regime)';
% %         fig_flag=flags.fig_switch(43);
% %         fig_name=flags.fig_switch(43,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_Each_EEZ_regime_ts_sp_ground_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             future_pollock_003_subroutine_043;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'each EEZ numegg(15d) time series (regime)';
% %         fig_flag=flags.fig_switch(44);
% %         fig_name=flags.fig_switch(44,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_Each_EEZ_regime_ts_egg_', num2str(tmp.checktime, '%02i'),'d_', ...
% %                     num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_044;
% %             end
% %             fig_flag=0;
% %         end
% %         
% % % % %       'each EEZ numegg(30d) time series (regime)';
% %         fig_flag=flags.fig_switch(45);
% %         fig_name=flags.fig_switch(45,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_Each_EEZ_regime_ts_egg_30d_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_045;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'EEZ plot'
% %         fig_flag=flags.fig_switch(46);
% %         fig_name=flags.fig_switch(46,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_EEZ', '.tif'); %% ~_year_month.jpg            
% %             pollock_paper01_002_subroutine_046;
% %             fig_flag=0;
% %         end    
% % 
% % % % %       'lon lat distance timeseries'
% %         fig_flag=flags.fig_switch(47);
% %         fig_name=flags.fig_switch(47,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_SEKB2_ts_egg_density_15d_nwv_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_047;
% %             fig_flag=0;
% %         end            
% % 
% % % % %       'egg probability (elapsed 30 days) diff plot between 80s and later'
% %         fig_flag=flags.fig_switch(48);
% %         fig_name=flags.fig_switch(48,1};
% %         while (fig_flag)
% %             jpgname=strcat(spdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_egg_30d_prob_diff_80s_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_048;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'Southern Korea EEZ spawning ground time series (regime)';
% %         fig_flag=flags.fig_switch(49);
% %         fig_name=flags.fig_switch(49,1};
% %         while (fig_flag)
% % %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_SKEEZ_regime_ts_sp_ground_', ...
% % %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% % %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_ekb2_regime_ts_sp_ground_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             future_pollock_003_subroutine_049;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'Southern Korea EEZ numegg(15d) time series (regime)';
% %         fig_flag=flags.fig_switch(50);
% %         fig_name=flags.fig_switch(50,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_SKEEZ_regime_ts_egg_',num2str(tmp.checktime, '%02i'), 'd_', ...
% %                     num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_050;
% %             end
% %             fig_flag=0;
% %         end    
% %         
% % % % %       'Southern Korea EEZ numegg(30d) time series (regime)';
% %         fig_flag=flags.fig_switch(51);
% %         fig_name=flags.fig_switch(51,1};
% %         while (fig_flag)
% %             jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_SKEEZ_regime_ts_egg_30d_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_051;
% %             fig_flag=0;
% %         end    
% %         
% % % % %       'probability plot for pollock location (elapsed 30 days)'
% %         fig_flag=flags.fig_switch(52);
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_loc_prob_30days', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_052;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'particle moved distance-lat time series (regime)';
% %         fig_flag=flags.fig_switch(53);
% %         fig_name=flags.fig_switch(53,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_regime_ts_mv_lat_', num2str(tmp.checktime, '%02i'), '_', ...
% %                     num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_053;
% %             end
% %             fig_flag=0;
% %         end
% %         
% % % % %       'latitudinal particle density (horizontal bar)';
% %         fig_flag=flags.fig_switch(54);
% %         fig_name=flags.fig_switch(54,1};
% %         while (fig_flag)
% % %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_lat_par_dens_', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 future_pollock_003_subroutine_054;
% % %             end
% %             fig_flag=0;
% %         end
% % 
% % % % %       'latitudinal particle density (elapsed ?? days, horizontal bar)'
% %         fig_flag=flags.fig_switch(55);
% %         while (fig_flag)
% %              for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_lat_par_dens_', num2str(tmp.checktime, '%02i'),'days', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_055;
% %              end
% %             fig_flag=0;
% %         end
% % 
% % % % %       'SST RMSEplot'
% %         fig_flag_num=56;
% %         fig_flag=flags.fig_switch(fig_flag_num);
% %         while (fig_flag)
% %             figvarname='SST_RMSE';
% %             figtypename = 'pcol';
% %             figtimename = [num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y_', ...
% %                             num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm'];
% %             extensionname = '.tif';
% %             scenname = 'historical';
% %             figname = Func_0002_get_figname(figrawdir, figvarname, figtypename, ...
% %                 tmp.regionname, scenname, tmp.testname, figtimename, extensionname);
% % %             jpgname=strcat(climoutfile, '_', tmp.testname,'_',tmp.regionname, '_SST_RMSE_', ...
% % %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% % %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_056;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'latitudinal particle density (horizontal bar)';
% %         fig_flag=flags.fig_switch(57);
% %         fig_name=flags.fig_switch(57,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(spdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_lat_par_dens_diff_', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_057;
% %             end
% %             fig_flag=0;
% %         end
% %         
% % % % %       'latitudinal particle density difference (elapsed ??days, horizontal bar)';
% %         fig_flag=flags.fig_switch(58);
% %         fig_name=flags.fig_switch(58,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(spdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_lat_par_dens_diff_', num2str(tmp.checktime, '%02i'),'d_', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_058;
% %             end
% %             fig_flag=0;
% %         end
% %         
% % % % %       'probability plot for pollock spwaning ground'
% %         fig_flag=flags.fig_switch(59);
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_spawning_ground_all_dot_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_059;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'Around Southern East Korean Bay2 easterly & AOI time series'
% %         fig_flag=flags.fig_switch(60);
% %         fig_name=flags.fig_switch(60,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_AOI_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_060;
% %             fig_flag=0;
% %         end  
% % 
% % % % %       'latitudinal particle density difference (elapsed ??days, horizontal bar, percentage)';
% %         fig_flag=flags.fig_switch(61);
% %         fig_name=flags.fig_switch(61,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(spdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_lat_par_dens_diff_percentage_', num2str(tmp.checktime, '%02i'),'d_', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_061;
% %             end
% %             fig_flag=0;
% %         end
% %         
% % % % %       'vec + isotherm (2,5) plot'
% %         fig_flag=flags.fig_switch(62);
% %         fig_name=flags.fig_switch(62,1};
% %         while (fig_flag)
% %             jpgname=strcat(climoutfile, '_', tmp.testname,'_',tmp.regionname, '_vec_isotherm_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_062;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'particle latitude boxplot (elapsed ?? days)';
% %         fig_flag=flags.fig_switch(63);
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_lat_par_box_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_063;
% %             fig_flag=0;
% %         end
% % 
% % % % %       'EAWMI time series'
% %         fig_flag=flags.fig_switch(64);
% %         fig_name=flags.fig_switch(64,1};
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_EAWMI_', ...
% %                 num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_064;
% %             fig_flag=0;
% %         end  
% %         
% %         
% % % % %       '# of individual plot at each grid (elapsed ?? days)'
% %         fig_flag=flags.fig_switch(65);
% %         while (fig_flag)
% %              for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 future_pollock_003_subroutine_065;
% %              end
% %             fig_flag=0;
% %         end
% %         
% %         
% % % % %       '# of individual (elapsed ?? days) diff plot between 80s and later';
% %         fig_flag=flags.fig_switch(66);
% %         fig_name=flags.fig_switch(66,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(spdiff_outfile, '_', tmp.testname,'_',tmp.regionname, '_egg_num_', num2str(tmp.checktime, '%02i'), 'd_diff_80s_', ...
% %                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 future_pollock_003_subroutine_066;
% %             end
% %             fig_flag=0;
% %         end
% %         
% % % % %       'probability plot for pollock spwaning ground (~41N)'
% %         fig_flag=flags.fig_switch(67);
% %         while (fig_flag)
% %             jpgname=strcat(outfile, '_', tmp.testname,'_',tmp.regionname, '_spawning_ground_probability_', ...
% %                 num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), 'y', ...
% %                 num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %             pollock_paper01_002_subroutine_067;
% %             fig_flag=0;
% %         end
% %         
% % % % %       'particle moved distance-lon time series (regime)';
% %         fig_flag=flags.fig_switch(68);
% %         fig_name=flags.fig_switch(68,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_regime_ts_mv_lon_', num2str(tmp.checktime, '%02i'), '_', ...
% %                     num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_068;
% %             end
% %             fig_flag=0;
% %         end   
% %         
% % % % %       'coastal area numegg(??d) time series (regime)';
% %         fig_flag=flags.fig_switch(69);
% %         fig_name=flags.fig_switch(69,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_coast_regime_ts_egg_',num2str(tmp.checktime, '%02i'), 'd_', ...
% %                     num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 future_pollock_003_subroutine_069;
% %             end
% %             fig_flag=0;
% %         end    
% %         
% % % % %       'Southern Korea coastal area numegg(??d) time series (regime)';
% %         fig_flag=flags.fig_switch(70);
% %         fig_name=flags.fig_switch(70,1};
% %         while (fig_flag)
% %             for checkti=1:length(RCM_info.checktime)
% %                 tmp.checktime=RCM_info.checktime(checkti);
% %                 jpgname=strcat(regime_outfile, '_', tmp.testname,'_',tmp.regionname, '_southern_coast_regime_ts_egg_',num2str(tmp.checktime, '%02i'), 'd_', ...
% %                     num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
% %                     num2str(min(RCM_info.months),'%02i'),'_',num2str(max(RCM_info.months),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
% %                 pollock_paper01_002_subroutine_070;
% %             end
% %             fig_flag=0;
% %         end    
        
    end
end



