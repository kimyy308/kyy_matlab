close all; clear all;  clc;
% %  get cmemsstructed SSH. compare model and reSSH. save. 

% RCM_info.name = {'prob_ens2203'}; 
% RCM_info.name = {'ens2203'}; 

% RCM_info.name = {'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
% RCM_info.name = {'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
% RCM_info.name = {'t17s_t27c', 't18s_t28c', 't19s_t29c', 't20s_t30c', 't21s_t31c'};
RCM_info.name = {'t27s_t17c', 't28s_t18c', 't29s_t19c', 't30s_t20c', 't31s_t21c'};

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
        tmp.dropboxpath='C:\users\user/Dropbox/';

        addpath(genpath('C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\MICT_pollack\2022_future_pollock\subroutine\'))
        addpath(genpath('C:\Users\User\Dropbox\source\matlab\function\'))
        [tmp.dropboxpath, error_status] = Func_0008_set_dropbox_path(computer);
        
        dl=1/10;
        
        [cmaps.byrmap3, tmp.error_status] = Func_0009_get_colormaps('byr3', tmp.dropboxpath);
        [cmaps.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr2', tmp.dropboxpath);        
        [cmaps.yrmap3, tmp.error_status] = Func_0009_get_colormaps('yr3', tmp.dropboxpath);        
        [cmaps.bymap3, tmp.error_status] = Func_0009_get_colormaps('by3', tmp.dropboxpath);        

        % for snu_desktop
        tmp.testname_ssp=RCM_info.name{testnameind};
        [tmp.testname_his, tmp.error_status] = Func_0023_RCM_CMIP6_testname_his(tmp.testname_ssp);
        disp(['testname_ssp is ', tmp.testname_ssp])
        disp(['testname_his is ', tmp.testname_his])
        
        RCM_info.years_ssp = [1995:2014]; % % put year which you want to plot [year year ...]
%         RCM_info.years = [1995:2014];
%         RCM_info.years = [2000:2009]; % % put year which you want to plot [year year ...]
%         RCM_info.years = [2010:2019]; % % put year which you want to plot [year year ...]
%         RCM_info.years_ssp = [0]; % % put year which you want to plot [year year ...]

%         RCM_info.years_his = [1995:2014];
        RCM_info.years_his = [2081:2100];
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

        tmp.param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/fig_param/fig_param2_kyy_', tmp.regionname, '.m'];
        dirs.figrawdir =strcat('D:\Research\Ph_D_course\2022_pollock_future\figure\', tmp.testname_ssp, '\'); % % where figure files will be saved
        dirs.filedir_ssp = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname_ssp, '\pollock\'); % % where data files are
        dirs.savedir_ssp = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname_ssp, '\pollock\');
        
        dirs.filedir_his = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname_his, '\pollock\'); % % where data files are
        dirs.savedir_his = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname_his, '\pollock\');
        
        dirs.figrawdir_allmean = strcat('D:\Research\Ph_D_course\2022_pollock_future\figure\', 'allmean', '\'); % % where figure files will be saved
        
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
            flags.fig_switch(7)=2;  %'plot intermodel std of individual plot of all test';
            flags.fig_switch(8)=0;  %'plot intermodel std of SSPR of all test';

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
% % %             flags.fig_switch(66)=1;  %'# of individual (elapsed ?? days) diff plot between 80s and later';
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

        
        dirs.figdir=[dirs.figrawdir,tmp.LTRANS_testname, '\', tmp.regionname, '\spawn\'];
        if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.figdir));
        end 
        
        dirs.spdiff_figdir=[dirs.figrawdir,tmp.LTRANS_testname, '\', tmp.regionname, '\spawn_diff\'];
        if (exist(strcat(dirs.spdiff_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.spdiff_figdir));
        end 
        
        dirs.climfigdir=[dirs.figrawdir,tmp.LTRANS_testname, '\', tmp.regionname, '\clim\'];
        if (exist(strcat(dirs.climfigdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.climfigdir));
        end 
        
        dirs.clim_atfigdir=[dirs.figrawdir,tmp.LTRANS_testname, '\', tmp.regionname, '\clim_at\'];
        if (exist(strcat(dirs.clim_atfigdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.clim_atfigdir));
        end 
        
        dirs.clim_diff_figdir=[dirs.figrawdir,tmp.LTRANS_testname, '\', tmp.regionname, '\clim_diff\'];
        if (exist(strcat(dirs.clim_diff_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.clim_diff_figdir));
        end 
        
        dirs.clim_atdiff_figdir=[dirs.figrawdir,tmp.LTRANS_testname, '\', tmp.regionname, '\clim_atdiff\'];
        if (exist(strcat(dirs.clim_atdiff_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.clim_atdiff_figdir));
        end 
        
        dirs.all_figdir=[dirs.figrawdir,tmp.LTRANS_testname, '\', tmp.regionname, '\all\'];
        if (exist(strcat(dirs.all_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.all_figdir));
        end 
        
        dirs.regime_figdir=[dirs.figrawdir,tmp.LTRANS_testname, '\', tmp.regionname, '\regime_shift\'];
        if (exist(strcat(dirs.regime_figdir) , 'dir') ~= 7)
            mkdir(strcat(dirs.regime_figdir));
        end         
        
        dirs.figdir_allmean=[dirs.figrawdir_allmean,tmp.LTRANS_testname, '\', tmp.regionname, '\spawn_t27s_t17c\'];
        if (exist(strcat(dirs.figdir_allmean) , 'dir') ~= 7)
            mkdir(strcat(dirs.figdir_allmean));
        end 
        
        
        
% % %         run
        for flagi=1:8
            if flags.fig_switch(flagi)>=1
                run(['future_pollock_003_subroutine_', num2str(flagi, '%03i')]);
            end
        end
        
      
        
    end
end



