close all; clear all;  clc;
% %  get cmemsstructed SSH. compare model and reSSH. save. 


all_testname = {'test06'};

all_region ={'pollock_egg3'};
% all_region ={'pollock_egg', 'pollock_egg2'};
% all_region ={'ES'};
% all_region ={'ES_KHOA','YS_KHOA', 'SS_KHOA'};


for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname
        % % % 
        dropboxpath='C:\users\user/Dropbox/';
%         addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
%         addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%         addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%         addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
%         addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/mat_tool']));
%         addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
        addpath(genpath('C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\MICT_pollack\paper\subroutine\'))
        addpath(genpath('C:\Users\User\Dropbox\source\matlab\function\'))
        [dropboxpath, erorr_status] = Func_0008_set_dropbox_path(computer);
        shadlev = [-2 2];
        % rms_shadlev = [0 4];
        trend_shadlev = [0 4];
        % bias_shadlev = [-4 4];
        conlev  = 0:5:35;
        dl=1/10;
        
        [byrmap3, error_status] = Func_0009_get_colormaps('byr3', dropboxpath);
        [byrmap, error_status] = Func_0009_get_colormaps('byr2', dropboxpath);        
        yrmap = byrmap(129:256,:);
        
        % for snu_desktop
        testname=all_testname{testnameind} 
%         inputyear = [1983];
%         inputyear = [1983:1987]; % % put year which you want to plot [year year ...]
        inputyear = [1983:1992]; % % put year which you want to plot [year year ...]
%         inputyear = [1987:1990]; % % put year which you want to plot [year year ...]
%         inputyear = [1988:1992]; % % put year which you want to plot [year year ...]
%         inputyear = [1987:1995]; % % put year which you want to plot [year year ...]
%         inputyear = [1987:2000]; % % put year which you want to plot [year year ...]
%         inputyear = [1987:2019]; % % put year which you want to plot [year year ...]
%         inputyear = [1990:1999]; % % put year which you want to plot [year year ...]
%         inputyear = [1993:2016];
%         inputyear = [2000:2009]; % % put year which you want to plot [year year ...]
%         inputyear = [2010:2019]; % % put year which you want to plot [year year ...]
        
%         allyear =[1983];
%         allyear=[1983:1990];
        allyear=[1983:1992];
%         allyear=[1983:1995];
%         allyear=[1983:2000];
%         allyear=[1983:2019];
        
%         refyear = [1983];
        refyear =[1983:1987];

%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1,2]; % % put month which you want to plot [month month ...]
%         inputmonth = [1,2,3,12]; % % put month which you want to plot [month month ...]
        checktime=[15,30];
        varname ='zeta';
        variable='zeta';
        regionname=all_region{regionind};
        run('nwp_polygon_point.m');
        
% % %         switch region
        [error_status, refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);

        param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_', regionname, '.m'];
        figrawdir =strcat('Z:\내 드라이브\research\Ph_D_course\Particle tracking of the walleye pollock\figure\', testname, '\DA\'); % % where figure files will be saved
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
        savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
        inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
        mkdir(savedir);
        LTRANS_testname='Pollock6';
        
        % % %         flag configuration
%         for folding=1:1
%             fig_flags{1,1}='probability plot for pollock spwaning ground';
%             fig_flags{2,1}='probability plot for pollock location (elapsed ?? days)';
%             fig_flags{3,1}='temp diff plot between 80s and later';
%             fig_flags{4,1}='vec diff plot between 80s and later';
%             fig_flags{5,1}='egg probability (elapsed ?? days) diff plot between 80s and later';
%             fig_flags{6,1}='spwaning ground diff plot between 80s and later';
%             fig_flags{7,1}='wind diff plot between 80s and later';
%             fig_flags{8,1}='tair diff plot between 80s and later';
%             fig_flags{9,1}='shflux diff plot between 80s and later';
%             fig_flags{10,1}='relative vorticity plot';
%             fig_flags{11,1}='wind stress curl plot';
%             fig_flags{12,1}='vorticity diff plot between 80s and later';
%             fig_flags{13,1}='wind stress curl diff plot between 80s and later';
%             fig_flags{14,1}='vec plot';
%             fig_flags{15,1}='probability plot for pollock spwaning ground (hovmoller)';
%             fig_flags{16,1}='wind plot';
%             fig_flags{17,1}='test';
%             fig_flags{18,1}='Southern East Korean Bay nwv & numegg(15d) time series';
%             fig_flags{19,1}='Around Vladivostok nev & numegg(15d) time series';
%             fig_flags{20,1}='Southern East Korean Bay nwv & easterly time series';
%             fig_flags{21,1}='Around Vladivostok nev & southerly time series';
%             fig_flags{22,1}='Southern East Korean Bay2 nwv & numegg(15d) time series';
%             fig_flags{23,1}='Southern East Korean Bay2 nwv & easterly time series';
%             fig_flags{24,1}='Southern East Korean Bay2 nwv & East Korean Bay curl time series';
%             fig_flags{25,1}='Southern East Korean Bay2 nwv & East Korean Bay curl(all) time series';
%             fig_flags{26,1}='Around Vladivostok2 nev & numegg(15d) time series';
%             fig_flags{27,1}='Around Vladivostok2 nev & southerly time series';
%             fig_flags{28,1}='Around Vladivostok southerly & AOI time series';
%             fig_flags{29,1}='Around Southern East Korean Bay2 easterly & AOI time series';
%             fig_flags{30,1}='pair diff plot between 80s and later';
%             fig_flags{31,1}='SST plot';
%             fig_flags{32,1}='OISST plot';
%             fig_flags{33,1}='CMEMS geostrophic vec plot';
%             fig_flags{34,1}='Southern East Korean Bay2 spawning ground time series (regime)';
%             fig_flags{35,1}='Around Vladivostok spawning ground time series (regime)';
%             fig_flags{36,1}='Southern East Korean Bay2 numegg(15d) time series (regime)';
%             fig_flags{37,1}='Around Vladivostok numegg(15d) time series (regime)';
%             fig_flags{38,1}='Around Southern East Korean Bay2 numegg2 & AOI time series';
%             fig_flags{39,1}='Nothern East Sea spawning ground time series (regime)';
%             fig_flags{40,1}='Nothern East Sea numegg(15d) time series (regime)';
%             fig_flags{41,1}='Russian EEZ spawning ground time series (regime)';
%             fig_flags{42,1}='Russian EEZ numegg(15d) time series (regime)';
%             fig_flags{43,1}='each EEZ spawning ground time series (regime)';
%             fig_flags{44,1}='each EEZ numegg(15d) time series (regime)';
%             fig_flags{45,1}='each EEZ numegg(30d) time series (regime)';
%             fig_flags{46,1}='EEZ plot';
%             fig_flags{47,1}='lon lat distance time series';
%             fig_flags{48,1}='egg probability (elapsed 30 days) diff plot between 80s and later';
%             fig_flags{49,1}='South Korea EEZ spawning ground time series (regime)';
%             fig_flags{50,1}='Southern Korea EEZ numegg(??d) time series (regime)';
%             fig_flags{51,1}='Southern Korea EEZ numegg(30d) time series (regime)';
%             fig_flags{52,1}='probability plot for pollock location (elapsed 30 days)';
%             fig_flags{53,1}='particle moved distance-lat time series (regime, EEZ)';
%             fig_flags{54,1}='latitudinal particle density (horizontal bar)';
%             fig_flags{55,1}='latitudinal particle density (elapsed ??days, horizontal bar)';
%             fig_flags{56,1}='SST RMSE plot';
%             fig_flags{57,1}='latitudinal particle density difference (horizontal bar)';
%             fig_flags{58,1}='latitudinal particle density difference (elapsed ??days, horizontal bar)';
%             fig_flags{59,1}='spawning ground, dot plot (all period)';
%             fig_flags{60,1}='AO Index plot (regime)';
%             fig_flags{61,1}='latitudinal particle density difference (elapsed ??days, horizontal bar, percentage)';
%             fig_flags{62,1}='vec + isotherm(2,5) plot';
%             fig_flags{63,1}='particle latitude boxplot (elapsed ?? days)';
%             fig_flags{64,1}='EAWMI plot (regime)';
%             fig_flags{65,1}='# of individual plot at each grid (elapsed ?? days)';
%             fig_flags{66,1}='# of individual (elapsed ?? days) diff plot between 80s and later';
%             fig_flags{67,1}='probability plot for pollock spwaning ground (~41N)';
%             
%             for flagi=1:100
%                 fig_flags{flagi,2}=0;
%             end
%             fig_flags{1,2}=1;  % 'probability plot for pollock spwaning ground';
%             fig_flags{2,2}=1;  % 'probability plot for pollock location (elapsed ?? days)';
%             fig_flags{3,2}=1;  % 'temp diff plot between 80s and later';
%             fig_flags{4,2}=1;  % 'vec diff plot between 80s and later';
%             fig_flags{5,2}=1;  % 'egg probability (elapsed ?? days) diff plot between 80s and later';
%             fig_flags{6,2}=1;  % 'spawning ground diff plot between 80s and later';
%             fig_flags{7,2}=1;  % 'wind diff plot between 80s and later';
%             fig_flags{8,2}=0;  % 'tair diff plot between 80s and later';
%             fig_flags{9,2}=0;  % 'shflux diff plot between 80s and later';
%             fig_flags{10,2}=0;  % 'relative vorticity plot';
%             fig_flags{11,2}=0;  % 'wind stress curl plot';
%             fig_flags{12,2}=0;  % 'vorticity diff plot between 80s and later';
%             fig_flags{13,2}=0;  % 'wind stress curl diff plot between 80s and later';
%             fig_flags{14,2}=1;  % 'vec plot';
% %             fig_flags{15,2}=1;  % 'probability plot for pollock spwaning ground (hovmoller)';
%             fig_flags{16,2}=1;  % 'wind plot';
%             fig_flags{17,2}=0;  %'test';
%             fig_flags{18,2}=0;  %'Southern East Korean Bay nwv & numegg(15d) time series';
%             fig_flags{19,2}=0;  %'Around Vladivostok nev & numegg(15d) time series'
%             fig_flags{20,2}=0;  %'Southern East Korean Bay nwv & easterly time series';
%             fig_flags{21,2}=0;  %'Around Vladivostok nev & southerly time series'
%             fig_flags{22,2}=0;  %'Southern East Korean Bay2 nwv & numegg(15d) time series';
%             fig_flags{23,2}=0;  %'Southern East Korean Bay2 nwv & easterly time series';
%             fig_flags{24,2}=0;  %'Southern East Korean Bay2 nwv & East Korean Bay curl time series';
%             fig_flags{25,2}=0;  %'Southern East Korean Bay2 nwv & East Korean Bay curl(all) time series'
%             fig_flags{26,2}=0;  %'Around Vladivostok2 nev & numegg(15d) time series';
%             fig_flags{27,2}=0;  %'Around Vladivostok2 nev & southerly time series';
%             fig_flags{28,2}=0;  %'Around Vladivostok southerly & AOI time series';
%             fig_flags{29,2}=0;  %'Around Southern East Korean Bay2 easterly & AOI time series';
%             fig_flags{30,2}=0;  %'pair diff plot between 80s and later';
%             fig_flags{31,2}=1;  %'SST plot';
%             fig_flags{32,2}=1;  %'OISST plot';
%             fig_flags{33,2}=0;  %'CMEMS geostrophic vec plot';
%             fig_flags{34,2}=0;  %'Southern East Korean Bay2 spawning ground time series (regime)';
%             fig_flags{35,2}=0;  %'Around Vladivostok spawning ground time series (regime)';
%             fig_flags{36,2}=0;  %'Southern East Korean Bay2 numegg(??d) time series (regime)';
%             fig_flags{37,2}=0;  %'Around Vladivostok numegg(15d) time series (regime)';
%             fig_flags{38,2}=0;  %'Around Southern East Korean Bay2 numegg(15d) & AOI time series';
%             fig_flags{39,2}=0;  %'Northern East Sea spawning ground time series (regime)';
%             fig_flags{40,2}=0;  %'Northern East Sea numegg(15d) time series (regime)';
%             fig_flags{41,2}=0;  %'Russian EEZ spawning ground time series (regime)';
%             fig_flags{42,2}=0;  %'Russian EEZ numegg(15d) time series (regime)';
%             fig_flags{43,2}=0;  %'each EEZ spawning ground time series (regime)';
%             fig_flags{44,2}=0;  %'each EEZ numegg(??d) time series (regime)';
%             fig_flags{45,2}=0;  %'each EEZ numegg(30d) time series (regime)';
%             fig_flags{46,2}=0;  %'EEZ plot';
%             fig_flags{47,2}=0;  %'lon lat distance time series (not completed)';
%             fig_flags{48,2}=0;  % 'egg probability (elapsed 30 days) diff plot between 80s and later';
%             fig_flags{49,2}=1;  %'Southern Korea EEZ spawning ground time series (regime)';
%             fig_flags{50,2}=1;  %'Southern Korea EEZ numegg(??d) time series (regime)';
%             fig_flags{51,2}=0;  %'Southern Korea EEZ numegg(30d) time series (regime)';
%             fig_flags{52,2}=0;  %'probability plot for pollock location (elapsed 30 days)';
%             fig_flags{53,2}=1;  %'particle moved distance-lat time series (regime)';
%             fig_flags{54,2}=1;  %'latitudinal particle density (horizontal bar)';
%             fig_flags{55,2}=1;  %'latitudinal particle density (elapsed ??days, horizontal bar)';
%             fig_flags{56,2}=1;  %'SST RMSE plot';
%             fig_flags{57,2}=1;  %'latitudinal particle density difference (horizontal bar)';
%             fig_flags{58,2}=1;  %'latitudinal particle density difference (elapsed ??days, horizontal bar)';
%             fig_flags{59,2}=1;  %'spawning ground, all period, dot plot';
%             fig_flags{60,2}=1;  %'AO Index plot (regime)';
%             fig_flags{61,2}=1;  %'latitudinal particle density difference (elapsed ??days, horizontal bar)';
%             fig_flags{62,2}=1;  %'vec + isotherm(2,5) plot';
%             fig_flags{63,2}=1;  %'particle latitude boxplot (elapsed ?? days)';
%             fig_flags{64,2}=1;  %'EAWMI plot (regime)';
%             fig_flags{65,2}=1;  %'# of individual plot at each grid (elapsed ?? days)';
%             fig_flags{66,2}=1;  %'# of individual (elapsed ?? days) diff plot between 80s and later';
%             fig_flags{67,2}=0;  %'probability plot for pollock spwaning ground (~41N)';
% 
%         end
%         
        for flagi=1:100
            fig_flags{flagi,2}=0;
        end
        fig_flags{31,2}=2; 
        fig_flags{32,2}=2; 
        
        
        figdir=[figrawdir,LTRANS_testname, '\', regionname, '\spawn\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        
        spdiff_figdir=[figrawdir,LTRANS_testname, '\', regionname, '\spawn_diff\'];
        if (exist(strcat(spdiff_figdir) , 'dir') ~= 7)
            mkdir(strcat(spdiff_figdir));
        end 
        spdiff_outfile = strcat(spdiff_figdir,regionname);
        
        climfigdir=[figrawdir,LTRANS_testname, '\', regionname, '\clim\'];
        if (exist(strcat(climfigdir) , 'dir') ~= 7)
            mkdir(strcat(climfigdir));
        end 
        climoutfile = strcat(climfigdir,regionname);
        
        clim_atfigdir=[figrawdir,LTRANS_testname, '\', regionname, '\clim_at\'];
        if (exist(strcat(clim_atfigdir) , 'dir') ~= 7)
            mkdir(strcat(clim_atfigdir));
        end 
        clim_at_outfile = strcat(clim_atfigdir,regionname);
        
        clim_diff_figdir=[figrawdir,LTRANS_testname, '\', regionname, '\clim_diff\'];
        if (exist(strcat(clim_diff_figdir) , 'dir') ~= 7)
            mkdir(strcat(clim_diff_figdir));
        end 
        clim_diff_outfile = strcat(clim_diff_figdir,regionname);
        
        clim_atdiff_figdir=[figrawdir,LTRANS_testname, '\', regionname, '\clim_atdiff\'];
        if (exist(strcat(clim_atdiff_figdir) , 'dir') ~= 7)
            mkdir(strcat(clim_atdiff_figdir));
        end 
        clim_atdiff_outfile = strcat(clim_atdiff_figdir,regionname);
        
        all_figdir=[figrawdir,LTRANS_testname, '\', regionname, '\all\'];
        if (exist(strcat(all_figdir) , 'dir') ~= 7)
            mkdir(strcat(all_figdir));
        end 
        all_outfile = strcat(all_figdir,regionname);
        
        regime_figdir=[figrawdir,LTRANS_testname, '\', regionname, '\regime_shift\'];
        if (exist(strcat(regime_figdir) , 'dir') ~= 7)
            mkdir(strcat(regime_figdir));
        end 
        regime_outfile = strcat(regime_figdir,regionname);
        
        
        
% % %       'probability plot for pollock spwaning ground'
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_spawning_ground_probability_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_001;
            fig_flag=0;
        end
        
% % %       'probability plot for pollock location (elapsed ?? days)'
        fig_flag=fig_flags{2,2};
        while (fig_flag)
             for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_loc_prob_', num2str(temp_checktime, '%02i'),'days', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_002;
             end
            fig_flag=0;
        end

% % %       'temp diff plot between 80s and later'
        fig_flag=fig_flags{3,2};
        while (fig_flag)
            jpgname=strcat(clim_diff_outfile, '_', testname,'_',regionname, '_temp_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_003;
            fig_flag=0;
        end
        
% % %       'vec diff plot between 80s and later'
        fig_flag=fig_flags{4,2};
        fig_name=fig_flags{4,1};
        while (fig_flag)
            jpgname=strcat(clim_diff_outfile, '_', testname,'_',regionname, '_vec_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_004;
            fig_flag=0;
        end

% % %       'egg probability (elapsed 15 days) diff plot between 80s and later'
        fig_flag=fig_flags{5,2};
        fig_name=fig_flags{5,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(spdiff_outfile, '_', testname,'_',regionname, '_egg_prob_', num2str(temp_checktime, '%02i'), 'd_diff_80s_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_005;
            end
            fig_flag=0;
        end

% % %       'spwaning ground diff plot between 80s and later'
        fig_flag=fig_flags{6,2};
        fig_name=fig_flags{6,1};
        while (fig_flag)
            jpgname=strcat(spdiff_outfile, '_', testname,'_',regionname, '_spawn_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_006;
            fig_flag=0;
        end

% % %       'wind diff plot between 80s and later'
        fig_flag=fig_flags{7,2};
        fig_name=fig_flags{7,1};
        while (fig_flag)
            jpgname=strcat(clim_atdiff_outfile, '_', testname,'_',regionname, '_wind_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_007;
            fig_flag=0;
        end

% % %       'tair diff plot between 80s and later'
        fig_flag=fig_flags{8,2};
        fig_name=fig_flags{8,1};
        while (fig_flag)
            jpgname=strcat(clim_atdiff_outfile, '_', testname,'_',regionname, '_tair_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_008;
            fig_flag=0;
        end

% % %       'shflux diff plot between 80s and later'
        fig_flag=fig_flags{9,2};
        fig_name=fig_flags{9,1};
        while (fig_flag)
            jpgname=strcat(clim_atdiff_outfile, '_', testname,'_',regionname, '_shflux_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_009;
            fig_flag=0;
        end
        
% % %       'vort plot'
        fig_flag=fig_flags{10,2};
        fig_name=fig_flags{10,1};
        while (fig_flag)
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_vort_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_010;
            fig_flag=0;
        end
        
% % %       'wind stress curl plot'
        fig_flag=fig_flags{11,2};
        fig_name=fig_flags{11,1};
        while (fig_flag)
            jpgname=strcat(clim_at_outfile, '_', testname,'_',regionname, '_wind_stress_curl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_011;
            fig_flag=0;
        end
        
% % %       'vort diff plot between 80s and later'
        fig_flag=fig_flags{12,2};
        fig_name=fig_flags{12,1};
        while (fig_flag)
            jpgname=strcat(clim_diff_outfile, '_', testname,'_',regionname, '_vort_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_012;
            fig_flag=0;
        end

% % %       'wstr curl diff plot between 80s and later'
        fig_flag=fig_flags{13,2};
        fig_name=fig_flags{13,1};
        while (fig_flag)
            jpgname=strcat(clim_atdiff_outfile, '_', testname,'_',regionname, '_wstr_curl_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_013;
            fig_flag=0;
        end
        
% % %       'vec plot'
        fig_flag=fig_flags{14,2};
        fig_name=fig_flags{14,1};
        while (fig_flag)
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_vec_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_014;
            fig_flag=0;
        end
        
% % %       'probability plot for pollock spwaning ground (hovmoller)'
        fig_flag=fig_flags{15,2};
        fig_name=fig_flags{15,1};
        while (fig_flag)
            jpgname=strcat(all_outfile, '_', testname,'_',regionname, '_sp_ground_prob_hovm_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_015;
            fig_flag=0;
        end

% % %       'wind plot'
        fig_flag=fig_flags{16,2};
        fig_name=fig_flags{16,1};
        while (fig_flag)
            jpgname=strcat(clim_at_outfile, '_', testname,'_',regionname, '_wind_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_016;
            fig_flag=0;
        end

% % %       'test'
        fig_flag=fig_flags{17,2};
        fig_name=fig_flags{17,1};
        while (fig_flag)
            jpgname=strcat(clim_at_outfile, '_', testname,'_',regionname, '_wsc_ts_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_017;
            fig_flag=0;
        end
        
% % %       'Southern East Korean Bay nwv & numegg(15d) time series'
        fig_flag=fig_flags{18,2};
        fig_name=fig_flags{18,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB_ts_egg_density_15d_nwv_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_018;
            fig_flag=0;
        end       
  
% % %       'Around Vladivostok nev & numegg(15d) time series'
        fig_flag=fig_flags{19,2};
        fig_name=fig_flags{19,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SVLA_ts_egg_density_15d_nev_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_019;
            fig_flag=0;
        end   

% % %       'Southern East Korean Bay nwv & easterly time series'
        fig_flag=fig_flags{20,2};
        fig_name=fig_flags{20,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB_ts_easterly_nwv_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_020;
            fig_flag=0;
        end    

% % %       'Around Vladivostok nev & southerly time series'
        fig_flag=fig_flags{21,2};
        fig_name=fig_flags{21,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SVLA_ts_southerly_nev_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_021;
            fig_flag=0;
        end  

% % %       'Southern East Korean Bay2 nwv & numegg(15d) time series'
        fig_flag=fig_flags{22,2};
        fig_name=fig_flags{22,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB2_ts_egg_density_15d_nwv_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_022;
            fig_flag=0;
        end    

% % %       'Southern East Korean Bay2 nwv & easterly time series'
        fig_flag=fig_flags{23,2};
        fig_name=fig_flags{23,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB2_ts_easterly_nwv_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_023;
            fig_flag=0;
        end   

% % %       'Southern East Korean Bay2 nwv & East Korean Bay curl time series'
        fig_flag=fig_flags{24,2};
        fig_name=fig_flags{24,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB2_ts_wstr_curl_nwv_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_024;
            fig_flag=0;
        end   

% % %       'Southern East Korean Bay2 nwv & East Korean Bay curl(all) time series'
        fig_flag=fig_flags{25,2};
        fig_name=fig_flags{25,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB2_ts_allwstr_curl_nwv_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_025;
            fig_flag=0;
        end 
        
% % %       'Around Vladivostok2 nev & numegg(15d) time series'
        fig_flag=fig_flags{26,2};
        fig_name=fig_flags{26,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SVLA2_ts_egg_density_15d_nev_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_026;
            fig_flag=0;
        end   
        
% % %       'Around Vladivostok2 nev & southerly time series'
        fig_flag=fig_flags{27,2};
        fig_name=fig_flags{27,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SVLA2_ts_southerly_nev_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_027;
            fig_flag=0;
        end  
 
% % %       'Around Vladivostok southerly & AOI time series'
        fig_flag=fig_flags{28,2};
        fig_name=fig_flags{28,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SVLA_ts_southerly_AOI_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_028;
            fig_flag=0;
        end  

% % %       'Around Southern East Korean Bay2 easterly & AOI time series'
        fig_flag=fig_flags{29,2};
        fig_name=fig_flags{29,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB2_ts_easterly_AOI_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_029;
            fig_flag=0;
        end  

% % %       'pair diff plot between 80s and later'
        fig_flag=fig_flags{30,2};
        fig_name=fig_flags{30,1};
        while (fig_flag)
            jpgname=strcat(clim_atdiff_outfile, '_', testname,'_',regionname, '_pair_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_030;
            fig_flag=0;
        end

% % %       'SST plot'
        fig_flag=fig_flags{31,2};
        fig_name=fig_flags{31,1};
        while (fig_flag)
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_SST_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_031;
            fig_flag=0;
        end

% % %       'OISST plot, from 003 code'
        fig_flag=fig_flags{32,2};
        fig_name=fig_flags{32,1};
        while (fig_flag)
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_OISST_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_032;
            fig_flag=0;
        end
        
% % %       'CMEMS geostrophic vec plot, from CMEMS data (DAMO)'
        fig_flag=fig_flags{33,2};
        fig_name=fig_flags{33,1};
        while (fig_flag)
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_vec_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_033;
            fig_flag=0;
        end
        
% % %       'Southern East Korean Bay2 spawning ground time series (regime)';
        fig_flag=fig_flags{34,2};
        fig_name=fig_flags{34,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_SEKB2_regime_ts_sp_ground_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_034;
            fig_flag=0;
        end

% % %       'Around Vladivostok spawning ground time series (regime)';
        fig_flag=fig_flags{35,2};
        fig_name=fig_flags{35,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_SVLA_regime_ts_sp_ground_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_035;
            fig_flag=0;
        end    

% % %       'Southern East Korean Bay2 numegg(15d) time series (regime)';
        fig_flag=fig_flags{36,2};
        fig_name=fig_flags{36,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_SEKB2_regime_ts_egg_', num2str(temp_checktime, '%02i'),'d_', ...
                    num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_036;
            end
            fig_flag=0;
        end    

% % %       'Around Vladivostok numegg(15d) time series (regime)';
        fig_flag=fig_flags{37,2};
        fig_name=fig_flags{37,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_SVLA_regime_ts_egg_15d_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_037;
            fig_flag=0;
        end   

% % %       'Around Southern East Korean Bay2 numegg(15d) & AOI time series'
        fig_flag=fig_flags{38,2};
        fig_name=fig_flags{38,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB2_ts_egg_15d_AOI_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_038;
            fig_flag=0;
        end  

% % %       'Northern East Sea spawning ground time series (regime)';
        fig_flag=fig_flags{39,2};
        fig_name=fig_flags{39,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_NES_regime_ts_sp_ground_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_039;
            fig_flag=0;
        end

% % %       'Northern East Sea numegg(15d) time series (regime)';
        fig_flag=fig_flags{40,2};
        fig_name=fig_flags{40,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_NES_regime_ts_egg_15d_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_040;
            fig_flag=0;
        end    

% % %       'Russian EEZ spawning ground time series (regime)';
        fig_flag=fig_flags{41,2};
        fig_name=fig_flags{41,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_REEZ_regime_ts_sp_ground_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_041;
            fig_flag=0;
        end

% % %       'Russian EEZ numegg(15d) time series (regime)';
        fig_flag=fig_flags{42,2};
        fig_name=fig_flags{42,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_REEZ_regime_ts_egg_15d_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_042;
            fig_flag=0;
        end    

% % %       'each EEZ spawning ground time series (regime)';
        fig_flag=fig_flags{43,2};
        fig_name=fig_flags{43,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_Each_EEZ_regime_ts_sp_ground_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_043;
            fig_flag=0;
        end

% % %       'each EEZ numegg(15d) time series (regime)';
        fig_flag=fig_flags{44,2};
        fig_name=fig_flags{44,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_Each_EEZ_regime_ts_egg_', num2str(temp_checktime, '%02i'),'d_', ...
                    num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_044;
            end
            fig_flag=0;
        end
        
% % %       'each EEZ numegg(30d) time series (regime)';
        fig_flag=fig_flags{45,2};
        fig_name=fig_flags{45,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_Each_EEZ_regime_ts_egg_30d_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_045;
            fig_flag=0;
        end

% % %       'EEZ plot'
        fig_flag=fig_flags{46,2};
        fig_name=fig_flags{46,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_EEZ', '.tif'); %% ~_year_month.jpg            
            pollock_paper01_002_subroutine_046;
            fig_flag=0;
        end    

% % %       'lon lat distance timeseries'
        fig_flag=fig_flags{47,2};
        fig_name=fig_flags{47,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SEKB2_ts_egg_density_15d_nwv_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_047;
            fig_flag=0;
        end            

% % %       'egg probability (elapsed 30 days) diff plot between 80s and later'
        fig_flag=fig_flags{48,2};
        fig_name=fig_flags{48,1};
        while (fig_flag)
            jpgname=strcat(spdiff_outfile, '_', testname,'_',regionname, '_egg_30d_prob_diff_80s_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_048;
            fig_flag=0;
        end
        
% % %       'Southern Korea EEZ spawning ground time series (regime)';
        fig_flag=fig_flags{49,2};
        fig_name=fig_flags{49,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_SKEEZ_regime_ts_sp_ground_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_049;
            fig_flag=0;
        end
        
% % %       'Southern Korea EEZ numegg(15d) time series (regime)';
        fig_flag=fig_flags{50,2};
        fig_name=fig_flags{50,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_SKEEZ_regime_ts_egg_',num2str(temp_checktime, '%02i'), 'd_', ...
                    num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_050;
            end
            fig_flag=0;
        end    
        
% % %       'Southern Korea EEZ numegg(30d) time series (regime)';
        fig_flag=fig_flags{51,2};
        fig_name=fig_flags{51,1};
        while (fig_flag)
            jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_SKEEZ_regime_ts_egg_30d_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_051;
            fig_flag=0;
        end    
        
% % %       'probability plot for pollock location (elapsed 30 days)'
        fig_flag=fig_flags{52,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_loc_prob_30days', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_052;
            fig_flag=0;
        end
        
% % %       'particle moved distance-lat time series (regime)';
        fig_flag=fig_flags{53,2};
        fig_name=fig_flags{53,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(regime_outfile, '_', testname,'_',regionname, '_regime_ts_mv_lat_', num2str(temp_checktime, '%02i'), '_', ...
                    num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_053;
            end
            fig_flag=0;
        end
        
% % %       'latitudinal particle density (horizontal bar)';
        fig_flag=fig_flags{54,2};
        fig_name=fig_flags{54,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_lat_par_dens_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_054;
            end
            fig_flag=0;
        end

% % %       'latitudinal particle density (elapsed ?? days, horizontal bar)'
        fig_flag=fig_flags{55,2};
        while (fig_flag)
             for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_lat_par_dens_', num2str(temp_checktime, '%02i'),'days', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_055;
             end
            fig_flag=0;
        end

% % %       'SST RMSEplot'
        fig_flag_num=56;
        fig_flag=fig_flags{fig_flag_num,2};
        while (fig_flag)
            figvarname='SST_RMSE';
            figtypename = 'pcol';
            figtimename = [num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y_', ...
                            num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm'];
            extensionname = '.tif';
            scenname = 'historical';
            figname = Func_0002_get_figname(figrawdir, figvarname, figtypename, ...
                regionname, scenname, testname, figtimename, extensionname);
%             jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_SST_RMSE_', ...
%                 num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
%                 num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_056;
            fig_flag=0;
        end
        
% % %       'latitudinal particle density (horizontal bar)';
        fig_flag=fig_flags{57,2};
        fig_name=fig_flags{57,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(spdiff_outfile, '_', testname,'_',regionname, '_lat_par_dens_diff_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_057;
            end
            fig_flag=0;
        end
        
% % %       'latitudinal particle density difference (elapsed ??days, horizontal bar)';
        fig_flag=fig_flags{58,2};
        fig_name=fig_flags{58,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(spdiff_outfile, '_', testname,'_',regionname, '_lat_par_dens_diff_', num2str(temp_checktime, '%02i'),'d_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_058;
            end
            fig_flag=0;
        end
        
% % %       'probability plot for pollock spwaning ground'
        fig_flag=fig_flags{59,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_spawning_ground_all_dot_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_059;
            fig_flag=0;
        end
        
% % %       'Around Southern East Korean Bay2 easterly & AOI time series'
        fig_flag=fig_flags{60,2};
        fig_name=fig_flags{60,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_AOI_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_060;
            fig_flag=0;
        end  

% % %       'latitudinal particle density difference (elapsed ??days, horizontal bar, percentage)';
        fig_flag=fig_flags{61,2};
        fig_name=fig_flags{61,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(spdiff_outfile, '_', testname,'_',regionname, '_lat_par_dens_diff_percentage_', num2str(temp_checktime, '%02i'),'d_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_061;
            end
            fig_flag=0;
        end
        
% % %       'vec + isotherm (2,5) plot'
        fig_flag=fig_flags{62,2};
        fig_name=fig_flags{62,1};
        while (fig_flag)
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_vec_isotherm_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_062;
            fig_flag=0;
        end
        
% % %       'particle latitude boxplot (elapsed ?? days)';
        fig_flag=fig_flags{63,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_lat_par_box_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_063;
            fig_flag=0;
        end

% % %       'EAWMI time series'
        fig_flag=fig_flags{64,2};
        fig_name=fig_flags{64,1};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_EAWMI_', ...
                num2str(min(allyear),'%04i'),'_',num2str(max(allyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_064;
            fig_flag=0;
        end  
        
        
% % %       '# of individual plot at each grid (elapsed ?? days)'
        fig_flag=fig_flags{65,2};
        while (fig_flag)
             for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_loc_num_', num2str(temp_checktime, '%02i'),'days', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_065;
             end
            fig_flag=0;
        end
        
        
% % %       '# of individual (elapsed ?? days) diff plot between 80s and later';
        fig_flag=fig_flags{66,2};
        fig_name=fig_flags{66,1};
        while (fig_flag)
            for checkti=1:length(checktime)
                temp_checktime=checktime(checkti);
                jpgname=strcat(spdiff_outfile, '_', testname,'_',regionname, '_egg_num_', num2str(temp_checktime, '%02i'), 'd_diff_80s_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                    num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
                pollock_paper01_002_subroutine_066;
            end
            fig_flag=0;
        end
        
% % %       'probability plot for pollock spwaning ground (~41N)'
        fig_flag=fig_flags{67,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_spawning_ground_probability_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), 'y', ...
                num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'), 'm', '.tif'); %% ~_year_month.jpg
            pollock_paper01_002_subroutine_067;
            fig_flag=0;
        end
        
        
    end
end



