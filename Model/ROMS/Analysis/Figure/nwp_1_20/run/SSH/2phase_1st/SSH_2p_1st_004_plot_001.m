close all; clear all;  clc;
warning off;

% all_testname2 = {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
% all_testname2 = {'test2107', 'test2108', 'test2109', 'test2110', 'test2111'};
all_testname2 = {'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
% all_testname2 = {'test2102'};

% all_region2 ={'NWP'};
all_region2 ={'NWP', 'AKP4'};
% all_region2 ={'AKP4'};

% all_region2 ={'TEST'};

for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2
        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIASRF
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\user\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
            addpath(genpath([dropboxpath, '\source\matlab\function']));
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end

        shadlev = [0 35];
        rms_shadlev = [0 4];
    %     trendlev = [-3 3];  %% trend lev
        trendlev = [-10 10];  %% trend lev
%         abstrendlev =[2 7]; % historical
        abstrendlev =[0 10]; % ssp585
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.15 0.15];
        trendplotlev = [0 7];
        trenddifflev = [-10 10];
        sshlev =[-0.3 0.3];
        sshdifflev = [0 20];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1985:2014]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        gcmtestname = Func_0004_get_GCMname_from_RCM(testname);
        scenname= Func_0013_RCM_CMIP6_scenname(testname);
        varname ='zeta';
        variable='SSH';
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        
% % %         switch region
        [refpolygon, lonlat, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(regionname);
        
% % %         flag configuration
        for folding=1:1
            fig_flags{1,1}='correlation_interped_plot';
            fig_flags{2,1}='correlation_interped_filtered_plot';
            fig_flags{3,1}='correlation_seasonal_plot';
            fig_flags{4,1}='correlation_interped_detrended_plot';
            fig_flags{5,1}='correlation_low_passed_plot';
            fig_flags{6,1}='correlation_corrected_low_passed_plot';
            fig_flags{7,1}='correlation_corrected_low_passed_difference_plot';
            fig_flags{8,1}='correlation_corrected_low_passed_difference_normalized_plot';
            fig_flags{9,1}='difference_mask_low_passed_corr_rms_plot';
            fig_flags{10,1}='difference_low_passed_corr_rms_plot';
            fig_flags{11,1}='correlation_climatological_plot';
            fig_flags{12,1}='trend_absolute_plot';
            fig_flags{13,1}='trend_absolute_bwrmap_plot';
            fig_flags{14,1}='trend_relative_plot';
            fig_flags{15,1}='trend_difference_plot';
            fig_flags{16,1}='trend_absolute_corrected_plot';
            fig_flags{17,1}='trend_climatological_plot';
            fig_flags{18,1}='time_series_ssh_seasonal_filtered';
            fig_flags{19,1}='time_series_ssh';
            fig_flags{20,1}='time_series_ssh_seasonal_filtered_corrected';
            fig_flags{21,1}='time_series_ssh_corrected';
            fig_flags{22,1}='time_series_ssh_low_passed';
            fig_flags{23,1}='time_series_ssh_detrended_low_passed';
            fig_flags{24,1}='time_series_ssh_climatological';
            fig_flags{25,1}='trend_climatological_line';
            fig_flags{26,1}='time_series_ssh_seasonal';
            fig_flags{27,1}='time_series_sla_cmemsstructed_seasonal';
            fig_flags{28,1}='amplitude_seasonal_ssh_cmemsstructed';
            fig_flags{29,1}='amplitude_seasonal_ssh_model';
            fig_flags{30,1}='amplitude_difference_seasonal_ssh_model';
            fig_flags{31,1}='time_series_ssh_yearly';
            fig_flags{32,1}='time_series_cmemsstructed_plot';
            fig_flags{33,1}='RMS_ssh_plot';
            fig_flags{34,1}='RMS_corrected_ssh_plot';
            fig_flags{35,1}='BIAS_ssh_plot';
            fig_flags{36,1}='variance_cmems_sla_plot';
            fig_flags{37,1}='std_cmems_sla_plot';
            fig_flags{38,1}='variance_cmems_sla_seasonal_plot';
            fig_flags{39,1}='variance_model_sla_plot';
            fig_flags{40,1}='variance_model_sla_seasonal_plot';
            fig_flags{41,1}='variance_cmemsstructed_sla_lowpass_plot';
            fig_flags{42,1}='variance_model_sla_lowpass_plot';
            fig_flags{43,1}='pattern_correlation';
            fig_flags{44,1}='correlation_movmean_plot';
            fig_flags{45,1}='correlation_corrected_movmean_plot';
            fig_flags{46,1}='time_series_ssh_movmean';
            fig_flags{47,1}='variance_cmems_sla_movmean_plot';
            fig_flags{48,1}='variance_model_sla_movmean_plot';
            fig_flags{49,1}='trend_cmems_relative_plot';
            fig_flags{50,1}='trend_cmems_absolute_bwrmap_plot';
            fig_flags{51,1}='trend_recon_relative_plot';
            fig_flags{52,1}='trend_recon_absolute_bwrmap_plot';
            fig_flags{53,1}='time_series_ssh_seasonal_cmems';
            fig_flags{54,1}='std_model_sla_movmean_plot';
            fig_flags{55,1}='std_cmems_sla_movmean_plot';
            fig_flags{56,1}='variance_model_sla_plot';
            fig_flags{57,1}='std_model_yearly_sla_detrended_plot';
            fig_flags{58,1}='std_cmems_yearly_sla_detrended_plot';
            fig_flags{59,1}='SLR plot from yearly exponential fit data';
            fig_flags{60,1}='SLR plot from yearly poly1 fit data';
            fig_flags{61,1}='SLR r plot from yearly exp fit';
            fig_flags{62,1}='SLR r plot from yearly poly1 fit';
            fig_flags{63,1}='SLR rmse plot from yearly exp fit';
            fig_flags{64,1}='SLR rmse plot from yearly poly1 fit';
            fig_flags{65,1}='std_cmems_yearly_sla_exp_fit_detrended_plot';
            fig_flags{66,1}='std_model_yearly_sla_exp_fit_detrended_plot';
            fig_flags{67,1}='cmems SLR plot from yearly exponential fit data';
            fig_flags{68,1}='cmems SLR plot from yearly poly1 fit data';
            fig_flags{69,1}='std_cmems_yearly_sla_detrended_plot';
            fig_flags{70,1}='std_model_yearly_sla_detrended_plot';
            fig_flags{71,1}='surf_u rms';
            fig_flags{76,1}='correlation_interped_spearman_plot';

            for flagi=1:100
                fig_flags{flagi,2}=0;
            end
            fig_flags{1,2}=0;
            fig_flags{2,2}=0;
            fig_flags{3,2}=0;
            fig_flags{4,2}=0;
            fig_flags{12,2}=2;
            fig_flags{13,2}=0;
            fig_flags{14,2}=0;
            fig_flags{18,2}=0;
            fig_flags{19,2}=0;
            fig_flags{26,2}=0;
            fig_flags{28,2}=0;
            fig_flags{29,2}=0;
            fig_flags{31,2}=0;
            fig_flags{32,2}=0;
            fig_flags{37,2}=0;
            fig_flags{49,2}=0;
            fig_flags{50,2}=0;
            fig_flags{51,2}=0;
            fig_flags{52,2}=0;
            fig_flags{53,2}=0;
            fig_flags{57,2}=0;
            fig_flags{58,2}=0;
            fig_flags{59,2}=0;
            fig_flags{60,2}=0;
            fig_flags{61,2}=0;
            fig_flags{62,2}=0;
            fig_flags{63,2}=0;
            fig_flags{64,2}=0;
            fig_flags{65,2}=0;
            fig_flags{66,2}=0;
            fig_flags{67,2}=0;
            fig_flags{68,2}=0;
            fig_flags{69,2}=0;
            fig_flags{70,2}=0;
            fig_flags{76,2}=0;
        end
%         for flagi=1:100
%             fig_flags{flagi,2}=0;
%         end
%         
%         fig_flags{43,2}=2;
        
%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
%         filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m';
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
            gcmfiledir = strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', gcmtestname, '\'); % % where data files are
            matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\analysis\');
            dirs.figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\'); % % where figure files will be saved
            tmp.variable=variable;
            tmp.fs=filesep;
            tmp.regionname=regionname;
            tmp.testname=testname;
            RCM_info.testname=testname;
            RCM_info.years=inputyear;
            RCM_info.regionname =regionname;
            RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
              'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
            RCM_info.savedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs, 'run', tmp.fs, ...
             'zeta', tmp.fs];
            
        elseif (strcmp(system_name,'GLNXA64'))
        end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

        modelfilename = strcat(filedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        corrfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        lpffilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        lpf_corrfilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        movfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        mov_corrfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        detrendfilename = strcat(filedir, testname,'_',regionname, 'detrended_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        cmemsvelfilename = strcat(filedir, testname,'_',regionname, 'cmems_vel_rms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        interpedvelfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_vel_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        avhrrfilename = strcat(filedir, testname, '_',regionname,'_rms_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        gcminterpedfilename = strcat(gcmfiledir, gcmtestname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        gcmcmemsvelfilename=strcat(cmip5dir,'\surface\', gcmtestname,'_',regionname, 'cmems_vel_rms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        
        valnum=0;
        run('C:\Users\user\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);
        
        

        run(param_script);
        
        figdir=[figrawdir,'Trend\SSH\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        
        outfile = strcat(figdir,regionname);
        trendoutfile = strcat(figdir,regionname);
        
%         cmems_trend=ncread(cmemsfilename, 'cmems_trend');
%         cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
% 
%         cmems_mask=ones(size(cmems_trend));
%         cmems_mask(isnan(cmems_trend))=NaN;
%         
%         lon_rho=ncread(modelfilename,'lon_rho');
%         lat_rho=ncread(modelfilename,'lat_rho');
%         trend_filtered = ncread(modelfilename,'trend_filtered');
% %         mean_trend_filtered = ncread(filename,'mean_trend_filtered');
%         cmems_trend_filtered = ncread(cmemsfilename,'cmems_trend_filtered');
%         interped_trend_filtered = ncread(interpedfilename,'interped_trend_filtered');
%         lon_cmems=ncread(cmemsfilename,'lon_cmems');
%         lat_cmems=ncread(cmemsfilename,'lat_cmems');
        
% start-------------------- make timedata(xData) for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
        

% start-------------------- absolute trend plot
        fig_flag=fig_flags{12,2};
        while (fig_flag)
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'ssh_abs_trend', tmp.fs, ...
                num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
            if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
                mkdir(strcat(dirs.figdir));
            end 
            tmp.tifname=strcat(dirs.figdir, testname, '_abs_ssh_trend_',num2str(min(RCM_info.years),'%04i'), ...
                '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
%                 matname = [matdir, testname, '_', regionname, '_', '012_SSH_abs_trend', ...
%                     '_', num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), '.mat'];
%                 
                RCM_info.matname_interped = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_cmems_interped_ssh_', ...
                    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
                load(RCM_info.matname_interped);
%                 CMEMS_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_CMEMS_ssh_', ...
%                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
%                 load(CMEMS_info.matname);
                RCM_info.matname_trends = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_ssh_trend_', ...
                    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
                load(RCM_info.matname_trends);
                m_proj(m_proj_name,'lon',[RCM_grid.domain(1)-0.5 RCM_grid.domain(2)+0.5],'lat',[RCM_grid.domain(3)-0.5 RCM_grid.domain(4)+0.5]);
                hold on;
                
%                 m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
                m_pcolor(double(RCM_grid.lon_rho)',RCM_grid.lat_rho',squeeze(RCM_data_trend.yearly_trend(:,:)'));
%                 mean_trend_filtered=mean(interped_trend_filtered(:),'omitnan');

%                 mean_trend_filtered=mean(trend_filtered(:),'omitnan');
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                [m_value, error_status] = Func_0011_get_area_weighted_mean(RCM_data_trend.yearly_trend, RCM_grid.lon_rho, RCM_grid.lat_rho);
%                 titlename = strcat('SSH trend(abs), ',RCM_info.abbs{testnameind2}, ',(',num2str(min(inputyear),'%04i'),'-', ...
%                     num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(m_value,2)), ' mm/y');  %% + glacier contribution
                titlename = strcat(RCM_info.abbs{testnameind2}, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(m_value,2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis(abstrendlev);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif'); RemoveWhiteSpace([], 'file', jpgname); 
                saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
                disp(' ')
                disp(['clim_', 'absolute_ssh_trend', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            tmp.tifname=strcat(dirs.figdir, gcmtestname, '_abs_ssh_trend_',num2str(min(RCM_info.years),'%04i'), ...
                '_',num2str(max(RCM_info.years),'%04i'), '.tif'); %% ~_year_month.jpg

            if (exist(tmp.tifname , 'file') ~= 2 || fig_flag==2)
%                 matname = [matdir, testname, '_', regionname, '_', '012_SSH_abs_trend', ...
%                     '_', num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), '.mat'];
%                 
                GCM_info.matname_interped = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_cmems_interped_ssh_', ...
                    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
                load(GCM_info.matname_interped);
%                 CMEMS_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_CMEMS_ssh_', ...
%                     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
%                 load(CMEMS_info.matname);
                GCM_info.matname_trends = [RCM_info.savedir, gcmtestname,'_',RCM_info.regionname, '_GCM_ssh_trend_', ...
                    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
                load(GCM_info.matname_trends);
                m_proj(m_proj_name,'lon',[RCM_grid.domain(1)-0.5 RCM_grid.domain(2)+0.5],'lat',[RCM_grid.domain(3)-0.5 RCM_grid.domain(4)+0.5]);
                hold on;
                
%                 m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
                m_pcolor(double(GCM_grid.lon)',GCM_grid.lat',squeeze(GCM_data_trend.yearly_trend(:,:)'));
%                 mean_trend_filtered=mean(interped_trend_filtered(:),'omitnan');

%                 mean_trend_filtered=mean(trend_filtered(:),'omitnan');
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                [m_value, error_status] = Func_0011_get_area_weighted_mean(GCM_data_trend.yearly_trend, GCM_grid.lon, GCM_grid.lat);
%                 titlename = strcat('SSH trend(abs), ',GCM_info.abbs{testnameind2}, ',(',num2str(min(inputyear),'%04i'),'-', ...
%                     num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(m_value,2)), ' mm/y');  %% + glacier contribution
                titlename = strcat(GCM_info.abbs{testnameind2}, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(m_value,2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis(abstrendlev);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
                disp(' ')
                disp(['clim_', 'absolute_ssh_trend', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            
            fig_flag=0;
        end
 
        
    end
end