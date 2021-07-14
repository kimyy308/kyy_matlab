close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'test53', 'test54', 'test55', 'test56', 'ens03'};
all_testname2 = {'test53', 'test54', 'test55', 'test56'};
% all_testname2 = {'test56'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP4', 'YS', 'YSECS', 'ECS2', 'NES', 'SES', 'ES'}

% all_region2 ={'NWP', 'AKP4'};
all_region2 ={'AKP4'};

% all_region2 ={'TEST'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2
        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIAS
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\user\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
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
        abstrendlev =[2 7];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.15 0.15];
        trendplotlev = [0 7];
        trenddifflev = [-10 10];
        sshlev =[-0.3 0.3];
        sshdifflev = [0 20];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1993:2005]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        gcmtestname = Func_0004_get_GCMname_from_RCM(testname);
        scenname= Func_0003_RCM_CMIP5_scenname(testname);
        varname ='zeta';
        variable='SSH';
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        
% % %         switch region
        [refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);
        
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
            fig_flags{69,1}='std_cmems_yearly_sla_exp_fit_detrended_plot';
            fig_flags{70,1}='std_model_yearly_sla_exp_fit_detrended_plot';
            fig_flags{71,1}='surf_u rms';

            for flagi=1:100
                fig_flags{flagi,2}=0;
            end
            fig_flags{1,2}=1;
%             fig_flags{2,2}=1;
%             fig_flags{3,2}=0;
%             fig_flags{4,2}=0;
%             fig_flags{5,2}=1;
%             fig_flags{6,2}=0;
%             fig_flags{7,2}=0;
%             fig_flags{8,2}=0;
%             fig_flags{9,2}=0;
%             fig_flags{10,2}=0;
%             fig_flags{11,2}=0;
            fig_flags{12,2}=1;
            fig_flags{13,2}=1;
            fig_flags{14,2}=1;
%             fig_flags{15,2}=1;
%             fig_flags{16,2}=0;
%             fig_flags{17,2}=0;
%             fig_flags{18,2}=1;
%             fig_flags{19,2}=1;
%             fig_flags{20,2}=0;
%             fig_flags{21,2}=0;
%             fig_flags{22,2}=2;
%             fig_flags{23,2}=0;
%             fig_flags{24,2}=0;
%             fig_flags{25,2}=1;
            fig_flags{26,2}=1;
%             fig_flags{27,2}=0;
%             fig_flags{28,2}=0;
%             fig_flags{29,2}=0;
%             fig_flags{30,2}=0;
%             fig_flags{31,2}=1;
            fig_flags{32,2}=1;
%             fig_flags{33,2}=0;
%             fig_flags{34,2}=0;
%             fig_flags{35,2}=0;
            fig_flags{36,2}=1;
            fig_flags{37,2}=1;
%             fig_flags{38,2}=0;
            fig_flags{39,2}=1;
%             fig_flags{40,2}=0;
%             fig_flags{41,2}=0;
%             fig_flags{42,2}=0;
%             fig_flags{43,2}=0;
            fig_flags{44,2}=1;
%             fig_flags{45,2}=0;
            fig_flags{46,2}=1;
            fig_flags{47,2}=1;
            fig_flags{48,2}=1;
            fig_flags{49,2}=1;
            fig_flags{50,2}=0;
            fig_flags{51,2}=0;
            fig_flags{52,2}=0;
            fig_flags{53,2}=1;
            fig_flags{54,2}=1;
            fig_flags{55,2}=1;
            fig_flags{56,2}=1;
            fig_flags{57,2}=1;
            fig_flags{58,2}=1;
            fig_flags{59,2}=1;
            fig_flags{60,2}=1;
            fig_flags{61,2}=1;
            fig_flags{62,2}=1;
            fig_flags{63,2}=1;
            fig_flags{64,2}=1;
            fig_flags{65,2}=1;
            fig_flags{66,2}=1;
            fig_flags{67,2}=1;
            fig_flags{68,2}=1;
            fig_flags{69,2}=1;
            fig_flags{70,2}=1;
        end
        for flagi=1:100
            fig_flags{flagi,2}=0;
        end
        
        fig_flags{43,2}=2;
        
%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
%         filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
            gcmfiledir = strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', gcmtestname, '\'); % % where data files are
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
        
        cmems_trend=ncread(cmemsfilename, 'cmems_trend');
        cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');

        cmems_mask=ones(size(cmems_trend));
        cmems_mask(isnan(cmems_trend))=NaN;
        
        lon_rho=ncread(modelfilename,'lon_rho');
        lat_rho=ncread(modelfilename,'lat_rho');
        trend_filtered = ncread(modelfilename,'trend_filtered');
%         mean_trend_filtered = ncread(filename,'mean_trend_filtered');
        cmems_trend_filtered = ncread(cmemsfilename,'cmems_trend_filtered');
        interped_trend_filtered = ncread(interpedfilename,'interped_trend_filtered');

% start-------------------- make timedata(xData) for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
        
% start-------------------- corr_interped plot
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_cmems_interped_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                corr_interped=ncread(corrfilename,'corr_interped');
                corr_interped=corr_interped.*cmems_mask;
                lon_cmems=ncread(cmemsfilename,'lon_cmems');
                lat_cmems=ncread(cmemsfilename,'lat_cmems');

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',squeeze(corr_interped(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('corr, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped(:), 'omitnan'),2)));  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([0.2 0.8]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- corr_interped_filtered plot
        fig_flag=fig_flags{2,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_cmems_interped_filtered_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                corr_interped_filtered=ncread(corrfilename,'corr_interped_filtered');
                lon_cmems=ncread(cmemsfilename,'lon_cmems');
                lat_cmems=ncread(cmemsfilename,'lat_cmems');

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',squeeze(corr_interped_filtered(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('corr clim filtered, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped_filtered(:), 'omitnan'),2)));  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([0.2 0.8]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            fig_flag=0;
        end

% % % % start-------------------- corr_seasonal plot
        fig_flag=fig_flags{3,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, 'cmems_interped_corr_spatial_mean_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                corr_spatial_mean=ncread(filename,'corr_spatial_mean');
                lon_cmems=ncread(filename,'lon_cmems');
                lat_cmems=ncread(filename,'lat_cmems');

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',squeeze(corr_spatial_mean(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('corr climatology, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_spatial_mean(:), 'omitnan'),2)));  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([0.2 0.8]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- corr_interped_detrended plot
        fig_flag=fig_flags{4,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_cmems_interped_detrended_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                corr_interped_detrended=ncread(filename,'corr_interped_detrended');
                lon_cmems=ncread(filename,'lon_cmems');
                lat_cmems=ncread(filename,'lat_cmems');

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',squeeze(corr_interped_detrended(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('corr filt detrend, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped_detrended(:), 'omitnan'),2)));  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([0.2 0.8]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- corr_ lowpassed plot
        fig_flag=fig_flags{5,2};
        while (fig_flag)
            nyears=[2,3,5];
    %         nc_varname_prefixes={'interped_', 'interped_detrended_'};
                    nc_varname_prefixes={'interped_'};

    %         nc_titlename_prefixes={'lowpass', 'lowpass-det'};
                    nc_titlename_prefixes={'lowpass'};

            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                    jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_cmems_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        eval(['corr_',nc_varname, '=ncread(lpf_corrfilename,','''','corr_', nc_varname,'''',');']);
                        eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);

                        lon_cmems=ncread(cmemsfilename,'lon_cmems');
                        lat_cmems=ncread(cmemsfilename,'lat_cmems');

                        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                        hold on;
                        eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                        shading(gca,m_pcolor_shading_method);
                        m_gshhs_i('color',m_gshhs_line_color);
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                        eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        tit_prefix=nc_titlename_prefixes{nc_varnameij};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr,2)));  

                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(wrmap);
                        set(h,'fontsize',colorbar_fontsize);
                        caxis([0.2 0.8]);

                        % set grid
                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                        saveas(gcf,jpgname,'tif');

                        disp(' ')
                        disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
                        disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- corr_ lowpassed plot (corrected)
        fig_flag=fig_flags{6,2};
        while (fig_flag)
            nyears=[2,5];
            nc_varname_prefixes={'corrected_interped_sla_'};
            nc_titlename_prefixes={'lowpass-c'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                    jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_corrected_lowpass_cmems_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        eval(['corr_',nc_varname, '=ncread(filename,','''','corr_', nc_varname,'''',');']);
                        eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);

                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

                        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                        hold on;
                        eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                        shading(gca,m_pcolor_shading_method);
                        m_gshhs_i('color',m_gshhs_line_color);
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                        eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        tit_prefix=nc_titlename_prefixes{nc_varnameij};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr,2)));  

                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(wrmap);
                        set(h,'fontsize',colorbar_fontsize);
                        caxis([0.2 0.8]);

                        % set grid
                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                        saveas(gcf,jpgname,'tif');

                        disp(' ')
                        disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
                        disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- corr_ lowpassed plot (corrected_diff)
        fig_flag=fig_flags{7,2};
        while (fig_flag)
            nyears=[2,5];
            nc_varname_prefixes={'corrected_interped_sla_'};
            nc_varname_prefixes2={'interped_'};
            nc_titlename_prefixes={'lowpass-c'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                    nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                    jpgname=strcat(outfile, '_', 'diff','_',regionname, '_corr_diff_cmems_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                                            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                                            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        eval(['corr_',nc_varname, '=ncread(diffname,','''','corr_', nc_varname2,'''',');']);
                        eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);
                        eval(['corr2_',nc_varname, '=ncread(diffname2,','''','corr_', nc_varname,'''',');']);
                        eval(['corr2_',nc_varname, '=corr2_', nc_varname, '.*cmems_mask;']);
                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

                        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                        hold on;
                        eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''','-','corr2_', nc_varname, '(:,:)', '''',  '));']);
                        shading(gca,m_pcolor_shading_method);
                        m_gshhs_i('color',m_gshhs_line_color);
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                        eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        eval(['mcorr2=mean(corr2_',nc_varname,'(:),', '''', 'omitnan', '''',');']);

                        tit_prefix=nc_titlename_prefixes{nc_varnameij};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr-mcorr2,2)));  

                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(bwrmap);
                        set(h,'fontsize',colorbar_fontsize);
                        caxis([-0.5 0.5]);
    % corr_corrected_interped_sla_2y_lowpass-corr2_corrected_interped_sla_2y_lowpass
                        % set grid
                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                        saveas(gcf,jpgname,'tif');

                        disp(' ')
                        disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
                        disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- corr_ lowpassed plot (corrected_diff, normalized)
        fig_flag=fig_flags{8,2};
        while (fig_flag)
            nyears=[2,5];
            nc_varname_prefixes={'corrected_interped_sla_'};
            nc_varname_prefixes2={'interped_'};
            nc_titlename_prefixes={'lowpass-c'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                    nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                    jpgname=strcat(outfile, '_', 'diff','_',regionname, '_corr_diff_norm_cmems_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                                            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                                            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        eval(['corr_',nc_varname, '=ncread(diffname,','''','corr_', nc_varname2,'''',');']);
                        eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);
                        eval(['corr2_',nc_varname, '=ncread(diffname2,','''','corr_', nc_varname,'''',');']);
                        eval(['corr2_',nc_varname, '=corr2_', nc_varname, '.*cmems_mask;']);
                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

                        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                        hold on;
                        eval(['mean_corr_diff=', 'squeeze(corr_', nc_varname, '(:,:)', '-', 'corr2_', nc_varname, '(:,:));'])
                        std_mean_corr_diff=std(mean_corr_diff(:), 'omitnan');
                        n_mean_corr_diff=mean_corr_diff/std_mean_corr_diff;
                        eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',n_mean_corr_diff','''',  ');']);

    %                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''','-','corr2_', nc_varname, '(:,:)', '''',  '));']);
                        shading(gca,m_pcolor_shading_method);
                        m_gshhs_i('color',m_gshhs_line_color);
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                        eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        eval(['mcorr2=mean(corr2_',nc_varname,'(:),', '''', 'omitnan', '''',');']);

                        tit_prefix=nc_titlename_prefixes{nc_varnameij};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(n_mean_corr_diff(:), 'omitnan'),2)));  

                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(bwrmap);
                        set(h,'fontsize',colorbar_fontsize);
                        caxis([-6 6]);
    % corr_corrected_interped_sla_2y_lowpass-corr2_corrected_interped_sla_2y_lowpass
                        % set grid
                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                        saveas(gcf,jpgname,'tif');

                        disp(' ')
                        disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
                        disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- diff mask between lowpassed corr and rms 
        fig_flag=fig_flags{9,2};
        while (fig_flag)
            nyears=[2,5];
            nc_varname_prefixes={'corrected_interped_sla_'};
            nc_varname_prefixes2={'interped_'};
            nc_titlename_prefixes={'lowpass-c'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                    nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                    jpgname=strcat(outfile, '_', 'diff','_',regionname, '_rms_corr_diff_cmems_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)


                        matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                        '_rms_diff_',nc_varname, '_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                        load(matname)
                        rms_mask=NaN(size(mean_rms_diff));
                        rms_mask(mean_rms_diff<0)=-1;
                        rms_mask(mean_rms_diff>0)=1;
                        diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                                            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                                            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        eval(['corr_',nc_varname, '=ncread(diffname,','''','corr_', nc_varname2,'''',');']);
                        eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);
                        eval(['corr2_',nc_varname, '=ncread(diffname2,','''','corr_', nc_varname,'''',');']);
                        eval(['corr2_',nc_varname, '=corr2_', nc_varname, '.*cmems_mask;']);
                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

                        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                        hold on;
                        eval(['mean_corr_diff=', 'squeeze(corr_', nc_varname, '(:,:)', '-', 'corr2_', nc_varname, '(:,:));'])

                        corr_mask=NaN(size(mean_corr_diff));
                        corr_mask(mean_corr_diff<0)=-1;
                        corr_mask(mean_corr_diff>0)=1;

                        diff_mask = rms_mask+corr_mask;
    %                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',mean_corr_diff','''',  ');']);
                        m_pcolor(lon_cmems', lat_cmems', diff_mask');
                        shading(gca,m_pcolor_shading_method);
                        m_gshhs_i('color',m_gshhs_line_color);
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                        eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        eval(['mcorr2=mean(corr2_',nc_varname,'(:),', '''', 'omitnan', '''',');']);

                        tit_prefix=nc_titlename_prefixes{nc_varnameij};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr-mcorr2,2)));  

                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(bwrmap);
                        set(h,'fontsize',colorbar_fontsize);
                        caxis([-2 2]);
    % corr_corrected_interped_sla_2y_lowpass-corr2_corrected_interped_sla_2y_lowpass
                        % set grid
                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                        saveas(gcf,jpgname,'tif');

                        disp(' ')
                        disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
                        disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- diff between lowpassed corr and rms 
        fig_flag=fig_flags{10,2};
        while (fig_flag)
            nyears=[2,5];
            nc_varname_prefixes={'corrected_interped_sla_'};
            nc_varname_prefixes2={'interped_'};
            nc_titlename_prefixes={'diff-val-norm'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname_prefix2=nc_varname_prefixes2{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                    nc_varname2=[nc_varname_prefix2, num2str(nyear),'y_lowpass'];
                    jpgname=strcat(outfile, '_', 'diff','_',regionname, '_rms_corr_diff_norm_cmems_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)


                        matname = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                        '_rms_diff_',nc_varname, '_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
                        load(matname)

                        std_mean_rms_diff=std(mean_rms_diff(:),'omitnan');
                        n_mean_rms_diff=mean_rms_diff/std_mean_rms_diff;

                        diffname = ['E:\Data\Model\ROMS\nwp_1_10\','test11','\run\','test11','_',regionname, ...
                                            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        diffname2 = ['E:\Data\Model\ROMS\nwp_1_10\','test12','\run\','test12','_',regionname, ...
                                            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
                        eval(['corr_',nc_varname, '=ncread(diffname,','''','corr_', nc_varname2,'''',');']);
                        eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);
                        eval(['corr2_',nc_varname, '=ncread(diffname2,','''','corr_', nc_varname,'''',');']);
                        eval(['corr2_',nc_varname, '=corr2_', nc_varname, '.*cmems_mask;']);
                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

                        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                        hold on;
                        eval(['mean_corr_diff=', 'squeeze(corr_', nc_varname, '(:,:)', '-', 'corr2_', nc_varname, '(:,:));'])

                        std_mean_corr_diff=std(mean_corr_diff(:),'omitnan');
                        n_mean_corr_diff=mean_corr_diff/std_mean_corr_diff;

    %                     diff_mask = rms_mask+corr_mask;
                        diff_value = n_mean_rms_diff + n_mean_corr_diff;
    %                     eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',mean_corr_diff','''',  ');']);
                        m_pcolor(lon_cmems', lat_cmems', diff_value');
                        shading(gca,m_pcolor_shading_method);
                        m_gshhs_i('color',m_gshhs_line_color);
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                        eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        eval(['mcorr2=mean(corr2_',nc_varname,'(:),', '''', 'omitnan', '''',');']);

                        tit_prefix=nc_titlename_prefixes{nc_varnameij};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(diff_value(:),'omitnan'),2)));  

                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(bwrmap);
                        set(h,'fontsize',colorbar_fontsize);
                        caxis([-4 4]);
    % corr_corrected_interped_sla_2y_lowpass-corr2_corrected_interped_sla_2y_lowpass
                        % set grid
                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                        saveas(gcf,jpgname,'tif');

                        disp(' ')
                        disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
                        disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- corr_clim plot
        fig_flag=fig_flags{11,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_corr_clim_cmems_interped_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    corr_clim=ncread(filename,'corr_clim');
                    lon_cmems=ncread(filename,'lon_cmems');
                    lat_cmems=ncread(filename,'lat_cmems');

                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(lon_cmems',lat_cmems',squeeze(corr_clim(:,:,monthij))');
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');

                    tempcorr= corr_clim(:,:,monthij);
                    titlename = strcat('corr clim, ', num2str(monthij, '%02i'), ', ', testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                        num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(tempcorr(:), 'omitnan'),2)));  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(wrmap);
                    set(h,'fontsize',colorbar_fontsize);
        %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                    caxis([0.2 0.8]);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
            end
            fig_flag=0;
        end

% start-------------------- absolute trend plot
        fig_flag=fig_flags{12,2};
        while (fig_flag)
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                 if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                 end
                
%                 m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(interped_trend_filtered(:,:)'));
                mean_trend_filtered=mean(interped_trend_filtered(:),'omitnan');

%                 mean_trend_filtered=mean(trend_filtered(:),'omitnan');
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

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

                saveas(gcf,jpgname,'tif');

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
        
% start-------------------- absolute trend plot (bwrmap)
        fig_flag=fig_flags{13,2};
        while (fig_flag)
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_bwrmap_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                mean_trend_filtered=mean(trend_filtered(:),'omitnan');

        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([-4 4]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

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

% start-------------------- relative trend plot (bwrmap)
        fig_flag=fig_flags{14,2};
        while (fig_flag)
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_relative_ssh_trend_bwrmap_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                mean_trend_filtered=mean(trend_filtered(:),'omitnan');
                m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'-mean_trend_filtered));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([-4 4]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

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

%  start-------------------- absolute trend difference plot
        fig_flag=fig_flags{15,2};
        while (fig_flag)
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_diff_cmems_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                lon_cmems=ncread(cmemsfilename,'lon_cmems');
                lat_cmems=ncread(cmemsfilename,'lat_cmems');
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
                trend_diff = interped_trend_filtered-cmems_trend_filtered;
                m_pcolor(double(lon_cmems)',lat_cmems',trend_diff(:,:)');
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(trend_diff(:), 'omitnan'),2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis(trenddifflev);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

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

% start-------------------- absolute trend_corrected plot
        fig_flag=fig_flags{16,2};
        while (fig_flag)
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_absolute_cmems_corrected_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                correction_trend=mean(cmems_trend_filtered(:), 'omitnan')-mean_trend_filtered;
                m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)')+correction_trend);
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered+correction_trend,2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([0 7]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

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

% start-------------------- climatological absolute trend plot
        fig_flag=fig_flags{17,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(clim_trend(:,:,monthij)'.*1000));
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                    titlename = strcat('SSH trend(abs), ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                        ', ', calendarname{monthij}(1:3), '), ','M=',num2str(round(mean_clim_trend(monthij)*1000,2)),'mm/y');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(wrmap);
                    set(h,'fontsize',colorbar_fontsize);
    %                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
                    caxis(abstrendlev);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['clim_', num2str(monthij), '_', 'relative_ssh_trend', ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
            end
            close all;
            fig_flag=0;
        end

% start-------------------- msl time series (seasonal filtered)
        fig_flag=fig_flags{18,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                for varind=1:length(inputyear)*12
                    msl_filt(varind)=squeeze(mean(mean(ncread(modelfilename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
                end
        %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
                msl_filt=msl_filt-mean(msl_filt);    
                p2=polyfit(1:length(msl_filt),msl_filt,1);
                p2=p2(1)*1000.0*12.0;

                p=polyfit(xData,msl_filt,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl_filt,'k')
                hold on
                mslplot2=plot(xData,msl2,'Color','r')
                xlabel('year')
                ylabel('Mean SSH (m)')
                title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
            end
            close all;
            fig_flag=0;
        end

% start-------------------- msl time series
        fig_flag=fig_flags{19,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                for varind=1:length(inputyear)*12
                    msl(varind)=squeeze(mean(mean(ncread(modelfilename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
                end
                msl=msl-mean(msl);     
                p2=polyfit(1:length(msl),msl,1);
                p2=p2(1)*1000.0*12.0;

                p=polyfit(xData,msl,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl,'k')
                hold on
                mslplot2=plot(xData,msl2,'Color','r')
                xlabel('year')
                ylabel('Mean SSH (m)')
    %             mean_trend=mean(trend(:), 'omitnan');
                title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- msl time series (seasonal filtered) (corrected, SROCC)
        fig_flag=fig_flags{20,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_cmems_corrected_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                for varind=1:length(inputyear)*12
                    msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
                end
        %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));

                for tind=1:length(msl_filt)
                    msl_filt(tind)=msl_filt(tind)+(0.56 + 0.46 + 0.29 + 0.09)/1000.0/12.0*(tind-1);
                end
                msl_filt=msl_filt-mean(msl_filt);    
                p2=polyfit(1:length(msl_filt),msl_filt,1);
                p2=p2(1)*1000.0*12.0;


                p=polyfit(xData,msl_filt,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl_filt,'k')
                hold on
                mslplot2=plot(xData,msl2,'Color','r')
                xlabel('year')
                ylabel('Mean SSH (m)')
                title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
            end
            close all;
        end
        
% start-------------------- msl time series (corrected, SROCC)
        fig_flag=fig_flags{21,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_cmems_corrected_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                for varind=1:length(inputyear)*12
                    msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
                end
                for tind=1:length(msl)
                    msl(tind)=msl(tind)+(0.56 + 0.46 + 0.29 + 0.09)/1000.0/12.0*(tind-1);
                end
                msl=msl-mean(msl);     
                p2=polyfit(1:length(msl),msl,1);
                p2=p2(1)*1000.0*12.0;

                p=polyfit(xData,msl,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl,'k')
                hold on
                mslplot2=plot(xData,msl2,'Color','r')
                xlabel('year')
                ylabel('Mean SSH (m)')
                mean_trend=ncread(filename,'mean_trend');
                title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- lowpassed msl time series
        fig_flag=fig_flags{22,2};
        while (fig_flag)
            nyears=[2];
%                         nyears=[2,3,5];
            nc_varname_prefixes={'cmems_', 'interped_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                    jpgname=strcat(outfile, '_', testname, '_',regionname, '_cmems_',nc_varname,'_msl_', ...
                        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        for varind=1:length(inputyear)*12
                            eval(['tempmsl=squeeze(ncread(lpffilename,', '''', nc_varname, '''', ',[1,1,varind], [inf inf 1]));']);
                            msl(varind)=mean(tempmsl(:),'omitnan');
                            if (strcmp(nc_varname_prefix, 'interped_')==1)
                                temp_raw_msl=ncread(modelfilename,'raw_ssh',[1 1 varind], [inf inf 1]);
                            elseif (strcmp(nc_varname_prefix, 'cmems_')==1)
                                temp_raw_msl=ncread(cmemsfilename,'cmems_sla',[1 1 varind], [inf inf 1]);
                            end
                            raw_msl(varind)=squeeze(mean(temp_raw_msl(:),'omitnan'));
                        end
                        msl=msl-mean(msl, 'omitnan');    
                        raw_msl=raw_msl-mean(raw_msl(:));

                        p=polyfit(xData,msl,1);
                        msl2=xData*p(1)+p(2);
                        mslplot=plot(xData,msl,'k')
                        hold on
                        mslplot2=plot(xData,msl2,'Color','r')
                        mslplot3=plot(xData,raw_msl,'Color', [0.8 0.8 0.8])
                        xlabel('year')
                        ylabel('Mean SSH (m)')
    %                     mean_trend=ncread(filename,'mean_trend');
                        title([regionname, ', ', num2str(nyear), 'y low-passed MSL(', ...
                            num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
                        datetick('x','yyyy','keepticks')
                        axis tight;
                        ylim(meanplotlev)
                        set(mslplot,'LineWidth',2);
                        set(gca,'FontSize',15);
                        grid on
                        hold off
                        saveas(gcf,jpgname,'jpg');
                        grid off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- detrended lowpassed msl time series
        fig_flag=fig_flags{23,2};
        while (fig_flag)
            nyears=[2,3,5];
            nc_varname_prefixes={'cmems_detrended_', 'interped_detrended_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                    jpgname=strcat(outfile, '_', testname, '_',regionname, '_cmems_',nc_varname,'_msl_', ...
                        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        for varind=1:length(inputyear)*12
                            eval(['tempmsl=squeeze(ncread(filename,', '''', nc_varname, '''', ',[1,1,varind], [inf inf 1]));']);
                            msl(varind)=mean(tempmsl(:),'omitnan');
                            if (strcmp(nc_varname_prefix, 'interped_detrended_')==1)
                                temp_raw_msl=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
                            elseif (strcmp(nc_varname_prefix, 'cmems_detrended_')==1)
                                temp_raw_msl=ncread(filename,'cmems_sla',[1 1 varind], [inf inf 1]);
                            end
                            raw_msl(varind)=squeeze(mean(temp_raw_msl(:),'omitnan'));
                        end
                        msl=msl-mean(msl, 'omitnan');    
                        raw_msl=raw_msl-mean(raw_msl(:));

                        p=polyfit(xData,raw_msl,1);
                        msl2=xData*p(1)+p(2);
                        mslplot=plot(xData,msl,'k')
                        hold on
    %                     mslplot2=plot(xData,msl2,'Color','r')
                        mslplot2=plot(xData,raw_msl-msl2,'Color', [0.8 0.8 0.8])
                        xlabel('year')
                        ylabel('Mean SSH (m)')
    %                     mean_trend=ncread(filename,'mean_trend');
                        title([regionname, ', ', num2str(nyear), 'y low-passed det MSL(', ...
                            num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
                        datetick('x','yyyy','keepticks')
                        axis tight;
                        ylim(meanplotlev)
                        set(mslplot,'LineWidth',2);
                        set(gca,'FontSize',15);
                        grid on
                        hold off
                        saveas(gcf,jpgname,'jpg');
                        grid off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- climatological msl time series
        fig_flag=fig_flags{24,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij=1:12
                jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    if (exist('msl')==0)
                        for varind=1:length(inputyear)*12
                            msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
                        end
                        climmsl=reshape(msl,[12,length(inputyear)]);
                    else
                        climmsl=reshape(msl,[12,length(inputyear)]);
                    end
                    tempmsl=squeeze(climmsl(monthij,:));
                    tempmsl=tempmsl-mean(tempmsl);
                    p=polyfit(inputyear,tempmsl,1);
                    msl2=inputyear*p(1)+p(2);
                    mslplot=plot(inputyear,tempmsl,'k')
                    hold on
                    mslplot2=plot(inputyear,msl2,'Color','r')
                    xlabel('year')
                    ylabel('Mean SSH (m)')
                    title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                        ', ', calendarname{monthij}(1:3), '), ',num2str(round(mean_clim_trend(monthij)*1000,2)), ' mm/y'])
    %                 datetick('x','yyyy','keepticks')
                    axis tight;
                    ylim(meanplotlev)
                    set(mslplot,'LineWidth',2);
                    set(gca,'FontSize',15);
                    grid on
                    hold off
                    saveas(gcf,jpgname,'jpg');
                    grid off
                    close all;
                end
            end
            fig_flag=0;
        end

% start-------------------- climatological msl trend (line)
        fig_flag=fig_flags{25,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            clim_trend=ncread(modelfilename,'clim_ssh_trend');
            clim_trend_reshaped=reshape(clim_trend, [size(clim_trend,1)*size(clim_trend,2), 12]);
            mean_clim_trend=mean(clim_trend_reshaped,1,'omitnan');
            if abs(mean_clim_trend(1))<0.01
                mean_clim_trend=mean_clim_trend*1000.0  %% m -> mm
            end
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                mslplot=plot(1:12,mean_clim_trend,'k')
                hold on
                xlabel('month')
                ylabel('trend (mm/yr)')
                title([regionname, ', climatological trend(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(trendplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- seasonal msl time series
        fig_flag=fig_flags{26,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

    %         clim_ssh=ncread(filename,'clim_ssh');
    %         for ijij=1:12
    %             temp_clim_ssh=clim_ssh(:,:,ijij);
    %             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %         end
    %         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                comb_interped_data=ncread(interpedfilename,'interped_ssh');
                lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                lat_cmems=ncread(cmemsfilename, 'lat_cmems');
                len_lon=size(lon_cmems,1);
                len_lat=size(lat_cmems,2);
                comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
                clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
                cmems_trend=ncread(cmemsfilename, 'cmems_trend');
                cmems_mask=ones(size(cmems_trend));
                cmems_mask(isnan(cmems_trend))=NaN;
                clim_ssh = clim_ssh .* cmems_mask;
                for ijij=1:12
                    temp_clim_ssh=clim_ssh(:,:,ijij);
                    mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
                end
                mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

                mslplot=plot(1:12,mean_clim_ssh,'k')
                hold on
                xlabel('month')
                ylabel('trend (mm/yr)')
                title([regionname, ', seasonal msl(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- cmems seasonal msl
        fig_flag=fig_flags{27,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_cmems_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg

    %         clim_ssh=ncread(filename,'clim_ssh');
    %         for ijij=1:12
    %             temp_clim_ssh=clim_ssh(:,:,ijij);
    %             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %         end
    %         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                comb_interped_data=ncread(filename,'interped_ssh');
                lon_cmems=ncread(filename, 'lon_cmems');
                lat_cmems=ncread(filename, 'lat_cmems');
                len_lon=length(lon_cmems);
                len_lat=length(lat_cmems);
                clim_ssh=ncread(filename, 'clim_cmems_ssh');
                clim_ssh(clim_ssh>1000000)=NaN;
                cmems_trend=ncread(filename, 'cmems_trend');
                cmems_mask=ones(size(cmems_trend));
                cmems_mask(isnan(cmems_trend))=NaN;
                clim_ssh = clim_ssh .* cmems_mask;
                for ijij=1:12
                    temp_clim_ssh=clim_ssh(:,:,ijij);
                    mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
                end
                mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

                mslplot=plot(1:12,mean_clim_ssh,'k')
                hold on
                xlabel('month')
                ylabel('trend (mm/yr)')
                title([regionname, ', seasonal msl(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- cmems seasonal msl amplitude
        fig_flag=fig_flags{28,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_amp_cmems_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg

    %         clim_ssh=ncread(filename,'clim_ssh');
    %         for ijij=1:12
    %             temp_clim_ssh=clim_ssh(:,:,ijij);
    %             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %         end
    %         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                comb_interped_data=ncread(filename,'interped_ssh');
                lon_cmems=ncread(filename, 'lon_cmems');
                lat_cmems=ncread(filename, 'lat_cmems');
                len_lon=length(lon_cmems);
                len_lat=length(lat_cmems);
                clim_ssh=ncread(filename, 'clim_cmems_ssh');
                clim_ssh(clim_ssh>1000000)=NaN;
                cmems_trend=ncread(filename, 'cmems_trend');
                cmems_mask=ones(size(cmems_trend));
                cmems_mask(isnan(cmems_trend))=NaN;
                clim_ssh = clim_ssh .* cmems_mask;
                for ijij=1:12
                    temp_clim_ssh=clim_ssh(:,:,ijij);
                    mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
                end
                mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
                for loni=1:size(clim_ssh,1)
                    for lati=1:size(clim_ssh,2)
                        amp_clim_ssh(loni,lati)=(max(clim_ssh(loni,lati,:))-min(clim_ssh(loni,lati,:)))/2.0;
                    end
                end
    %             pcolor(amp_clim_ssh'*100)

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',squeeze(amp_clim_ssh(:,:)'*100.0));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('seasonal amp, ','cmems',',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(amp_clim_ssh(:)*100.0, 'omitnan'),2)));  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([2 16]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'cmems seasonal msl amplitude', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- seasonal msl amplitude
        fig_flag=fig_flags{29,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_amp_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg

    %         clim_ssh=ncread(filename,'clim_ssh');
    %         for ijij=1:12
    %             temp_clim_ssh=clim_ssh(:,:,ijij);
    %             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %         end
    %         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)

                lon_cmems=ncread(filename, 'lon_cmems');
                lat_cmems=ncread(filename, 'lat_cmems');
                len_lon=length(lon_cmems);
                len_lat=length(lat_cmems);
                comb_interped_data=ncread(filename,'interped_ssh');
                comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
                clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
                clim_ssh(clim_ssh>1000000)=NaN;
                cmems_trend=ncread(filename, 'cmems_trend');
                cmems_mask=ones(size(cmems_trend));
                cmems_mask(isnan(cmems_trend))=NaN;
                clim_ssh = clim_ssh .* cmems_mask;
                for ijij=1:12
                    temp_clim_ssh=clim_ssh(:,:,ijij);
                    mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
                end
                mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
                for loni=1:size(clim_ssh,1)
                    for lati=1:size(clim_ssh,2)
                        amp_clim_ssh(loni,lati)=(max(clim_ssh(loni,lati,:))-min(clim_ssh(loni,lati,:)))/2.0;
                    end
                end
    %             pcolor(amp_clim_ssh'*100)

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',squeeze(amp_clim_ssh(:,:)'*100.0));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('seasonal amp, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(amp_clim_ssh(:)*100.0, 'omitnan'),2)));  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([2 16]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'cmems seasonal msl amplitude', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- seasonal msl amplitude diff
        fig_flag=fig_flags{30,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_amp_diff', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)

                lon_cmems=ncread(filename, 'lon_cmems');
                lat_cmems=ncread(filename, 'lat_cmems');
                len_lon=length(lon_cmems);
                len_lat=length(lat_cmems);
                clim_cmems_ssh=ncread(filename, 'clim_cmems_ssh');
                clim_cmems_ssh(clim_cmems_ssh>1000000)=NaN;
                comb_interped_data=ncread(filename,'interped_ssh');
                comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
                clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
                clim_ssh(clim_ssh>1000000)=NaN;
                cmems_trend=ncread(filename, 'cmems_trend');
                cmems_mask=ones(size(cmems_trend));
                cmems_mask(isnan(cmems_trend))=NaN;
                clim_ssh = clim_ssh .* cmems_mask;
                clim_cmems_ssh=clim_cmems_ssh .* cmems_mask;

    %             for ijij=1:12
    %                 temp_clim_ssh=clim_ssh(:,:,ijij);
    %                 mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %             end
    %             mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
                for loni=1:size(clim_ssh,1)
                    for lati=1:size(clim_ssh,2)
                        amp_clim_ssh(loni,lati)=(max(clim_ssh(loni,lati,:))-min(clim_ssh(loni,lati,:)))/2.0;
                        amp_clim_cmems_ssh(loni,lati)=(max(clim_cmems_ssh(loni,lati,:))-min(clim_cmems_ssh(loni,lati,:)))/2.0;
                    end
                end
    %             pcolor(amp_clim_ssh'*100)

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',amp_clim_cmems_ssh(:,:)'*100-amp_clim_ssh(:,:)'*100);
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('seasonal amp diff, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(amp_clim_cmems_ssh(:)*100.0, 'omitnan')-mean(amp_clim_ssh(:)*100.0, 'omitnan'),2)));  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
    %             caxis([2 16]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'cmems seasonal msl amplitude', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- yearly msl time series
        fig_flag=fig_flags{31,2};
        while (fig_flag)
            for i =1:length(inputyear) 
                tempyear=inputyear(i);
        %         for month=1:12
                    xData2(i) = datenum([num2str(6,'%02i'),'-30-',num2str(tempyear)]);
        %         end
            end
            %     plot(squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan')))
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                comb_data_filtered=ncread(interpedfilename, 'interped_sla_filtered');
                isize = size(comb_data_filtered,1);
                jsize = size(comb_data_filtered,2);
                lsize = size(comb_data_filtered,3);
    %             comb_yearly_data_filtered=reshape(comb_data_filtered,[isize, jsize, 12, lsize/12]);
                comb_yearly_data_filtered=reshape(comb_data_filtered,[isize*jsize*12, lsize/12]);

    %             mean_yearly_data_filtered=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
                mean_yearly_data_filtered=squeeze(mean(comb_yearly_data_filtered,1,'omitnan'));

                trendtime=1:length(xData2);
                p=polyfit(trendtime,mean_yearly_data_filtered(1:length(xData2)),1);
                yearly_interped_trend=p(1);
                yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y

%                 yearly_msl=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
                yearly_msl=mean_yearly_data_filtered;

                yearly_msl=yearly_msl-mean(yearly_msl);    
                p=polyfit(xData2,yearly_msl,1);
                yearly_msl2=xData2*p(1)+p(2);
                yearly_mslplot=plot(xData2,yearly_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
                hold on
                yearly_mslplot2=plot(xData2,yearly_msl2,'Color','r')
                xlabel('Year')
                ylabel('Mean SSH (m)')
                title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(yearly_interped_trend,2)), ' mm/y'])
                ylim(meanplotlev)
                datetick('x','yymmm','keepticks')
                axis tight;
                ylim(meanplotlev)
                set(yearly_mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;            
        end

% start-------------------- cmems MSL plot
        fig_flag=fig_flags{32,2};
        while (fig_flag)
%             figdir2=[figrawdir,'CLIM\'];
%             if (exist(strcat(figdir2) , 'dir') ~= 7)
%                 mkdir(strcat(figdir2));
%             end 
%             outfile = strcat(figdir2,regionname);

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij;
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                            lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(cmemsfilename,'cmems_sla',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                        end
                        mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
                    end
                end
                cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                mean_data= mean_data .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' mean, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(m)','fontsize',colorbar_title_fontsize);
                caxis(colorbar_lev);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
            fig_flag=0;            
        end

% start-------------------- RMS plot
        fig_flag=fig_flags{33,2};
        while (fig_flag)

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_rms_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                            mean_model_data=mean_data;
                        end
                        mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
                        mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
                    end
                end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
                mean_data= mean_data .* mask_model;
                mean_model_data = mean_model_data .* mask_model;
                cmems_msl=mean(mean_data(:), 'omitnan');
                model_msl=mean(mean_model_data(:), 'omitnan');
%                 mean_data=mean_data-mean(mean_data(:),'omitnan');
%                 mean_model_data=mean_model_data-mean(mean_model_data(:),'omitnan');
%                 mean_rms = sqrt((mean_model_data - mean_data).^2) * 100.0;

                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('sq_diff' , 'var') ~= 1)
                            sq_diff=zeros(size(data));
                        end
                            sq_diff=sq_diff + ((model_data-model_msl)*100-(data-cmems_msl)*100).^2;
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_rms=sqrt(sq_diff./(length(inputyear) * length(inputmonth)));
                mean_rms=mean_rms.*mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([5 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data sq_diff
            end
            fig_flag=0;
        end

% start-------------------- RMS_corrected plot
        fig_flag=fig_flags{34,2};
        while (fig_flag)

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_rms_cmems_corrected_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'corrected_interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                            mean_model_data=mean_data;
                        end
                        mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
                        mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
                    end
                end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
                mean_data= mean_data .* mask_model;
                mean_model_data = mean_model_data .* mask_model;
                cmems_msl=mean(mean_data(:), 'omitnan');
                model_msl=mean(mean_model_data(:), 'omitnan');
%                 mean_data=mean_data-mean(mean_data(:),'omitnan');
%                 mean_model_data=mean_model_data-mean(mean_model_data(:),'omitnan');
%                 mean_rms = sqrt((mean_model_data - mean_data).^2) * 100.0;

                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'corrected_interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('sq_diff' , 'var') ~= 1)
                            sq_diff=zeros(size(data));
                        end
                            sq_diff=sq_diff + ((model_data-model_msl)*100-(data-cmems_msl)*100).^2;
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_rms=sqrt(sq_diff./(length(inputyear) * length(inputmonth)));
                mean_rms=mean_rms.*mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' rms',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_rms(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([5 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data sq_diff
            end
        end

% start-------------------- BIAS plot
        fig_flag=fig_flags{35,2};
        while (fig_flag)

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_bias_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_data' , 'var') ~= 1)
                            mean_data=zeros(size(data));
                            mean_model_data=mean_data;
                        end
                        mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
                        mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
                    end
                end
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
                mean_data= mean_data .* mask_model;
                mean_model_data = mean_model_data .* mask_model;
                cmems_msl=mean(mean_data(:), 'omitnan');
                model_msl=mean(mean_model_data(:), 'omitnan');
%                 mean_data=mean_data-mean(mean_data(:),'omitnan');
%                 mean_model_data=mean_model_data-mean(mean_model_data(:),'omitnan');
%                 mean_rms = sqrt((mean_model_data - mean_data).^2) * 100.0;

                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);
                    yearstr=num2str(tempyear, '%04i');
                    for monthij=1:length(inputmonth)
                        tempmonth=inputmonth(monthij);
                        monthstr=num2str(tempmonth, '%02i');
                        varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
                        data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        if (exist('mean_bias' , 'var') ~= 1)
                            mean_bias=zeros(size(data));
                        end
                            mean_bias=mean_bias + ((model_data-model_msl)*100-(data-cmems_msl)*100)/(length(inputyear) * length(inputmonth));
                    end
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;    

                mean_bias=mean_bias.*mask_model;
               
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_bias');
%                 m_pcolor(cut_lon_rho(1:30,1:20)',cut_lat_rho(1:30,1:20)',mean_bias(1:30,1:20)');


                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' bias',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'M=', num2str(round(mean(mean_bias(:),'omitnan'),1)));  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis(colorbar_lev);
                caxis([-15 15]);

            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data mean_bias
            end
        end

% start-------------------- cmems variance plot
        fig_flag=fig_flags{36,2};
        while (fig_flag)

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);


                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                    lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                              
                cmems_sla=ncread(cmemsfilename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                for loni=1:size(cmems_sla,1)
                    for lati=1:size(cmems_sla,2)
                        cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:).*100);
                    end
                end
%                 cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
%                 mean_data= mean_data .* mask_model;
                cmems_sla_var=cmems_sla_var .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',cmems_sla_var');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' VAR, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 200]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
            fig_flag=0;
        end
        
% start-------------------- cmems std plot
        fig_flag=fig_flags{37,2};
        while (fig_flag)

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_', 'STD','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
%                 for yearij=1:length(inputyear)
%                     tempyear=inputyear(yearij);
%                     yearstr=num2str(tempyear, '%04i');
%                     for monthij=1:length(inputmonth)
%                         tempmonth=inputmonth(monthij);
%                         monthstr=num2str(tempmonth, '%02i');
%                         varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                            lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
%                         data(:,:)=ncread(filename,'cmems_sla',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         if (exist('mean_data' , 'var') ~= 1)
%                             mean_data=zeros(size(data));
%                         end
%                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
%                     end
%                 end
                cmems_sla=ncread(cmemsfilename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                for loni=1:size(cmems_sla,1)
                    for lati=1:size(cmems_sla,2)
                        cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:).*100);
                    end
                end
%                 cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
%                 mean_data= mean_data .* mask_model;
                cmems_sla_var=cmems_sla_var .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',sqrt(cmems_sla_var)');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' STD, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 15]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
            fig_flag=0;
        end

% start-------------------- cmems variance_seasonal plot
        fig_flag=fig_flags{38,2};
        while (fig_flag)
           
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '_', num2str(tempmonth,'%02i'),'.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                    run(param_script);
    %                 for yearij=1:length(inputyear)
    %                     tempyear=inputyear(yearij);
    %                     yearstr=num2str(tempyear, '%04i');
    %                     for monthij=1:length(inputmonth)
    %                         tempmonth=inputmonth(monthij);
    %                         monthstr=num2str(tempmonth, '%02i');
    %                         varind=((yearij-1)*12)+monthij
                            if (exist('lon_min' , 'var') ~= 1)
                                lon_cmems=ncread(filename, 'lon_cmems');
                                lat_cmems=ncread(filename, 'lat_cmems');
                                [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            end
    %                         data_info = ncinfo(filename, varname); 
    %                         data(:,:)=ncread(filename,'cmems_sla',[lon_min(1) lat_min(1) varind], ...
    %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
    %                         if (exist('mean_data' , 'var') ~= 1)
    %                             mean_data=zeros(size(data));
    %                         end
    %                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
    %                     end
    %                 end
                    for yearij=1:length(inputyear)
                        cmems_sla(:,:,yearij)=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) (yearij-1)*12+tempmonth], ...
                                    [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                    end

                    for loni=1:size(cmems_sla,1)
                        for lati=1:size(cmems_sla,2)
                            cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:).*100);
                        end
                    end
    %                 cmems_msl=mean(mean_data(:), 'omitnan');
                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
    %                 mean_data= mean_data .* mask_model;
                    cmems_sla_var=cmems_sla_var .* mask_model;
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                    m_pcolor(cut_lon_rho',cut_lat_rho',cmems_sla_var');

                    shading(gca,m_pcolor_shading_method);   

                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat(variable, ' VAR, ','cmems,', calendarname{tempmonth}(1:3), ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
                    caxis([0 1]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    close all;
                    clear lon_rho mean_data
                end
            end
        end

% start-------------------- model variance plot
        fig_flag=fig_flags{39,2};
        while (fig_flag)
        

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_interped_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                    lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end

                model_sla=ncread(interpedfilename,'interped_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_var(loni,lati)=var(model_sla(loni,lati,:).*100); %% (m -> cm)
                    end
                end
%                 cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
%                 mean_data= mean_data .* mask_model;
                model_sla_var=model_sla_var .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_var');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' VAR, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 200]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
            fig_flag=0;
        end

% start-------------------- model variance_seasonal plot
        fig_flag=fig_flags{40,2};
        while (fig_flag)
        
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_interped_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'),'_', num2str(tempmonth,'%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                    run(param_script);
                    if (exist('lon_min' , 'var') ~= 1)
                        lon_cmems=ncread(filename, 'lon_cmems');
                        lat_cmems=ncread(filename, 'lat_cmems');
                        [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                        cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end

    %                 model_sla=ncread(filename,'interped_sla', [lon_min(1) lat_min(1) 1], ...
    %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
    %                         
                    for yearij=1:length(inputyear)
                        model_sla(:,:,yearij)=ncread(filename,'interped_sla', [lon_min(1) lat_min(1) (yearij-1)*12+tempmonth], ...
                                    [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                    end        

                    for loni=1:size(model_sla,1)
                        for lati=1:size(model_sla,2)
                            model_sla_var(loni,lati)=var(model_sla(loni,lati,:).*100);
                        end
                    end
    %                 cmems_msl=mean(mean_data(:), 'omitnan');
                    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                    hold on;

                    mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                    mask_model(mask_model==0)=NaN;
    %                 mean_data= mean_data .* mask_model;
                    model_sla_var=model_sla_var .* mask_model;
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
    %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                    m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_var');

                    shading(gca,m_pcolor_shading_method);   

                    m_gshhs_i('color',m_gshhs_line_color)  
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                    titlename = strcat(variable, ' VAR, ',testname,',', calendarname{tempmonth}(1:3), ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'(cm)','fontsize',colorbar_title_fontsize);
                    caxis([0 1]);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                    saveas(gcf,pngname,'tif');
                    close all;
                    clear lon_rho mean_data
                end
            end
        end

% start-------------------- cmems variance_lowpass plot
        fig_flag=fig_flags{41,2};
        while (fig_flag)
           
            nyears=[1:5];
            nc_varname_prefixes={'cmems_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                
                    pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_', nc_varname,'_','VAR_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['cmems_sla', '=ncread(filename,', '''', nc_varname,'''',');']);
%                         cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
%                                     [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                        cmems_sla=cmems_sla(:,:,nyear*6+1:end-(nyear*6)-1);
                        for loni=1:size(cmems_sla,1)
                            for lati=1:size(cmems_sla,2)
                                cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:));
                            end
                        end
        %                 cmems_msl=mean(mean_data(:), 'omitnan');
                        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                        hold on;

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
        %                 mean_data= mean_data .* mask_model;
                        cmems_sla_var=cmems_sla_var .* mask_model;
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                        m_pcolor(cut_lon_rho',cut_lat_rho',cmems_sla_var'.*100);

                        shading(gca,m_pcolor_shading_method);   

                        m_gshhs_i('color',m_gshhs_line_color)  
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                        titlename = strcat(nc_varname, ' VAR, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(jet);
                        set(h,'fontsize',colorbar_fontsize);
                        title(h,'(cm)','fontsize',colorbar_title_fontsize);
                        caxis([0 1]);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                        saveas(gcf,pngname,'tif');
                        close all;
                        clear lon_rho mean_data
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- model variance_lowpass plot
        fig_flag=fig_flags{42,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            nyears=[1:5];
            nc_varname_prefixes={'interped_sla_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                
                    pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', nc_varname,'_','VAR_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                            lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['model_sla', '=ncread(lpffilename,', '''', nc_varname,'''',');']);
%                         cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
%                                     [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                        model_sla=model_sla(:,:,nyear*6+1:end-(nyear*6)-1);
                        for loni=1:size(model_sla,1)
                            for lati=1:size(model_sla,2)
                                model_sla_var(loni,lati)=var(model_sla(loni,lati,:));
                            end
                        end
        %                 cmems_msl=mean(mean_data(:), 'omitnan');
                        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                        hold on;

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
        %                 mean_data= mean_data .* mask_model;
                        model_sla_var=model_sla_var .* mask_model;
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                        m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_var'.*100);

                        shading(gca,m_pcolor_shading_method);   

                        m_gshhs_i('color',m_gshhs_line_color)  
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                        titlename = strcat(nc_varname, ' VAR, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(jet);
                        set(h,'fontsize',colorbar_fontsize);
                        title(h,'(cm)','fontsize',colorbar_title_fontsize);
                        caxis([0 1]);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                        saveas(gcf,pngname,'tif');
                        close all;
                        clear lon_rho mean_data
                    end
                end
            end
            fig_flag=0;
        end
        
% start-------------------- model pattern correlation    
        fig_flag=fig_flags{43,2};
        while (fig_flag)
            cmems_adt=ncread(cmemsfilename,'cmems_adt');
            interped_ssh=ncread(interpedfilename,'interped_ssh');
            cmems_sla=ncread(cmemsfilename,'cmems_sla');
            interped_sla=ncread(interpedfilename,'interped_sla');
            gcm_interped_ssh = ncread(gcminterpedfilename, 'interped_ssh');
            m_gcm_interped_ssh = mean(gcm_interped_ssh,3,'omitnan');
            cut_lon_rho = ncread(cmemsfilename, 'lon_cmems');
            cut_lat_rho = ncread(cmemsfilename, 'lat_cmems');

            
%             month_season = [6,7,8];
%             for meani = 1:length(month_season)
%                 tempind{meani}= (inputyear-1982)*12 + month_season(meani);
%             end
%             ind_season = sort([tempind{:}]);
%             season_interped_sst = comb_interped_sst(:,:,ind_season);
%             comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
%             season_avhrr_sst = comb_avhrr_sst(:,:,ind_season);
%             mean_avhrr_sst = mean(season_avhrr_sst, 3);
            
            
            gcm_mask = m_gcm_interped_ssh;
            gcm_mask(isfinite(gcm_mask))=1;
            m_cmems_adt=mean(cmems_adt,3,'omitnan');
%             m_interped_ssh=mean(interped_ssh,3,'omitnan').*gcm_mask-0.25;
            m_interped_ssh=mean(interped_ssh,3,'omitnan');
            m_cmems_sla=mean(cmems_sla,3,'omitnan');
            m_interped_sla=mean(interped_sla,3,'omitnan');
            m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
            m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
            
            ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecs_khoapolygon);
            ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
            esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, es_khoapolygon);
            
            ysecs_m_interped_ssh = m_interped_ssh.*ysecsmask;
            ysecs_m_cmems_adt = m_cmems_adt.*ysecsmask;
            ysecs_m_gcm_interped_ssh = m_gcm_interped_ssh.*ysecsmask;

            es_m_interped_ssh = m_interped_ssh.*esmask;
            es_m_gcm_interped_ssh = m_gcm_interped_ssh.*esmask;
            es_m_cmems_adt = m_cmems_adt.*esmask;

% % % % % %             check 2d area for correlation
% % % % % %             m_interped_ssh_mask=m_interped_ssh;
% % % % % %             m_interped_ssh_mask(isfinite(m_interped_ssh_mask))=1;
% % % % % %             m_cmems_adt_mask=m_cmems_adt;
% % % % % %             m_cmems_adt_mask(isfinite(m_cmems_adt_mask))=1;
% % % % % %             pcolor((m_interped_ssh'-mean(m_interped_ssh', 'omitnan')).*m_interped_ssh_mask'.*m_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             m_gcm_interped_ssh_mask=m_gcm_interped_ssh;
% % % % % %             m_gcm_interped_ssh_mask(isfinite(m_gcm_interped_ssh_mask))=1;
% % % % % %             m_gcm_cmems_adt_mask=m_cmems_adt;
% % % % % %             m_gcm_cmems_adt_mask(isfinite(m_gcm_cmems_adt_mask))=1;
% % % % % %             pcolor((m_gcm_interped_ssh'-mean(m_gcm_interped_ssh', 'omitnan')).*m_gcm_interped_ssh_mask'.*m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             ysecs_m_interped_ssh_mask=ysecs_m_interped_ssh;
% % % % % %             ysecs_m_interped_ssh_mask(isfinite(ysecs_m_interped_ssh_mask))=1;
% % % % % %             ysecs_m_cmems_adt_mask=ysecs_m_cmems_adt;
% % % % % %             ysecs_m_cmems_adt_mask(isfinite(ysecs_m_cmems_adt_mask))=1;
% % % % % %             pcolor((ysecs_m_interped_ssh'-mean(ysecs_m_interped_ssh', 'omitnan')).*ysecs_m_interped_ssh_mask'.*ysecs_m_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             ysecs_m_gcm_interped_ssh_mask=ysecs_m_gcm_interped_ssh;
% % % % % %             ysecs_m_gcm_interped_ssh_mask(isfinite(ysecs_m_gcm_interped_ssh_mask))=1;
% % % % % %             ysecs_m_gcm_cmems_adt_mask=ysecs_m_cmems_adt;
% % % % % %             ysecs_m_gcm_cmems_adt_mask(isfinite(ysecs_m_gcm_cmems_adt_mask))=1;
% % % % % %             pcolor((ysecs_m_gcm_interped_ssh'-mean(ysecs_m_gcm_interped_ssh', 'omitnan')).*ysecs_m_gcm_interped_ssh_mask'.*ysecs_m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             pcolor((ysecs_m_cmems_adt'-mean(ysecs_m_cmems_adt', 'omitnan')).*ysecs_m_gcm_interped_ssh_mask'.*ysecs_m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % % 
% % % % % %             es_m_interped_ssh_mask=es_m_interped_ssh;
% % % % % %             es_m_interped_ssh_mask(isfinite(es_m_interped_ssh_mask))=1;
% % % % % %             es_m_cmems_adt_mask=es_m_cmems_adt;
% % % % % %             es_m_cmems_adt_mask(isfinite(es_m_cmems_adt_mask))=1;
% % % % % %             pcolor((es_m_interped_ssh'-mean(es_m_interped_ssh', 'omitnan')).*es_m_interped_ssh_mask'.*es_m_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             es_m_gcm_interped_ssh_mask=es_m_gcm_interped_ssh;
% % % % % %             es_m_gcm_interped_ssh_mask(isfinite(es_m_gcm_interped_ssh_mask))=1;
% % % % % %             es_m_gcm_cmems_adt_mask=es_m_cmems_adt;
% % % % % %             es_m_gcm_cmems_adt_mask(isfinite(es_m_gcm_cmems_adt_mask))=1;
% % % % % %             pcolor((es_m_gcm_interped_ssh'-mean(es_m_gcm_interped_ssh', 'omitnan')).*es_m_gcm_interped_ssh_mask'.*es_m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % % 
% % % % % %             es_m_gcm_cmems_adt_mask=es_m_cmems_adt;
% % % % % %             es_m_gcm_cmems_adt_mask(isfinite(es_m_gcm_cmems_adt_mask))=1;
% % % % % %             pcolor((es_m_cmems_adt'-mean(es_m_cmems_adt', 'omitnan')).*es_m_gcm_interped_ssh_mask'.*es_m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);


            
            
%              pcolor(m_interped_ssh'-0.25); shading flat; colorbar; caxis([-0.2 1.0]); colormap jet
             aaa= movmean(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))),1);
             bbb= movmean(m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))),1);
             corrcoef(aaa,bbb);
             ccc= movmean(m_gcm_interped_ssh(logical(~isnan(m_gcm_interped_ssh).*~isnan(m_cmems_adt))),1);
             ddd= movmean(m_cmems_adt(logical(~isnan(m_gcm_interped_ssh).*~isnan(m_cmems_adt))),1);
             corrcoef(ccc,ddd);
             
             eee= ysecs_m_interped_ssh(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
             fff= ysecs_m_cmems_adt(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
             ysecscorr=corrcoef(eee,fff)
             ysecscorr_rank = corr(eee,fff, 'Type', 'Spearman')
             ysecscorr_tau = corr(eee,fff, 'Type', 'Kendall')
             
             eeeg= ysecs_m_interped_ssh(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt).*~isnan(ysecs_m_gcm_interped_ssh)));
             fffg= ysecs_m_cmems_adt(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt).*~isnan(ysecs_m_gcm_interped_ssh)));
             rcmg_ysecscorr=corrcoef(eeeg,fffg)
             rcmg_ysecscorr_rank = corr(eeeg,fffg, 'Type', 'Spearman')
             rcmg_ysecscorr_tau = corr(eeeg,fffg, 'Type', 'Kendall')
             
             
             gg= es_m_interped_ssh(logical(~isnan(es_m_interped_ssh).*~isnan(es_m_cmems_adt)));
             hh= es_m_cmems_adt(logical(~isnan(es_m_interped_ssh).*~isnan(es_m_cmems_adt)));
             escorr=corrcoef(gg,hh)
             escorr_rank = corr(gg,hh, 'Type', 'Spearman')
%              ysecsescorr=corrcoef([eee; gg], [fff; hh])
                
             ggg= es_m_interped_ssh(logical(~isnan(es_m_interped_ssh).*~isnan(es_m_cmems_adt).*~isnan(es_m_gcm_interped_ssh)));
             hhg= es_m_cmems_adt(logical(~isnan(es_m_interped_ssh).*~isnan(es_m_cmems_adt).*~isnan(es_m_gcm_interped_ssh)));
             rcmg_escorr=corrcoef(ggg,hhg)
             [rcmg_escorr_tau, pval]=corr(ggg,hhg, 'Type', 'Kendall')

             ii= ysecs_m_gcm_interped_ssh(logical(~isnan(ysecs_m_gcm_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
             jj= ysecs_m_cmems_adt(logical(~isnan(ysecs_m_gcm_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
             gcmysecscorr=corrcoef(ii,jj)
             gcmysecscorr_rank = corr(ii,jj, 'Type', 'Spearman')
             gcmysecscorr_tau = corr(ii,jj, 'Type', 'Kendall')
             
             
             kk= es_m_gcm_interped_ssh(logical(~isnan(es_m_gcm_interped_ssh).*~isnan(es_m_cmems_adt)));
             ll= es_m_cmems_adt(logical(~isnan(es_m_gcm_interped_ssh).*~isnan(es_m_cmems_adt)));
             gcmescorr=corrcoef(kk,ll)
             gcmescorr_rank = corr(kk,ll, 'Type', 'Spearman')
             [gcmescorr_tau, pval] = corr(kk,ll, 'Type', 'Kendall')
%              gcmysecsescorr=corrcoef([ii; kk], [jj; ll])
           
            plot(aaa-mean(aaa)); hold on; plot(bbb-mean(bbb)); plot(ccc-mean(ccc)-0.5); plot(ddd-mean(ddd)-0.5); legend('rcm', 'sat', 'gcm', 'sat-gcm'); hold off
            plot(eee-mean(eee)); hold on; plot(fff-mean(fff)); plot(ii-mean(ii)-0.5); plot(jj-mean(jj)-0.5); legend('ysecs-rcm', 'ysecs-sat', 'ysecs-gcm', 'ysecs-sat-gcm'); hold off
            plot(gg-mean(gg)); hold on; plot(hh-mean(hh)); plot(kk-mean(kk)-0.5); plot(ll-mean(ll)-0.5); legend('es-rcm', 'es-sat', 'es-gcm', 'es-sat-gcm'); hold off
            
%             corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))))
            corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))))
            corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))))
%             corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
%             cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)

        % end-------------------- model pattern correlation


%         % start-------------------- model pattern correlation (yearly)
%             pngname=strcat(outfile, '_', testname,'_',regionname, '_', nc_varname,'_','yearly_pattern_corr_',num2str(min(inputyear),'%04i'), ...
%                                 '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%             if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
%                 for yearij=1:length(inputyear)
%                     tempyear=inputyear(yearij);
% 
%                     yearly_cmems_adt=ncread(filename,'cmems_adt',[1 1 (yearij-1)*12+1], [inf inf 12]);
%                     yearly_interped_ssh=ncread(filename,'interped_ssh',[1 1 (yearij-1)*12+1], [inf inf 12]);
%                     yearly_cmems_sla=ncread(filename,'cmems_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);
%                     yearly_interped_sla=ncread(filename,'interped_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);
% 
%                     m_cmems_adt=mean(yearly_cmems_adt,3,'omitnan');
%                     m_interped_ssh=mean(yearly_interped_ssh,3,'omitnan');
%                     m_cmems_sla=mean(yearly_cmems_sla,3,'omitnan');
%                     m_interped_sla=mean(yearly_interped_sla,3,'omitnan');
%         %             m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
%         %             m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
% 
%                     ssh_corr(yearij) = corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))));
%                     sla_corr(yearij) = corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))));
%         %             corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
%         %                 cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)
%                 end
%                 plot(ssh_corr);
%                 hold on
%                 plot(sla_corr);
%                 hold off
%                 close all
%             end
            fig_flag=0;
        end

% start-------------------- corr_ movmean plot
        fig_flag=fig_flags{44,2};
        while (fig_flag)
            nyears=[1,2,5];
    %         nc_varname_prefixes={'interped_', 'interped_detrended_'};
                    nc_varname_prefixes={'interped_'};

    %         nc_titlename_prefixes={'movmean', 'movmean-det'};
                    nc_titlename_prefixes={'movmean'};

            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];
                    jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_cmems_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        eval(['corr_',nc_varname, '=ncread(mov_corrfilename,','''','corr_', nc_varname,'''',');']);
                        eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);

                        lon_cmems=ncread(cmemsfilename,'lon_cmems');
                        lat_cmems=ncread(cmemsfilename,'lat_cmems');

                        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                        hold on;
                        eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                        shading(gca,m_pcolor_shading_method);
                        m_gshhs_i('color',m_gshhs_line_color);
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                        eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        tit_prefix=nc_titlename_prefixes{nc_varnameij};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr,2)));  

                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(wrmap);
                        set(h,'fontsize',colorbar_fontsize);
                        caxis([0.2 0.8]);

                        % set grid
                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                        saveas(gcf,jpgname,'tif');

                        disp(' ')
                        disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
                        disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end
        
% start-------------------- corr_ movmean plot (corrected)
        fig_flag=fig_flags{45,2};
        while (fig_flag)
            nyears=[1,2,5];
            nc_varname_prefixes={'corrected_interped_sla_'};
            nc_titlename_prefixes={'movmean-c'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];
                    jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_cmems_movmean_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        eval(['corr_',nc_varname, '=ncread(mov_corrfilename,','''','corr_', nc_varname,'''',');']);
                        eval(['corr_',nc_varname, '=corr_', nc_varname, '.*cmems_mask;']);

                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

                        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                        hold on;
                        eval(['m_pcolor(lon_cmems', '''', ',lat_cmems', '''',',squeeze(corr_', nc_varname, '(:,:)', '''', '));']);
                        shading(gca,m_pcolor_shading_method);
                        m_gshhs_i('color',m_gshhs_line_color);
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                        eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        tit_prefix=nc_titlename_prefixes{nc_varnameij};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mcorr,2)));  

                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(wrmap);
                        set(h,'fontsize',colorbar_fontsize);
                        caxis([0.2 0.8]);

                        % set grid
                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                        saveas(gcf,jpgname,'tif');

                        disp(' ')
                        disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
                        disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end
       
% start-------------------- movmean msl time series
        fig_flag=fig_flags{46,2};
        while (fig_flag)
%             nyears=[2];
                        nyears=[1,2,3,5];
            nc_varname_prefixes={'cmems_', 'interped_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];
                    jpgname=strcat(outfile, '_', testname, '_',regionname, '_cmems_',nc_varname,'_msl_', ...
                        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        for varind=1:length(inputyear)*12
                            eval(['tempmsl=squeeze(ncread(movfilename,', '''', nc_varname, '''', ',[1,1,varind], [inf inf 1]));']);
                            msl(varind)=mean(tempmsl(:),'omitnan');
                            if (strcmp(nc_varname_prefix, 'interped_')==1)
                                temp_raw_msl=ncread(modelfilename,'raw_ssh',[1 1 varind], [inf inf 1]);
                            elseif (strcmp(nc_varname_prefix, 'cmems_')==1)
                                temp_raw_msl=ncread(cmemsfilename,'cmems_sla',[1 1 varind], [inf inf 1]);
                            end
                            raw_msl(varind)=squeeze(mean(temp_raw_msl(:),'omitnan'));
                        end
                        msl=msl-mean(msl, 'omitnan');    
                        raw_msl=raw_msl-mean(raw_msl(:));

                        p=polyfit(xData,msl,1);
                        msl2=xData*p(1)+p(2);
                        mslplot=plot(xData,msl,'k')
                        hold on
                        mslplot2=plot(xData,msl2,'Color','r')
                        mslplot3=plot(xData,raw_msl,'Color', [0.8 0.8 0.8])
                        xlabel('year')
                        ylabel('Mean SSH (m)')
    %                     mean_trend=ncread(filename,'mean_trend');
                        title([regionname, ', ', num2str(nyear), 'y movmean MSL(', ...
                            num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
                        datetick('x','yyyy','keepticks')
                        axis tight;
                        ylim(meanplotlev)
                        set(mslplot,'LineWidth',2);
                        set(gca,'FontSize',15);
                        grid on
                        hold off
                        saveas(gcf,jpgname,'jpg');
                        grid off
                        close all;
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- cmems variance_movmean plot
        fig_flag=fig_flags{47,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            nyears=[1,2,5];
            nc_varname_prefixes={'cmems_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];
                
                    pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', nc_varname,'_','VAR_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                            lat_cmems=ncread(cmemsfilename, 'lat_cmems');
        %                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['cmems_sla', '=ncread(movfilename,', '''', nc_varname,'''',');']);
%                         cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
%                                     [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                        cmems_sla=cmems_sla(:,:,nyear*6+1:end-(nyear*6)-1);
                        for loni=1:size(cmems_sla,1)
                            for lati=1:size(cmems_sla,2)
                                cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:).*100.0);
                            end
                        end
        %                 cmems_msl=mean(mean_data(:), 'omitnan');
                        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                        hold on;

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
        %                 mean_data= mean_data .* mask_model;
                        cmems_sla_var=cmems_sla_var .* mask_model;
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                        m_pcolor(cut_lon_rho',cut_lat_rho',cmems_sla_var');

                        shading(gca,m_pcolor_shading_method);   

                        m_gshhs_i('color',m_gshhs_line_color)  
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                        titlename = strcat(nc_varname, ' VAR, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(jet);
                        set(h,'fontsize',colorbar_fontsize);
                        title(h,'(cm)','fontsize',colorbar_title_fontsize);
                        caxis([0 60]);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                        saveas(gcf,pngname,'tif');
                        close all;
                        clear lon_rho mean_data
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- model variance_movmean plot
        fig_flag=fig_flags{48,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            nyears=[1,2,5];
            nc_varname_prefixes={'interped_sla_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];
                
                    pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', nc_varname,'_','VAR_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                            lat_cmems=ncread(cmemsfilename, 'lat_cmems');
        %                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['model_sla', '=ncread(movfilename,', '''', nc_varname,'''',');']);
%                         cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
%                                     [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                        model_sla=model_sla(:,:,nyear*6+1:end-(nyear*6)-1);
                        for loni=1:size(model_sla,1)
                            for lati=1:size(model_sla,2)
                                model_sla_var(loni,lati)=var(model_sla(loni,lati,:).*100.0);
                            end
                        end
        %                 cmems_msl=mean(mean_data(:), 'omitnan');
                        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                        hold on;

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
        %                 mean_data= mean_data .* mask_model;
                        model_sla_var=model_sla_var .* mask_model;
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                        m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_var');

                        shading(gca,m_pcolor_shading_method);   

                        m_gshhs_i('color',m_gshhs_line_color)  
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                        titlename = strcat(nc_varname, ' VAR, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(jet);
                        set(h,'fontsize',colorbar_fontsize);
                        title(h,'(cm)','fontsize',colorbar_title_fontsize);
                        caxis([0 60]);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                        saveas(gcf,pngname,'tif');
                        close all;
                        clear lon_rho mean_data
                    end
                end
            end
            fig_flag=0;
        end
        
% start-------------------- cmems_relative trend plot (bwrmap)
        fig_flag=fig_flags{49,2};
        while (fig_flag)
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_cmems_relative_ssh_trend_bwrmap_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                    lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                        
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                mean_trend_filtered=mean(cmems_trend_filtered(:),'omitnan');
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(cmems_trend_filtered(:,:)'-mean_trend_filtered));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([-4 4]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

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

% start-------------------- cmems absolute trend plot (bwrmap)
        fig_flag=fig_flags{50,2};
        while (fig_flag)
            outfile = strcat(figdir,regionname);
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_cmems_absolute_ssh_trend_bwrmap_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                    lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(cmems_trend_filtered(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                mean_trend_filtered=mean(cmems_trend_filtered(:),'omitnan');

        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ','CMEMS', ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([-4 4]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

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
        
% start-------------------- cmems seasonal msl time series
        fig_flag=fig_flags{53,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', 'cmems', '_',regionname, '_seasonal_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

    %         clim_ssh=ncread(filename,'clim_ssh');
    %         for ijij=1:12
    %             temp_clim_ssh=clim_ssh(:,:,ijij);
    %             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %         end
    %         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                comb_interped_data=ncread(cmemsfilename,'cmems_sla');
                lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                lat_cmems=ncread(cmemsfilename, 'lat_cmems');
                len_lon=size(lon_cmems,1);
                len_lat=size(lat_cmems,2);
                comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
                clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
                cmems_trend=ncread(cmemsfilename, 'cmems_trend');
                cmems_mask=ones(size(cmems_trend));
                cmems_mask(isnan(cmems_trend))=NaN;
                clim_ssh = clim_ssh .* cmems_mask;
                for ijij=1:12
                    temp_clim_ssh=clim_ssh(:,:,ijij);
                    mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
                end
                mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

                mslplot=plot(1:12,mean_clim_ssh,'k')
                hold on
                xlabel('month')
                ylabel('SSH (m)')
                title([regionname, ', seasonal msl(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end
        
% start-------------------- model std_movmean plot
        fig_flag=fig_flags{54,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            nyears=[1,2,5];
            nc_varname_prefixes={'interped_sla_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];
                
                    pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_interped_', nc_varname,'_','std_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(interpedfilename, 'lon_cmems');
                            lat_cmems=ncread(interpedfilename, 'lat_cmems');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['model_sla', '=ncread(movfilename,', '''', nc_varname,'''',');']);

                        model_sla=model_sla(:, :, nyear*(12*nyear/2)+1 : end-(12*nyear/2)-1);
                        for loni=1:size(model_sla,1)
                            for lati=1:size(model_sla,2)
                                model_sla_std(loni,lati)=std(model_sla(loni,lati,:).*100.0);
                            end
                        end
                        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                        hold on;

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        model_sla_std=model_sla_std .* mask_model;
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                        m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_std');

                        shading(gca,m_pcolor_shading_method);   

                        m_gshhs_i('color',m_gshhs_line_color)  
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                        titlename = strcat(nc_varname, ' VAR, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(jet);
                        set(h,'fontsize',colorbar_fontsize);
                        title(h,'(cm)','fontsize',colorbar_title_fontsize);
                        caxis([0 15]);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                        saveas(gcf,pngname,'tif');
                        close all;
                        clear lon_rho model_sla model_sla_std
                    end
                end
            end
            fig_flag=0;
        end
        
% start-------------------- cmems std_movmean plot
        fig_flag=fig_flags{55,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            nyears=[1,2,5];
            nc_varname_prefixes={'cmems_'};
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];
                
                    pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', nc_varname,'_','std_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                            lat_cmems=ncread(cmemsfilename, 'lat_cmems');
        %                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['cmems_sla', '=ncread(movfilename,', '''', nc_varname,'''',');']);
%                         cmems_sla=ncread(filename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
%                                     [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                        cmems_sla=cmems_sla(:,:,nyear*6+1:end-(nyear*6)-1);
                        for loni=1:size(cmems_sla,1)
                            for lati=1:size(cmems_sla,2)
                                cmems_sla_var(loni,lati)=std(cmems_sla(loni,lati,:).*100.0);
                            end
                        end
        %                 cmems_msl=mean(mean_data(:), 'omitnan');
                        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                        hold on;

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
        %                 mean_data= mean_data .* mask_model;
                        cmems_sla_var=cmems_sla_var .* mask_model;
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
        %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                        m_pcolor(cut_lon_rho',cut_lat_rho',cmems_sla_var');

                        shading(gca,m_pcolor_shading_method);   

                        m_gshhs_i('color',m_gshhs_line_color)  
                        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                        titlename = strcat(nc_varname, ' std, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        % set colorbar 
                        h = colorbar;
                        colormap(jet);
                        set(h,'fontsize',colorbar_fontsize);
                        title(h,'(cm)','fontsize',colorbar_title_fontsize);
                        caxis([0 15]);

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                        saveas(gcf,pngname,'tif');
                        close all;
                        clear lon_rho mean_data
                    end
                end
            end
            fig_flag=0;
        end

% start-------------------- model std plot
        fig_flag=fig_flags{56,2};
        while (fig_flag)

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_interped_', 'std','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                    lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end

                model_sla=ncread(interpedfilename,'interped_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_std(loni,lati)=std(model_sla(loni,lati,:).*100); %% (m -> cm)
                    end
                end
%                 cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
%                 mean_data= mean_data .* mask_model;
                model_sla_std=model_sla_std .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_std');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' std, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 15]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
            fig_flag=0;
        end
 
% start-------------------- model std_yearly_detrended plot
        fig_flag=fig_flags{57,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_interped_', 'yearly_detrended','_','std_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                model_sla=ncread(interpedfilename,'interped_sla_yearly_detrended');

                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_std(loni,lati)=std(model_sla(loni,lati,:));
                    end
                end
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                model_sla_std=model_sla_std .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_std');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('yearly sla detrened', ' std, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho model_sla model_sla_std
            end

            fig_flag=0;
        end
 
% start-------------------- cmems std_yearly_detrended plot
        fig_flag=fig_flags{58,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'yearly_detrended','_','std_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                model_sla=ncread(cmemsfilename,'cmems_sla_yearly_detrended');

                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_std(loni,lati)=std(model_sla(loni,lati,:));
                    end
                end
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                model_sla_std=model_sla_std .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_std');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('yearly sla detrened', ' std, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho model_sla model_sla_std
            end

            fig_flag=0;
        end
        
% start-------------------- SLR plot from yearly exponential fit data
        fig_flag=fig_flags{59,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, 'yearly_SLR_exp_fit_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=ncread(interpedfilename,'interped_sla_yearly_exp_fit');

                slr=sl(:,:,end)-sl(:,:,1);
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

                slr=slr.*mask_model;
%                 pcolor(sl');
                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('SLR-exp fit , ',testname, ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)), ') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{59,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end
        
% start-------------------- SLR plot from yearly poly1 fit data
        fig_flag=fig_flags{60,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, 'yearly_SLR_poly1_fit_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=ncread(interpedfilename,'interped_sla_yearly_poly1_fit');

                slr=sl(:,:,end)-sl(:,:,1);
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

                slr=slr.*mask_model;
%                 pcolor(sl');
                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('SLR-poly1 fit , ',testname, ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{60,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end
        
% start-------------------- SLR r plot from yearly exp fit
        fig_flag=fig_flags{61,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, 'yearly_SLR_exp_fit_r_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=sqrt(ncread(interpedfilename,'interped_sla_yearly_exp_fit_rsquare'));

                slr=sl;
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

                slr=slr.*mask_model;
%                 pcolor(sl');
                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('r -exp fit , ',testname, ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'( )','fontsize',colorbar_title_fontsize);
                caxis([0 1]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{61,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end
   
% start-------------------- SLR r plot from yearly poly1 fit
        fig_flag=fig_flags{62,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, 'yearly_SLR_poly1_fit_r_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=sqrt(ncread(interpedfilename,'interped_sla_yearly_poly1_fit_rsquare'));

                slr=sl;
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

                slr=slr.*mask_model;
%                 pcolor(sl');
                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('r -poly1 fit , ',testname, ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'( )','fontsize',colorbar_title_fontsize);
                caxis([0 1]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{62,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end
        
% start-------------------- SLR rmse plot from yearly exp fit
        fig_flag=fig_flags{63,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, 'yearly_SLR_exp_fit_rmse_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=ncread(interpedfilename,'interped_sla_yearly_exp_fit_rmse');

                slr=sl;
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

                slr=slr.*mask_model;
%                 pcolor(sl');
                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('RMSE-exp fit , ',testname, ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis([0 15]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{63,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end       
        
% start-------------------- SLR rmse plot from yearly poly1 fit
        fig_flag=fig_flags{64,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, 'yearly_SLR_poly1_fit_rmse_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=ncread(interpedfilename,'interped_sla_yearly_poly1_fit_rmse');

                slr=sl;
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

                slr=slr.*mask_model;
%                 pcolor(sl');
                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('RMSE-poly1 fit , ',testname, ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
%                 caxis([0 15]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{64,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end               

% start-------------------- cmems std_yearly_exp_fit_detrended plot
        fig_flag=fig_flags{65,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'yearly_exp_fit_detrended','_','std_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                model_sla=ncread(cmemsfilename,'cmems_sla_yearly_exp_fit_detrended');

                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_std(loni,lati)=std(model_sla(loni,lati,:));
                    end
                end
                cmems_land=ones(size(model_sla));
                cmems_land(isnan(model_sla))=1;
                cmems_land(isfinite(model_sla))=NaN;
                save([filedir,regionname, '_', testname, '_cmems_std_land','.mat'], 'cmems_land');
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                model_sla_std=model_sla_std .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));

                yearstr_min=num2str(inputyear(1));
                yearstr_max=num2str(inputyear(end));
                save([filedir,regionname, '_cmems_', 'det_std','_',yearstr_min, '_', yearstr_max, '.mat'], 'cut_lon_rho', 'cut_lat_rho', 'model_sla_std');
                
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_std');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('y-sla-detrened', ' std, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(model_sla_std(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho model_sla model_sla_std
            end

            fig_flag=0;
        end

% start-------------------- model std_yearly_exp_fit_detrended plot
        fig_flag=fig_flags{66,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_interped_', 'yearly_exp_fit_detrended','_','std_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                model_sla=ncread(interpedfilename,'interped_sla_yearly_exp_fit_detrended');

                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_std(loni,lati)=std(model_sla(loni,lati,:));
                    end
                end
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                model_sla_std=model_sla_std .* mask_model;
                
                yearstr_min=num2str(inputyear(1));
                yearstr_max=num2str(inputyear(end));
                save([filedir,regionname,'_',testname, '_model_', 'det_std','_',yearstr_min, '_', yearstr_max, '.mat'], 'cut_lon_rho', 'cut_lat_rho', 'model_sla_std');
   
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_std');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('y-sla-exp det', ' std, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(model_sla_std(:),'omitnan'),2)),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho model_sla model_sla_std
            end

            fig_flag=0;
        end
    
% start-------------------- cmems SLR plot from yearly exponential fit data
        fig_flag=fig_flags{67,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', 'cmems','_',regionname, 'yearly_SLR_exp_fit_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=ncread(cmemsfilename,'cmems_sla_yearly_exp_fit');

                slr=sl(:,:,end)-sl(:,:,1);
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

                slr=slr.*mask_model;
%                 pcolor(sl');
                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('SLR-exp fit , ','cmems', ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)), ') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{67,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end
        
% start-------------------- cmems SLR plot from yearly poly1 fit data
        fig_flag=fig_flags{68,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', 'cmems','_',regionname, 'yearly_SLR_poly1_fit_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=ncread(cmemsfilename,'cmems_sla_yearly_poly1_fit');

                slr=sl(:,:,end)-sl(:,:,1);
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;

                slr=slr.*mask_model;
%                 pcolor(sl');
                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('SLR-poly1 fit , ','cmems', ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{68,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end
        
% start-------------------- cmems std_yearly_detrended plot
        fig_flag=fig_flags{69,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'yearly_detrended','_','std_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                model_sla=ncread(cmemsfilename,'cmems_sla_yearly_detrended');

                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_std(loni,lati)=std(model_sla(loni,lati,:));
                    end
                end
                cmems_land=ones(size(model_sla));
                cmems_land(isnan(model_sla))=1;
                cmems_land(isfinite(model_sla))=NaN;
                save([filedir,regionname, '_', testname, '_cmems_std_land','.mat'], 'cmems_land');
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                model_sla_std=model_sla_std .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));

                yearstr_min=num2str(inputyear(1));
                yearstr_max=num2str(inputyear(end));
                save([filedir,regionname, '_cmems_linear_', 'det_std','_',yearstr_min, '_', yearstr_max, '.mat'], 'cut_lon_rho', 'cut_lat_rho', 'model_sla_std');
                
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_std');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('y-sla-detrened', ' std, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(model_sla_std(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho model_sla model_sla_std
            end

            fig_flag=0;
        end
        
% start-------------------- model std_yearly_detrended plot
        fig_flag=fig_flags{70,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_interped_', 'yearly_detrended','_','std_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                model_sla=ncread(interpedfilename,'interped_sla_yearly_detrended');

                for loni=1:size(model_sla,1)
                    for lati=1:size(model_sla,2)
                        model_sla_std(loni,lati)=std(model_sla(loni,lati,:));
                    end
                end
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                model_sla_std=model_sla_std .* mask_model;
                
                yearstr_min=num2str(inputyear(1));
                yearstr_max=num2str(inputyear(end));
                save([filedir,regionname,'_',testname, '_model_linear_', 'det_std','_',yearstr_min, '_', yearstr_max, '.mat'], 'cut_lon_rho', 'cut_lat_rho', 'model_sla_std');
   
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',model_sla_std');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('y-sla-exp det', ' std, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(model_sla_std(:),'omitnan'),2)),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 10]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho model_sla model_sla_std
            end

            fig_flag=0;
        end

% start-------------------- u_rms
        fig_flag=fig_flags{71,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'u_','rms_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                rms_u_rho=ncread(cmemsvelfilename,'rms_u_rho');
                comb_interped_u_rho = ncread(cmemsvelfilename, 'interped_u_rho');
                mean_interped_u_rho = mean(comb_interped_u_rho, 3);
                comb_cmems_u_rho = ncread(cmemsvelfilename, 'cmems_u_rho');
                mean_cmems_u_rho = mean(comb_cmems_u_rho, 3);
                
                
                gcm_rms_u_rho = ncread(gcmcmemsvelfilename,'rms_u_rho');
                gcm_rms_mask = gcm_rms_u_rho;
                gcm_rms_mask(isfinite(gcm_rms_u_rho))=1;
                
                u_rho_corr = corrcoef(mean_interped_u_rho(logical(~isnan(mean_interped_u_rho).*~isnan(mean_cmems_u_rho))), mean_cmems_u_rho(logical(~isnan(mean_interped_u_rho).*~isnan(mean_cmems_u_rho))))
                u_rho_gcmcorr = corrcoef(mean_interped_u_rho(logical(~isnan(mean_interped_u_rho).*~isnan(mean_cmems_u_rho).*~isnan(gcm_rms_mask))), mean_cmems_u_rho(logical(~isnan(mean_interped_u_rho).*~isnan(mean_cmems_u_rho).*~isnan(gcm_rms_mask))))
                
                
                m_interped_u_rho = mean_interped_u_rho;
                m_cmems_u_rho=mean_cmems_u_rho;
                comb_gcm_interped_u_rho = ncread(gcmcmemsvelfilename, 'interped_u_rho');
                m_gcm_interped_u_rho = mean(comb_gcm_interped_u_rho,3);
                
                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecspolygon);
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, espolygon);

                ysecs_m_interped_u_rho = m_interped_u_rho.*ysecsmask;
                ysecs_m_cmems_u_rho = m_cmems_u_rho.*ysecsmask;
                ysecs_m_gcm_interped_u_rho = m_gcm_interped_u_rho.*ysecsmask;

                es_m_interped_u_rho = m_interped_u_rho.*esmask;
                es_m_gcm_interped_u_rho = m_gcm_interped_u_rho.*esmask;
                es_m_cmems_u_rho = m_cmems_u_rho.*esmask;
                
                eeeg= ysecs_m_interped_u_rho(logical(~isnan(ysecs_m_interped_u_rho).*~isnan(ysecs_m_cmems_u_rho).*~isnan(ysecs_m_gcm_interped_u_rho)));
                fffg= ysecs_m_cmems_u_rho(logical(~isnan(ysecs_m_interped_u_rho).*~isnan(ysecs_m_cmems_u_rho).*~isnan(ysecs_m_gcm_interped_u_rho)));
                rcmg_ysecscorr=corrcoef(eeeg,fffg)

                ggg= es_m_interped_u_rho(logical(~isnan(es_m_interped_u_rho).*~isnan(es_m_cmems_u_rho).*~isnan(es_m_gcm_interped_u_rho)));
                hhg= es_m_cmems_u_rho(logical(~isnan(es_m_interped_u_rho).*~isnan(es_m_cmems_u_rho).*~isnan(es_m_gcm_interped_u_rho)));
                rcmg_escorr=corrcoef(ggg,hhg)

                ii= ysecs_m_gcm_interped_u_rho(logical(~isnan(ysecs_m_gcm_interped_u_rho).*~isnan(ysecs_m_cmems_u_rho)));
                jj= ysecs_m_cmems_u_rho(logical(~isnan(ysecs_m_gcm_interped_u_rho).*~isnan(ysecs_m_cmems_u_rho)));
                gcmysecscorr=corrcoef(ii,jj)

                kk= es_m_gcm_interped_u_rho(logical(~isnan(es_m_gcm_interped_u_rho).*~isnan(es_m_cmems_u_rho)));
                ll= es_m_cmems_u_rho(logical(~isnan(es_m_gcm_interped_u_rho).*~isnan(es_m_cmems_u_rho)));
                gcmescorr=corrcoef(kk,ll)
                
                
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',rms_u_rho');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(u_rho_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('rms-u,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(rms_u_rho(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 0.7]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_u_rho
            end

            fig_flag=0;
        end

% start-------------------- v_rms
        fig_flag=fig_flags{72,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'v_','rms_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                rms_v_rho=ncread(cmemsvelfilename,'rms_v_rho');
                comb_interped_v_rho = ncread(cmemsvelfilename, 'interped_v_rho');
                mean_interped_v_rho = mean(comb_interped_v_rho, 3);
                comb_cmems_v_rho = ncread(cmemsvelfilename, 'cmems_v_rho');
                mean_cmems_v_rho = mean(comb_cmems_v_rho, 3);
                
                gcm_rms_v_rho = ncread(gcmcmemsvelfilename,'rms_v_rho');
                gcm_rms_mask = gcm_rms_v_rho;
                gcm_rms_mask(isfinite(gcm_rms_v_rho))=1;
                
                v_rho_corr = corrcoef(mean_interped_v_rho(logical(~isnan(mean_interped_v_rho).*~isnan(mean_cmems_v_rho))), mean_cmems_v_rho(logical(~isnan(mean_interped_v_rho).*~isnan(mean_cmems_v_rho))));
                v_rho_gcmcorr = corrcoef(mean_interped_v_rho(logical(~isnan(mean_interped_v_rho).*~isnan(mean_cmems_v_rho).*~isnan(gcm_rms_mask))), mean_cmems_v_rho(logical(~isnan(mean_interped_v_rho).*~isnan(mean_cmems_v_rho).*~isnan(gcm_rms_mask))))

                
                m_interped_v_rho = mean_interped_v_rho;
                m_cmems_v_rho=mean_cmems_v_rho;
                comb_gcm_interped_v_rho = ncread(gcmcmemsvelfilename, 'interped_v_rho');
                m_gcm_interped_v_rho = mean(comb_gcm_interped_v_rho,3);
                
                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecspolygon);
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, espolygon);

                ysecs_m_interped_v_rho = m_interped_v_rho.*ysecsmask;
                ysecs_m_cmems_v_rho = m_cmems_v_rho.*ysecsmask;
                ysecs_m_gcm_interped_v_rho = m_gcm_interped_v_rho.*ysecsmask;

                es_m_interped_v_rho = m_interped_v_rho.*esmask;
                es_m_gcm_interped_v_rho = m_gcm_interped_v_rho.*esmask;
                es_m_cmems_v_rho = m_cmems_v_rho.*esmask;
                
                eeeg= ysecs_m_interped_v_rho(logical(~isnan(ysecs_m_interped_v_rho).*~isnan(ysecs_m_cmems_v_rho).*~isnan(ysecs_m_gcm_interped_v_rho)));
                fffg= ysecs_m_cmems_v_rho(logical(~isnan(ysecs_m_interped_v_rho).*~isnan(ysecs_m_cmems_v_rho).*~isnan(ysecs_m_gcm_interped_v_rho)));
                rcmg_ysecscorr=corrcoef(eeeg,fffg)

                ggg= es_m_interped_v_rho(logical(~isnan(es_m_interped_v_rho).*~isnan(es_m_cmems_v_rho).*~isnan(es_m_gcm_interped_v_rho)));
                hhg= es_m_cmems_v_rho(logical(~isnan(es_m_interped_v_rho).*~isnan(es_m_cmems_v_rho).*~isnan(es_m_gcm_interped_v_rho)));
                rcmg_escorr=corrcoef(ggg,hhg)

                ii= ysecs_m_gcm_interped_v_rho(logical(~isnan(ysecs_m_gcm_interped_v_rho).*~isnan(ysecs_m_cmems_v_rho)));
                jj= ysecs_m_cmems_v_rho(logical(~isnan(ysecs_m_gcm_interped_v_rho).*~isnan(ysecs_m_cmems_v_rho)));
                gcmysecscorr=corrcoef(ii,jj)

                kk= es_m_gcm_interped_v_rho(logical(~isnan(es_m_gcm_interped_v_rho).*~isnan(es_m_cmems_v_rho)));
                ll= es_m_cmems_v_rho(logical(~isnan(es_m_gcm_interped_v_rho).*~isnan(es_m_cmems_v_rho)));
                gcmescorr=corrcoef(kk,ll)
% 
%                 for loni=1:size(model_sla,1)
%                     for lati=1:size(model_sla,2)
%                         model_sla_std(loni,lati)=std(model_sla(loni,lati,:));
%                     end
%                 end
%                 cmems_land=ones(size(model_sla));
%                 cmems_land(isnan(model_sla))=1;
%                 cmems_land(isfinite(model_sla))=NaN;
%                 save([filedir,regionname, '_', testname, '_cmems_std_land','.mat'], 'cmems_land');
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;

%                 mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%                 mask_model(mask_model==0)=NaN;
%                 rms_u_rho=rms_u_rho .* mask_model;

%                 yearstr_min=num2str(inputyear(1));
%                 yearstr_max=num2str(inputyear(end));
%                 save([filedir,regionname, '_cmems_linear_', 'det_std','_',yearstr_min, '_', yearstr_max, '.mat'], 'cut_lon_rho', 'cut_lat_rho', 'model_sla_std');
                
                m_pcolor(cut_lon_rho',cut_lat_rho',rms_v_rho');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(v_rho_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('rms-v,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(rms_v_rho(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 0.7]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_v_rho
            end

            fig_flag=0;
        end

% start-------------------- speed_std
        fig_flag=fig_flags{73,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_interped_', 'speed_','std_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                rms_u_rho=ncread(cmemsvelfilename,'rms_u_rho');
                comb_interped_u_rho = ncread(cmemsvelfilename, 'interped_u_rho');
                comb_interped_v_rho = ncread(cmemsvelfilename, 'interped_v_rho');
                comb_interped_speed = sqrt ( comb_interped_u_rho.^2 + comb_interped_v_rho.^2 );
%                 mean_interped_speed = mean(comb_interped_speed, 3);
                for loni = 1:size(cut_lon_rho,1)
                    for lati = 1:size(cut_lat_rho,2)
                        std_interped_speed(loni,lati) = std(comb_interped_speed(loni,lati,:));
                    end
                end
                pcolor(std_interped_speed'); shading flat; colorbar;
                
                
%                 comb_cmems_u_rho = ncread(cmemsvelfilename, 'cmems_u_rho');
%                 comb_cmems_v_rho = ncread(cmemsvelfilename, 'cmems_v_rho');
%                 comb_cmems_speed = sqrt ( comb_cmems_u_rho.^2 + comb_cmems_v_rho.^2 );
% %                 mean_cmems_speed = mean(comb_cmems_speed, 3);
%                 for loni = 1:size(cut_lon_rho,1)
%                     for lati = 1:size(cut_lat_rho,2)
%                         std_cmems_speed(loni,lati) = std(comb_cmems_speed(loni,lati,:));
%                     end
%                 end
%                 pcolor(std_cmems_speed'); shading flat; colorbar;
                
                

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',std_interped_speed');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%                 pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(u_rho_corr(1,2),2))]);
%                 set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('std-speed,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(std_interped_speed(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 0.2]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_u_rho
            end

            fig_flag=0;
        end

% start-------------------- cmems_ speed_std
        fig_flag=fig_flags{74,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'speed_','std_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
%                 rms_u_rho=ncread(cmemsvelfilename,'rms_u_rho');
%                 comb_interped_u_rho = ncread(cmemsvelfilename, 'interped_u_rho');
%                 comb_interped_v_rho = ncread(cmemsvelfilename, 'interped_v_rho');
%                 comb_interped_speed = sqrt ( comb_interped_u_rho.^2 + comb_interped_v_rho.^2 );
% %                 mean_interped_speed = mean(comb_interped_speed, 3);
%                 for loni = 1:size(cut_lon_rho,1)
%                     for lati = 1:size(cut_lat_rho,2)
%                         std_interped_speed(loni,lati) = std(comb_interped_speed(loni,lati,:));
%                     end
%                 end
%                 pcolor(std_interped_speed'); shading flat; colorbar;
                
                
                comb_cmems_u_rho = ncread(cmemsvelfilename, 'cmems_u_rho');
                comb_cmems_v_rho = ncread(cmemsvelfilename, 'cmems_v_rho');
                comb_cmems_speed = sqrt ( comb_cmems_u_rho.^2 + comb_cmems_v_rho.^2 );
%                 mean_cmems_speed = mean(comb_cmems_speed, 3);
                for loni = 1:size(cut_lon_rho,1)
                    for lati = 1:size(cut_lat_rho,2)
                        std_cmems_speed(loni,lati) = std(comb_cmems_speed(loni,lati,:));
                    end
                end
                std_cmems_speed = std_cmems_speed.*cmems_mask;
%                 pcolor(std_cmems_speed'); shading flat; colorbar;
                
                

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',std_cmems_speed');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%                 pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(u_rho_corr(1,2),2))]);
%                 set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('std-speed,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(std_cmems_speed(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 0.2]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_u_rho
            end

            fig_flag=0;
        end

% start-------------------- sst_rms
        fig_flag=fig_flags{75,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'sst_','rms_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_avhrr=ncread(avhrrfilename, 'lon_avhrr');
                    lat_avhrr=ncread(avhrrfilename, 'lat_avhrr');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_avhrr, lat_avhrr);
                    cut_lon_rho = lon_avhrr(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_avhrr(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                rms_sst=ncread(avhrrfilename,'rms_sst');
                comb_interped_sst = ncread(avhrrfilename, 'interped_sst');
                mean_interped_sst = mean(comb_interped_sst, 3);
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                mean_avhrr_sst = mean(comb_avhrr_sst, 3);
                
                sst_corr = corrcoef(mean_interped_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))), mean_avhrr_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))));

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',rms_sst');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(sst_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('rms-u,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(rms_sst(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 0.7]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_sst
            end

            fig_flag=0;
        end
        
    end
end