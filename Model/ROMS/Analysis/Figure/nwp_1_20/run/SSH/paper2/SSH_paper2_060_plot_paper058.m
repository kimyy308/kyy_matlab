close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR', 'ens03'};
all_testname2 = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
% all_testname2 = {'test53'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP4', 'YS', 'YSECS', 'ECS2', 'NES', 'SES', 'ES'}

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
            dropboxpath='C:\Users\USER\Dropbox';
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
        inputyear = [1982:2005]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]

        varname ='zeta'
        variable='SSH'
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
            fig_flags{67,1}='std_model_yearly_sla_detrended_plot';
            
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
%             fig_flags{32,2}=0;
%             fig_flags{33,2}=0;
%             fig_flags{34,2}=0;
%             fig_flags{35,2}=0;
            fig_flags{36,2}=1;
%             fig_flags{37,2}=1;
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
            fig_flags{53,2}=1;
            fig_flags{54,2}=1;
            fig_flags{55,2}=1;
            fig_flags{56,2}=1;
            fig_flags{57,2}=1;
            fig_flags{58,2}=0;
            fig_flags{59,2}=1;
            fig_flags{60,2}=1;
            fig_flags{61,2}=1;
            fig_flags{62,2}=1;
            fig_flags{63,2}=1;
            fig_flags{64,2}=1;
            fig_flags{65,2}=1;
            fig_flags{66,2}=1;
            fig_flags{67,2}=2;
        end
        for flagi=1:100
            fig_flags{flagi,2}=0;
        end
        fig_flags{75,2}=2;

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
%         filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
        
        scenname='historical';
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\CMIP5\',testname,'\',regionname,'\', scenname, '\'); % % where figure files will be saved
%             param_script ='C:\Users\USER\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\USER\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
            filedir = strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\'); % % where data files are
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
        cmemsvelfilename = strcat(cmip5dir,'\surface\', testname,'_',regionname, 'cmems_vel_rms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        avhrrfilename = strcat(cmip5dir,'\surface\', scenname, '_', testname, '_',regionname,'_rms_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

        valnum=0;
        run('C:\Users\USER\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);

        run(param_script);
        
        figdir=[figrawdir,'Trend\SSH\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);


% start-------------------- make timedata(xData) for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end

        if (exist('lon_min' , 'var') ~= 1)
            lon=ncread(avhrrfilename, 'lon');
            lat=ncread(avhrrfilename, 'lat');
            [lat_avhrr, lon_avhrr]= meshgrid(lat, lon);
            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_avhrr, lat_avhrr);
            cut_lon_rho = lon_avhrr(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            cut_lat_rho = lat_avhrr(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
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
                
                rms_sst=ncread(avhrrfilename,'rms');
                comb_interped_sst = ncread(avhrrfilename, 'sst');
                mean_interped_sst = mean(comb_interped_sst, 3);
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                mean_avhrr_sst = mean(comb_avhrr_sst, 3);
                
                sst_corr = corrcoef(mean_interped_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))), mean_avhrr_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))));

                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecspolygon);
                rms_ysecs= rms_sst .* ysecsmask;
                mean_rms_ysecs=mean(rms_ysecs(:),'omitnan')
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, espolygon);
                rms_es= rms_sst .* esmask;
                mean_rms_es=mean(rms_es(:),'omitnan')
                
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',rms_sst');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(sst_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('rms-sst,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(rms_sst(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([1.5 6]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_sst
            end

            fig_flag=0;
        end

% start-------------------- sst_bias
        fig_flag=fig_flags{76,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'sst_','bias_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                
                comb_bias=ncread(avhrrfilename,'bias');
                bias_sst=mean(comb_bias,3);
                comb_interped_sst = ncread(avhrrfilename, 'sst');
                mean_interped_sst = mean(comb_interped_sst, 3);
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                mean_avhrr_sst = mean(comb_avhrr_sst, 3);
                
                sst_corr = corrcoef(mean_interped_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))), mean_avhrr_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))));

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',bias_sst');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(sst_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('bias-sst,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(bias_sst(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([-4 4]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_sst
            end
            fig_flag=0;
        end

% start-------------------- summer sst_rms
        fig_flag=fig_flags{77,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'summer_sst_','rms_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
%                 rms_sst=ncread(avhrrfilename,'rms');
                month_season = [6,7,8];
                for meani = 1:length(month_season)
                    tempind{meani}= (inputyear-1982)*12 + month_season(meani);
                end
                ind_season = sort([tempind{:}]);
                comb_interped_sst = ncread(avhrrfilename, 'sst');
                season_interped_sst = comb_interped_sst(:,:,ind_season);
                mean_interped_sst = mean(season_interped_sst, 3);
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                season_avhrr_sst = comb_avhrr_sst(:,:,ind_season);
                comb_mse = ncread(avhrrfilename, 'mse');
                season_mse = comb_mse(:,:,ind_season);
                season_rms = sqrt(mean(season_mse,3));
                mean_avhrr_sst = mean(season_avhrr_sst, 3);

                ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
                season_rms_ys= season_rms .* ysmask;
                mean_rms_ys=mean(season_rms_ys(:),'omitnan');

%                 sst_corr = corrcoef(season_interped_sst(logical(~isnan(season_interped_sst).*~isnan(season_avhrr_sst))), season_avhrr_sst(logical(~isnan(season_interped_sst).*~isnan(season_avhrr_sst))));
                sst_corr = corrcoef(mean_interped_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))), mean_avhrr_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))));

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',season_rms');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(sst_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                ys_rms_text=m_text(120, 46, ['ys rms = ', num2str(round(mean_rms_ys,2))]);
                set(ys_rms_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('rms-summer,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(season_rms(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(^oC)','fontsize',colorbar_title_fontsize);
                caxis([1.5 6]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_sst
            end

            fig_flag=0;
        end
        
% start-------------------- winter sst_rms
        fig_flag=fig_flags{78,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'winter_sst_','rms_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
%                 rms_sst=ncread(avhrrfilename,'rms');
                month_season = [1,2,12];
                for meani = 1:length(month_season)
                    tempind{meani}= (inputyear-1982)*12 + month_season(meani);
                end
                ind_season = sort([tempind{:}]);
                comb_interped_sst = ncread(avhrrfilename, 'sst');
                season_interped_sst = comb_interped_sst(:,:,ind_season);
%                 mean_interped_sst = mean(comb_interped_sst, 3);
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                season_avhrr_sst = comb_avhrr_sst(:,:,ind_season);
                comb_mse = ncread(avhrrfilename, 'mse');
                season_mse = comb_mse(:,:,ind_season);
                season_rms = sqrt(mean(season_mse,3));
%                 mean_avhrr_sst = mean(comb_avhrr_sst, 3);
                
                sst_corr = corrcoef(season_interped_sst(logical(~isnan(season_interped_sst).*~isnan(season_avhrr_sst))), season_avhrr_sst(logical(~isnan(season_interped_sst).*~isnan(season_avhrr_sst))));

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',season_rms');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(sst_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('rms-winter,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(season_rms(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([1.5 6]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_sst
            end

            fig_flag=0;
        end

% start-------------------- summer sst_bias
        fig_flag=fig_flags{79,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'summer_sst_','bias_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
%                 rms_sst=ncread(avhrrfilename,'rms');
                month_season = [6,7,8];
                for meani = 1:length(month_season)
                    tempind{meani}= (inputyear-1982)*12 + month_season(meani);
                end
                ind_season = sort([tempind{:}]);
                comb_interped_sst = ncread(avhrrfilename, 'sst');
                season_interped_sst = comb_interped_sst(:,:,ind_season);
                mean_interped_sst = mean(season_interped_sst, 3);
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                season_avhrr_sst = comb_avhrr_sst(:,:,ind_season);
%                 comb_mse = ncread(avhrrfilename, 'mse');
%                 season_mse = comb_mse(:,:,ind_season);
                comb_season_bias=season_interped_sst - season_avhrr_sst;
%                 season_rms = sqrt(mean(season_mse,3));
                season_bias = mean(comb_season_bias,3);
                mean_avhrr_sst = mean(season_avhrr_sst, 3);

                ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
                season_bias_ys= season_bias .* ysmask;
                mean_bias_ys=mean(season_bias_ys(:),'omitnan');

%                 sst_corr = corrcoef(season_interped_sst(logical(~isnan(season_interped_sst).*~isnan(season_avhrr_sst))), season_avhrr_sst(logical(~isnan(season_interped_sst).*~isnan(season_avhrr_sst))));
                sst_corr = corrcoef(mean_interped_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))), mean_avhrr_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))));

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',season_bias');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(sst_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                ys_bias_text=m_text(120, 46, ['ys bias = ', num2str(round(mean_bias_ys,2))]);
                set(ys_bias_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('bias-summer,','gcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(season_bias(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(^oC)','fontsize',colorbar_title_fontsize);
                caxis([-4 4]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho rms_sst
            end

            fig_flag=0;
        end

% start-------------------- winter sst_bias
        fig_flag=fig_flags{80,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'winter_sst_','bias_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
%                 rms_sst=ncread(avhrrfilename,'rms');
                month_season = [1,2,12];
                for meani = 1:length(month_season)
                    tempind{meani}= (inputyear-1982)*12 + month_season(meani);
                end
                ind_season = sort([tempind{:}]);
                comb_interped_sst = ncread(avhrrfilename, 'sst');
                season_interped_sst = comb_interped_sst(:,:,ind_season);
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                season_avhrr_sst = comb_avhrr_sst(:,:,ind_season);
%                 comb_mse = ncread(avhrrfilename, 'mse');
%                 season_mse = comb_mse(:,:,ind_season);
                comb_season_bias=season_interped_sst - season_avhrr_sst;
%                 season_rms = sqrt(mean(season_mse,3));
                season_bias = mean(comb_season_bias,3);
                mean_interped_sst = mean(season_interped_sst, 3);
                mean_avhrr_sst = mean(season_avhrr_sst, 3);

                ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
                season_bias_ys= season_bias .* ysmask;
                mean_bias_ys=mean(season_bias_ys(:),'omitnan');

%                 sst_corr = corrcoef(season_interped_sst(logical(~isnan(season_interped_sst).*~isnan(season_avhrr_sst))), season_avhrr_sst(logical(~isnan(season_interped_sst).*~isnan(season_avhrr_sst))));
                sst_corr = corrcoef(mean_interped_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))), mean_avhrr_sst(logical(~isnan(mean_interped_sst).*~isnan(mean_avhrr_sst))));

                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
        
                m_pcolor(cut_lon_rho',cut_lat_rho',season_bias');

                shading(gca,m_pcolor_shading_method);   

                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                pt_corr_text=m_text(120, 48, ['ptrn corr = ', num2str(round(sst_corr(1,2),2))]);
                set(pt_corr_text, 'fontsize', colorbar_title_fontsize);
                ys_bias_text=m_text(120, 46, ['ys bias = ', num2str(round(mean_bias_ys,2))]);
                set(ys_bias_text, 'fontsize', colorbar_title_fontsize);
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat('bias-winter,','gcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(season_bias(:),'omitnan'),2)), ') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(bwrmap);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(^oC)','fontsize',colorbar_title_fontsize);
                caxis([-4 4]);

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