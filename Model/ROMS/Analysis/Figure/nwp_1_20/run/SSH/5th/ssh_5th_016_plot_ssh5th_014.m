close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'test57', 'test58', 'test59', 'test60'};

% all_testname2 = {'test65', 'test66', 'test67', 'test68'};
% all_testname2 = {'test57', 'test58', 'test59', 'test60', 'test65', 'test66', 'test67', 'test68'};

% all_testname2 = {'test61', 'test62', 'test63', 'test64', 'ens09'};
% all_testname2 = {'test61', 'test62', 'test63', 'test64', 'test65', 'test66', 'test67', 'test68', 'ens09', 'ens10'};
% all_testname2 = {'test57', 'test58', 'test59', 'test60', 'ens08'};

all_testname2 = {'test57', 'test58', 'test59', 'test60', 'ens08', 'test61', 'test62', 'test63', 'test64', 'test65', 'test66', 'test67', 'test68', 'ens09', 'ens10'};

% all_testname2 = {'ens09', 'ens08', 'ens10'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP4', 'YS', 'YSECS', 'ECS2', 'NES', 'SES', 'ES'}

% all_region2 ={'NWP', 'AKP4'};
% all_region2 ={'ES_KHOA','YS_KHOA', 'SS_KHOA'};
% all_region2 ={'SS_KHOA'};

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
            dropboxpath='C:\Users\User\Dropbox';
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
        abstrendlev =[3 9];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.15 0.15];
        trendplotlev = [0 7];
        trenddifflev = [-10 10];
        sshlev =[-0.3 0.3];
        sshdifflev = [0 20];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]

        varname ='zeta';
        variable='SSH';
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        
        
table=1;
fccur_list = {'103', '253', '302'};
for i=1:3
    eval(['sort_', fccur_list{i}, ' = table']);
end
        
% % %         switch region
        for folding=1:1
            switch(regionname)
                case('NWP') %% North western Pacific
                    lonlat = [115, 164, 15, 52];  %% whole data area
                    refpolygon(1,1)=lonlat(1);
                    refpolygon(2,1)=lonlat(2);
                    refpolygon(1,2)=lonlat(3);
                    refpolygon(2,2)=lonlat(4);
                case('NWP2') %% North western Pacific
                    lonlat = [115, 145, 25, 52];  %% whole data area
                    refpolygon(1,1)=lonlat(1);
                    refpolygon(2,1)=lonlat(2);
                    refpolygon(1,2)=lonlat(3);
                    refpolygon(2,2)=lonlat(4);
                case('ES') %% East Sea
                    refpolygon=espolygon;
                case('NES') %% Northern East Sea
                    refpolygon=nespolygon;
                case('SES') %% Southern East Sea
                    refpolygon=sespolygon;
                case('SS') %% South Sea
                    refpolygon=sspolygon;
                case('YS') %% Yellow Sea
                    refpolygon=yspolygon;
                case('ECS') %% East China Sea
                    refpolygon=ecspolygon;
                 case('ECS2') %% East China Sea
                    refpolygon=ecs2polygon;
                case('YSECS') %% East China Sea
                    refpolygon=ysecspolygon;
                case('AKP') %% Around Korea Peninsula
                    refpolygon=akppolygon;
                case('AKP2') %% Around Korea Peninsula
                    refpolygon=akp2polygon;
                case('AKP3') %% Around Korea Peninsula
                    refpolygon=akp3polygon;
                case('AKP4') %% Around Korea Peninsula
                    refpolygon=akp4polygon;
                case('CA') %% Around Korea Peninsula
                    refpolygon=capolygon;
                case('EKB') %% Around Korea Peninsula
                    refpolygon=akp2polygon;
                case('BOH') %% Around Korea Peninsula
                    refpolygon=bohpolygon;
                case('ES_KHOA') %% East Sea
                    refpolygon=es_khoapolygon;
                case('YS_KHOA') %% Yellow sea
                    refpolygon=ys_khoapolygon;
                case('SS_KHOA') %% South sea
                    refpolygon=ss_khoapolygon;
                case('TEST') %% for debugging
                    refpolygon=testpolygon;
                otherwise
                    ('?')
            end
            lonlat(1)=min(refpolygon(:,1));
            lonlat(2)=max(refpolygon(:,1));
            lonlat(3)=min(refpolygon(:,2));
            lonlat(4)=max(refpolygon(:,2));
        end
        
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
            fig_flags{27,1}='time_series_sla_cmems_seasonal';
            fig_flags{28,1}='amplitude_seasonal_ssh_cmems';
            fig_flags{29,1}='amplitude_seasonal_ssh_model';
            fig_flags{30,1}='amplitude_difference_seasonal_ssh_model';
            fig_flags{31,1}='time_series_ssh_yearly';
            fig_flags{32,1}='time_series_cmems_plot';
            fig_flags{33,1}='RMS_ssh_plot';
            fig_flags{34,1}='RMS_corrected_ssh_plot';
            fig_flags{35,1}='BIAS_ssh_plot';
            fig_flags{36,1}='variance_cmems_sla_plot';
            fig_flags{37,1}='std_cmems_sla_plot';
            fig_flags{38,1}='variance_cmems_sla_seasonal_plot';
            fig_flags{39,1}='variance_model_sla_plot';
            fig_flags{40,1}='variance_model_sla_seasonal_plot';
            fig_flags{41,1}='variance_cmems_sla_lowpass_plot';
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
            fig_flags{55,1}='slr-ratio_plot';
            fig_flags{56,1}='slr-ratio_without_tide_plot';
            fig_flags{57,1}='SLR_plot from yearly movmean data(linear fit)';
            fig_flags{58,1}='SLR plot from yearly exponential fit data';
            fig_flags{59,1}='SLR plot from yearly poly1 fit data';
            fig_flags{60,1}='SLR r plot from yearly exp fit';
            fig_flags{61,1}='SLR r plot from yearly poly1 fit';
            fig_flags{62,1}='SLR rmse plot from yearly exp fit';
            fig_flags{63,1}='SLR rmse plot from yearly poly1 fit';
            fig_flags{64,1}='slr-ratio_plot with yearly exp fit SLR';
            fig_flags{65,1}='slr-ratio_plot with yearly exp fit SLR + cmems yearly detrended std';
            fig_flags{66,1}='slr-ratio_plot with yearly exp fit SLR, all test';
            fig_flags{67,1}='tgtyear SLR plot from yearly exponential fit data';

            for flagi=1:67
                fig_flags{flagi,2}=0;
            end
%             fig_flags{1,2}=0;
%             fig_flags{2,2}=0;
%             fig_flags{3,2}=0;
%             fig_flags{4,2}=0;
%             fig_flags{5,2}=1;
%             fig_flags{6,2}=0;
%             fig_flags{7,2}=0;
%             fig_flags{8,2}=0;
%             fig_flags{9,2}=0;
%             fig_flags{10,2}=0;
%             fig_flags{11,2}=0;
            fig_flags{12,2}=2;
            fig_flags{13,2}=1;
            fig_flags{14,2}=1;
%             fig_flags{15,2}=1;
%             fig_flags{16,2}=0;
%             fig_flags{17,2}=0;
            fig_flags{18,2}=1;
            fig_flags{19,2}=1;
%             fig_flags{20,2}=0;
%             fig_flags{21,2}=0;
%             fig_flags{22,2}=2;
%             fig_flags{23,2}=0;
%             fig_flags{24,2}=0;
%             fig_flags{25,2}=1;0
            fig_flags{26,2}=1;
%             fig_flags{27,2}=0;
%             fig_flags{28,2}=0;
            fig_flags{29,2}=1;
%             fig_flags{30,2}=0;
            fig_flags{31,2}=1;
%             fig_flags{32,2}=1;
%             fig_flags{33,2}=0;
%             fig_flags{34,2}=0;
%             fig_flags{35,2}=0;
            fig_flags{36,2}=1;
%             fig_flags{37,2}=1;
%             fig_flags{38,2}=0;
            fig_flags{39,2}=1;
%             fig_flags{40,2}=1;
%             fig_flags{41,2}=0;
%             fig_flags{42,2}=0;
%             fig_flags{43,2}=0;
%             fig_flags{44,2}=1;
%             fig_flags{45,2}=0;
%             fig_flags{46,2}=1;
%             fig_flags{47,2}=1;
%             fig_flags{48,2}=1;
%             fig_flags{49,2}=0;
%             fig_flags{50,2}=0;
%             fig_flags{51,2}=0;
%             fig_flags{52,2}=0;
%             fig_flags{53,2}=1;
%             fig_flags{54,2}=1;
            fig_flags{55,2}=0;
            fig_flags{56,2}=0;
            fig_flags{57,2}=0;
            fig_flags{58,2}=1;
            fig_flags{59,2}=1;
            fig_flags{60,2}=1;
            fig_flags{61,2}=1;
            fig_flags{62,2}=0;
            fig_flags{63,2}=0;
            fig_flags{64,2}=0;
            fig_flags{65,2}=0;
            fig_flags{66,2}=0;
            fig_flags{67,2}=0;
        end

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
%         filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            drivename='D:\';
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
            filedir = strcat(drivename, 'Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
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



        valnum=0;
        run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);
        parula512=parula(512);
        parula_bg=parula512(1:256,:);
        parula_gy=parula512(257:512,:);
        
        cmems_trend=ncread(interpedfilename, 'interped_trend');
        cmems_trend_filtered=ncread(interpedfilename, 'interped_trend_filtered');

        cmems_mask=ones(size(cmems_trend));
        cmems_mask(isnan(cmems_trend))=NaN;

        run(param_script);
        
        figdir=[figrawdir,'Trend\SSH\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        
%         lon_rho=ncread(modelfilename,'lon_rho');
%         lat_rho=ncread(modelfilename,'lat_rho');
        lon_interped=ncread(interpedfilename, 'lon_cmems');
        lat_interped=ncread(interpedfilename, 'lat_cmems');
%         trend_filtered = ncread(modelfilename,'trend_filtered');
%         mean_trend_filtered = ncread(filename,'mean_trend_filtered');
%         cmems_trend_filtered = ncread(cmemsfilename,'cmems_trend_filtered');
        interped_trend_filtered = ncread(interpedfilename,'interped_trend_filtered');
                
        
        for folding=1:1
        if (exist('lon_min' , 'var') ~= 1)
            lon_cmems=ncread(interpedfilename, 'lon_cmems');
            lat_cmems=ncread(interpedfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
        end
        end
        
        
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
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(lon_interped)',lat_interped',squeeze(interped_trend_filtered(:,:)'));
                mean_trend_filtered=mean(interped_trend_filtered(:),'omitnan');
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                
                switch testname
                    case {'test61', 'test62', 'test63', 'test64','ens09'}
                        colormap(winter);
                        caxis([3 6]);
                    case {'test57', 'test58', 'test59', 'test60','ens08'}
                        colormap(parula_gy);
                        caxis([3 7]);
                    case {'test65', 'test66', 'test67', 'test68','ens10'}
                        colormap(flip(autumn));
                        caxis([6 9]);
                end
                
                % set colorbar 
                h = colorbar;
%                 colormap(wrmap);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%                 caxis(abstrendlev);

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
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_bwrmap_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(lon_interped)',lat_interped',squeeze(interped_trend_filtered(:,:)'));
                mean_trend_filtered=mean(interped_trend_filtered(:),'omitnan');
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

% start-------------------- relative trend plot (bwrmap)
        fig_flag=fig_flags{14,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_relative_ssh_trend_bwrmap_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                mean_trend_filtered=mean(interped_trend_filtered(:),'omitnan');
                m_pcolor(double(lon_interped)',lat_interped',squeeze(interped_trend_filtered(:,:)'-mean_trend_filtered));
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
                    msl_filt(varind)=squeeze(mean(mean(ncread(interpedfilename,'interped_sla_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
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
                    msl(varind)=squeeze(mean(mean(ncread(interpedfilename,'interped_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
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
            clim_trend=ncread(interpedfilename,'clim_interped_ssh_trend');
            clim_trend_reshaped=reshape(clim_trend, [size(clim_trend,1)*size(clim_trend,2), 12]);
            mean_clim_trend=mean(clim_trend_reshaped,1,'omitnan');
            if abs(mean_clim_trend(1))<0.01
                mean_clim_trend=mean_clim_trend*1000.0;  %% m -> mm
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
                lon_cmems=ncread(interpedfilename, 'lon_cmems');
                lat_cmems=ncread(interpedfilename, 'lat_cmems');
                len_lon=size(lon_cmems,1);
                len_lat=size(lat_cmems,2);
                comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
                clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
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

% start-------------------- seasonal msl amplitude
        fig_flag=fig_flags{29,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_amp_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

    %         clim_ssh=ncread(filename,'clim_ssh');
    %         for ijij=1:12
    %             temp_clim_ssh=clim_ssh(:,:,ijij);
    %             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %         end
    %         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)

                lon_cmems=ncread(interpedfilename, 'lon_cmems');
                lat_cmems=ncread(interpedfilename, 'lat_cmems');
                len_lon=size(lon_cmems,1);
                len_lat=size(lat_cmems,2);
                comb_interped_data=ncread(interpedfilename,'interped_ssh');
                comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
                clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
                clim_ssh(clim_ssh>1000000)=NaN;
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

% start-------------------- model variance plot
        fig_flag=fig_flags{39,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_interped_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
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

% start-------------------- model variance_climatological plot
        fig_flag=fig_flags{40,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_interped_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'),'_', num2str(tempmonth,'%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                    run(param_script);
                    if (exist('lon_min' , 'var') ~= 1)
                        lon_cmems=ncread(interpedfilename, 'lon_cmems');
                        lat_cmems=ncread(interpedfilename, 'lat_cmems');
%                         [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                        [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                        cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    end

    %                 model_sla=ncread(filename,'interped_sla', [lon_min(1) lat_min(1) 1], ...
    %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
    %                         
                    for yearij=1:length(inputyear)
                        model_sla(:,:,yearij)=ncread(interpedfilename,'interped_sla', [lon_min(1) lat_min(1) (yearij-1)*12+tempmonth], ...
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
                    caxis([0 200]);

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
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            nyears=[1:5];
            nc_varname_prefixes={'cmems_'};
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
            cmems_adt=ncread(filename,'cmems_adt');
            interped_ssh=ncread(filename,'interped_ssh');
            cmems_sla=ncread(filename,'cmems_sla');
            interped_sla=ncread(filename,'interped_sla');
            m_cmems_adt=mean(cmems_adt,3,'omitnan');
            m_interped_ssh=mean(interped_ssh,3,'omitnan');
            m_cmems_sla=mean(cmems_sla,3,'omitnan');
            m_interped_sla=mean(interped_sla,3,'omitnan');
            m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
            m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));

            corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))))
            corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))))
            corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
                cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)

        % end-------------------- model pattern correlation


        % start-------------------- model pattern correlation (yearly)
            pngname=strcat(outfile, '_', testname,'_',regionname, '_', nc_varname,'_','yearly_pattern_corr_',num2str(min(inputyear),'%04i'), ...
                                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);

                    yearly_cmems_adt=ncread(filename,'cmems_adt',[1 1 (yearij-1)*12+1], [inf inf 12]);
                    yearly_interped_ssh=ncread(filename,'interped_ssh',[1 1 (yearij-1)*12+1], [inf inf 12]);
                    yearly_cmems_sla=ncread(filename,'cmems_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);
                    yearly_interped_sla=ncread(filename,'interped_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);

                    m_cmems_adt=mean(yearly_cmems_adt,3,'omitnan');
                    m_interped_ssh=mean(yearly_interped_ssh,3,'omitnan');
                    m_cmems_sla=mean(yearly_cmems_sla,3,'omitnan');
                    m_interped_sla=mean(yearly_interped_sla,3,'omitnan');
        %             m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
        %             m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));

                    ssh_corr(yearij) = corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))));
                    sla_corr(yearij) = corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))));
        %             corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
        %                 cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)
                end
                plot(ssh_corr);
                hold on
                plot(sla_corr);
                hold off
                close all
            end
            fig_flag=0;
        end
        
       
% start-------------------- movmean msl time series
        fig_flag=fig_flags{46,2};
        while (fig_flag)
%             nyears=[2];
                        nyears=[2,3,5];
            nc_varname_prefixes={'interped_'};
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
                                temp_raw_msl=ncread(interpedfilename,'interped_ssh',[1 1 varind], [inf inf 1]);
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
                            lon_cmems=ncread(interpedfilename, 'lon_cmems');
                            lat_cmems=ncread(interpedfilename, 'lat_cmems');
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
        
% start-------------------- model slr-ratio index plot (from linear trend, monthly movmean data)
        fig_flag=fig_flags{55,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_slr-ratio_ind_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)

                switch(testname)
                    case {'test57', 'test61', 'test65'}
                        testname_his='test53';
                    case {'test58', 'test62', 'test66'}
                        testname_his='test54';
                    case {'test59', 'test63', 'test67'}
                        testname_his='test55';
                    case {'test60', 'test64', 'test68'}
                        testname_his='test56';
                end
                his_period = [1993 2005];
                hisdir = strcat(drivename, 'Data\Model\ROMS\nwp_1_20\', testname_his, '\run\'); % % where data files are
                movfilename_his = strcat(hisdir, testname_his,'_',regionname, ...
                    'moving_averaged_cmems_interped_ssh_trend_',num2str(min(his_period),'%04i'),'_',num2str(max(his_period),'%04i'), '.nc');

                nyears=[1];
                nc_varname_prefixes={'interped_sla_'};
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for nc_varnameij=1:length(nc_varname_prefixes)
                        nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                        nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];  
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(interpedfilename, 'lon_cmems');
                            lat_cmems=ncread(interpedfilename, 'lat_cmems');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['model_sla_his', '=ncread(movfilename_his,', '''', nc_varname,'''',');']);

                        model_sla_his=model_sla_his(:, :, nyear*(12*nyear/2)+1 : end-(12*nyear/2)-1);
                        for loni=1:size(model_sla_his,1)
                            for lati=1:size(model_sla_his,2)
                                model_sla_his_std(loni,lati)=std(model_sla_his(loni,lati,:).*100.0);
                            end
                        end

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        model_sla_his_std=model_sla_his_std .* mask_model; 
    %                     trend_1y_movmean=ncread(movfilename, 'trend_1y_movmean');
                        eval(['trend_',num2str(nyear),'y_movmean=ncread(movfilename,', '''', 'trend_', num2str(nyear), 'y_movmean', '''', ');']);
                    end
                end
                slr=trend_1y_movmean/10.0 * (max(inputyear)-min(inputyear));

                tide_info.name{1}='M2  ';
                tide_info.name{2}='S2  ';
                tide_info.name{3}='K1  ';
                tide_info.name{4}='O1  ';
                earlytidename=strcat(filedir, '\', num2str(min(inputyear),'%04i'), '\', testname,'_harmonic_analysis_zeta_',regionname, ...
                    '_', num2str(min(inputyear),'%04i'), '.nc');
    %             ncinfo(earlytidename)
                lon_rho = ncread(earlytidename, 'lon_rho');
                lat_rho = ncread(earlytidename, 'lat_rho');
                tname=ncread(earlytidename, 'tname')';
                num_tide_all=size(tname,1);
                num_tide_tgt=length(tide_info.name);
                for coni=1:num_tide_all
                    for tide_namei=1:num_tide_tgt
                        if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                            tide_info.index(tide_namei)=coni;
                        end
                    end
                end

                tcon=ncread(earlytidename, 'tcon');
                m2_amp=tcon(:,:,tide_info.index(1),1).*100;
                s2_amp=tcon(:,:,tide_info.index(2),1).*100;
                k1_amp=tcon(:,:,tide_info.index(3),1).*100;
                o1_amp=tcon(:,:,tide_info.index(4),1).*100;
                early_sum4=m2_amp+s2_amp+k1_amp+o1_amp;

                latetidename=strcat(filedir, '\', num2str(max(inputyear),'%04i'), '\', testname,'_harmonic_analysis_zeta_',regionname, ...
                    '_', num2str(max(inputyear),'%04i'), '.nc');
                tname=ncread(latetidename, 'tname')';
                num_tide_all=size(tname,1);
                num_tide_tgt=length(tide_info.name);
                for coni=1:num_tide_all
                    for tide_namei=1:num_tide_tgt
                        if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                            tide_info.index(tide_namei)=coni;
                        end
                    end
                end

                tcon=ncread(latetidename, 'tcon');
                m2_amp=tcon(:,:,tide_info.index(1),1).*100;
                s2_amp=tcon(:,:,tide_info.index(2),1).*100;
                k1_amp=tcon(:,:,tide_info.index(3),1).*100;
                o1_amp=tcon(:,:,tide_info.index(4),1).*100;
                late_sum4=m2_amp+s2_amp+k1_amp+o1_amp;
                diff_sum4=late_sum4 - early_sum4;

                diff_sum4_interped = griddata(double(lon_rho), double(lat_rho), diff_sum4,double(cut_lon_rho),double(cut_lat_rho));   
                diff_sum4_interped=diff_sum4_interped .* mask_model; 

                slr_ind=(slr+diff_sum4_interped) ./model_sla_his_std;

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('slr-ratio, ',testname, ',(', ...
                    num2str(max(inputyear),'%04i'),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                caxis([5 60]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{55,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            
            fig_flag=0;
        end

% start-------------------- model slr-ratio index plot without tide (from linear trend, monthly movmean data)
        fig_flag=fig_flags{56,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_slr-ratio_without_tide_ind_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)

                switch(testname)
                    case {'test57', 'test61', 'test65'}
                        testname_his='test53';
                    case {'test58', 'test62', 'test66'}
                        testname_his='test54';
                    case {'test59', 'test63', 'test67'}
                        testname_his='test55';
                    case {'test60', 'test64', 'test68'}
                        testname_his='test56';
                end
                his_period = [1993 2005];
                hisdir = strcat(drivename, 'Data\Model\ROMS\nwp_1_20\', testname_his, '\run\'); % % where data files are
                movfilename_his = strcat(hisdir, testname_his,'_',regionname, ...
                    'moving_averaged_cmems_interped_ssh_trend_',num2str(min(his_period),'%04i'),'_',num2str(max(his_period),'%04i'), '.nc');

                nyears=[1];
                nc_varname_prefixes={'interped_sla_'};
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for nc_varnameij=1:length(nc_varname_prefixes)
                        nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                        nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];  
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(interpedfilename, 'lon_cmems');
                            lat_cmems=ncread(interpedfilename, 'lat_cmems');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['model_sla_his', '=ncread(movfilename_his,', '''', nc_varname,'''',');']);

                        model_sla_his=model_sla_his(:, :, nyear*(12*nyear/2)+1 : end-(12*nyear/2)-1);
                        for loni=1:size(model_sla_his,1)
                            for lati=1:size(model_sla_his,2)
                                model_sla_his_std(loni,lati)=std(model_sla_his(loni,lati,:).*100.0);
                            end
                        end

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                        model_sla_his_std=model_sla_his_std .* mask_model; 
    %                     trend_1y_movmean=ncread(movfilename, 'trend_1y_movmean');
                        eval(['trend_',num2str(nyear),'y_movmean=ncread(movfilename,', '''', 'trend_', num2str(nyear), 'y_movmean', '''', ');']);
                    end
                end
                slr=trend_1y_movmean/10.0 * (max(inputyear)-min(inputyear));

                slr_ind=(slr) ./model_sla_his_std;

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('slr-ratio-notide, ',testname, ',(', ...
                    num2str(max(inputyear),'%04i'),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                caxis([5 60]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{56,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end
        
% start-------------------- SLR plot (from linear trend, monthly movmean data)
        fig_flag=fig_flags{57,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_SLR_', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                
                nyears=[1];
                nc_varname_prefixes={'interped_sla_'};
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for nc_varnameij=1:length(nc_varname_prefixes)
                        nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                        nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];  
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(interpedfilename, 'lon_cmems');
                            lat_cmems=ncread(interpedfilename, 'lat_cmems');
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end

                        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
    %                     trend_1y_movmean=ncread(movfilename, 'trend_1y_movmean');
                        eval(['trend_',num2str(nyear),'y_movmean=ncread(movfilename,', '''', 'trend_', num2str(nyear), 'y_movmean', '''', ');']);
                    end
                end
                trend_1y_movmean(trend_1y_movmean==0)=NaN;
                slr=trend_1y_movmean/10.0 * (max(inputyear)-min(inputyear));
                slr=slr.*mask_model;

                slr_ind=(slr);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('Sea Level Rise, ',testname, ',(', ...
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([40 100]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{57,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end

% start-------------------- SLR plot from yearly exponential fit data
        fig_flag=fig_flags{58,2};
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
                
%                 for i=1:size(sl,3)
%                     sl=sl.*mask_model;
%                 end
%                 sly=squeeze(mean(mean(sl,1,'omitnan'),2,'omitnan'))
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
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                switch inputyear(end)
                    case 2030
                        cbar_ratio=0.18;
                    case 2050
                        cbar_ratio=0.37;
                    case 2070
                        cbar_ratio=0.59;
                    otherwise
                        cbar_ratio=1;
                end

%                  2030 ens09=0.1980 ens08=0.1958 ens10=0.1823
%                  2050 ens09=0.3896 ens08=0.3864 ens10=0.3668
%                  2070 ens09=0.6091 ens08=0.6060 ens10=0.5871

                switch testname
                    case {'test61', 'test62', 'test63', 'test64','ens09'}
%                         colormap(parula_bg);
                        colormap(winter);
                        switch inputyear(end)
                            case 2030
                                caxis([7 17]);
                            case 2050
                                caxis([14 27]);
                            case 2070
                                caxis([22 33]);
                            case 2100
                                caxis([32 44]);
                        end
%                         caxis([30*cbar_ratio *1.15 30*cbar_ratio *1.15+20]);
                    case {'test57', 'test58', 'test59', 'test60','ens08'}
                        colormap(parula_gy);
                        switch inputyear(end)
                            case 2030
                                caxis([-3 18]);
                            case 2050
                                caxis([10 32]);
                            case 2070
                                caxis([18 50]);
                            case 2100
                                caxis([40 64]);
                        end
%                         caxis([30*cbar_ratio 30*cbar_ratio+40]);
                    case {'test65', 'test66', 'test67', 'test68','ens10'}
                        colormap(flip(autumn));
                        switch inputyear(end)
                            case 2030
                                caxis([7 19]);
                            case 2050
                                caxis([13 38]);
                            case 2070
                                caxis([30 53]);
                            case 2100
                                caxis([60 85]);
                        end
%                         caxis([60.*cbar_ratio *0.9 60.*cbar_ratio *0.9+40]);
                end

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{58,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end
        
% start-------------------- SLR plot from yearly poly1 fit data
        fig_flag=fig_flags{59,2};
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
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                switch testname
                    case {'test61', 'test62', 'test63', 'test64','ens09'}
                    case {'test57', 'test58', 'test59', 'test60','ens08'}
                        colormap(parula);
                        caxis([30 70]);
                    case {'test65', 'test66', 'test67', 'test68','ens10'}
                        colormap(flip(autumn));
                        caxis([60 100]);
                end

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
        
% start-------------------- SLR r plot from yearly exp fit
        fig_flag=fig_flags{60,2};
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
                caxis([0.7 1]);

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
   
% start-------------------- SLR r plot from yearly poly1 fit
        fig_flag=fig_flags{61,2};
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
                caxis([0.7 1]);

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
        
% start-------------------- SLR rmse plot from yearly exp fit
        fig_flag=fig_flags{62,2};
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
                disp([fig_flags{62,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end       
        
% start-------------------- SLR rmse plot from yearly poly1 fit
        fig_flag=fig_flags{63,2};
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
                disp([fig_flags{63,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            fig_flag=0;
        end               
 
% start-------------------- model slr-ratio index plot (from exp fit yearly SLR)
        fig_flag=fig_flags{64,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_slr-ratio_ind_exp_SLR', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)

                switch(testname)
                    case {'test57', 'test61', 'test65'}
                        testname_his='test53';
                    case {'test58', 'test62', 'test66'}
                        testname_his='test54';
                    case {'test59', 'test63', 'test67'}
                        testname_his='test55';
                    case {'test60', 'test64', 'test68'}
                        testname_his='test56';
                end
                his_period = [1993 2005];
                hisdir = strcat(drivename, 'Data\Model\ROMS\nwp_1_20\', testname_his, '\run\'); % % where data files are
                interpedfilename_his = strcat(hisdir, testname_his,'_',regionname, ...
                    'cmems_interped_ssh_trend_',num2str(min(his_period),'%04i'),'_',num2str(max(his_period),'%04i'), '.nc');
                
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                model_sla_his_yearly_detrended=ncread(interpedfilename_his,'interped_sla_yearly_detrended');
                for loni=1:size(model_sla_his_yearly_detrended,1)
                    for lati=1:size(model_sla_his_yearly_detrended,2)
                        model_sla_his_std(loni,lati)=std(model_sla_his_yearly_detrended(loni,lati,:));
                    end
                end

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                model_sla_his_std=model_sla_his_std .* mask_model;    
                
                sl=ncread(interpedfilename,'interped_sla_yearly_exp_fit');
                slr=sl(:,:,end)-sl(:,:,1);
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                slr=slr.*mask_model;

                tide_info.name{1}='M2  ';
                tide_info.name{2}='S2  ';
                tide_info.name{3}='K1  ';
                tide_info.name{4}='O1  ';
                earlytidename=strcat(filedir, '\', num2str(min(inputyear),'%04i'), '\', testname,'_harmonic_analysis_zeta_',regionname, ...
                    '_', num2str(min(inputyear),'%04i'), '.nc');
    %             ncinfo(earlytidename)
                lon_rho = ncread(earlytidename, 'lon_rho');
                lat_rho = ncread(earlytidename, 'lat_rho');
                tname=ncread(earlytidename, 'tname')';
                num_tide_all=size(tname,1);
                num_tide_tgt=length(tide_info.name);
                for coni=1:num_tide_all
                    for tide_namei=1:num_tide_tgt
                        if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                            tide_info.index(tide_namei)=coni;
                        end
                    end
                end

                tcon=ncread(earlytidename, 'tcon');
                m2_amp=tcon(:,:,tide_info.index(1),1).*100;
                s2_amp=tcon(:,:,tide_info.index(2),1).*100;
                k1_amp=tcon(:,:,tide_info.index(3),1).*100;
                o1_amp=tcon(:,:,tide_info.index(4),1).*100;
                early_sum4=m2_amp+s2_amp+k1_amp+o1_amp;

                latetidename=strcat(filedir, '\', num2str(max(inputyear),'%04i'), '\', testname,'_harmonic_analysis_zeta_',regionname, ...
                    '_', num2str(max(inputyear),'%04i'), '.nc');
                tname=ncread(latetidename, 'tname')';
                num_tide_all=size(tname,1);
                num_tide_tgt=length(tide_info.name);
                for coni=1:num_tide_all
                    for tide_namei=1:num_tide_tgt
                        if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                            tide_info.index(tide_namei)=coni;
                        end
                    end
                end

                tcon=ncread(latetidename, 'tcon');
                m2_amp=tcon(:,:,tide_info.index(1),1).*100;
                s2_amp=tcon(:,:,tide_info.index(2),1).*100;
                k1_amp=tcon(:,:,tide_info.index(3),1).*100;
                o1_amp=tcon(:,:,tide_info.index(4),1).*100;
                late_sum4=m2_amp+s2_amp+k1_amp+o1_amp;
                diff_sum4=late_sum4 - early_sum4;

                diff_sum4_interped = griddata(double(lon_rho), double(lat_rho), diff_sum4,double(cut_lon_rho),double(cut_lat_rho));   
                diff_sum4_interped=diff_sum4_interped .* mask_model; 

                slr_ind=(slr+diff_sum4_interped) ./model_sla_his_std;

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('slr-ratio, ',testname, ',(', ...
                    num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                caxis([0 70]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');
                save([filedir,testname,'_slr_ind_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'], 'slr_ind');

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

% start-------------------- model slr-ratio index plot (from exp fit yearly SLR + cmems std)
        fig_flag=fig_flags{65,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_slr-ratio_cmems_ind_exp_SLR', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)

                switch(testname)
                    case {'test57', 'test61', 'test65'}
                        testname_his='test53';
                    case {'test58', 'test62', 'test66'}
                        testname_his='test54';
                    case {'test59', 'test63', 'test67'}
                        testname_his='test55';
                    case {'test60', 'test64', 'test68'}
                        testname_his='test56';
                end
                his_period = [1993 2005];
                hisdir = strcat(drivename, 'Data\Model\ROMS\nwp_1_20\', testname_his, '\run\'); % % where data files are
                interpedfilename_his = strcat(hisdir, testname_his,'_',regionname, ...
                    'cmems_ssh_trend_',num2str(min(his_period),'%04i'),'_',num2str(max(his_period),'%04i'), '.nc');
                
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                model_sla_his_yearly_detrended=ncread(interpedfilename_his,'cmems_sla_yearly_detrended');
                for loni=1:size(model_sla_his_yearly_detrended,1)
                    for lati=1:size(model_sla_his_yearly_detrended,2)
                        model_sla_his_std(loni,lati)=std(model_sla_his_yearly_detrended(loni,lati,:));
                    end
                end

                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                model_sla_his_std=model_sla_his_std .* mask_model;    
                
                sl=ncread(interpedfilename,'interped_sla_yearly_exp_fit');
                slr=sl(:,:,end)-sl(:,:,1);
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                slr=slr.*mask_model;

                tide_info.name{1}='M2  ';
                tide_info.name{2}='S2  ';
                tide_info.name{3}='K1  ';
                tide_info.name{4}='O1  ';
                earlytidename=strcat(filedir, '\', num2str(min(inputyear),'%04i'), '\', testname,'_harmonic_analysis_zeta_',regionname, ...
                    '_', num2str(min(inputyear),'%04i'), '.nc');
    %             ncinfo(earlytidename)
                lon_rho = ncread(earlytidename, 'lon_rho');
                lat_rho = ncread(earlytidename, 'lat_rho');
                tname=ncread(earlytidename, 'tname')';
                num_tide_all=size(tname,1);
                num_tide_tgt=length(tide_info.name);
                for coni=1:num_tide_all
                    for tide_namei=1:num_tide_tgt
                        if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                            tide_info.index(tide_namei)=coni;
                        end
                    end
                end

                tcon=ncread(earlytidename, 'tcon');
                m2_amp=tcon(:,:,tide_info.index(1),1).*100;
                s2_amp=tcon(:,:,tide_info.index(2),1).*100;
                k1_amp=tcon(:,:,tide_info.index(3),1).*100;
                o1_amp=tcon(:,:,tide_info.index(4),1).*100;
                early_sum4=m2_amp+s2_amp+k1_amp+o1_amp;

                latetidename=strcat(filedir, '\', num2str(max(inputyear),'%04i'), '\', testname,'_harmonic_analysis_zeta_',regionname, ...
                    '_', num2str(max(inputyear),'%04i'), '.nc');
                tname=ncread(latetidename, 'tname')';
                num_tide_all=size(tname,1);
                num_tide_tgt=length(tide_info.name);
                for coni=1:num_tide_all
                    for tide_namei=1:num_tide_tgt
                        if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                            tide_info.index(tide_namei)=coni;
                        end
                    end
                end

                tcon=ncread(latetidename, 'tcon');
                m2_amp=tcon(:,:,tide_info.index(1),1).*100;
                s2_amp=tcon(:,:,tide_info.index(2),1).*100;
                k1_amp=tcon(:,:,tide_info.index(3),1).*100;
                o1_amp=tcon(:,:,tide_info.index(4),1).*100;
                late_sum4=m2_amp+s2_amp+k1_amp+o1_amp;
                diff_sum4=late_sum4 - early_sum4;

                diff_sum4_interped = griddata(double(lon_rho), double(lat_rho), diff_sum4,double(cut_lon_rho),double(cut_lat_rho));   
                diff_sum4_interped=diff_sum4_interped .* mask_model; 

                slr_ind=(slr+diff_sum4_interped) ./model_sla_his_std;

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('slr-ratio exp fit cmems, ',testname, ',(', ...
                    num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                caxis([0 70]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_flags{65,1}, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            
            
            fig_flag=0;
        end
        
% start-------------------- model slr-ratio index plot (from exp fit yearly SLR) (all test)
        fig_flag=fig_flags{66,2};
        figalldir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\','all','\',regionname,'\slr-ratio\'); % % where figure files will be saved
        while (fig_flag)
            outfile = strcat(figalldir,regionname);
            
            jpgname=strcat(outfile, '_', all_testname2{1},'_',all_testname2{end},'_',regionname, '_slr-ratio_ind_exp_SLR', ...
                num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                if (exist(strcat(figalldir) , 'dir') ~= 7)
                    mkdir(strcat(figalldir));
                end
                for temp_testind=1:length(all_testname2)
                    temp_testname=all_testname2{temp_testind};
                    temp_filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', temp_testname, '\run\'); % % where data files are
%                         temp_filename=[temp_filedir, num2str(tempyear,'%04i'), '\', temp_testname, '_harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
%                         temp_tcon=ncread(temp_filename, 'tcon');
%                         all_m2_amp(1:size(temp_tcon,1), 1:size(temp_tcon,2), temp_testind)=temp_tcon(:,:,tide_info.index(1),1).*100;
                    temp_diffname=[temp_filedir, temp_testname, '_slr_ind_', num2str(max(inputyear),'%04i'), '-', num2str(min(inputyear),'%04i'), '.mat'];
                    clear slr_ind
                    load(temp_diffname)
                    all_slr_ind(:,:,temp_testind)=slr_ind;
                end
                mean_slr_ind=mean(all_slr_ind,3);

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(mean_slr_ind(:,:)'));
                
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('slr-ratio, ',testname, ',(', ...
                    num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(mean_slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                caxis([0 70]);

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
                clear all_slr_ind
                hold off
                close all;
            end
            
            
            fig_flag=0;
        end
        
% start-------------------- SLR plot from yearly exponential fit(2006-2100) data
        fig_flag=fig_flags{67,2};
        while (fig_flag)
            figdir2=[figrawdir,'slr-ratio\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            
            tgtyear_slr=2070;
            
            jpgname=strcat(outfile, '_', testname,'_',regionname, 'yearly_SLR_', num2str(tgtyear_slr, '%04i'), 'exp_fit_', ...
                num2str(min(inputyear), '%04i'), '_', num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
              
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                sl=ncread(interpedfilename,'interped_sla_yearly_exp_fit');
                
                endind=tgtyear_slr-min(inputyear)+1;
                maxslr=sl(:,:,end)-sl(:,:,1);
                slr=sl(:,:,endind)-sl(:,:,1);
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
                
%                 for i=1:size(sl,3)
%                     sl=sl.*mask_model;
%                 end
%                 sly=squeeze(mean(mean(sl,1,'omitnan'),2,'omitnan'))
                maxslr=maxslr.*mask_model;
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
                    num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), ', ', num2str(round(mean(slr_ind(:),'omitnan'),2)),') ' );  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                
                cbar_ratio=mean(slr(:), 'omitnan')/mean(maxslr(:), 'omitnan')
%                  2030 ens09=0.1980 ens08=0.1958 ens10=0.1823
%                  2050 ens09=0.3896 ens08=0.3864 ens10=0.3668
%                  2070 ens09=0.6091 ens08=0.6060 ens10=0.5871

                switch testname
                    case {'test61', 'test62', 'test63', 'test64','ens09'}
                        colormap(parula);
                        caxis([30 50].*cbar_ratio);
                    case {'test57', 'test58', 'test59', 'test60','ens08'}
                        colormap(parula);
                        caxis([30 70].*cbar_ratio);
                    case {'test65', 'test66', 'test67', 'test68','ens10'}
                        colormap(flip(autumn));
                        caxis([60 100].*cbar_ratio);
                end

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
    end
end