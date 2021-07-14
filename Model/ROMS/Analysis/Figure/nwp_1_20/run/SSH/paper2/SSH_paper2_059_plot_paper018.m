close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'test53', 'test54', 'test55', 'test56', 'ens03'};
all_testname2 = {'test53', 'test54', 'test55', 'test56'};
% all_testname2 = {'test54'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP4', 'YS', 'YSECS', 'ECS2', 'NES', 'SES', 'ES'}

% all_region2 ={'NWP', 'AKP4'};
all_region2 ={'AKP4'};

% all_region2 ={'TEST'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 mean_rms_ysecs mean_rms_es ...
             mean_gcm_rms_es  mean_gcm_rms_ysecs w_mean_rms_ysecs w_mean_rms_es w_mean_gcm_rms_es  w_mean_gcm_rms_ysecs
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
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\SSH\function']));
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
        gcmtestname = Func_0004_get_GCMname_from_RCM(testname);
        scenname= Func_0003_RCM_CMIP5_scenname(testname);
        inputyear = [1982:2005]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]

        varname ='zeta';
        variable='SSH';
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        
% % %         switch region
        [refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);
        
% % %         flag configuration        
        fig_flags{71,1}='surf_u rms';
        for flagi=1:100
            fig_flags{flagi,2}=0;
        end
        fig_flags{1,2}=1;
        fig_flags{70,2}=1;
        
        for flagi=1:100
            fig_flags{flagi,2}=0;
        end
        fig_flags{79,2}=2;
        
%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
%         filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
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
        gcmavhrrfilename = strcat(cmip5dir,'\surface\', scenname, '_', gcmtestname, '_',regionname,'_rms_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

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
                
                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecspolygon);
                rms_ysecs= rms_sst .* ysecsmask;
                mean_rms_ysecs=mean(rms_ysecs(:),'omitnan')
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, espolygon);
                rms_es= rms_sst .* esmask;
                mean_rms_es=mean(rms_es(:),'omitnan')
                
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
                titlename = strcat('rms-summer,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(season_rms(:),'omitnan'),2)), ') ');  %% + glacier contribution
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
        
% start-------------------- sst_rms (gcm region)
        fig_flag=fig_flags{79,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'sst_gcmr_','rms_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                
                rms_sst=ncread(avhrrfilename,'rms');
                gcm_rms_sst = ncread(gcmavhrrfilename,'rms');
                gcm_rms_mask = gcm_rms_sst;
                gcm_rms_mask(isfinite(gcm_rms_sst))=1;
                rms_sst = rms_sst.*gcm_rms_mask;
                comb_interped_sst = ncread(avhrrfilename, 'sst');
                mean_interped_sst = mean(comb_interped_sst, 3).*gcm_rms_mask;
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                mean_avhrr_sst = mean(comb_avhrr_sst, 3);
                
                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecspolygon);
                rms_ysecs= rms_sst .* ysecsmask;
                mean_rms_ysecs=mean(rms_ysecs(:),'omitnan')
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, espolygon);
                rms_es= rms_sst .* esmask;
                mean_rms_es=mean(rms_es(:),'omitnan')
                
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
                titlename = strcat('rms-sst,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(rms_sst(:),'omitnan'),2)), ') ');  %% + glacier contribution
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
        
% start-------------------- summer sst_rms
        fig_flag=fig_flags{80,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'summer_sst_gcmr_','rms_',num2str(min(inputyear),'%04i'), ...
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
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                season_avhrr_sst = comb_avhrr_sst(:,:,ind_season);
                mean_avhrr_sst = mean(season_avhrr_sst, 3);
                comb_mse = ncread(avhrrfilename, 'mse');
                season_mse = comb_mse(:,:,ind_season);
                season_rms = sqrt(mean(season_mse,3));
                gcm_rms_sst = ncread(gcmavhrrfilename,'rms');
                gcm_rms_mask = gcm_rms_sst;
                gcm_rms_mask(isfinite(gcm_rms_sst))=1;
                season_rms = season_rms.*gcm_rms_mask;
                mean_interped_sst = mean(season_interped_sst, 3).*gcm_rms_mask;
                mean_avhrr_sst = mean(season_avhrr_sst, 3).*gcm_rms_mask;
                
                comb_gcm_mse = ncread(gcmavhrrfilename, 'mse');
                season_gcm_mse = comb_gcm_mse(:,:,ind_season);
                season_gcm_rms = sqrt(mean(season_gcm_mse,3));
                
                
                ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
                season_rms_ys= season_rms .* ysmask;
                mean_rms_ys=mean(season_rms_ys(:),'omitnan');
                
                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecs_khoapolygon);
                season_rms_ysecs= season_rms .* ysecsmask;
                mean_rms_ysecs(testnameind2)=mean(season_rms_ysecs(:),'omitnan');
                season_gcm_rms_ysecs = season_gcm_rms .* ysecsmask;
                mean_gcm_rms_ysecs(testnameind2)=mean(season_gcm_rms_ysecs(:), 'omitnan');
                
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, es_khoapolygon);
                season_rms_es= season_rms .* esmask;
                mean_rms_es(testnameind2)=mean(season_rms_es(:),'omitnan');
                season_gcm_rms_es = season_gcm_rms .* esmask;
                mean_gcm_rms_es(testnameind2)=mean(season_gcm_rms_es(:), 'omitnan');
                
%                 pcolor(ysecsmask'); shading flat; colorbar; hold on; pcolor(esmask'); hold off;
                
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
                clear lon_rho season_rms
            end

            fig_flag=0;
        end

% start-------------------- winter sst_rms
        fig_flag=fig_flags{81,2};
        while (fig_flag)
           figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_', 'winter_sst_gcmr_','rms_',num2str(min(inputyear),'%04i'), ...
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
                mean_avhrr_sst = mean(season_avhrr_sst, 3);
                comb_mse = ncread(avhrrfilename, 'mse');
                season_mse = comb_mse(:,:,ind_season);
                season_rms = sqrt(mean(season_mse,3));
                gcm_rms_sst = ncread(gcmavhrrfilename,'rms');
                gcm_rms_mask = gcm_rms_sst;
                gcm_rms_mask(isfinite(gcm_rms_sst))=1;
                season_rms = season_rms.*gcm_rms_mask;
                mean_interped_sst = mean(season_interped_sst, 3).*gcm_rms_mask;
                mean_avhrr_sst = mean(season_avhrr_sst, 3).*gcm_rms_mask;
                
                comb_gcm_mse = ncread(gcmavhrrfilename, 'mse');
                season_gcm_mse = comb_gcm_mse(:,:,ind_season);
                season_gcm_rms = sqrt(mean(season_gcm_mse,3));
                
                
                ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
                season_rms_ys= season_rms .* ysmask;
                mean_rms_ys=mean(season_rms_ys(:),'omitnan')
                season_gcm_rms_ys = season_gcm_rms .* ysmask;
                mean_gcm_rms_ys=mean(season_gcm_rms_ys(:), 'omitnan')
                
                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecs_khoapolygon);
                season_rms_ysecs= season_rms .* ysecsmask;
                w_mean_rms_ysecs(testnameind2)=mean(season_rms_ysecs(:),'omitnan');
                season_gcm_rms_ysecs = season_gcm_rms .* ysecsmask;
                w_mean_gcm_rms_ysecs(testnameind2)=mean(season_gcm_rms_ysecs(:), 'omitnan');
                
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, es_khoapolygon);
                season_rms_es= season_rms .* esmask;
                w_mean_rms_es(testnameind2)=mean(season_rms_es(:),'omitnan');
                season_gcm_rms_es = season_gcm_rms .* esmask;
                w_mean_gcm_rms_es(testnameind2)=mean(season_gcm_rms_es(:), 'omitnan');
                
%                 pcolor(ysecsmask'); shading flat; colorbar; hold on; pcolor(esmask'); hold off;
                
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
                titlename = strcat('rms-winter,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(season_rms(:),'omitnan'),2)), ') ');  %% + glacier contribution
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
                clear lon_rho season_rms
            end

            fig_flag=0;
        end
        
% start-------------------- summer sst_bias
        fig_flag=fig_flags{82,2};
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
                comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
                season_avhrr_sst = comb_avhrr_sst(:,:,ind_season);
%                 comb_mse = ncread(avhrrfilename, 'mse');
%                 season_mse = comb_mse(:,:,ind_season);
                comb_season_bias=season_interped_sst - season_avhrr_sst;
%                 season_rms = sqrt(mean(season_mse,3));
                gcm_rms_sst = ncread(gcmavhrrfilename,'rms');
                gcm_rms_mask = gcm_rms_sst;
                gcm_rms_mask(isfinite(gcm_rms_sst))=1;
                season_bias = mean(comb_season_bias,3).*gcm_rms_mask;
                mean_interped_sst = mean(season_interped_sst, 3).*gcm_rms_mask;
                mean_avhrr_sst = mean(season_avhrr_sst, 3).*gcm_rms_mask;

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
                titlename = strcat('bias-summer,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(season_bias(:),'omitnan'),2)), ') ');  %% + glacier contribution
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
        fig_flag=fig_flags{83,2};
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
                gcm_rms_sst = ncread(gcmavhrrfilename,'rms');
                gcm_rms_mask = gcm_rms_sst;
                gcm_rms_mask(isfinite(gcm_rms_sst))=1;
                season_bias = mean(comb_season_bias,3).*gcm_rms_mask;
                mean_interped_sst = mean(season_interped_sst, 3).*gcm_rms_mask;
                mean_avhrr_sst = mean(season_avhrr_sst, 3).*gcm_rms_mask;

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
                titlename = strcat('bias-winter,','rcm',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),  ', ', num2str(round(mean(season_bias(:),'omitnan'),2)), ') ');  %% + glacier contribution
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