close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'GODAS'};
all_testname2 = {'ORAS5'};

% all_testname2 = {'test60'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}
% all_region2 ={'SES'}

all_region2 ={'AKP4'};
% all_region2 ={'ES'}
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
            dropboxpath='C:\Users\KYY\Dropbox';
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
        abstrendlev =[-4 4];
        absssttrendlev=[-0.02 0.02];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.15 0.15];
        meansstfilteredlev =[-2 2];
        meansstlev = [10 25]; 
        trendplotlev = [0 7];
        trenddifflev = [-10 10];
        sshlev =[-0.3 0.3];
        sshdifflev = [0 20];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1994:2017]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [2 8]; % % put month which you want to plot [month month ...]

        varname ='zeta'
        variable='SSH'
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
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
            otherwise
                ('?')
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));

        % % % for EKB
        % regionname='EKB';
        % lonlat = [127, 129.5, 38, 40.5];

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        filename = ['E:\Data\Reanalysis\', testname, '\',testname,'_',regionname, ...
                    '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);
        
       cmems_trend=ncread(filename, 'cmems_trend');
        cmems_mask=ones(size(cmems_trend));
        cmems_mask(isnan(cmems_trend))=NaN;
        
    % %     valid cell number
    %      for vi=1:size(comb_spatial_meanressh,1)
    %          for vj=1:size(comb_spatial_meanressh,2)
    %              if (isnan(comb_spatial_meanressh(vi,vj,1))~=1)
    %                 valnum=valnum+1;
    %              end
    %          end
    %      end

    %     plot(squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan')))
    %     isize = size(comb_data_filtered,1)
    %     jsize = size(comb_data_filtered,2)
    %     lsize = size(comb_data_filtered,3)
    %     comb_yearly_data_filtered=reshape(comb_data_filtered,[isize, jsize, 12, lsize/12]);
    %     mean_yearly_data_filtered=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
    %     trendtime=14:29;
    %     p=polyfit(trendtime,mean_yearly_data_filtered(14:29)',1);
    %     yearly_interped_trend=p(1);
    %     yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y


        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_soda_', regionname, '.m']
            filedir = strcat('E:\Data\Reanalysis\', testname, '\'); % % where data files are
            cmemsdir='E:\Data\Observation\CMEMS\';
        elseif (strcmp(system_name,'GLNXA64'))
        end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))
        
        run(param_script);
        
        figdir=[figrawdir,'Trend\SSH\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        
        
        lon_rho=ncread(filename,'lon_rho');
        lat_rho=ncread(filename,'lat_rho');
        trend_filtered = ncread(filename,'trend_filtered');
        mean_trend_filtered = ncread(filename,'mean_trend_filtered');
        cmems_trend_filtered = ncread(filename,'cmems_trend_filtered');
        interped_trend_filtered = ncread(filename,'interped_trend_filtered');

        
% start-------------------- corr_interped plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_interped_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            corr_interped=ncread(filename,'corr_interped');
            corr_interped=corr_interped.*cmems_mask;
            lon_cmems=ncread(filename,'lon_cmems');
            lat_cmems=ncread(filename,'lat_cmems');
            
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
% end-------------------- corr_interped plot 

% % start-------------------- corr_interped_sst plot
%         jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_interped_sst_',num2str(min(inputyear),'%04i'), ...
%             '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             corr_interped_sst=ncread(filename,'corr_interped_sst');
%             lon_cmems=ncread(filename,'lon_cmems');
%             lat_cmems=ncread(filename,'lat_cmems');
%             
%             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%             hold on;
%             m_pcolor(lon_cmems',lat_cmems',squeeze(corr_interped_sst(:,:)'));
%             shading(gca,m_pcolor_shading_method);
%             m_gshhs_i('color',m_gshhs_line_color);
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%     %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
%             titlename = strcat('corr, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
%                 num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped_sst(:), 'omitnan'),2)));  %% + glacier contribution
% 
%             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%             % set colorbar 
%             h = colorbar;
%             colormap(wrmap);
%             set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis([0.2 0.8]);
% 
%             % set grid
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%             saveas(gcf,jpgname,'tif');
% 
%             disp(' ')
%             disp(['corr_', 'interped_sst_', ' plot is created.'])
%             disp(' ')
%             disp([' File path is : ',jpgname])
%             disp(' ')
% 
%             hold off
%             close all;
%         end
% % end-------------------- corr_interped_sst plot 


% start-------------------- corr_interped_filtered plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_interped_filtered_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            corr_interped_filtered=ncread(filename,'corr_interped_filtered');
            lon_cmems=ncread(filename,'lon_cmems');
            lat_cmems=ncread(filename,'lat_cmems');
            
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
% end-------------------- corr_interped_filtered plot plot

% start-------------------- corr_spatial mean plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_spatial_mean_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
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
% end-------------------- corr_spatial mean plot


% % start-------------------- corr_interped_sst_filtered plot
%         jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_interped_sst_filtered_',num2str(min(inputyear),'%04i'), ...
%             '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             corr_interped_sst_filtered=ncread(filename,'corr_interped_sst_filtered');
%             lon_cmems=ncread(filename,'lon_cmems');
%             lat_cmems=ncread(filename,'lat_cmems');
%             
%             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%             hold on;
%             m_pcolor(lon_cmems',lat_cmems',squeeze(corr_interped_sst_filtered(:,:)'));
%             shading(gca,m_pcolor_shading_method);
%             m_gshhs_i('color',m_gshhs_line_color);
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%     %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
%             titlename = strcat('corr sst filtered, ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
%                 num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(corr_interped_sst_filtered(:), 'omitnan'),2)));  %% + glacier contribution
% 
%             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%             % set colorbar 
%             h = colorbar;
%             colormap(wrmap);
%             set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis([0.2 0.8]);
% 
%             % set grid
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%             saveas(gcf,jpgname,'tif');
% 
%             disp(' ')
%             disp(['interped_sst_filtered', ' plot is created.'])
%             disp(' ')
%             disp([' File path is : ',jpgname])
%             disp(' ')
% 
%             hold off
%             close all;
%         end
% % end-------------------- corr_interped_sst_filtered plot plot



% start-------------------- corr_interped_detrended plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_interped_detrended_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
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
% end-------------------- corr_interped_detrended plot plot


% start-------------------- corr_ lowpassed plot
        nyears=[2,3,5];
        nc_varname_prefixes={'interped_', 'interped_detrended_'};
        nc_titlename_prefixes={'lowpass', 'lowpass-det'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
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
% end-------------------- corr_ lowpassed plot


% start-------------------- corr_clim plot
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        for monthij = 1:12
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_corr_clim_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
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
% end-------------------- corr_clim plot plot

        clim_trend = ncread(filename,'clim_ssh_trend');
        mean_clim_trend = squeeze(mean(mean(clim_trend, 1, 'omitnan'), 2, 'omitnan'));


% start-------------------- absolute trend plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
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
% end-------------------- absolute trend plot

% start-------------------- absolute trend plot (wr)
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_wr_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
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
            caxis([2 7]);

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
% end-------------------- absolute trend plot (wr)

%  start-------------------- absolute trend difference plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_diff_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            lon_cmems=ncread(filename,'lon_cmems');
            lat_cmems=ncread(filename,'lat_cmems');
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
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
% end-------------------- absolute trend difference plot

% % start-------------------- absolute sst trend filtered plot
%         jpgname=strcat(outfile, '_', testname,'_',regionname, '_sst_trend_', ...
%             num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             interped_sst_trend_filtered=ncread(filename,'interped_sst_trend_filtered')/1000.0;
%             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%             hold on;
%             lon_cmems=ncread(filename,'lon_cmems');
%             lat_cmems=ncread(filename,'lat_cmems');
%             m_pcolor(lon_cmems',lat_cmems',squeeze(interped_sst_trend_filtered(:,:)'));
%             shading(gca,m_pcolor_shading_method);
%             m_gshhs_i('color',m_gshhs_line_color);
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%             mean_sst_trend_filtered=mean(interped_sst_trend_filtered(:),'omitnan');
%     %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
%             titlename = strcat('SST trend, ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
%                 num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_sst_trend_filtered,4)), ' ^o/y');  %% + glacier contribution
% 
%             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%             % set colorbar 
%             h = colorbar;
%             colormap(bwrmap);
%             set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis(absssttrendlev);
% 
%             % set grid
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%             saveas(gcf,jpgname,'tif');
% 
%             disp(' ')
%             disp(['clim_', 'absolute_sst_trend', ' plot is created.'])
%             disp(' ')
%             disp([' File path is : ',jpgname])
%             disp(' ')
% 
%             hold off
%             close all;
%         end
% % end-------------------- absolute sst trend filtered plot

% start-------------------- climatological absolute trend plot
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        for monthij = 1:12
            jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
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
                colormap(bwrmap);
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
% end-------------------- climatological absolute trend plot
   
% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
% end-------------------- make timedata for time series  

% start-------------------- msl time series (seasonal filtered)
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            for varind=1:length(inputyear)*12
                msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
            end
    %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
            msl_filt=msl_filt-mean(msl_filt);    
            
            
            p=polyfit(xData,msl_filt,1);
            msl2=xData*p(1)+p(2);
            mslplot=plot(xData,msl_filt,'k')
            hold on
            mslplot2=plot(xData,msl2,'Color','r')
            xlabel('year')
            ylabel('Mean SSH (m)')
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend_filtered,2)), ' mm/y'])
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
% end-------------------- msl time series (seasonal filtered)

% % start-------------------- msst time series (seasonal filtered)
%         jpgname=strcat(outfile, '_', testname, '_',regionname, '_msst_', ...
%             num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             for varind=1:length(inputyear)*12
%                 msst_filt(varind)=squeeze(mean(mean(ncread(filename,'sst_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%             end
%     %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
% %             msst_filt=msst_filt-mean(msst_filt);    
%             
%             sst_trend_filtered=ncread(filename,'interped_sst_trend_filtered')/1000.0;
%             mean_sst_trend_filtered=mean(sst_trend_filtered(:), 'omitnan');
%             p=polyfit(xData,msst_filt,1);
%             msst2=xData*p(1)+p(2);
%             msstplot=plot(xData,msst_filt,'k')
%             hold on
%             msstplot2=plot(xData,msst2,'Color','r')
%             xlabel('year')
%             ylabel('Mean SSH (m)')
%             title([regionname, ', Mean SST(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_sst_trend_filtered,3)), ' ^oC/y'])
%             datetick('x','yyyy','keepticks')
%             axis tight;
%             ylim(meansstfilteredlev)
%             set(msstplot,'LineWidth',2);
%             set(gca,'FontSize',15);
%             grid on
%             hold off
%             saveas(gcf,jpgname,'jpg');
%             grid off
%         end
%         close all;
% % end-------------------- msst time series (seasonal filtered)


% start-------------------- msl time series
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            for varind=1:length(inputyear)*12
                msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
            end
            msl=msl-mean(msl);     
            
            p=polyfit(xData,msl,1);
            msl2=xData*p(1)+p(2);
            mslplot=plot(xData,msl,'k')
            hold on
            mslplot2=plot(xData,msl2,'Color','r')
            xlabel('year')
            ylabel('Mean SSH (m)')
            mean_trend=ncread(filename,'mean_trend');
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend,2)), ' mm/y'])
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
% end-------------------- msl time series

% % start-------------------- msst time series
%         jpgname=strcat(outfile, '_', testname, '_',regionname, '_msst_nonfilt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             for varind=1:length(inputyear)*12
%                 msst(varind)=squeeze(mean(mean(ncread(filename,'raw_sst',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%             end
% %             msst=msst-mean(msst);     
%             sst_trend=ncread(filename,'interped_sst_trend')/1000.0;
%             mean_sst_trend=mean(sst_trend(:), 'omitnan');
%             p=polyfit(xData,msst,1);
%             msst2=xData*p(1)+p(2);
%             msstplot=plot(xData,msst,'k')
%             hold on
%             msstplot2=plot(xData,msst2,'Color','r')
%             xlabel('year')
%             ylabel('Mean SSH (m)')
%             mean_trend=ncread(filename,'mean_trend');
%             title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_sst_trend,3)), ' ^oC/y'])
%             datetick('x','yyyy','keepticks')
%             axis tight;
%             ylim(meansstlev)
%             set(msstplot,'LineWidth',2);
%             set(gca,'FontSize',15);
%             grid on
%             hold off
%             saveas(gcf,jpgname,'jpg');
%             grid off
%             close all;
%         end
% % end-------------------- msst time series


% start-------------------- lowpassed msl time series
        nyears=[2,3,5];
        nc_varname_prefixes={'cmems_', 'interped_'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', testname, '_',regionname, '_',nc_varname,'_msl_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
                    for varind=1:length(inputyear)*12
                        eval(['tempmsl=squeeze(ncread(filename,', '''', nc_varname, '''', ',[1,1,varind], [inf inf 1]));']);
                        msl(varind)=mean(tempmsl(:),'omitnan');
                        if (strcmp(nc_varname_prefix, 'interped_')==1)
                            temp_raw_msl=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
                        elseif (strcmp(nc_varname_prefix, 'cmems_')==1)
                            temp_raw_msl=ncread(filename,'cmems_sla',[1 1 varind], [inf inf 1]);
                        end
                        raw_msl(varind)=squeeze(mean(temp_raw_msl(:),'omitnan'));
                    end
                    msl=msl-mean(msl,'omitnan');    
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
% end-------------------- lowpassed msl time series

% start-------------------- detrended lowpassed msl time series
        nyears=[2,3,5];
        nc_varname_prefixes={'cmems_detrended_', 'interped_detrended_'};
        for nyeari=1:length(nyears)
            nyear=nyears(nyeari);
            for nc_varnameij=1:length(nc_varname_prefixes)
                nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                jpgname=strcat(outfile, '_', testname, '_',regionname, '_',nc_varname,'_msl_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
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
                    msl=msl-mean(msl,'omitnan');    
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
% end-------------------- detrended lowpassed msl time series



% start-------------------- climatological msl time series
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        for monthij=1:12
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
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
% end-------------------- climatological msl time series

% start-------------------- climatological msl trend
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_trend_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
        if abs(mean_clim_trend(1))<0.01
            mean_clim_trend=mean_clim_trend*1000.0  %% m -> mm
        end
        if (exist(jpgname , 'file') ~= 2)
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
% end-------------------- climatological msl trend

% start-------------------- seasonal msl
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_msl_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
        
%         clim_ssh=ncread(filename,'clim_ssh');
%         for ijij=1:12
%             temp_clim_ssh=clim_ssh(:,:,ijij);
%             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
%         end
%         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
        
        if (exist(jpgname , 'file') ~= 2)
            comb_interped_data=ncread(filename,'interped_ssh');
            lon_cmems=ncread(filename, 'lon_cmems');
            lat_cmems=ncread(filename, 'lat_cmems');
            len_lon=length(lon_cmems);
            len_lat=length(lat_cmems);
            comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
            clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
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
% end-------------------- seasonal msl


% start-------------------- cmems seasonal msl
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
        
        if (exist(jpgname , 'file') ~= 2)
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
% end-------------------- cmems seasonal msl



% start-------------------- get distance weighted msl
%         for i= 1:length(recon_lon)
%             for j=1:length(recon_lat)
%                 for l=1:size(comb_data_filtered,3)
%                     comb_data_filtered_weight(i,j,l)=comb_data_filtered(i,j,l).*cos(recon_lat(j)/180.0*pi);
%                 end
%             end
%         end
%         
%         temp_varvar=comb_data_filtered_weight(:,:,1);
%         for i= 1:length(recon_lon)
%             for j=1:length(recon_lat)
%                 temp_varvar(i,j)=cos(recon_lat(j)/180.0*pi).*temp_varvar(i,j)./temp_varvar(i,j);
%             end
%         end
%         for i= 1:length(recon_lon)
%             for j=1:length(recon_lat)
%                 temp_varvar(isnan(temp_varvar))=NaN;
%             end
%         end
%         m_cidfw=squeeze(sum(sum(comb_data_filtered_weight,1,'omitnan'),2,'omitnan'));
%         m_cidfw=m_cidfw./sum(sum(temp_varvar,1,'omitnan'),2,'omitnan');
% end-------------------- get distance weighted msl        
        



% start-------------------- climatological msl time series

        for i =1:length(inputyear) 
            tempyear=inputyear(i);
    %         for month=1:12
                xData2(i) = datenum([num2str(6,'%02i'),'-30-',num2str(tempyear)]);
    %         end
        end
        %     plot(squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan')))
%         isize = size(comb_data_filtered,1)
%         jsize = size(comb_data_filtered,2)
%         lsize = size(comb_data_filtered,3)
%         comb_yearly_data_filtered=reshape(comb_data_filtered,[isize, jsize, 12, lsize/12]);
%         mean_yearly_data_filtered=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
        
%         trendtime=1:length(xData2);
%         p=polyfit(trendtime,mean_yearly_data_filtered(1:length(xData2))',1);
%         yearly_interped_trend=p(1);
%         yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y
% 
%         yearly_msl=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         yearly_msl=yearly_msl-mean(yearly_msl);    
%         p=polyfit(xData2,yearly_msl',1);
%         yearly_msl2=xData2*p(1)+p(2);
%         yearly_mslplot=plot(xData2,yearly_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
%         hold on
%         yearly_mslplot2=plot(xData2,yearly_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('Year')
%         ylabel('Mean SSH (m)')
%         title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(yearly_interped_trend,2)), ' mm/y'])
%         ylim(meanplotlev)
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(yearly_mslplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
%         close all;


% start-------------------- RMS plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_rms_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
% %             if (exist(pngname , 'file') ~= 2)        
% %                 run(param_script);
% %                 for yearij=1:length(inputyear)
% %                     tempyear=inputyear(yearij);
% %                     yearstr=num2str(tempyear, '%04i');
% %                     for monthij=1:length(inputmonth)
% %                         tempmonth=inputmonth(monthij);
% %                         monthstr=num2str(tempmonth, '%02i');
% %                         varind=((yearij-1)*12)+monthij
% %                         if (exist('lon_min' , 'var') ~= 1)
% %                             lon_cmems=ncread(filename, 'lon_cmems');
% %                             lat_cmems=ncread(filename, 'lat_cmems');
% %                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
% %                             [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
% %                             cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% %                             cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% %                         end
% % %                         data_info = ncinfo(filename, varname); 
% %                         data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
% %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% %                         model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
% %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% %                         if (exist('mean_data' , 'var') ~= 1)
% %                             mean_data=zeros(size(data));
% %                             mean_model_data=mean_data;
% %                         end
% %                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
% %                         mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
% %                     end
% %                 end
% %                 cmems_msl=mean(mean_data(:), 'omitnan');
% %                 model_msl=mean(mean_model_data(:), 'omitnan');
% %                 m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% %                 hold on;
% %                 
% %                 mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
% %                 mask_model(mask_model==0)=NaN;
% %                 mean_data= mean_data .* mask_model;
% %                 mean_model_data = mean_model_data .* mask_model;
% %                 mean_data=mean_data-mean(mean_data(:),'omitnan');
% %                 mean_model_data=mean_model_data-mean(mean_model_data(:),'omitnan');
% %                 mean_rms = sqrt((mean_model_data - mean_data).^2) * 100.0;
% % %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
% % %                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%             if (exist(pngname , 'file') ~= 2)        
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


                m_pcolor(cut_lon_rho',cut_lat_rho',mean_rms');
%                 m_pcolor(cut_lon_rho(1:30,1:20)',cut_lat_rho(1:30,1:20)',mean_rms(1:30,1:20)');


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
                caxis([0 15]);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
%             end
% end-------------------- RMS plot

% start-------------------- BIAS plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_bias_', 'SSH','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
           if (exist(pngname , 'file') ~= 2)        
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
                clear lon_rho mean_data
            end
% end-------------------- BIAS plot

% start-------------------- model variance plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);

            pngname=strcat(outfile, '_', testname,'_',regionname, '_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
           if (exist(pngname , 'file') ~= 2)        
                run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end

                model_sla=ncread(filename,'interped_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
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
                titlename = strcat(variable, ' VAR, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 2]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
% end-------------------- model variance plot

% start-------------------- model variance_seasonal plot
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                pngname=strcat(outfile, '_', testname,'_',regionname, '_', 'VAR','_',num2str(min(inputyear),'%04i'), ...
                    '_',num2str(max(inputyear),'%04i'),'_', num2str(tempmonth,'%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(pngname , 'file') ~= 2)        
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
% end-------------------- model variance_seasonal plot


% start-------------------- model variance_lowpass plot
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
                
                    pngname=strcat(outfile, '_', testname,'_',regionname, '_', nc_varname,'_','VAR_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(pngname , 'file') ~= 2)        
                        run(param_script);
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(filename, 'lon_cmems');
                            lat_cmems=ncread(filename, 'lat_cmems');
                            [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                            cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
                        eval(['model_sla', '=ncread(filename,', '''', nc_varname,'''',');']);
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
% end-------------------- model variance_lowpass plot


% % % % % start-------------------- YS_ msl plot
% % % %             figdir2=[figrawdir,'CLIM\'];
% % % %             if (exist(strcat(figdir2) , 'dir') ~= 7)
% % % %                 mkdir(strcat(figdir2));
% % % %             end 
% % % %             outfile = strcat(figdir2,regionname);
% % % % 
% % % %             pngname=strcat(outfile, '_', testname,'_',regionname, '_ys_msl','_',num2str(min(inputyear),'%04i'), ...
% % % %                 '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
% % % %             if (exist(pngname , 'file') ~= 2)        
% % % %                 run(param_script);
% % % %                 for yearij=1:length(inputyear)
% % % %                     tempyear=inputyear(yearij);
% % % %                     yearstr=num2str(tempyear, '%04i');
% % % %                     for monthij=1:length(inputmonth)
% % % %                         tempmonth=inputmonth(monthij);
% % % %                         monthstr=num2str(tempmonth, '%02i');
% % % %                         varind=((yearij-1)*12)+monthij
% % % %                         if (exist('lon_min' , 'var') ~= 1)
% % % %                             lon_cmems=ncread(filename, 'lon_cmems');
% % % %                             lat_cmems=ncread(filename, 'lat_cmems');
% % % %                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
% % % %                             [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
% % % %                             cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % % %                             cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % % %                         end
% % % % %                         data_info = ncinfo(filename, varname); 
% % % %                         data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
% % % %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% % % %                         model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
% % % %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% % % %                         if (exist('mean_data' , 'var') ~= 1)
% % % %                             mean_data=zeros(size(data));
% % % %                             mean_model_data=mean_data;
% % % %                         end
% % % %                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
% % % %                         mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
% % % %                     end
% % % %                 end
% % % %                 cmems_msl=mean(mean_data(:), 'omitnan');
% % % %                 model_msl=mean(mean_model_data(:), 'omitnan');
% % % % %                 hold on;
% % % %                 
% % % % %                 ys_data=ncread(filename,'cmems_adt', [1 1 1], [30 20 inf]);
% % % % %                 ys_model_data=ncread(filename,'interped_ssh', [1 1 1], [30 20 inf]);
% % % % %                 for l=1:288
% % % % %                     tmp_ys_data=ys_data(:,:,l);
% % % % %                     tmp_ys_model_data=ys_model_data(:,:,l);
% % % % %                     ys_mean_data(l)=mean(tmp_ys_data(:),'omitnan') - cmems_msl;
% % % % %                     ys_mean_model_data(l)=mean(tmp_ys_model_data(:),'omitnan') - model_msl;
% % % % %                 end
% % % % 
% % % % %                 plot(ys_mean_data - 0.6291);
% % % % %                 hold on
% % % % %                 plot(ys_mean_model_data - 0.0840);
% % % % %                 hold off
% % % %                 
% % % %                 mslplot=plot(ys_mean_data .*100.0,'k')
% % % %                 hold on
% % % %                 msl2plot=plot(ys_mean_model_data.*100.0,'r')
% % % %                 xlabel('month')
% % % %                 ylabel('msl(cm)')
% % % %                 title(['YS-SW', ', msl(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
% % % %                 axis tight;
% % % %                 ylim([-20 40])
% % % %                 set(mslplot,'LineWidth',2);
% % % %                 set(msl2plot,'LineWidth',2);
% % % % 
% % % %                 set(gca,'FontSize',20);
% % % %                 grid on
% % % %                 hold off
% % % %                 grid off
% % % % %                 lgd=legend('CMEMS','SODA');                
% % % %                 
% % % %                 saveas(gcf,pngname,'tif');
% % % %                 close all;
% % % %                 clear lon_rho mean_data
% % % %             end
% % % % % end-------------------- YS_ msl plot
% % % % 
% % % % % start-------------------- YS_ msl_diff plot
% % % %             figdir2=[figrawdir,'CLIM\'];
% % % %             if (exist(strcat(figdir2) , 'dir') ~= 7)
% % % %                 mkdir(strcat(figdir2));
% % % %             end 
% % % %             outfile = strcat(figdir2,regionname);
% % % % 
% % % %             pngname=strcat(outfile, '_', testname,'_',regionname, '_ys_msl_diff','_',num2str(min(inputyear),'%04i'), ...
% % % %                 '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
% % % %             if (exist(pngname , 'file') ~= 2)        
% % % %                 run(param_script);
% % % %                 for yearij=1:length(inputyear)
% % % %                     tempyear=inputyear(yearij);
% % % %                     yearstr=num2str(tempyear, '%04i');
% % % %                     for monthij=1:length(inputmonth)
% % % %                         tempmonth=inputmonth(monthij);
% % % %                         monthstr=num2str(tempmonth, '%02i');
% % % %                         varind=((yearij-1)*12)+monthij
% % % %                         if (exist('lon_min' , 'var') ~= 1)
% % % %                             lon_cmems=ncread(filename, 'lon_cmems');
% % % %                             lat_cmems=ncread(filename, 'lat_cmems');
% % % %                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
% % % %                             [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
% % % %                             cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % % %                             cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
% % % %                         end
% % % % %                         data_info = ncinfo(filename, varname); 
% % % %                         data(:,:)=ncread(filename,'cmems_adt',[lon_min(1) lat_min(1) varind], ...
% % % %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% % % %                         model_data(:,:)=ncread(filename,'interped_ssh',[lon_min(1) lat_min(1) varind], ...
% % % %                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% % % %                         if (exist('mean_data' , 'var') ~= 1)
% % % %                             mean_data=zeros(size(data));
% % % %                             mean_model_data=mean_data;
% % % %                         end
% % % %                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
% % % %                         mean_model_data=mean_model_data + (model_data / (length(inputyear) * length(inputmonth)));
% % % %                     end
% % % %                 end
% % % %                 cmems_msl=mean(mean_data(:), 'omitnan');
% % % %                 model_msl=mean(mean_model_data(:), 'omitnan');
% % % % %                 hold on;
% % % %                 
% % % %                 ys_data=ncread(filename,'cmems_adt', [1 1 1], [30 20 inf]);
% % % %                 ys_model_data=ncread(filename,'interped_ssh', [1 1 1], [30 20 inf]);
% % % %                 for l=1:288
% % % %                     tmp_ys_data=ys_data(:,:,l);
% % % %                     tmp_ys_model_data=ys_model_data(:,:,l);
% % % %                     ys_mean_data(l)=mean(tmp_ys_data(:),'omitnan') - cmems_msl;
% % % %                     ys_mean_model_data(l)=mean(tmp_ys_model_data(:),'omitnan') - model_msl;
% % % %                 end
% % % % %                 plot(ys_mean_data - 0.6291);
% % % % %                 hold on
% % % % %                 plot(ys_mean_model_data - 0.0840);
% % % % %                 hold off
% % % %                 
% % % %                 mslplot=plot(ys_mean_model_data .*100.0 - ys_mean_data .*100,'k')
% % % %                 hold on
% % % %                 xlabel('month')
% % % %                 ylabel('msl(cm)')
% % % %                 title(['YS-SW', ', msl diff(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
% % % %                 axis tight;
% % % %                 ylim([0 25])
% % % %                 set(mslplot,'LineWidth',2);
% % % % %                 set(msl2plot,'LineWidth',2);
% % % % 
% % % %                 set(gca,'FontSize',20);
% % % %                 grid on
% % % %                 hold off
% % % %                 
% % % %                 saveas(gcf,pngname,'tif');
% % % %                 close all;
% % % %                 clear lon_rho mean_data
% % % %             end
% % % % % end-------------------- YS_ msl_diff plot


    end
end