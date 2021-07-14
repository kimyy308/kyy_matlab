close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
all_testname2 = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
all_region2 ={'NWP', 'AKP2'}

% all_region2 ={'AKP2'};

% all_region2 ={'NWP'}
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
        abstrendlev =[4 7];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.3 0.3];
        trendplotlev = [4 6.5];
        sshlev =[-0.7 1.3];
        sshdifflev = [30 70];
        
        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
        inputyear1 = [2006:2015];
        inputyear2 = [2091:2100];
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

        varname ='zeta'
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
            case('CA') %% Around Korea Peninsula
                refpolygon=capolygon;
            case('EKB') %% Around Korea Peninsula
                refpolygon=akp2polygon;
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

%         load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
          load(['G:\Data\Model\CMIP5\zos\rcp45\Omon\',testname,'\',testname,'_',regionname, ...
                    'ssh_trend_rcp45_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
          filename = ['G:\Data\Model\CMIP5\zos\rcp45\Omon\',testname,'\',testname,'_',regionname, ...
                    '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
        valnum=0;
        wrmap = bwrmap(51:100,:);
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


%         if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\CMIP5\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m']
            filedir = strcat('G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
%         elseif (strcmp(system_name,'GLNXA64'))
%         end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))
        
        var='SSH';
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
        
% start-------------------- 1 year msl plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_1y_ssh_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            run(param_script);
            for varind=1:12
                msl_month(:,:,varind)=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
            end
            msl_1y=squeeze(mean(msl_month,3,'omitnan'));
            m_msl_1y=mean(mean(msl_1y,'omitnan'),'omitnan');
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(msl_1y')-m_msl_1y);
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH, ',testname,',(',num2str(min(inputyear),'%04i'),')');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'(m)','fontsize',colorbar_title_fontsize);
            caxis(colorbar_lev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['1y_', 'msl', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- 1 year msl plot

% start-------------------- 100 year msl plot (1y ssh + trend*year)
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_100y_ssh_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            for varind=1:12
                msl_month(:,:,varind)=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
            end
            msl_1y=squeeze(mean(msl_month,3,'omitnan'));
            m_msl_1y=mean(mean(msl_1y,'omitnan'),'omitnan');
            diff_2100=(trend_filtered/1000.0)*(2100-2006);
            diff_cm =mean_trend_filtered*10.0;
           
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(msl_1y')-m_msl_1y+diff_2100');
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH, ',testname,' (',num2str(2100,'%04i'),')');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'(m)','fontsize',colorbar_title_fontsize);
            caxis(colorbar_lev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['1y_', 'msl', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- 100 year msl plot

% start-------------------- 2100 - 2005 year diff plot (1y ssh + trend*year)
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_diff_ssh_2005_2100', '.tif'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
            for varind=1:12
                msl_month(:,:,varind)=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
            end
            msl_1y=squeeze(mean(msl_month,3,'omitnan'));
            m_msl_1y=mean(mean(msl_1y,'omitnan'),'omitnan');
            diff_2100=(trend_filtered/1000.0)*(2100-2006);
            diff_cm =mean_trend_filtered*(2100-2006)/10.0;
           
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',diff_2100'*(100));  %% (m -> cm)
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH, ',testname,' (',num2str(2100,'%04i'),') rise =', num2str(round(diff_cm)), 'cm');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
            caxis(sshdifflev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['1y_', 'msl', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
%         end
% end-------------------- 2100 - 2005 year diff plot (1y ssh + trend*year)


% start-------------------- first year msl plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_1st_y_ssh_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            for varind=1:12
                msl_month(:,:,varind)=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
            end
            msl_1y=squeeze(mean(msl_month,3,'omitnan'));
            m_msl_1y=mean(mean(msl_1y,'omitnan'),'omitnan');
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(msl_1y')-m_msl_1y);
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH, ',testname,',(',num2str(min(inputyear),'%04i'),')');  %% + glacier contribution
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'(m)','fontsize',colorbar_title_fontsize);
            caxis(colorbar_lev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['1y_', 'msl', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
% end-------------------- first year msl plot


% start-------------------- last year & diff msl plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_last_y_ssh_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        varind_s=(max(inputyear)-min(inputyear))*12 +1;
        varind_e=(max(inputyear)-min(inputyear))*12 +12;
%         if (exist(jpgname , 'file') ~= 2)
            for varind=1:12
                msl_month_1y(:,:,varind)=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
            end
            msl_1y=squeeze(mean(msl_month_1y,3,'omitnan'));
            for varind=varind_s:varind_e
                msl_month(:,:,varind-varind_s+1)=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
            end
            msl_ly=squeeze(mean(msl_month,3,'omitnan'));
            m_msl_1y=mean(mean(msl_1y,'omitnan'),'omitnan');
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',squeeze(msl_ly')-m_msl_1y);
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH, ',testname,',(',num2str(max(inputyear),'%04i'),')');  %% + glacier contribution
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'(m)','fontsize',colorbar_title_fontsize);
            caxis(colorbar_lev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['1y_', 'msl', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
  
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_diff_1st_last_y_ssh_',num2str(min(inputyear),'%04i'), ...
                        '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg

            diff_cm=(msl_ly - msl_1y)*100.0;
            
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(lon_rho)',lat_rho',diff_cm');
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('dSSH, ',testname,' (',num2str(max(inputyear),'%04i'),') rise =', num2str(round(mean(diff_cm(:),'omitnan'))), 'cm');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'(cm)','fontsize',colorbar_title_fontsize);
            caxis(sshdifflev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['1y_', 'msl', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
%         end

% end-------------------- last year & diff msl plot
        
% % start-------------------- relative trend plot
%         jpgname=strcat(outfile, '_', testname,'_',regionname, '_relative_ssh_trend_',num2str(min(inputyear),'%04i'), ...
%             '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%             hold on;
%             m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'-mean_trend_filtered));
%     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
%     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
%             shading(gca,m_pcolor_shading_method);
%             m_gshhs_i('color',m_gshhs_line_color);
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%     %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
%             titlename = strcat('SSH trend(rel), ',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
%                 num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution
% 
%             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%             % set colorbar 
%             h = colorbar;
%             colormap(bwrmap);
%             set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis(reltrendlev);
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
%             disp(['clim_', 'relative_ssh_trend', ' plot is created.'])
%             disp(' ')
%             disp([' File path is : ',jpgname])
%             disp(' ')
% 
%             hold off
%             close all;
%         end
% % end-------------------- relative trend plot

        clim_trend = ncread(filename,'clim_ssh_trend');
        mean_clim_trend = squeeze(mean(mean(clim_trend, 1, 'omitnan'), 2, 'omitnan'));
% % start-------------------- climatological relative trend plot
%         climdir = [figdir,'\CLIM\'];
%         if (exist(strcat(climdir) , 'dir') ~= 7)
%             mkdir(strcat(climdir));
%         end 
%         climoutfile = strcat(climdir,regionname);
%         for monthij = 1:12
%             jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_relative_ssh_trend_', ...
%                 num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
% 
%             if (exist(jpgname , 'file') ~= 2)
%                 m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%                 hold on;
%                 m_pcolor(double(lon_rho)',lat_rho',squeeze(clim_trend(:,:,monthij)')-mean_clim_trend(monthij));
%         %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
%         %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
%                 shading(gca,m_pcolor_shading_method);
%                 m_gshhs_i('color',m_gshhs_line_color);
%                 m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
%                 titlename = strcat('SSH trend(rel), ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
%                     ', ', calendarname{monthij}(1:3), '), ','M=',num2str(round(mean_clim_trend(monthij),2)), 'mm/y');  %% + glacier contribution
% 
%                 title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%                 % set colorbar 
%                 h = colorbar;
%                 colormap(bwrmap);
%                 set(h,'fontsize',colorbar_fontsize);
% %                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
%                 caxis(reltrendlev);
% 
%                 % set grid
%                 m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%                 set(gcf, 'PaperUnits', 'points');
%                 set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%                 set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%                 saveas(gcf,jpgname,'tif');
% 
%                 disp(' ')
%                 disp(['clim_', num2str(monthij), '_', 'relative_ssh_trend', ' plot is created.'])
%                 disp(' ')
%                 disp([' File path is : ',jpgname])
%                 disp(' ')
% 
%                 hold off
%                 close all;
%             end
%         end
% % end-------------------- climatological relative trend plot



% start-------------------- absolute trend plot
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'-mean_trend_filtered));
    %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
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
 
% start-------------------- earlier decadal SSH plot
        figdir=[figrawdir,'CLIM\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);

        var='SSH';
        pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', var,'_',num2str(min(inputyear1),'%04i'), ...
            '_',num2str(max(inputyear1),'%04i'), '.tif'); %% ~_year_month.jpg
%         if (exist(pngname , 'file') ~= 2)     
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m']

            run(param_script);
            for yearij=1:length(inputyear1)
                tempyear=inputyear1(yearij);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');

                    tind=(tempyear-min(inputyear1))*12+tempmonth;

%                     if (exist('lon_rho' , 'var') ~= 1)
%                         lon_rho=ncread(filename, 'lon');
%                         lat_rho=ncread(filename, 'lat');
%                         [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
%                         cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                         cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                     end
%                     data_info = ncinfo(filename, varname); 

                    data = comb_data(:,:,tind);

                    if (exist('mean_data' , 'var') ~= 1)
                        mean_data=zeros(size(data));
                    end
                    mean_data=mean_data + (data / (length(inputyear1) * length(inputmonth)));
                end
            end
            
            mean_SSH_former =mean(mean(mean_data,'omitnan'),'omitnan');
            mean_data1 = mean_data - mean_SSH_former*2/3;
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            hold on;

            m_pcolor(lon_rho',lat_rho',mean_data1');
            shading(gca,m_pcolor_shading_method);   

            m_gshhs_i('color',m_gshhs_line_color)  
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
            titlename = strcat(var, ' mean, ',testname,',(',num2str(min(inputyear1),'%04i'),'-',num2str(max(inputyear1),'%04i'),') ');  %% + glacier contribution
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
            caxis(colorbar_lev);


            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
            saveas(gcf,pngname,'tif');
            close all;
            clear mean_data     
%         end
% end-------------------- eariler decadal SSH plot

% start-------------------- later decadal SSH plot
        figdir=[figrawdir,'CLIM\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);

        var='SSH';
        pngname=strcat(outfile, '_', testname,'_',regionname, '_clim_', var,'_',num2str(min(inputyear2),'%04i'), ...
            '_',num2str(max(inputyear2),'%04i'), '.tif'); %% ~_year_month.jpg
%         if (exist(pngname , 'file') ~= 2)     
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m']

            run(param_script);
            for yearij=1:length(inputyear2)
                tempyear=inputyear2(yearij);
                yearstr=num2str(tempyear, '%04i');
                for monthij=1:length(inputmonth)
                    tempmonth=inputmonth(monthij);
                    monthstr=num2str(tempmonth, '%02i');

                    tind=(inputyear2(yearij)-min(inputyear1))*12+tempmonth;

%                     if (exist('lon_rho' , 'var') ~= 1)
%                         lon_rho=ncread(filename, 'lon');
%                         lat_rho=ncread(filename, 'lat');
%                         [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
%                         cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                         cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                     end
%                     data_info = ncinfo(filename, varname); 

                    data = comb_data(:,:,tind);

                    if (exist('mean_data' , 'var') ~= 1)
                        mean_data=zeros(size(data));
                    end
                    mean_data=mean_data + (data / (length(inputyear2) * length(inputmonth)));
                end
            end
            
            mean_data2 = mean_data - mean_SSH_former*2/3;
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            hold on;

            m_pcolor(lon_rho',lat_rho',mean_data2');
            shading(gca,m_pcolor_shading_method);   

            m_gshhs_i('color',m_gshhs_line_color)  
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
            titlename = strcat(var, ' mean, ',testname,',(',num2str(min(inputyear2),'%04i'),'-',num2str(max(inputyear2),'%04i'),') ');  %% + glacier contribution
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
            caxis(colorbar_lev);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
            saveas(gcf,pngname,'tif');
            close all;
            clear mean_data     
%         end
% end-------------------- later decadal SSH plot




% % start-------------------- climatological absolute trend plot
%         climdir = [figdir,'\CLIM\'];
%         if (exist(strcat(climdir) , 'dir') ~= 7)
%             mkdir(strcat(climdir));
%         end 
%         climoutfile = strcat(climdir,regionname);
%         for monthij = 1:12
%             jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
%                 num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
%             if (exist(jpgname , 'file') ~= 2)
%                 m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%                 hold on;
%                 m_pcolor(double(lon_rho)',lat_rho',squeeze(clim_trend(:,:,monthij)'));
%         %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
%         %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
%                 shading(gca,m_pcolor_shading_method);
%                 m_gshhs_i('color',m_gshhs_line_color);
%                 m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
%                 titlename = strcat('SSH trend(abs), ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
%                     ', ', calendarname{monthij}(1:3), '), ','M=',num2str(round(mean_clim_trend(monthij),2)),'mm/y');  %% + glacier contribution
% 
%                 title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%                 % set colorbar 
%                 h = colorbar;
%                 colormap(wrmap);
%                 set(h,'fontsize',colorbar_fontsize);
% %                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
%                 caxis(abstrendlev);
% 
%                 % set grid
%                 m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%                 set(gcf, 'PaperUnits', 'points');
%                 set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%                 set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%                 saveas(gcf,jpgname,'tif');
% 
%                 disp(' ')
%                 disp(['clim_', num2str(monthij), '_', 'relative_ssh_trend', ' plot is created.'])
%                 disp(' ')
%                 disp([' File path is : ',jpgname])
%                 disp(' ')
% 
%                 hold off
%                 close all;
%             end
%         end
%         close all;
% % end-------------------- climatological absolute trend plot



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
        

% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
% end-------------------- make timedata for time series  

% start-------------------- msl time series (seasonal filtered)

        figdir=[figrawdir,'Trend\SSH\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        
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
            set(gca,'FontSize',20);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
        end
        close all;
% end-------------------- msl time series (seasonal filtered)

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
            set(gca,'FontSize',20);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
% end-------------------- msl time series

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
                    ', ', calendarname{monthij}(1:3), '), ',num2str(round(mean_clim_trend(monthij),2)), ' mm/y'])
%                 datetick('x','yyyy','keepticks')
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
        if (exist(jpgname , 'file') ~= 2)
            mslplot=plot(1:12,mean_clim_trend,'k')
            hold on
            xlabel('month')
            ylabel('trend (mm/yr)')
            title([regionname, ', climatological trend(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
            axis tight;
            ylim(trendplotlev)
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',20);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
% end-------------------- climatological msl trend

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

    end
end