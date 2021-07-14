close all; clear all;  clc;
warning off;
% all_region2 ={'AKP2'}
% all_region2 ={'SS', 'YS'}

all_region2 ={'AKP4'};
% all_region2 ={'BOH'};

% all_region2 ={'NWP', 'AKP2', 'NES', 'SES', 'YS'}
for regionind2=1:length(all_region2)
    close all;
    clearvars '*' -except regionind2 all_region2
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
    sshlev =[-0.3 0.3];
    sshdifflev = [0 20];

    % for snu_desktop
%     'cmems'='test49'   % % need to change
    inputyear = [1994:2014]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

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
            msl_abs_lev=[-500 1500];
            msllev= [-1000 1000];
        case('NWP2') %% North western Pacific
            lonlat = [115, 145, 25, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('ES') %% East Sea
            refpolygon=espolygon;
        case('NES') %% East Sea
            refpolygon=nespolygon;
        case('SES') %% East Sea
            refpolygon=sespolygon;
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        case('AKP2')
            refpolygon=akp2polygon;
        case('AKP4') %% Around Korea Peninsula
            refpolygon=akp4polygon;
        case('BOH') %% Bohai bay
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

%     load(['E:\Data\Model\ROMS\nwp_1_20\SODA\',regionname, ...
%         'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
    load(['E:\Data\Observation\CMEMS\',regionname, ...
        'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
    
    testname='cmems';

    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\CMEMS\',regionname,'\'); % % where figure files will be saved
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_', regionname, '.m']
%         filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\SODA\'); % % where data files are
        filedir = strcat('E:\Data\Observation\CMEMS\'); % % where data files are
%         cmemsdir='E:\Data\Model\ROMS\nwp_1_20\SODA\';
        cmemsdir='E:\Data\Observation\CMEMS\';
    elseif (strcmp(system_name,'GLNXA64'))
    end
    
    run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
    wrmap = bwrmap(51:100,:);
    
    varname ='zeta'
%     var='SSH'
    run(param_script);

    figdir=[figrawdir,'Trend\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
%     outfile = strcat(figdir,regionname);
    outfile = strcat(figdir);
    
    
    m_err=mean(comb_cmems_err,3);
    mean_1y=comb_cmems_data(:,:,1:12);
    mean_1y=mean(mean_1y(:), 'omitnan');
    mean_ly=comb_cmems_data(:,:,end-12:end);
    mean_ly=mean(mean_ly(:), 'omitnan');

    for y_ind=1:length(inputyear)
        yearly_sealevel(:,:,y_ind)=mean(comb_cmems_data(:,:,(y_ind-1)*12+1:y_ind*12),3,'omitnan');
    end
    clear var
    for xind=1:size(yearly_sealevel,1)
        for yind=1:size(yearly_sealevel,2)                
            var_yearly_sealevel(xind,yind)=var(yearly_sealevel(xind,yind,:));
        end
    end
    diff_all_y=mean(sqrt(var_yearly_sealevel(:)).*100, 'omitnan');
    m_err(m_err.*100<diff_all_y.*2)=NaN;
    err_mask = mask_cmems;
    err_mask(isfinite(m_err))=NaN;
    
    
% % % %     % start-------------------- 1 year msl plot
% % % %         jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_1y_ssh_',num2str(min(inputyear),'%04i'), ...
% % % %             '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
% % % %         if (exist(jpgname , 'file') ~= 2)
% % % %             msl_1y=squeeze(mean(comb_cmems_adt(:,:,1:12),3,'omitnan')).*mask_cmems;
% % % %             m_msl_1y=mean(msl_1y(:),'omitnan');
% % % %             
% % % %             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% % % %             hold on;
% % % %             m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(msl_1y')-m_msl_1y);
% % % %     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
% % % %     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
% % % %             shading(gca,m_pcolor_shading_method);
% % % %             m_gshhs_i('color',m_gshhs_line_color);
% % % %             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % %     %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
% % % %             titlename = strcat('SSH, ','CMEMS',',(',num2str(min(inputyear),'%04i'),')');  %% + glacier contribution
% % % %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % % 
% % % %             % set colorbar 
% % % %             h = colorbar;
% % % %             colormap(jet);
% % % %             set(h,'fontsize',colorbar_fontsize);
% % % %             title(h,'(m)','fontsize',colorbar_title_fontsize);
% % % %             caxis(sshlev);
% % % % 
% % % %             % set grid
% % % %             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % % % 
% % % %             set(gcf, 'PaperUnits', 'points');
% % % %             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % % %             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % % 
% % % %             saveas(gcf,jpgname,'tif');
% % % % 
% % % %             disp(' ')
% % % %             disp(['1y_', 'msl', ' plot is created.'])
% % % %             disp(' ')
% % % %             disp([' File path is : ',jpgname])
% % % %             disp(' ')
% % % % 
% % % %             hold off
% % % %             close all;
% % % %         end
% % % % % end-------------------- 1 year msl plot
% % % % 
% % % % % start-------------------- last year & diff msl plot
% % % %         jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_last_y_ssh_',num2str(min(inputyear),'%04i'), ...
% % % %             '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
% % % % %         varind_s=(max(inputyear)-min(inputyear))*12 +1;
% % % % %         varind_e=(max(inputyear)-min(inputyear))*12 +12;
% % % %         if (exist(jpgname , 'file') ~= 2)
% % % % % %             for varind=1:12
% % % % % %                 msl_month_1y(:,:,varind)=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
% % % % % %             end
% % % % % %             for varind=varind_s:varind_e
% % % % % %                 msl_month(:,:,varind-varind_s+1)=ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]);
% % % % % %             end
% % % %             msl_1y=squeeze(mean(comb_cmems_adt(:,:,1:12),3,'omitnan')).*mask_cmems;
% % % %             msl_ly=squeeze(mean(comb_cmems_adt(:,:,end-11:end),3,'omitnan')).*mask_cmems;
% % % %             m_msl_1y=mean(mean(msl_1y,'omitnan'),'omitnan');
% % % %             m_msl_ly=mean(mean(msl_ly,'omitnan'),'omitnan');
% % % % 
% % % %             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% % % %             hold on;
% % % %             m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(msl_ly')-m_msl_1y);
% % % %     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
% % % %     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
% % % %             shading(gca,m_pcolor_shading_method);
% % % %             m_gshhs_i('color',m_gshhs_line_color);
% % % %             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % %     %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
% % % %             titlename = strcat('SSH, ','cmems',',(',num2str(max(inputyear),'%04i'),')');  %% + glacier contribution
% % % %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % % 
% % % %             % set colorbar 
% % % %             h = colorbar;
% % % %             colormap(jet);
% % % %             set(h,'fontsize',colorbar_fontsize);
% % % %             title(h,'(m)','fontsize',colorbar_title_fontsize);
% % % %             caxis(colorbar_lev);
% % % % 
% % % %             % set grid
% % % %             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % % % 
% % % %             set(gcf, 'PaperUnits', 'points');
% % % %             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % % %             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % % 
% % % %             saveas(gcf,jpgname,'tif');
% % % % 
% % % %             disp(' ')
% % % %             disp(['1y_', 'msl', ' plot is created.'])
% % % %             disp(' ')
% % % %             disp([' File path is : ',jpgname])
% % % %             disp(' ')
% % % % 
% % % %             hold off
% % % %             close all;
% % % %   
% % % %             jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_diff_1st_last_y_ssh_',num2str(min(inputyear),'%04i'), ...
% % % %                         '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
% % % % 
% % % %             diff_cm=(msl_ly - msl_1y)*100.0;
% % % %             
% % % %             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% % % %             hold on;
% % % %             m_pcolor(double(cmems_lon2)',cmems_lat2',diff_cm');
% % % %     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
% % % %     %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
% % % %             shading(gca,m_pcolor_shading_method);
% % % %             m_gshhs_i('color',m_gshhs_line_color);
% % % %             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % %     %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
% % % %             titlename = strcat('dSSH, ','cmems',' (',num2str(max(inputyear),'%04i'),') rise =', num2str(round(mean(diff_cm(:),'omitnan'))), 'cm');  %% + glacier contribution
% % % % 
% % % % %             mean(mean(diff_cm,'omitnan'),'omitnan')
% % % % %             mean(mean(trend_filtered,'omitnan'),'omitnan')
% % % %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % % 
% % % %             % set colorbar 
% % % %             h = colorbar;
% % % %             colormap(jet);
% % % %             set(h,'fontsize',colorbar_fontsize);
% % % %             title(h,'(cm)','fontsize',colorbar_title_fontsize);
% % % %             caxis(sshdifflev);
% % % % 
% % % %             % set grid
% % % %             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % % % 
% % % %             set(gcf, 'PaperUnits', 'points');
% % % %             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % % %             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % % 
% % % %             saveas(gcf,jpgname,'tif');
% % % % 
% % % %             disp(' ')
% % % %             disp(['1y_', 'msl', ' plot is created.'])
% % % %             disp(' ')
% % % %             disp([' File path is : ',jpgname])
% % % %             disp(' ')
% % % % 
% % % %             hold off
% % % %             close all;
% % % %         end
% % % % 
% % % % % end-------------------- last year & diff msl plot

% start-------------------- absolute trend plot
        jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_absolute_ssh_trend_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(cmems_trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ','cmems', ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(cmems_trend_filtered(:),'omitnan'),2)), ' mm/y');  %% + glacier contribution

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

% start-------------------- absolute trend plot_set colorbar 
        jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_absolute_corrected_ssh_trend_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(cmems_trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ','cmems', ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(cmems_trend_filtered(:),'omitnan'),2)), ' mm/y');  %% + glacier contribution

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
% end-------------------- absolute trend plot_set colorbar 

    
% start-------------------- absolute trend plot bwrmap
        jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_absolute_ssh_trend_bwrmap_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(cmems_trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ','cmems', ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(cmems_trend_filtered(:),'omitnan'),2)), ' mm/y');  %% + glacier contribution

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
% end-------------------- absolute trend plot bwrmap



% start-------------------- relative trend plot bwrmap
        jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_relative_ssh_trend_bwrmap_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(cmems_trend_filtered(:,:)'- mean(cmems_trend_filtered(:),'omitnan')));
            shading(gca,m_pcolor_shading_method)
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH trend(abs), ','cmems', ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(cmems_trend_filtered(:),'omitnan'),2)), ' mm/y');  %% + glacier contribution

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
% end-------------------- relative trend plot bwrmap




% start-------------------- err plot bwrmap
        jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_err_wrmap_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(mean(comb_cmems_err,3)'*100.0));
            shading(gca,m_pcolor_shading_method)
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
            titlename = strcat('SSH err(abs), ','cmems', ',(',num2str(min(inputyear),'%04i'),'-', ...
                num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(comb_cmems_err(:)*100.0,'omitnan'),2)), ' cm');  %% + glacier contribution

            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(wrmap);
            set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis([-4 4]);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', 'err', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
%         end
% end-------------------- err plot bwrmap





% %      % cmems trend plot
% %      jpgname=strcat(outfile,regionname, '_cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
% %      if (exist(jpgname , 'file') ~= 2)
% % 
% %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% %         hold on;
% %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
% %         shading(gca,m_pcolor_shading_method);
% %         m_gshhs_i('color',m_gshhs_line_color);
% %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % %         titlename = strcat(regionname,', SODA SSH trend','(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'Mean=',num2str(round(mean_cmems_trend_filtered,2)),'mm/y');
% %         titlename = strcat(regionname,', CMEMS SLA trend','(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'Mean=',num2str(round(mean_cmems_trend_filtered,2)));
% %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % 
% %         % set colorbar 
% %         h = colorbar;
% %         colormap(bwrmap);
% %         set(h,'fontsize',colorbar_fontsize);
% %         title(h,'mm/y','fontsize',colorbar_title_fontsize);
% %         caxis(trendlev);
% % 
% %         % set grid
% %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % 
% %         set(gcf, 'PaperUnits', 'points');
% %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % 
% %         saveas(gcf,jpgname,'tif');
% % 
% %         disp(' ')
% %         disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
% %         disp(' ')
% %         disp([' File path is : ',jpgname])
% %         disp(' ')
% % 
% %         hold off
% %         close all;
% %      end
% %      
     figdir2=[figrawdir,'MSL\'];
        if (exist(strcat(figdir2) , 'dir') ~= 7)
            mkdir(strcat(figdir2));
        end 
    %     outfile = strcat(figdir,regionname);
        outfile2 = strcat(figdir2);
% %         
% %         isize = size(comb_cmems_data_filtered,1);
% %         jsize = size(comb_cmems_data_filtered,2);
% %         lsize = size(comb_cmems_data_filtered,3);
% %         comb_monthly_cmems_data=reshape(comb_cmems_adt,[isize, jsize, 12, lsize/12]);
% %         comb_yearly_cmems_data=squeeze(mean(comb_monthly_cmems_data,3,'omitnan'));  
        
% %          % % % %             % msl pcolor (anomaly)  
% %         jpgname=strcat(outfile2, '_cmems_',regionname, '_msl_ano_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg 
% %         if (exist(jpgname , 'file') ~= 2)
% %             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% %             hold on;
% %             msl_alltime=squeeze(mean(comb_yearly_cmems_data,3,'omitnan'))'.*1000;
% %             msl_alltime=msl_alltime-mean(mean(msl_alltime,1,'omitnan'),2,'omitnan');
% %             m_pcolor(double(cmems_lon),cmems_lat',msl_alltime);
% %     %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
% %             shading(gca,m_pcolor_shading_method);
% %             m_gshhs_i('color',m_gshhs_line_color);
% %             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %             titlename = strcat('MSL',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
% %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % 
% %             % set colorbar 
% %             h = colorbar;
% %             colormap(bwrmap);
% %             set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm','fontsize',colorbar_title_fontsize);
% %             caxis(msllev);
% % 
% %             % set grid
% %             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % 
% %             set(gcf, 'PaperUnits', 'points');
% %             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% %             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % 
% %             saveas(gcf,jpgname,'tif');
% % 
% %             disp(' ')
% %             disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
% %             disp(' ')
% %             disp([' File path is : ',jpgname])
% %             disp(' ')
% % 
% %             hold off
% %             close all;
% %         end  
        
% %         jpgname=strcat(outfile2, '_cmems_',regionname, '_msl_abs_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
% %         if (exist(jpgname , 'file') ~= 2)        
% % % % % %             % msl pcolor (absolute, total year mean)
% %             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% %             hold on;
% %             msl_alltime=squeeze(mean(comb_yearly_cmems_data,3,'omitnan'))'.*1000;
% %             m_msl=mean(mean(msl_alltime,1,'omitnan'),2,'omitnan');
% %             m_pcolor(double(cmems_lon),cmems_lat',msl_alltime);
% %     %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
% %             shading(gca,m_pcolor_shading_method);
% %             m_gshhs_i('color',m_gshhs_line_color);
% %             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %             titlename = strcat('MSL',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ','M=',num2str(round(m_msl,1)));
% %             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % 
% %             % set colorbar 
% %             h = colorbar;
% %             colormap(jet);
% %             set(h,'fontsize',colorbar_fontsize);
% %             title(h,'mm','fontsize',colorbar_title_fontsize);
% % %             caxis(msllev);
% % 
% %             % set grid
% %             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % 
% %             set(gcf, 'PaperUnits', 'points');
% %             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% %             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % 
% %             saveas(gcf,jpgname,'tif');
% % 
% %             disp(' ')
% %             disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
% %             disp(' ')
% %             disp([' File path is : ',jpgname])
% %             disp(' ')
% % 
% %             hold off
% %             close all;          
% %         end
        
        
          for yearij = 1:length(inputyear)
          % % % %             % msl pcolor (absolute)
            
            tempyear= inputyear(yearij);
            
            
% %             jpgname=strcat(outfile2, '_cmems_',regionname, '_msl_abs_',num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
% %             if (exist(jpgname , 'file') ~= 2)               
% %                 m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% %                 hold on;
% %                 msl_alltime=squeeze(comb_yearly_cmems_data(:,:,yearij))'.*1000;
% %                 m_msl=mean(mean(msl_alltime,1,'omitnan'),2,'omitnan');
% %                 m_pcolor(double(cmems_lon),cmems_lat',msl_alltime);
% %         %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
% %                 shading(gca,m_pcolor_shading_method);
% %                 m_gshhs_i('color',m_gshhs_line_color);
% %                 m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %                 titlename = strcat('MSL',',(',num2str(tempyear,'%04i'),') ','M=',num2str(round(m_msl,1)));
% %                 title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % 
% %                 % set colorbar 
% %                 h = colorbar;
% %                 colormap(jet);
% %                 set(h,'fontsize',colorbar_fontsize);
% %                 title(h,'mm','fontsize',colorbar_title_fontsize);
% %                 caxis(msl_abs_lev);
% % 
% %                 % set grid
% %                 m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % 
% %                 set(gcf, 'PaperUnits', 'points');
% %                 set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% %                 set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % 
% %                 saveas(gcf,jpgname,'tif');
% % 
% %                 disp(' ')
% %                 disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
% %                 disp(' ')
% %                 disp([' File path is : ',jpgname])
% %                 disp(' ')
% % 
% %                 hold off
% %                 close all;   
% %             end
            
% %             if (yearij~=length(inputyear))    
% %                 jpgname=strcat(outfile2, '_cmems_',regionname, '_msl_diff_',num2str(inputyear(yearij+1),'%04i'),'_',num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
% %                 if (exist(jpgname , 'file') ~= 2)
% %                      % % % %             % msl pcolor (diff)
% %                     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% %                     hold on;
% %                     msl_diff=squeeze(comb_yearly_cmems_data(:,:,yearij+1))'.*1000 - squeeze(comb_yearly_cmems_data(:,:,yearij))'.*1000;
% %                     m_msl=mean(mean(msl_diff,1,'omitnan'),2,'omitnan');
% %                     m_pcolor(double(cmems_lon),cmems_lat',msl_diff);
% %                     shading(gca,m_pcolor_shading_method);
% %                     m_gshhs_i('color',m_gshhs_line_color);
% %                     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %                     titlename = strcat('MSL',',(',num2str(inputyear(yearij+1),'%04i'),'-',num2str(tempyear,'%04i'),') ','M=',num2str(round(m_msl,1)));
% %                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % 
% %                     % set colorbar 
% %                     h = colorbar;
% %                     colormap(bwrmap);
% %                     set(h,'fontsize',colorbar_fontsize);
% %                     title(h,'mm','fontsize',colorbar_title_fontsize);
% %                     caxis(msl_diff_lev);
% % 
% %                     % set grid
% %                     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % 
% %                     set(gcf, 'PaperUnits', 'points');
% %                     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% %                     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % 
% %                     saveas(gcf,jpgname,'tif');
% % 
% %                     disp(' ')
% %                     disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
% %                     disp(' ')
% %                     disp([' File path is : ',jpgname])
% %                     disp(' ')
% % 
% %                     hold off
% %                     close all;
% %                 end
% %             end
             
            
            
            figdir3=[figrawdir,'MSL\monthly\'];
            outfile3 = strcat(figdir3,regionname);
            if (exist(strcat(figdir3) , 'dir') ~= 7)
                mkdir(strcat(figdir3));
            end
            
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
% %                 jpgname=strcat(outfile3, '_cmems_',regionname, '_msl_abs_',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'), '.tif'); %% ~_year_month.jpg
% %                 if (exist(jpgname , 'file') ~= 2)
% %                     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% %                     hold on;
% %                     msl_alltime=squeeze(comb_monthly_cmems_data(:,:,monthij,yearij))'.*1000;
% %                     m_msl=mean(mean(msl_alltime,1,'omitnan'),2,'omitnan');
% %                     m_pcolor(double(cmems_lon),cmems_lat',msl_alltime);
% %             %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
% %                     shading(gca,m_pcolor_shading_method);
% %                     m_gshhs_i('color',m_gshhs_line_color);
% %                     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %                     titlename = strcat('MSL',',(',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'),') ','M=',num2str(round(m_msl,1)));
% %                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % 
% %                     % set colorbar 
% %                     h = colorbar;
% %                     colormap(jet);
% %                     set(h,'fontsize',colorbar_fontsize);
% %                     title(h,'mm','fontsize',colorbar_title_fontsize);
% %                     caxis(msl_abs_lev);
% % 
% %                     % set grid
% %                     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % 
% %                     set(gcf, 'PaperUnits', 'points');
% %                     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% %                     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % 
% %                     saveas(gcf,jpgname,'tif');
% % 
% %                     disp(' ')
% %                     disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
% %                     disp(' ')
% %                     disp([' File path is : ',jpgname])
% %                     disp(' ')
% % 
% %                     hold off
% %                     close all; 
% %                 end


            end
          end
        
        
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
                m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(clim_cmems_trend_divided(:,:,monthij)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                    ', ', calendarname{monthij}(1:3), '), ','M=',num2str(round(mean_clim_cmems_trend_divided(monthij),2)),'mm/y');  %% + glacier contribution

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
% end-------------------- climatological absolute trend plot
          
        
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
            xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',num2str(tempyear)]);
        end
    end

    cmems_msl=squeeze(mean(mean(comb_cmems_data_filtered,1,'omitnan'),2,'omitnan'));
    cmems_msl=cmems_msl-mean(cmems_msl);  
    p2=polyfit(1:length(cmems_msl),cmems_msl',1);
    p2=p2(1)*1000.0*12.0;
            
    p=polyfit(xData,cmems_msl',1);
    cmems_msl2=xData*p(1)+p(2);
    cmems_mslplot=plot(xData,cmems_msl,'k')
    hold on
    cmems_mslplot2=plot(xData,cmems_msl2,'Color','r')
    jpgname=strcat(outfile2, regionname,'_cmems_msl_', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('Year')
    ylabel('Mean SSH (m)')
    title([regionname ', Mean CMEMS SLA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
    datetick('x','yyyy','keepticks')
    axis tight;
    ylim(meanplotlev);
    set(cmems_mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
    
    cmems_adt=squeeze(mean(mean(comb_cmems_adt,1,'omitnan'),2,'omitnan'));
    cmems_adt=cmems_adt-mean(cmems_adt);     
    p2=polyfit(1:length(cmems_adt),cmems_adt',1);
    p2=p2(1)*1000.0*12.0;
    
    p=polyfit(xData,cmems_adt',1);
    cmems_adt2=xData*p(1)+p(2);
    cmems_adtplot=plot(xData,cmems_adt,'k')
    hold on
    cmems_adtplot2=plot(xData,cmems_adt2,'Color','r')
    jpgname=strcat(outfile2, regionname,'_cmems_adt_', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('Year')
    ylabel('Mean SSH (m)')
    title([regionname ', Mean CMEMS SLA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
    datetick('x','yyyy','keepticks')
    axis tight;
    ylim(meanplotlev);
    set(cmems_adtplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
    
    
%     
%     for lati=1:size(cmems_lat2,2)
%         weight_msl(1:size(cmems_lon2,1),lati)=m_lldist(cmems_lon2(1:2,lati), cmems_lat2(1:2,lati),1);
%     end
%     weight_msl_masked=weight_msl;
%     weight_msl_masked(isnan(comb_cmems_data_filtered(:,:,1))==1)=NaN;
%     weighted_sum_sl=sum(sum(comb_cmems_data_filtered,1,'omitnan'),2,'omitnan')
%     
%     cmems_msl=squeeze(mean(mean(comb_cmems_data_filtered,1,'omitnan'),2,'omitnan'));
%     cmems_msl=cmems_msl-mean(cmems_msl);        
%     p=polyfit(xData,cmems_msl',1);
%     cmems_msl2=xData*p(1)+p(2);
%     cmems_mslplot=plot(xData,cmems_msl,'k')
%     hold on
%     cmems_mslplot2=plot(xData,cmems_msl2,'Color','r')
%     jpgname=strcat(outfile2, regionname,'_cmems_msl_', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%     xlabel('Year')
%     ylabel('Mean SSH (m)')
%     title([regionname ', Mean CMEMS SLA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_cmems_trend_filtered,2)), ' mm/y'])
%     datetick('x','yyyy','keepticks')
%     axis tight;
%     ylim(meanplotlev);
%     set(cmems_mslplot,'LineWidth',2);
%     set(gca,'FontSize',15);
%     grid on
%     hold off
%     saveas(gcf,jpgname,'jpg');
%     grid off
    
    
    % start-------------------- msl abs time series

    jpgname=strcat(outfile2, 'cmems_',regionname, '_msl_abs_',num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%     var='SSH';
    param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_', regionname, '.m']
    run(param_script);

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    hold on;
    
    madt = mean(comb_cmems_adt,3, 'omitnan');
    m_pcolor(cmems_lon,cmems_lat,madt'-mean(madt(:),'omitnan')+0.6);
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
    title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
    caxis(colorbar_lev);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    saveas(gcf,jpgname,'tif');
    close all;
    % end-------------------- msl abs time series

    
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
                        tempdata=squeeze(comb_cmems_data(:,:,varind));
                        msl(varind)=mean(tempdata(:),'omitnan');
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
                    ', ', calendarname{monthij}(1:3), '), ',num2str(round(mean_clim_cmems_trend_divided(monthij),2)), ' mm/y'])
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
        if abs(mean_clim_cmems_trend_divided(1))<0.01
            mean_clim_trend=mean_clim_trend*1000.0  %% m -> mm
        end
        if (exist(jpgname , 'file') ~= 2)
            mslplot=plot(1:12,squeeze(mean_clim_cmems_trend_divided),'k')
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
    
end


% start-------------------- absolute err > std plot


% pcolor(judge_err)
% % pcolor(m_err')
% shading flat
% colorbar
% colormap(jet)
% % caxis([-0.15 0.15])

    jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_err_over_std_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    if (exist(jpgname , 'file') ~= 2)
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
%         m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(cmems_trend_filtered(:,:)'));
        m_err=mean(comb_cmems_err,3);
        m_var=mean(sqrt(cmems_sla_var(:)),'omitnan').*100;
        judge_err=sqrt(cmems_sla_var')-m_err';
%         judge_err=cmems_sla_var'-m_err';
        judge_err(judge_err>0)=NaN;
        m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(abs(judge_err)*100.0));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
        titlename = strcat('SSH err(>std), ','cmems', ',(',num2str(min(inputyear),'%04i'),'-', ...
            num2str(max(inputyear),'%04i'),')');  %% + glacier contribution

        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
%         colormap(wrmap);
        colormap([0 0 0]);

        set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%         caxis(abstrendlev);

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
% end-------------------- absolute err > std plot



% start-------------------- absolute err > mean plot


% pcolor(judge_err)
% % pcolor(m_err')
% shading flat
% colorbar
% colormap(jet)
% % caxis([-0.15 0.15])

    jpgname=strcat(outfile, '_', 'cmems','_',regionname, '_err_over_mean_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%     if (exist(jpgname , 'file') ~= 2)
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
%         m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(cmems_trend_filtered(:,:)'));
        m_err=mean(comb_cmems_err,3);
        mean_1y=comb_cmems_data(:,:,1:12);
        mean_1y=mean(mean_1y(:), 'omitnan');
        mean_ly=comb_cmems_data(:,:,end-12:end);
        mean_ly=mean(mean_ly(:), 'omitnan');
%         diff_all_y=100.0*(mean_ly-mean_1y);
%         diff_all_y=mean(m_err(:),'omitnan').*100;
        
        for y_ind=1:length(inputyear)
            yearly_sealevel(:,:,y_ind)=mean(comb_cmems_data(:,:,(y_ind-1)*12+1:y_ind*12),3,'omitnan');
        end
        clear var
        for xind=1:size(yearly_sealevel,1)
            for yind=1:size(yearly_sealevel,2)                
                var_yearly_sealevel(xind,yind)=var(yearly_sealevel(xind,yind,:));
            end
        end
        diff_all_y=mean(sqrt(var_yearly_sealevel(:)).*100, 'omitnan');

%         judge_err=sqrt(cmems_sla_var')-m_err';
%         judge_err=cmems_sla_var'-m_err';
        m_err(m_err.*100<diff_all_y.*2)=NaN;
        m_pcolor(double(cmems_lon2)',cmems_lat2',squeeze(abs(m_err').*100.0));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
        titlename = strcat('SSH err(>2std-yearly), ','cmems', ',(',num2str(min(inputyear),'%04i'),'-', ...
            num2str(max(inputyear),'%04i'),')');  %% + glacier contribution

        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
%         colormap(wrmap);
        colormap([0.4 0 0]);

        set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%         caxis(abstrendlev);

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
        err_mask = mask_cmems;
        err_mask(isfinite(m_err))=NaN;
%     end
% end-------------------- absolute err > mean plot


% m_lldist([115,164], [52,52], 1)
