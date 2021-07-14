close all; clear all;  clc;
warning off;
all_region2 ={'NWP','ES', 'SS', 'YS', 'ECS'}
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'SS', 'YS'}

% all_region2 ={'ECS'};

% all_region2 ={'NWP'}
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
    trendlev = [-5 5];  %% trend lev
    conlev  = 0:5:35;
    meanplotlev =[-0.2 0.2];

    % for snu_desktopd
    testname='test03'   % % need to change
    inputyear = [1980:2008]; % % put year which you want to plot [year year ...]
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

    load(['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
        'ssh_recon_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\nwp_1_10\',testname,'\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
        recondir='E:\Data\Observation\OISST\monthly\';
    elseif (strcmp(system_name,'GLNXA64'))
    end

    figdir=[figrawdir,'Trend\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);
        % relative trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'-mean_trend_filtered));
%         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
%         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
        titlename = strcat('SSH trend(rel), ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)));  %% + glacier contribution

        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
        colormap(bwrmap);
        set(h,'fontsize',colorbar_fontsize);
        title(h,'mm/y','fontsize',colorbar_title_fontsize);
        caxis(trendlev);

        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

        jpgname=strcat(outfile, '_', testname,'_',regionname, '_relative_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'tif');

        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'relative_ssh_trend', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
        
     % absolute trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
%         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'-mean_trend_filtered));
%         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
        m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
        titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)));  %% + glacier contribution

        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
        colormap(bwrmap);
        set(h,'fontsize',colorbar_fontsize);
        title(h,'mm/y','fontsize',colorbar_title_fontsize);
        caxis(trendlev);

        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'tif');

        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'absolute_ssh_trend', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
        
%         % corrected absolute trend plot
%         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%         hold on;
% %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'-mean_trend_filtered));
% %         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));  %% + glacier contribution
%         m_pcolor(double(recon_lon),recon_lat,squeeze(trend_filtered(:,:)'+0.92));
%         shading(gca,m_pcolor_shading_method);
%         m_gshhs_i('color',m_gshhs_line_color);
%         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
%         titlename = strcat('SSH trend(cor)',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','Mean=',num2str(round(mean_trend_filtered+0.92,2)),'mm/y');  %% + glacier contribution
% 
%         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%         % set colorbar 
%         h = colorbar;
%         colormap(bwrmap);
%         set(h,'fontsize',colorbar_fontsize);
%         title(h,'mm/y','fontsize',colorbar_title_fontsize);
%         caxis(trendlev);
% 
%         % set grid
%         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%         set(gcf, 'PaperUnits', 'points');
%         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%         jpgname=strcat(outfile, '_', testname,'_',regionname, '_corrected_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%         saveas(gcf,jpgname,'tif');
% 
%         disp(' ')
%         disp(['clim_', num2str(tempmonth), '_', 'corrected_ssh_trend', ' plot is created.'])
%         disp(' ')
%         disp([' File path is : ',jpgname])
%         disp(' ')
% 
%         hold off
%         close all;
        
        
     % recon rel trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'-mean_recon_trend_filtered));
%         m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        titlename = strcat('re SSH trend(rel)',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_recon_trend_filtered,2)));
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
        colormap(bwrmap);
        set(h,'fontsize',colorbar_fontsize);
        title(h,'mm/y','fontsize',colorbar_title_fontsize);
        caxis(trendlev);

        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

        jpgname=strcat(outfile, '_', testname,'_',regionname, '_relative_recon_ssh_trend',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'tif');

        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'relative_recon_ssh_trend', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
        
        % recon absolute trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'));
%         m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        titlename = strcat('re SSH trend(abs)',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_recon_trend_filtered,2)));
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
        colormap(bwrmap);
        set(h,'fontsize',colorbar_fontsize);
        title(h,'mm/y','fontsize',colorbar_title_fontsize);
        caxis(trendlev);

        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_recon_ssh_trend',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'tif');

        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'absolute_recon_ssh_trend', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
        
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
            xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
        end
    end

    msl=squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan'));
    msl=msl-mean(msl);    
    p=polyfit(xData,msl',1);
    msl2=xData*p(1)+p(2);
    mslplot=plot(xData,msl,'k')
    hold on
    mslplot2=plot(xData,msl2,'Color','r')
    jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('year')
    ylabel('Mean SSH (m)')
    title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend_filtered,2)), ' mm/y'])
    datetick('x','yymmm','keepticks')
    axis tight;
    ylim(meanplotlev)
    set(mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off

    recon_msl=squeeze(mean(mean(comb_recon_data_filtered,1,'omitnan'),2,'omitnan'));
    recon_msl=recon_msl-mean(recon_msl);        
    p=polyfit(xData,recon_msl',1);
    recon_msl2=xData*p(1)+p(2);
    recon_mslplot=plot(xData,recon_msl,'k')
    hold on
    recon_mslplot2=plot(xData,recon_msl2,'Color','r')
    jpgname=strcat(outfile, '_', testname, '_',regionname,'_recon_msl_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('year')
    ylabel('Mean SSH (m)')
    title([regionname ', Mean recon SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_recon_trend_filtered,2)), ' mm/y'])
    datetick('x','yymmm','keepticks')
    axis tight;
    ylim(meanplotlev)
    set(recon_mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
    close all;
end