close all; clear all;  clc;
warning off;
all_region2 ={'NWP','ES', 'SS', 'YS', 'ECS'}
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
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end

    shadlev = [0 35];
    rms_shadlev = [0 4];
    trendlev = [-3 3];  %% trend lev
    conlev  = 0:5:35;
    meanplotlev =[-0.1 0.1];

    % for snu_desktop
    testname='test03'   % % need to change
    inputyear = [1993:2009]; % % put year which you want to plot [year year ...]
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

    load(['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname,'ssh_aviso_trend_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end),'%04i'),'.mat']);

    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\nwp_1_10\',testname,'\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
        avisodir='E:\Data\Observation\OISST\monthly\';
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
        m_pcolor(double(aviso_lon),aviso_lat,squeeze(trend_filtered(:,:)'-mean_trend_filtered));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
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
        disp(['clim_', num2str(tempmonth), '_', 'realtive_ssh_trend', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
        
        
        % absolute trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        m_pcolor(double(aviso_lon),aviso_lat,squeeze(trend_filtered(:,:)'));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         titlename = strcat(regionname,', SSH trend(abs), ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
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
        
     % aviso relative trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        m_pcolor(double(aviso_lon),aviso_lat',squeeze(aviso_trend_filtered(:,:)'-mean_aviso_trend_filtered));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         titlename = strcat(regionname,', AVISO SSH trend(rel), ','Mean=',num2str(round(mean_aviso_trend_filtered,2)),'mm/y');
        titlename = strcat('AVISO SSH trend(rel)',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_aviso_trend_filtered,2)));
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

        jpgname=strcat(outfile, '_', testname,'_',regionname, '_relative_aviso_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'tif');

        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'relative_aviso_ssh_trend', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
        
        % aviso absolute trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        m_pcolor(double(aviso_lon),aviso_lat',squeeze(aviso_trend_filtered(:,:)'));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        titlename = strcat('AVISO SSH trend(abs)',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_aviso_trend_filtered,2)));
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

        jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_aviso_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'tif');

        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'absolute_aviso_ssh_trend', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
        
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
            xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',num2str(tempyear)]);
        end
    end

    msl=squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan'));
    msl=msl-mean(msl);    
    p=polyfit(xData,msl',1);
    msl2=xData*p(1)+p(2);
    mslplot=plot(xData,msl,'k')
    hold on
    mslplot2=plot(xData,msl2,'Color','r')
    jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl', '.jpg'); %% ~_year_month.jpg
    xlabel('month')
    ylabel('Mean SSH (m)')
    title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend_filtered,2)), ' mm/y'])
    ylim(meanplotlev)
    datetick('x','yyyy')
    set(mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
    close all;



    aviso_msl=squeeze(mean(mean(comb_aviso_data_filtered,1,'omitnan'),2,'omitnan'));
    aviso_msl=aviso_msl-mean(aviso_msl);        
    p=polyfit(xData,aviso_msl',1);
    aviso_msl2=xData*p(1)+p(2);
    aviso_mslplot=plot(xData,aviso_msl,'k')
    hold on
    aviso_mslplot2=plot(xData,aviso_msl2,'Color','r')
    jpgname=strcat(outfile, '_', testname, '_',regionname,'_aviso_msl', '.jpg'); %% ~_year_month.jpg
    xlabel('month')
    ylabel('Mean SSH (m)')
    title([regionname ', Mean AVISO SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_aviso_trend_filtered,2)), ' mm/y'])
    ylim(meanplotlev)
    datetick('x','yyyy')
    set(aviso_mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
    close all;
end