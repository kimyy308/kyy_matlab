close all; clear all;  clc;
warning off;
all_region2 ={'NWP','NWP2','ES', 'SS', 'YS', 'ECS'}
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
    testname='test42'   % % need to change
    inputyear = [1980:1994]; % % put year which you want to plot [year year ...]
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

    load(['E:\Data\Model\ROMS\nwp_1_20\SODA\',regionname, ...
        'soda_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
    
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\SODA\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\SODA\'); % % where data files are
        sodadir='E:\Data\Model\ROMS\nwp_1_20\SODA\';
    elseif (strcmp(system_name,'GLNXA64'))
    end

    figdir=[figrawdir,'Trend\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);

     % soda trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        m_pcolor(double(soda_lon),soda_lat',squeeze(soda_trend_filtered(:,:)'-mean_soda_trend_filtered));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        titlename = strcat(regionname,', SODA SSH trend, ','Mean=',num2str(round(mean_soda_trend_filtered,2)),'mm/y');
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

        jpgname=strcat(outfile, '_', testname,'_',regionname, '_soda_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'tif');

        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
            xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',num2str(tempyear)]);
        end
    end

    soda_msl=squeeze(mean(mean(comb_soda_data_filtered,1,'omitnan'),2,'omitnan'));
    soda_msl=soda_msl-mean(soda_msl);        
    p=polyfit(xData,soda_msl',1);
    soda_msl2=xData*p(1)+p(2);
    soda_mslplot=plot(xData,soda_msl,'k')
    hold on
    soda_mslplot2=plot(xData,soda_msl2,'Color','r')
    jpgname=strcat(outfile, '_', testname, '_',regionname,'_soda_msl', '.jpg'); %% ~_year_month.jpg
    xlabel('month')
    ylabel('Mean SSH (m)')
    title([regionname ', Mean SODA SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_soda_trend_filtered,2)), ' mm/y'])
    ylim(meanplotlev)
    datetick('x','yyyy')
    set(soda_mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
end