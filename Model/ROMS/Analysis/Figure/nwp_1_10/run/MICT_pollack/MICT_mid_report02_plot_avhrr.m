close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region2 ={'EKB'};

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
trendlev = [-0.1 0.1];
conlev  = 0:5:35;
meanplotlev =[-2 2];


% for snu_desktop
testname='test06'   % % need to change
inputyear = [1989:2017]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

    varname ='temp'
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
    case('EKB') %% North western Pacific
        lonlat = [127, 131, 37, 42];  %% East Korea Bay
        refpolygon(1,1)=lonlat(1);
        refpolygon(2,1)=lonlat(2);
        refpolygon(1,2)=lonlat(3);
        refpolygon(2,2)=lonlat(4);
    case('AKP2')
        refpolygon=akp2polygon;
    otherwise
        ('?')
    end
    lonlat(1)=min(refpolygon(:,1));
    lonlat(2)=max(refpolygon(:,1));
    lonlat(3)=min(refpolygon(:,2));
    lonlat(4)=max(refpolygon(:,2));

    load(['E:\Data\Observation\OISST\monthly_kimyy\',regionname,'sst_trend_',num2str(inputyear(1)),'_',num2str(inputyear(end)),'.mat']);
    
    
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\AVHRR\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat('E:\Data\Observation\OISST\monthly\'); % % where data files are
        avhrrdir='E:\Data\Observation\OISST\monthly_kimyy\';
    elseif (strcmp(system_name,'GLNXA64'))
    end


    trendlev = [-0.1 0.1];

    figdir=[figrawdir,'Trend\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);

     % avhrr trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        [avhrr_lat2 avhrr_lon2]=meshgrid(avhrr_lat, avhrr_lon);
        m_pcolor(double(avhrr_lon2),avhrr_lat2,squeeze(avhrr_trend_filtered(:,:)));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        titlename = strcat(regionname,', AVHRR SST trend, ','Mean=',num2str(round(mean_avhrr_trend_filtered,3)),'^oC/y');
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    
        % set colorbar 
        h = colorbar;
        colormap(bwrmap);
        set(h,'fontsize',colorbar_fontsize);
        title(h,'^oC/y','fontsize',colorbar_title_fontsize);
        caxis(trendlev);
        
        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    
        tifname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.tif
        drawnow;
        saveas(gcf,tifname,'tif');
    
        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',tifname])
        disp(' ')
    
        hold off
    for k=1:12
        close all;
        % avhrr climatological trend plot
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        
        m_pcolor(double(avhrr_lon2),avhrr_lat2,squeeze(clim_avhrr_trend_divided(:,:,k)));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        
        titlename = strcat(regionname,', ',calendarname{k},', SST trend, ','Mean=',num2str(round(mean_clim_avhrr_trend_divided(k),3)),'^oC/y');
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    
        % set colorbar 
        h = colorbar;
        colormap(bwrmap);
        set(h,'fontsize',colorbar_fontsize);
        title(h,'^oC/y','fontsize',colorbar_title_fontsize);
        caxis(trendlev);
        
        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    
        tifname=strcat(outfile,regionname, '_avhrr_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'_',num2str(k,'%02i'),'.tif'); %% ~_year_month.tif
        drawnow;
        saveas(gcf,tifname,'tif');
    
        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',tifname])
        disp(' ')
    
        hold off
    end
        
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
            xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',num2str(tempyear)]);
        end
    end
    
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
%             clim_xData(i,month) = datenum([num2str(month,'%02i'),'-01-',num2str(tempyear)]);
            clim_xData(i,month) = datenum([num2str(1,'%02i'),'-01-',num2str(tempyear)]);
        end
    end
    
    avhrr_msst=squeeze(mean(mean(comb_avhrr_data_filtered,1,'omitnan'),2,'omitnan'));
    % avhrr_msst=avhrr_msst-mean(avhrr_msst);        
    p=polyfit(xData,avhrr_msst',1);
    avhrr_msst2=xData*p(1)+p(2);
    avhrr_msstplot=plot(xData,avhrr_msst,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
    hold on
    avhrr_msstplot2=plot(xData,avhrr_msst2,'Color','r')
    tifname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_msst', '.tif'); %% ~_year_month.tif
    xlabel('year')
    ylabel('Mean SSTA (^oC)')
    title([regionname,', Mean AVHRR SSTA(', ...
        num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'_', ...
        num2str(k,'%02i'),'), ',num2str(round(mean_avhrr_trend_filtered,3)), '^oC/y'])
    set(avhrr_msstplot,'LineWidth',2);
    datetick('x','yy')
    set(gca,'FontSize',15);
    axis tight
    ylim(meanplotlev);
    grid on
    grid minor
    hold off
    drawnow;
    saveas(gcf,tifname,'tif');
    grid off
    comb_clim_avhrr_data=reshape(comb_avhrr_data,[len_lon,len_lat,12,length(inputyear)]);
    
    for k=1:12
        close all
        clim_avhrr_msst=squeeze(mean(mean(comb_clim_avhrr_data(:,:,k,:),1,'omitnan'),2,'omitnan'));
        % clim_avhrr_msst=clim_avhrr_msst-mean(clim_avhrr_msst);        
        p=polyfit(clim_xData(:,k)',clim_avhrr_msst',1);
        clim_avhrr_msst2=clim_xData(:,k)'*p(1)+p(2);
        clim_avhrr_msstplot=plot(clim_xData(:,k)',clim_avhrr_msst','k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
        hold on
        clim_avhrr_msstplot2=plot(clim_xData(:,k)',clim_avhrr_msst2','Color','r')
        tifname=strcat(outfile, '_', testname,'_',regionname, '_clim_avhrr_msst, ',num2str(k,'%02i'), '.tif'); %% ~_year_month.tif
        xlabel('year')
        ylabel('Mean SST (^oC)')
        title([regionname,', ',calendarname{k},', Mean SST(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_clim_avhrr_trend_divided(:,:,k),3)), '^oC/y'])
        set(clim_avhrr_msstplot,'LineWidth',2);
%         datetick('x','yy','keepticks')
        datetick('x','yy')
        set(gca,'FontSize',15);
        axis tight
%         ylim(meanplotlev);
%         ylim([5 25]);
        grid on
        grid minor
        hold off
        drawnow;
        saveas(gcf,tifname,'tif');
        grid off
    end
end