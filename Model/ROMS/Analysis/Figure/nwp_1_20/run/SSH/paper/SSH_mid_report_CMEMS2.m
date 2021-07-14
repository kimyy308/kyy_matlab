close all; clear all;  clc;
warning off;
all_region2 ={'AKP2'}
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
    trendlev = [-10 10];  %% trend lev
    conlev  = 0:5:35;
    meanplotlev =[-0.15 0.15];
    msl_diff_lev=[-300 300];

    % for snu_desktop
%     testname='test49'   % % need to change
    inputyear = [1994:2017]; % % put year which you want to plot [year year ...]
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
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        case('AKP2')
            refpolygon=akp2polygon;
            msl_abs_lev=[-100 500];
            msllev= [-300 300];
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
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\CMEMS\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
%         filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\SODA\'); % % where data files are
        filedir = strcat('E:\Data\Observation\CMEMS\'); % % where data files are
%         cmemsdir='E:\Data\Model\ROMS\nwp_1_20\SODA\';
        cmemsdir='E:\Data\Observation\CMEMS\';
    elseif (strcmp(system_name,'GLNXA64'))
    end

    figdir=[figrawdir,'Trend\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
%     outfile = strcat(figdir,regionname);
    outfile = strcat(figdir);
    
     % cmems trend plot
     jpgname=strcat(outfile,regionname, '_cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
     if (exist(jpgname , 'file') ~= 2)

        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        hold on;
        m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%         titlename = strcat(regionname,', SODA SSH trend','(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'Mean=',num2str(round(mean_cmems_trend_filtered,2)),'mm/y');
        titlename = strcat(regionname,', CMEMS SLA trend','(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ,', 'Mean=',num2str(round(mean_cmems_trend_filtered,2)));
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

        saveas(gcf,jpgname,'tif');

        disp(' ')
        disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
     end
     
     figdir2=[figrawdir,'MSL\'];
        if (exist(strcat(figdir2) , 'dir') ~= 7)
            mkdir(strcat(figdir2));
        end 
    %     outfile = strcat(figdir,regionname);
        outfile2 = strcat(figdir2);
        
        isize = size(comb_cmems_data_filtered,1);
        jsize = size(comb_cmems_data_filtered,2);
        lsize = size(comb_cmems_data_filtered,3);
        comb_monthly_cmems_data=reshape(comb_cmems_adt,[isize, jsize, 12, lsize/12]);
        comb_yearly_cmems_data=squeeze(mean(comb_monthly_cmems_data,3,'omitnan'));  
        
         % % % %             % msl pcolor (anomaly)  
        jpgname=strcat(outfile2, '_cmems_',regionname, '_msl_ano_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg 
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            msl_alltime=squeeze(mean(comb_yearly_cmems_data,3,'omitnan'))'.*1000;
            msl_alltime=msl_alltime-mean(mean(msl_alltime,1,'omitnan'),2,'omitnan');
            m_pcolor(double(cmems_lon),cmems_lat',msl_alltime);
    %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            titlename = strcat('MSL',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(bwrmap);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'mm','fontsize',colorbar_title_fontsize);
            caxis(msllev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end  
        
        jpgname=strcat(outfile2, '_cmems_',regionname, '_msl_abs_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)        
% % % %             % msl pcolor (absolute, total year mean)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            msl_alltime=squeeze(mean(comb_yearly_cmems_data,3,'omitnan'))'.*1000;
            m_msl=mean(mean(msl_alltime,1,'omitnan'),2,'omitnan');
            m_pcolor(double(cmems_lon),cmems_lat',msl_alltime);
    %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            titlename = strcat('MSL',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ','M=',num2str(round(m_msl,1)));
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'mm','fontsize',colorbar_title_fontsize);
%             caxis(msllev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;          
        end
%           for yearij = 1:length(inputyear)
%           % % % %             % msl pcolor (absolute)
%             
%             tempyear= inputyear(yearij);
%             
%             
%             jpgname=strcat(outfile2, '_cmems_',regionname, '_msl_abs_',num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
%             if (exist(jpgname , 'file') ~= 2)               
%                 m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%                 hold on;
%                 msl_alltime=squeeze(comb_yearly_cmems_data(:,:,yearij))'.*1000;
%                 m_msl=mean(mean(msl_alltime,1,'omitnan'),2,'omitnan');
%                 m_pcolor(double(cmems_lon),cmems_lat',msl_alltime);
%         %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
%                 shading(gca,m_pcolor_shading_method);
%                 m_gshhs_i('color',m_gshhs_line_color);
%                 m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%                 titlename = strcat('MSL',',(',num2str(tempyear,'%04i'),') ','M=',num2str(round(m_msl,1)));
%                 title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%                 % set colorbar 
%                 h = colorbar;
%                 colormap(jet);
%                 set(h,'fontsize',colorbar_fontsize);
%                 title(h,'mm','fontsize',colorbar_title_fontsize);
%                 caxis(msl_abs_lev);
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
%                 disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
%                 disp(' ')
%                 disp([' File path is : ',jpgname])
%                 disp(' ')
% 
%                 hold off
%                 close all;   
%             end
%             
%             if (yearij~=length(inputyear))    
%                 jpgname=strcat(outfile2, '_cmems_',regionname, '_msl_diff_',num2str(inputyear(yearij+1),'%04i'),'_',num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
%                 if (exist(jpgname , 'file') ~= 2)
%                      % % % %             % msl pcolor (diff)
%                     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%                     hold on;
%                     msl_diff=squeeze(comb_yearly_cmems_data(:,:,yearij+1))'.*1000 - squeeze(comb_yearly_cmems_data(:,:,yearij))'.*1000;
%                     m_msl=mean(mean(msl_diff,1,'omitnan'),2,'omitnan');
%                     m_pcolor(double(cmems_lon),cmems_lat',msl_diff);
%                     shading(gca,m_pcolor_shading_method);
%                     m_gshhs_i('color',m_gshhs_line_color);
%                     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%                     titlename = strcat('MSL',',(',num2str(inputyear(yearij+1),'%04i'),'-',num2str(tempyear,'%04i'),') ','M=',num2str(round(m_msl,1)));
%                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%                     % set colorbar 
%                     h = colorbar;
%                     colormap(bwrmap);
%                     set(h,'fontsize',colorbar_fontsize);
%                     title(h,'mm','fontsize',colorbar_title_fontsize);
%                     caxis(msl_diff_lev);
% 
%                     % set grid
%                     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%                     set(gcf, 'PaperUnits', 'points');
%                     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%                     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%                     saveas(gcf,jpgname,'tif');
% 
%                     disp(' ')
%                     disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
%                     disp(' ')
%                     disp([' File path is : ',jpgname])
%                     disp(' ')
% 
%                     hold off
%                     close all;
%                 end
%             end
%              
%             
%             
%             figdir3=[figrawdir,'MSL\monthly\'];
%             outfile3 = strcat(figdir3,regionname);
%             if (exist(strcat(figdir3) , 'dir') ~= 7)
%                 mkdir(strcat(figdir3));
%             end
%             
%             for monthij=1:length(inputmonth)
%                 tempmonth=inputmonth(monthij);
%                 jpgname=strcat(outfile3, '_cmems_',regionname, '_msl_abs_',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'), '.tif'); %% ~_year_month.jpg
%                 if (exist(jpgname , 'file') ~= 2)
%                     m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%                     hold on;
%                     msl_alltime=squeeze(comb_monthly_cmems_data(:,:,monthij,yearij))'.*1000;
%                     m_msl=mean(mean(msl_alltime,1,'omitnan'),2,'omitnan');
%                     m_pcolor(double(cmems_lon),cmems_lat',msl_alltime);
%             %         m_pcolor(double(cmems_lon),cmems_lat',squeeze(cmems_trend_filtered(:,:)'));
%                     shading(gca,m_pcolor_shading_method);
%                     m_gshhs_i('color',m_gshhs_line_color);
%                     m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%                     titlename = strcat('MSL',',(',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'),') ','M=',num2str(round(m_msl,1)));
%                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%                     % set colorbar 
%                     h = colorbar;
%                     colormap(jet);
%                     set(h,'fontsize',colorbar_fontsize);
%                     title(h,'mm','fontsize',colorbar_title_fontsize);
%                     caxis(msl_abs_lev);
% 
%                     % set grid
%                     m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%                     set(gcf, 'PaperUnits', 'points');
%                     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%                     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%                     saveas(gcf,jpgname,'tif');
% 
%                     disp(' ')
%                     disp(['clim_', num2str(tempmonth), '_', 'absolute_cmems_ssh_trend', ' plot is created.'])
%                     disp(' ')
%                     disp([' File path is : ',jpgname])
%                     disp(' ')
% 
%                     hold off
%                     close all; 
%                 end
%             end
%           end
        
        
        
        
        
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
            xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',num2str(tempyear)]);
        end
    end

    cmems_msl=squeeze(mean(mean(comb_cmems_data_filtered,1,'omitnan'),2,'omitnan'));
    cmems_msl=cmems_msl-mean(cmems_msl);        
    p=polyfit(xData,cmems_msl',1);
    cmems_msl2=xData*p(1)+p(2);
    cmems_mslplot=plot(xData,cmems_msl,'k')
    hold on
    cmems_mslplot2=plot(xData,cmems_msl2,'Color','r')
    jpgname=strcat(outfile2, regionname,'_cmems_msl_', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('Year')
    ylabel('Mean SSH (m)')
    title([regionname ', Mean CMEMS SLA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_cmems_trend_filtered,2)), ' mm/y'])
    datetick('x','yyyy','keepticks')
    axis tight;
    ylim(meanplotlev);
    set(cmems_mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
    
    
    cmems_msl=squeeze(mean(mean(comb_cmems_data,1,'omitnan'),2,'omitnan'));
    cmems_msl=cmems_msl-mean(cmems_msl);        
    p=polyfit(xData,cmems_msl',1);
    cmems_msl2=xData*p(1)+p(2);
    cmems_mslplot=plot(xData,cmems_msl,'k')
    hold on
    cmems_mslplot2=plot(xData,cmems_msl2,'Color','r')
    jpgname=strcat(outfile2, regionname,'_cmems_adt_', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('Year')
    ylabel('Mean SSH (m)')
    title([regionname ', Mean CMEMS ADT(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_cmems_trend_filtered,2)), ' mm/y'])
    datetick('x','yyyy','keepticks')
    axis tight;
    ylim(meanplotlev);
    set(cmems_mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
    
    
    for lati=1:size(cmems_lat2,2)
        weight_msl(1:size(cmems_lon2,1),lati)=m_lldist(cmems_lon2(1:2,lati), cmems_lat2(1:2,lati),1);
    end
    weight_msl_masked=weight_msl;
    weight_msl_masked(isnan(comb_cmems_data_filtered(:,:,1))==1)=NaN;
    weighted_sum_sl=sum(sum(comb_cmems_data_filtered,1,'omitnan'),2,'omitnan')
    
    cmems_msl=squeeze(mean(mean(comb_cmems_data_filtered,1,'omitnan'),2,'omitnan'));
    cmems_msl=cmems_msl-mean(cmems_msl);        
    p=polyfit(xData,cmems_msl',1);
    cmems_msl2=xData*p(1)+p(2);
    cmems_mslplot=plot(xData,cmems_msl,'k')
    hold on
    cmems_mslplot2=plot(xData,cmems_msl2,'Color','r')
    jpgname=strcat(outfile2, regionname,'_cmems_msl_', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('Year')
    ylabel('Mean SSH (m)')
    title([regionname ', Mean CMEMS SLA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_cmems_trend_filtered,2)), ' mm/y'])
    datetick('x','yyyy','keepticks')
    axis tight;
    ylim(meanplotlev);
    set(cmems_mslplot,'LineWidth',2);
    set(gca,'FontSize',15);
    grid on
    hold off
    saveas(gcf,jpgname,'jpg');
    grid off
    
    
    jpgname=strcat(outfile2, 'cmems_',regionname, '_msl_abs_',num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    var='SSH';
    param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
    run(param_script);

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    hold on;
    
    madt = mean(comb_cmems_adt,3, 'omitnan');
    mask_model = double(inpolygon(cmems_lon2,cmems_lat2,refpolygon(:,1),refpolygon(:,2)));
    mask_model(mask_model==0)=NaN;
    madt= (madt .* mask_model);
    
    m_pcolor(cmems_lon,cmems_lat,madt'-mean(madt(:),'omitnan')+0.1);
    shading(gca,m_pcolor_shading_method);   

    m_gshhs_i('color',m_gshhs_line_color)  
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    titlename = strcat(var, ' mean, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(jet);
    set(h,'fontsize',colorbar_fontsize);
    title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
    caxis([-0.2 0.6]);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    saveas(gcf,jpgname,'tif');
    close all;
end

% m_lldist([115,164], [52,52], 1)
