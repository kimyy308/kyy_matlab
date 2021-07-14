close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
all_testname2 = {'test57'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
all_region2 ={'NWP','ES', 'AKP2', 'SS', 'YS'}

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
        conlev  = 0:5:35;
        meanplotlev =[-0.3 0.3];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [2006:2085]; % % put year which you want to plot [year year ...]
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

        load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
            'ssh_recon_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

        valnum=0;
    % %     valid cell number
    %      for vi=1:size(comb_spatial_meanressh,1)
    %          for vj=1:size(comb_spatial_meanressh,2)
    %              if (isnan(comb_spatial_meanressh(vi,vj,1))~=1)
    %                 valnum=valnum+1;
    %              end
    %          end
    %      end

    %     plot(squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan')))
    %     isize = size(comb_interped_data_filtered,1)
    %     jsize = size(comb_interped_data_filtered,2)
    %     lsize = size(comb_interped_data_filtered,3)
    %     comb_yearly_interped_data_filtered=reshape(comb_interped_data_filtered,[isize, jsize, 12, lsize/12]);
    %     mean_yearly_interped_data_filtered=squeeze(mean(mean(mean(comb_yearly_interped_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
    %     trendtime=14:29;
    %     p=polyfit(trendtime,mean_yearly_interped_data_filtered(14:29)',1);
    %     yearly_interped_trend=p(1);
    %     yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y


        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            recondir='E:\Data\Observation\OISST\monthly\';
        elseif (strcmp(system_name,'GLNXA64'))
        end
        std(mean(mean(comb_interped_data_filtered,'omitnan'),'omitnan'))
        std(mean(mean(comb_recon_data_filtered,'omitnan'),'omitnan'))


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


%          % recon rel trend plot
%             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%             hold on;
%             m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'-mean_recon_trend_filtered));
%     %         m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'));
%             shading(gca,m_pcolor_shading_method);
%             m_gshhs_i('color',m_gshhs_line_color);
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%             titlename = strcat('re SSH trend(rel)',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_recon_trend_filtered,2)));
%             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%             % set colorbar 
%             h = colorbar;
%             colormap(bwrmap);
%             set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis(trendlev);
% 
%             % set grid
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%             jpgname=strcat(outfile, '_', testname,'_',regionname, '_relative_recon_ssh_trend',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%             saveas(gcf,jpgname,'tif');
% 
%             disp(' ')
%             disp(['clim_', num2str(tempmonth), '_', 'relative_recon_ssh_trend', ' plot is created.'])
%             disp(' ')
%             disp([' File path is : ',jpgname])
%             disp(' ')
% 
%             hold off
%             close all;

%             % recon absolute trend plot
%             m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%             hold on;
%             m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'));
%     %         m_pcolor(double(recon_lon),recon_lat',squeeze(recon_trend_filtered(:,:)'));
%             shading(gca,m_pcolor_shading_method);
%             m_gshhs_i('color',m_gshhs_line_color);
%             m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
%             titlename = strcat('re SSH trend(abs)',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_recon_trend_filtered,2)));
%             title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% 
%             % set colorbar 
%             h = colorbar;
%             colormap(bwrmap);
%             set(h,'fontsize',colorbar_fontsize);
%             title(h,'mm/y','fontsize',colorbar_title_fontsize);
%             caxis(trendlev);
% 
%             % set grid
%             m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% 
%             set(gcf, 'PaperUnits', 'points');
%             set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%             set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% 
%             jpgname=strcat(outfile, '_', testname,'_',regionname, '_absolute_recon_ssh_trend',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%             saveas(gcf,jpgname,'tif');
% 
%             disp(' ')
%             disp(['clim_', num2str(tempmonth), '_', 'absolute_recon_ssh_trend', ' plot is created.'])
%             disp(' ')
%             disp([' File path is : ',jpgname])
%             disp(' ')
% 
%             hold off
%             close all;

        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
        for i= 1:length(recon_lon)
            for j=1:length(recon_lat)
                for l=1:size(comb_interped_data_filtered,3)
                    comb_interped_data_filtered_weight(i,j,l)=comb_interped_data_filtered(i,j,l).*cos(recon_lat(j)/180.0*pi);
                end
            end
        end
        
        temp_varvar=comb_interped_data_filtered_weight(:,:,1);
        for i= 1:length(recon_lon)
            for j=1:length(recon_lat)
                temp_varvar(i,j)=cos(recon_lat(j)/180.0*pi).*temp_varvar(i,j)./temp_varvar(i,j);
            end
        end
        for i= 1:length(recon_lon)
            for j=1:length(recon_lat)
                temp_varvar(isnan(temp_varvar))=NaN;
            end
        end
        m_cidfw=squeeze(sum(sum(comb_interped_data_filtered_weight,1,'omitnan'),2,'omitnan'));
        m_cidfw=m_cidfw./sum(sum(temp_varvar,1,'omitnan'),2,'omitnan');
        
        
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
        
        msl=squeeze(mean(mean(comb_interped_data,1,'omitnan'),2,'omitnan'));
        msl=msl-mean(msl);    
        p=polyfit(xData,msl',1);
        msl2=xData*p(1)+p(2);
        mslplot=plot(xData,msl,'k')
        hold on
        mslplot2=plot(xData,msl2,'Color','r')
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        xlabel('year')
        ylabel('Mean SSH (m)')
        title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend,2)), ' mm/y'])
        datetick('x','yymmm','keepticks')
        axis tight;
        ylim(meanplotlev)
        set(mslplot,'LineWidth',2);
        set(gca,'FontSize',15);
        grid on
        hold off
        saveas(gcf,jpgname,'jpg');
        grid off
        

%         recon_msl=squeeze(mean(mean(comb_recon_data_filtered,1,'omitnan'),2,'omitnan'));
%         recon_msl=recon_msl-mean(recon_msl);        
%         p=polyfit(xData,recon_msl',1);
%         recon_msl2=xData*p(1)+p(2);
%         recon_mslplot=plot(xData,recon_msl,'k')
%         hold on
%         recon_mslplot2=plot(xData,recon_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname,'_recon_msl_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('year')
%         ylabel('Mean SSH (m)')
%         title([regionname ', Mean recon SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_recon_trend_filtered,2)), ' mm/y'])
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(recon_mslplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
        close all;



        for i =1:length(inputyear) 
            tempyear=inputyear(i);
    %         for month=1:12
                xData2(i) = datenum([num2str(6,'%02i'),'-30-',num2str(tempyear)]);
    %         end
        end
        %     plot(squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan')))
        isize = size(comb_interped_data_filtered,1)
        jsize = size(comb_interped_data_filtered,2)
        lsize = size(comb_interped_data_filtered,3)
        comb_yearly_interped_data_filtered=reshape(comb_interped_data_filtered,[isize, jsize, 12, lsize/12]);
        mean_yearly_interped_data_filtered=squeeze(mean(mean(mean(comb_yearly_interped_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
        trendtime=1:length(xData2);
        p=polyfit(trendtime,mean_yearly_interped_data_filtered(1:length(xData2))',1);
        yearly_interped_trend=p(1);
        yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y

        yearly_msl=squeeze(mean(mean(mean(comb_yearly_interped_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
        yearly_msl=yearly_msl-mean(yearly_msl);    
        p=polyfit(xData2,yearly_msl',1);
        yearly_msl2=xData2*p(1)+p(2);
        yearly_mslplot=plot(xData2,yearly_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
        hold on
        yearly_mslplot2=plot(xData2,yearly_msl2,'Color','r')
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
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


        yearly_corrected_msl=squeeze(mean(mean(mean(comb_yearly_interped_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
        for i=1:length(xData2)
            yearly_corrected_msl(i)=yearly_corrected_msl(i)+(i-1)*0.0017;
        end
        yearly_corrected_msl=yearly_corrected_msl-mean(yearly_corrected_msl);    
        p=polyfit(xData2,yearly_corrected_msl',1);
        yearly_corrected_msl2=xData2*p(1)+p(2);
        yearly_corrected_mslplot=plot(xData2,yearly_corrected_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
        hold on
        yearly_corrected_mslplot2=plot(xData2,yearly_corrected_msl2,'Color','r')
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_corrected_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        xlabel('Year')
        ylabel('Mean SSH (m)')
        title([regionname, ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
            '), ',num2str(round(yearly_interped_trend+1.7,2)),' (',num2str(round(yearly_interped_trend,2)),'+',num2str(round(1.7,2)),')', ' mm/y'])
        ylim(meanplotlev)
        datetick('x','yymmm','keepticks')
        axis tight;
        ylim(meanplotlev)
        set(yearly_corrected_mslplot,'LineWidth',2);
        set(gca,'FontSize',15);
        grid on
        hold off
        saveas(gcf,jpgname,'jpg');
        grid off
        close all;

        comb_yearly_recon_data_filtered=reshape(comb_recon_data_filtered,[isize, jsize, 12, lsize/12]);
        mean_yearly_recon_data_filtered=squeeze(mean(mean(mean(comb_yearly_recon_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
        trendtime=1:length(xData2);
        p=polyfit(trendtime,mean_yearly_recon_data_filtered(1:length(xData2))',1);
        yearly_recon_trend=p(1);
        yearly_recon_trend = yearly_recon_trend * 1000.0; %% m/y -> mm/y

        yearly_recon_msl=squeeze(mean(mean(mean(comb_yearly_recon_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
        yearly_recon_msl=yearly_recon_msl-mean(yearly_recon_msl);        
        p=polyfit(xData2,yearly_recon_msl',1);
        yearly_recon_msl2=xData2*p(1)+p(2);
        yearly_recon_mslplot=plot(xData2,yearly_recon_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
        hold on
        yearly_recon_mslplot2=plot(xData2,yearly_recon_msl2,'Color','r')
        jpgname=strcat(outfile, '_', testname, '_',regionname,'_yearly_recon_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        xlabel('Year')
        ylabel('Mean SSH (m)')
        title([regionname ', Mean recon SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(yearly_recon_trend,2)), ' mm/y'])
        ylim(meanplotlev)
        datetick('x','yymmm','keepticks')
        axis tight;
        ylim(meanplotlev)
        set(yearly_recon_mslplot,'LineWidth',2);
        set(gca,'FontSize',15);
        grid on
        hold off
        saveas(gcf,jpgname,'jpg');
        grid off
        close all;
    end
end