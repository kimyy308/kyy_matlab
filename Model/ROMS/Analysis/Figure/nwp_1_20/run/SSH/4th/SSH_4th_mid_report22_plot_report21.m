close all; clear all;  clc;
% % horizontal SST trend plot (connected with SSH_mid_report1)
warning off;
% all_region2 ={'AKP','NWP','ES', 'SS', 'YS'}
all_region2 ={'NWP', 'AKP2'};

all_testname2 = {'test49', 'test52'};

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
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end
   
    % for snu_desktop
    testname=all_testname2{testnameind2}
    inputyear = [1980:2008]; % % put year which you want to plot [year year ...]
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
        case('AKP') %% Around Korea Peninsula
            refpolygon=akppolygon;
        case('AKP2') %% Around Korea Peninsula
            refpolygon=akp2polygon;
        case('EKB')
            refpolygon=ekbpolygon;    
        otherwise
            ('?')
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));

        load(['E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_model_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
        
        shadlev = [-2 33];
        rms_shadlev = [0 4];
        trendlev = [-0.1 0.1];
        conlev  = -2:5:33;
        meanplotlev =[-4 4];
        meanplotlev2 =[-2 33];
        msst_lev= [-2 33];
        msst_abs_lev=[-2 33];
        msst_diff_lev=[-10 10];
    
            valnum=0;
    % %     valid cell number
         for vi=1:size(comb_spatial_meanavhrr,1)
             for vj=1:size(comb_spatial_meanavhrr,2)
                 if (isnan(comb_spatial_meanavhrr(vi,vj,1))~=1)
                    valnum=valnum+1;
                 end
             end
         end


        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('H:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            avhrrdir='E:\Data\Observation\OISST\monthly\';
        elseif (strcmp(system_name,'GLNXA64'))
        end


        trendlev = [-0.1 0.1];

        figdir=[figrawdir,'Trend\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        
            % trend plot
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(trend_filtered(:,:)));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname, ', SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)));  %% + glacier contribution
                titlename = strcat('SST trend, ',testname, ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ', 'M=',num2str(round(mean_trend_filtered,3)));
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

                saveas(gcf,jpgname,'jpg');

                disp(' ')
                disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
        
        
        isize = size(comb_interped_data_filtered,1)
        jsize = size(comb_interped_data_filtered,2)
        lsize = size(comb_interped_data_filtered,3)
        comb_monthly_interped_data=reshape(comb_interped_data,[isize, jsize, 12, lsize/12]);
        comb_yearly_interped_data=squeeze(mean(comb_monthly_interped_data,3,'omitnan'));
        figdir2=[figrawdir,'MSST\'];
        outfile2 = strcat(figdir2,regionname);
        if (exist(strcat(figdir2) , 'dir') ~= 7)
            mkdir(strcat(figdir2));
        end
            
          % % % %             % msst pcolor (anomaly)  
        jpgname=strcat(outfile2, '_', testname,'_',regionname, '_msst_ano_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg 
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            msst_alltime=squeeze(mean(comb_yearly_interped_data,3,'omitnan'));
            msst_alltime=msst_alltime-mean(mean(msst_alltime,1,'omitnan'),2,'omitnan');
            m_pcolor(double(avhrr_lon),avhrr_lat,msst_alltime);
    %         m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(avhrr_trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            titlename = strcat('MSST',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
        
            % set colorbar 
            h = colorbar;
            colormap(bwrmap);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'^oC','fontsize',colorbar_title_fontsize);
            caxis(msst_diff_lev);

            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', num2str(tempmonth), '_', 'sst_ano', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end  
        
        jpgname=strcat(outfile2, '_', testname,'_',regionname, '_msst_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)        
% % % %             % msst pcolor (total year mean)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            msst_alltime=squeeze(mean(comb_yearly_interped_data,3,'omitnan'));
            m_msst=mean(mean(msst_alltime,1,'omitnan'),2,'omitnan');
            m_pcolor(double(avhrr_lon),avhrr_lat,msst_alltime);
    %         m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(avhrr_trend_filtered(:,:)'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            titlename = strcat('MSST',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ','M=',num2str(round(m_msst,1)));
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);
            set(h,'fontsize',colorbar_fontsize);
            title(h,'^oC','fontsize',colorbar_title_fontsize);
%             caxis(msstlev);

            % contour
            [C,h2]=m_contour(double(avhrr_lon),avhrr_lat, msst_alltime, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
            clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
            
            % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'tif');

            disp(' ')
            disp(['clim_', num2str(tempmonth), '_', 'msst', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;          
        end
        
        for yearij = 1:length(inputyear)
        % % % %             % yearly msst pcolor (absolute)
            tempyear= inputyear(yearij);
            jpgname=strcat(outfile2, '_', testname,'_',regionname, '_msst_',num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)               
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                msst_alltime=squeeze(comb_yearly_interped_data(:,:,yearij));
                m_msst=mean(mean(msst_alltime,1,'omitnan'),2,'omitnan');
                m_pcolor(double(avhrr_lon),avhrr_lat,msst_alltime);
            %         m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(avhrr_trend_filtered(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                titlename = strcat('MSST',',(',num2str(tempyear,'%04i'),') ','M=',num2str(round(m_msst,1)));
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'^oC','fontsize',colorbar_title_fontsize);
                caxis(msst_lev);
                
                % contour
                [C,h2]=m_contour(double(avhrr_lon),avhrr_lat, msst_alltime, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
                clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                    'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
                
                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', num2str(tempmonth), '_', 'yearly_msst', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;   
            end

            if (yearij~=length(inputyear))    
                jpgname=strcat(outfile2, '_', testname,'_',regionname, '_msst_diff_',num2str(inputyear(yearij+1),'%04i'),'_',num2str(tempyear,'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
                     % % % %             % msst pcolor (diff)
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    msst_diff=squeeze(comb_yearly_interped_data(:,:,yearij+1)) - squeeze(comb_yearly_interped_data(:,:,yearij));
                    m_msst=mean(mean(msst_diff,1,'omitnan'),2,'omitnan');
                    m_pcolor(double(avhrr_lon),avhrr_lat,msst_diff);
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('MSST',',(',num2str(inputyear(yearij+1),'%04i'),'-',num2str(tempyear,'%04i'),') ','M=',num2str(round(m_msst,1)));
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(bwrmap);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'^oC','fontsize',colorbar_title_fontsize);
                    caxis(msst_diff_lev);

                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['clim_', num2str(tempmonth), '_', 'yearly_sst_diff', ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all;
                end
            end

            figdir3=[figrawdir,'MSST\monthly\'];
            outfile3 = strcat(figdir3,regionname);
            if (exist(strcat(figdir3) , 'dir') ~= 7)
                mkdir(strcat(figdir3));
            end

            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                jpgname=strcat(outfile3, '_', testname,'_',regionname, '_msst_',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2)
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    msst_alltime=squeeze(comb_monthly_interped_data(:,:,monthij,yearij));
                    m_msst=mean(mean(msst_alltime,1,'omitnan'),2,'omitnan');
                    m_pcolor(double(avhrr_lon),avhrr_lat,msst_alltime);
            %         m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(avhrr_trend_filtered(:,:)'));
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                    titlename = strcat('MSST',',(',num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'),') ','M=',num2str(round(m_msst,1)));
                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(jet);
                    set(h,'fontsize',colorbar_fontsize);
                    title(h,'^oC','fontsize',colorbar_title_fontsize);
                    caxis(msst_lev);
                    
                    % contour
                    [C,h2]=m_contour(double(avhrr_lon),avhrr_lat, msst_alltime, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
                    clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
                        'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
            
                    % set grid
                    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                    set(gcf, 'PaperUnits', 'points');
                    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                    saveas(gcf,jpgname,'tif');

                    disp(' ')
                    disp(['clim_', num2str(tempmonth), '_', 'monthly_sst', ' plot is created.'])
                    disp(' ')
                    disp([' File path is : ',jpgname])
                    disp(' ')

                    hold off
                    close all; 
                end
            end
        end    
            
            
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',num2str(tempyear)]);
            end
        end

        jpgname=strcat(outfile, '_', testname,'_',regionname, '_msst', '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            msst=squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan'));
            % msst=msst-mean(msst);    
            p=polyfit(xData,msst',1);
            msst2=xData*p(1)+p(2);
            msstplot=plot(xData,msst,'k')
            hold on
            msstplot2=plot(xData,msst2,'Color','r')
            xlabel('month')
            ylabel('Mean SSTA (^oC)')
            title([regionname, ', Mean SSTA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend_filtered,3)), ' ^oC/y'])
            set(msstplot,'LineWidth',2);
            datetick('x','yyyy','keepticks')
            set(gca,'FontSize',15);
            axis tight;
            ylim(meanplotlev);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
        
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_msst_nonfilt', '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            msst=squeeze(mean(mean(comb_interped_data,1,'omitnan'),2,'omitnan'));
            % msst=msst-mean(msst);    
            p=polyfit(xData,msst',1);
            msst2=xData*p(1)+p(2);
            msstplot=plot(xData,msst,'k')
            hold on
            msstplot2=plot(xData,msst2,'Color','r')
            xlabel('month')
            ylabel('Mean SSTA (^oC)')
            title([regionname, ', Mean SSTA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend,3)), ' ^oC/y'])
            set(msstplot,'LineWidth',2);
            datetick('x','yyyy','keepticks')
            set(gca,'FontSize',15);
            axis tight;
            ylim(meanplotlev2);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
    end
end

% SSH_4th_mid_report6