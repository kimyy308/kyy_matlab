close all; clear all;  clc;
% % horizontal SST trend plot (connected with SSH_mid_report1)
warning off;
% all_region2 ={'AKP','NWP','ES', 'SS', 'YS'}
all_region2 ={'AKP4'};

all_testname2 = {'SODA_3_4_2'};

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

    shadlev = [0 35];
    rms_lev = [0 3];
    bias_lev= [-4 4];
    trendlev = [-0.1 0.1];
    trenddifflev= [-0.1 0.1];
    conlev  = 0:5:35;
    meanplotlev =[-4 4];
    meanplotlev2 =[-2 33];

    % for snu_desktop
    testname=all_testname2{testnameind2}
    inputyear = [1994:2017]; % % put year which you want to plot [year year ...]
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
        case('AKP4') %% Around Korea Peninsula
            refpolygon=akp4polygon;
        case('EKB')
            refpolygon=ekbpolygon;    
        otherwise
            ('?')
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));

        load(['E:\Data\Reanalysis\SODA\',testname,'\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
        
        run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);

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
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\SODA\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
%             filedir = strcat('H:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
            avhrrdir='E:\Data\Observation\OISST\monthly\';
        elseif (strcmp(system_name,'GLNXA64'))
        end


        trendlev = [-0.1 0.1];

        figdir=[figrawdir,'Trend\SST\'];
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

         % avhrr trend plot
        
         jpgname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
         if (exist(jpgname , 'file') ~= 2)        

            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(avhrr_trend_filtered(:,:)));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            titlename = strcat('AVHRR SST trend, ','(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ', 'M=',num2str(round(mean_avhrr_trend_filtered,3)));
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
        
         %             % rms plot
        jpgname=strcat(outfile, '_', testname, '_rms_', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'),'.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)        

            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            
            m_pcolor(double(avhrr_lon),avhrr_lat,mean(comb_rms_data,3,'omitnan'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            titlename = strcat('RMS', ' (', num2str(min(inputyear)),'-', num2str(max(inputyear)),')','Mrms=',num2str(mean(comb_rms_data(:), 'omitnan')));
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

% %             set colorbar 
            h = colorbar;
            colormap(wrmap);
            set(h,'fontsize',colorbar_fontsize);
            caxis(rms_lev);
% %             set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'jpg');


            disp(' ')
            disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'RMS', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;
        end
        
% % %                     % bias plot
        jpgname=strcat(outfile, '_', testname, '_bias_', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'),'.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_pcolor(double(avhrr_lon),avhrr_lat,mean(comb_bias_data,3,'omitnan'));
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            titlename = strcat('BIAS', ' (', num2str(min(inputyear)),'-', num2str(max(inputyear)),')','Mbias=',num2str(mean(comb_bias_data(:), 'omitnan')));
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

% %             % set colorbar 
            h = colorbar;
            colormap(bwrmap);
            set(h,'fontsize',colorbar_fontsize);
            caxis(bias_lev);
            

% %             % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

            saveas(gcf,jpgname,'jpg');


            disp(' ')
            disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')
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
        
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_msst', '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            avhrr_msst=squeeze(mean(mean(comb_avhrr_data_filtered,1,'omitnan'),2,'omitnan'));
            % avhrr_msst=avhrr_msst-mean(avhrr_msst);        
            p=polyfit(xData,avhrr_msst',1);
            avhrr_msst2=xData*p(1)+p(2);
            avhrr_msstplot=plot(xData,avhrr_msst,'k')
            hold on
            avhrr_msstplot2=plot(xData,avhrr_msst2,'Color','r')
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_msst', '.jpg'); %% ~_year_month.jpg
            xlabel('month')
            ylabel('Mean SSTA (^oC)')
            title([regionname,', Mean AVHRR SSTA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_avhrr_trend_filtered,3)), '^oC/y'])
            set(avhrr_msstplot,'LineWidth',2);
            datetick('x','yyyy','keepticks')
            set(gca,'FontSize',15);
            axis tight
            ylim(meanplotlev);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
        
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_msst_nonfilt_', '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            avhrr_msst=squeeze(mean(mean(comb_avhrr_data,1,'omitnan'),2,'omitnan'));
            p=polyfit(xData,avhrr_msst',1);
            avhrr_msst2=xData*p(1)+p(2);
            avhrr_msstplot=plot(xData,avhrr_msst,'k')
            hold on
            avhrr_msstplot2=plot(xData,avhrr_msst2,'Color','r')
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_msst_nonfilt_', '.jpg'); %% ~_year_month.jpg
            xlabel('month')
            ylabel('Mean SSTA (^oC)')
            title([regionname,', Mean AVHRR SSTA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_avhrr_trend,3)), '^oC/y'])
            set(avhrr_msstplot,'LineWidth',2);
            datetick('x','yyyy','keepticks')
            set(gca,'FontSize',15);
            axis tight
            ylim(meanplotlev2);
            grid on
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
        end
        
        jpgname=strcat(outfile, '_', testname,'_',regionname,'_climrms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)
            mrms=mean(comb_rms_data(:),'omitnan')
            rmsplot=plot(mean(comb_meanrms,1),'k')
            jpgname=strcat(outfile, '_', testname,'_',regionname,'_climrms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
            xlabel('month')
            ylabel('RMS(^oC)')
            title(['RMS, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), M=', num2str(round(mrms,2))])
            ylim([0 4])
            set(rmsplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            saveas(gcf,jpgname,'jpg');
            grid off
        end
        
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_climbias_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
        if (exist(jpgname , 'file') ~= 2)        
            mbias=mean(comb_bias_data(:),'omitnan')
            biasplot=plot(mean(comb_meanbias,1) ,'k')
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_climbias_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
            xlabel('month')
            ylabel('bias(^o)')
            title(['BIAS, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), M=', num2str(round(mbias,2))])
            ylim([-4 4])
            set(biasplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            saveas(gcf,jpgname,'jpg');
            grid off
        end
        
        % %         trend difference
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_trend_diff_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
            hold on;
            m_trend_diff=trend_filtered-avhrr_trend_filtered;
            m_pcolor(double(avhrr_lon),avhrr_lat,m_trend_diff);
            shading(gca,m_pcolor_shading_method);
            m_gshhs_i('color',m_gshhs_line_color);
            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            titlename = strcat('Trend Diff', ' (', num2str(min(inputyear)),'-', num2str(max(inputyear)),')','M=',num2str(round(mean(m_trend_diff(:), 'omitnan'),4)));
            title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

% %             % set colorbar 
            h = colorbar;
            colormap(bwrmap);
            set(h,'fontsize',colorbar_fontsize);
            caxis(trenddifflev);
            

% %             % set grid
            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
            
            hold off
            saveas(gcf,jpgname,'jpg');


            disp(' ')
            disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')
%         end
    end
end

% SSH_4th_mid_report6