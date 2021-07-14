close all; clear all;  clc;
% % horizontal SST trend plot (connected with SSH_mid_report1)
warning off;
% all_region2 ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region2 ={'AKP'};

all_testname = {'test55', 'test56'};

% all_region2 ={'NWP'}
for testnameind=1:length(all_testname)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind testnameind all_region all_testname
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
    testname=all_testname{testnameind}
    inputyear = [1982:2005]; % % put year which you want to plot [year year ...]
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
        otherwise
            ('?')
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));

        load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);

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
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
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

            jpgname=strcat(outfile, '_', testname, '_',regionname, '_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            saveas(gcf,jpgname,'jpg');

            disp(' ')
            disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
            disp(' ')
            disp([' File path is : ',jpgname])
            disp(' ')

            hold off
            close all;

         % avhrr trend plot
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

            jpgname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            saveas(gcf,jpgname,'jpg');

            disp(' ')
            disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
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

        msst=squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan'));
        % msst=msst-mean(msst);    
        p=polyfit(xData,msst',1);
        msst2=xData*p(1)+p(2);
        msstplot=plot(xData,msst,'k')
        hold on
        msstplot2=plot(xData,msst2,'Color','r')
        jpgname=strcat(outfile, '_', testname,'_',regionname, '_msst', '.jpg'); %% ~_year_month.jpg
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
end