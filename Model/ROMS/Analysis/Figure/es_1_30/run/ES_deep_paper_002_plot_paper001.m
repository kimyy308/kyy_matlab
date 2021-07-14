close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
all_testname2 = {'v52'};
% all_testname2 = {'test54'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP4', 'YS', 'YSECS', 'ECS2', 'NES', 'SES', 'ES'}

% all_region2 ={'NWP', 'AKP4'};
% all_region2 ={'ES_deep', 'SES_deep', 'NES'};
all_region2 ={'NWES', 'NEES','SES_deep', 'ES_deep'};

% all_region2 ={'TEST'}
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
            dropboxpath='C:\Users\user\Dropbox';
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
        abstrendlev =[2 7];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev_wind =[0 15];
        meanplotlev_shflux =[-400 200];
        meanplotlev_yearly_wind =[3 5];
        meanplotlev_yearly_shflux =[-100 0];
        meanplotlev_yearly_latent =[-150 -50];
        meanplotlev_monthly_wind =[0 15];
        meanplotlev_monthly_shflux =[-400 200];
        stdspeedlev = [0.3 2.7];
        stdshfluxlev = [0 90];
       
        
%         ana_stddepth=-[1, 500, 1000, 2000];
        stddepth=[-1,-10,-20,-30,-50,-75,-100,-125,-150,-200,-300,-400,-500,-1000,-1500,-2000,-2500,-3000];
        ana_stddepth=stddepth;
        variable_energy={'MKE', 'MPE', 'EKE', 'EPE'};
        variable_atm={'speed', 'shflux'};
        lag_month=[-3:-1, 1:3];
                            
                                        
        byrmap = customcolormap_preset('red-yellow-blue');
        yrmap = byrmap(129:256,:);

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1994:2012]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]

        varname ='wind';
        variable='wind';
        run('es_polygon_point.m');
        regionname=all_region2{regionind2};
        
% % %         switch region
        for folding=1:1
            switch(regionname)
                case('ES') %% East Sea
                    refpolygon=espolygon;
                case('ES_deep') %% East Sea
                    refpolygon=es_deeppolygon;
                case('ES_KHOA') %% East Sea
                    refpolygon=es_khoapolygon;
                case('NES') %% Northern East Sea
                    refpolygon=nespolygon;
                case('NWES') %% North Western East Sea
                    refpolygon=nwespolygon;
                case('NEES') %% North Eastern East Sea
                    refpolygon=neespolygon;
                case('CES_deep') %% Northern East Sea (Center)
                    refpolygon=ces_deeppolygon;
                case('NES_deep') %% Northern East Sea
                    refpolygon=nes_deeppolygon;
                case('SES_deep') %% Southern East Sea
                    refpolygon=ses_deeppolygon;
                case('SES') %% Southern East Sea
                    refpolygon=sespolygon;
                case('CA') %% Coastal Area around korea peninsula
                    refpolygon=capolygon;
                case('EKB') %% Coastal Area around korea peninsula
                    refpolygon=ekbpolygon;
                case('TEST') %% for debugging
                    refpolygon=testpolygon;
                otherwise
                    ('?')
            end
            lonlat(1)=min(refpolygon(:,1));
            lonlat(2)=max(refpolygon(:,1));
            lonlat(3)=min(refpolygon(:,2));
            lonlat(4)=max(refpolygon(:,2));
        end
        
% % %         flag configuration
        for folding=1:1
            fig_flags{1,1}='trend_absolute_bwrmap_plot';
            fig_flags{2,1}='trend_climatological_plot';
            fig_flags{3,1}='time_series_ssh_seasonal_filtered';
            fig_flags{4,1}='time_series_ssh';
            fig_flags{5,1}='time_series_ssh_climatological';
            fig_flags{6,1}='trend_climatological_line';
            fig_flags{7,1}='time_series_ssh_seasonal';
            fig_flags{8,1}='amplitude_seasonal_ssh_cmemsstructed';
            fig_flags{9,1}='time_series_ssh_yearly';
            fig_flags{10,1}='std_cmems_sla_plot';
            fig_flags{11,1}='pattern_correlation';
            fig_flags{12,1}='time_series_wind_speed';
            fig_flags{13,1}='time_series_shflux';
            fig_flags{14,1}='time_series_latent';
            fig_flags{15,1}='time_series_sensible';

            
            for flagi=1:68
                fig_flags{flagi,2}=0;
            end
            fig_flags{12,2}=1;
            fig_flags{13,2}=1;
            fig_flags{14,2}=1;
            fig_flags{15,2}=0;
            fig_flags{16,2}=1;
            fig_flags{17,2}=1;
            fig_flags{18,2}=1;
            fig_flags{19,2}=1;
            fig_flags{20,2}=1;
            fig_flags{21,2}=1;
            fig_flags{22,2}=1;
            fig_flags{23,2}=1;
            fig_flags{24,2}=1;
            fig_flags{25,2}=1;
            fig_flags{26,2}=1;
            fig_flags{27,2}=1;
            fig_flags{28,2}=1;
            fig_flags{29,2}=1;
            fig_flags{30,2}=1;
            fig_flags{31,2}=1;
            fig_flags{32,2}=0;
            fig_flags{33,2}=0;
            fig_flags{34,2}=0;
            fig_flags{35,2}=1;
            fig_flags{36,2}=1;
            fig_flags{37,2}=1;
            fig_flags{38,2}=1;
        end

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
%         filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\research\Ph_D_course\EJS_deep_circulation\figure\es_1_30\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\es_1_30\run\fig_param\fig_param_kyy_', regionname, '.m'];
            filedir = strcat('D:\Data\Model\ROMS\es_1_30\', testname, '\output\run\short_monthly\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        rawfilename = [filedir,testname,'_',regionname,'model_atm_variable_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

%         modelfilename = strcat(filedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         corrfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         lpffilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         lpf_corrfilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         movfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         mov_corrfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         detrendfilename = strcat(filedir, testname,'_',regionname, 'detrended_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

        valnum=0;
        run('C:\Users\user\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);

        run(param_script);
        
        figdir=[figrawdir,'ATM\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
        trendoutfile = strcat(figdir,regionname);
        
%         lon_rho=ncread(modelfilename,'lon_rho');
%         lat_rho=ncread(modelfilename,'lat_rho');

% start-------------------- make timedata(xData) for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
            xData2(i) = datenum([num2str(6,'%02i'),'-30-',num2str(tempyear)]);
        end
        
        load(rawfilename)
% start-------------------- absolute trend plot (bwrmap)
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            jpgname=strcat(trendoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_bwrmap_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(double(lon_rho)',lat_rho',squeeze(trend_filtered(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                mean_trend_filtered=mean(trend_filtered(:),'omitnan');

        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('SSH trend(abs), ',testname, ',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend_filtered,2)), ' mm/y');  %% + glacier contribution

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
            fig_flag=0;
        end

% start-------------------- climatological absolute trend plot
        fig_flag=fig_flags{2,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_absolute_ssh_trend_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(clim_trend(:,:,monthij)'.*1000));
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                    titlename = strcat('SSH trend(abs), ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                        ', ', calendarname{monthij}(1:3), '), ','M=',num2str(round(mean_clim_trend(monthij)*1000,2)),'mm/y');  %% + glacier contribution

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
            fig_flag=0;
        end

% start-------------------- msl time series (seasonal filtered)
        fig_flag=fig_flags{3,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                for varind=1:length(inputyear)*12
                    msl_filt(varind)=squeeze(mean(mean(ncread(modelfilename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
                end
        %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
                msl_filt=msl_filt-mean(msl_filt);    
                p2=polyfit(1:length(msl_filt),msl_filt,1);
                p2=p2(1)*1000.0*12.0;

                p=polyfit(xData,msl_filt,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl_filt,'k')
                hold on
                mslplot2=plot(xData,msl2,'Color','r')
                xlabel('year')
                ylabel('Mean SSH (m)')
                title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
            end
            close all;
            fig_flag=0;
        end

% start-------------------- msl time series
        fig_flag=fig_flags{4,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                for varind=1:length(inputyear)*12
                    msl(varind)=squeeze(mean(mean(ncread(modelfilename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
                end
                msl=msl-mean(msl);     
                p2=polyfit(1:length(msl),msl,1);
                p2=p2(1)*1000.0*12.0;

                p=polyfit(xData,msl,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl,'k')
                hold on
                mslplot2=plot(xData,msl2,'Color','r')
                xlabel('year')
                ylabel('Mean SSH (m)')
    %             mean_trend=mean(trend(:), 'omitnan');
                title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' mm/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- climatological msl time series
        fig_flag=fig_flags{5,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij=1:12
                jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    if (exist('msl')==0)
                        for varind=1:length(inputyear)*12
                            msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
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
                        ', ', calendarname{monthij}(1:3), '), ',num2str(round(mean_clim_trend(monthij)*1000,2)), ' mm/y'])
    %                 datetick('x','yyyy','keepticks')
                    axis tight;
                    ylim(meanplotlev)
                    set(mslplot,'LineWidth',2);
                    set(gca,'FontSize',15);
                    grid on
                    hold off
                    saveas(gcf,jpgname,'jpg');
                    grid off
                    close all;
                end
            end
            fig_flag=0;
        end

% start-------------------- climatological msl trend (line)
        fig_flag=fig_flags{6,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            clim_trend=ncread(modelfilename,'clim_ssh_trend');
            clim_trend_reshaped=reshape(clim_trend, [size(clim_trend,1)*size(clim_trend,2), 12]);
            mean_clim_trend=mean(clim_trend_reshaped,1,'omitnan');
            if abs(mean_clim_trend(1))<0.01
                mean_clim_trend=mean_clim_trend*1000.0  %% m -> mm
            end
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                mslplot=plot(1:12,mean_clim_trend,'k')
                hold on
                xlabel('month')
                ylabel('trend (mm/yr)')
                title([regionname, ', climatological trend(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(trendplotlev)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- seasonal msl time series
        fig_flag=fig_flags{7,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

    %         clim_ssh=ncread(filename,'clim_ssh');
    %         for ijij=1:12
    %             temp_clim_ssh=clim_ssh(:,:,ijij);
    %             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %         end
    %         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                comb_interped_data=ncread(interpedfilename,'interped_ssh');
                lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                lat_cmems=ncread(cmemsfilename, 'lat_cmems');
                len_lon=size(lon_cmems,1);
                len_lat=size(lat_cmems,2);
                comb_interped_clim_divided=reshape(comb_interped_data,[len_lon, len_lat, 12, length(inputyear)]);
                clim_ssh=mean(comb_interped_clim_divided,4,'omitnan');
                cmems_trend=ncread(cmemsfilename, 'cmems_trend');
                cmems_mask=ones(size(cmems_trend));
                cmems_mask(isnan(cmems_trend))=NaN;
                clim_ssh = clim_ssh .* cmems_mask;
                for ijij=1:12
                    temp_clim_ssh=clim_ssh(:,:,ijij);
                    mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
                end
                mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

                mslplot=plot(1:12,mean_clim_ssh,'k')
                hold on
                xlabel('month')
                ylabel('trend (mm/yr)')
                title([regionname, ', seasonal msl(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
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
            fig_flag=0;
        end

% start-------------------- cmems seasonal msl amplitude
        fig_flag=fig_flags{8,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_amp_cmems_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg

    %         clim_ssh=ncread(filename,'clim_ssh');
    %         for ijij=1:12
    %             temp_clim_ssh=clim_ssh(:,:,ijij);
    %             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
    %         end
    %         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                comb_interped_data=ncread(filename,'interped_ssh');
                lon_cmems=ncread(filename, 'lon_cmems');
                lat_cmems=ncread(filename, 'lat_cmems');
                len_lon=length(lon_cmems);
                len_lat=length(lat_cmems);
                clim_ssh=ncread(filename, 'clim_cmems_ssh');
                clim_ssh(clim_ssh>1000000)=NaN;
                cmems_trend=ncread(filename, 'cmems_trend');
                cmems_mask=ones(size(cmems_trend));
                cmems_mask(isnan(cmems_trend))=NaN;
                clim_ssh = clim_ssh .* cmems_mask;
                for ijij=1:12
                    temp_clim_ssh=clim_ssh(:,:,ijij);
                    mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
                end
                mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
                for loni=1:size(clim_ssh,1)
                    for lati=1:size(clim_ssh,2)
                        amp_clim_ssh(loni,lati)=(max(clim_ssh(loni,lati,:))-min(clim_ssh(loni,lati,:)))/2.0;
                    end
                end
    %             pcolor(amp_clim_ssh'*100)

                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                m_pcolor(lon_cmems',lat_cmems',squeeze(amp_clim_ssh(:,:)'*100.0));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('seasonal amp, ','cmems',',(',num2str(min(inputyear),'%04i'),'-', ...
                    num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean(amp_clim_ssh(:)*100.0, 'omitnan'),2)));  

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                caxis([2 16]);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['clim_', 'cmems seasonal msl amplitude', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- yearly msl time series
        fig_flag=fig_flags{9,2};
        while (fig_flag)
            for i =1:length(inputyear) 
                tempyear=inputyear(i);
        %         for month=1:12
                    xData2(i) = datenum([num2str(6,'%02i'),'-30-',num2str(tempyear)]);
        %         end
            end
            %     plot(squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan')))
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                comb_data_filtered=ncread(interpedfilename, 'interped_sla_filtered');
                isize = size(comb_data_filtered,1);
                jsize = size(comb_data_filtered,2);
                lsize = size(comb_data_filtered,3);
    %             comb_yearly_data_filtered=reshape(comb_data_filtered,[isize, jsize, 12, lsize/12]);
                comb_yearly_data_filtered=reshape(comb_data_filtered,[isize*jsize*12, lsize/12]);

    %             mean_yearly_data_filtered=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
                mean_yearly_data_filtered=squeeze(mean(comb_yearly_data_filtered,1,'omitnan'));

                trendtime=1:length(xData2);
                p=polyfit(trendtime,mean_yearly_data_filtered(1:length(xData2)),1);
                yearly_interped_trend=p(1);
                yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y

%                 yearly_msl=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
                yearly_msl=mean_yearly_data_filtered;

                yearly_msl=yearly_msl-mean(yearly_msl);    
                p=polyfit(xData2,yearly_msl,1);
                yearly_msl2=xData2*p(1)+p(2);
                yearly_mslplot=plot(xData2,yearly_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
                hold on
                yearly_mslplot2=plot(xData2,yearly_msl2,'Color','r')
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
            end
            fig_flag=0;            
        end

% start-------------------- cmems std plot
        fig_flag=fig_flags{10,2};
        while (fig_flag)

            pngname=strcat(climoutfile, '_', testname,'_',regionname, '_cmems_', 'STD','_',num2str(min(inputyear),'%04i'), ...
                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
%                 for yearij=1:length(inputyear)
%                     tempyear=inputyear(yearij);
%                     yearstr=num2str(tempyear, '%04i');
%                     for monthij=1:length(inputmonth)
%                         tempmonth=inputmonth(monthij);
%                         monthstr=num2str(tempmonth, '%02i');
%                         varind=((yearij-1)*12)+monthij
                        if (exist('lon_min' , 'var') ~= 1)
                            lon_cmems=ncread(cmemsfilename, 'lon_cmems');
                            lat_cmems=ncread(cmemsfilename, 'lat_cmems');
%                             [lat_rho, lon_rho]= meshgrid(lat_cmems, lon_cmems);
                            [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                            cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                            cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                        end
%                         data_info = ncinfo(filename, varname); 
%                         data(:,:)=ncread(filename,'cmems_sla',[lon_min(1) lat_min(1) varind], ...
%                             [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                         if (exist('mean_data' , 'var') ~= 1)
%                             mean_data=zeros(size(data));
%                         end
%                         mean_data=mean_data + (data / (length(inputyear) * length(inputmonth)));
%                     end
%                 end
                cmems_sla=ncread(cmemsfilename,'cmems_sla', [lon_min(1) lat_min(1) 1], ...
                            [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf]);
                for loni=1:size(cmems_sla,1)
                    for lati=1:size(cmems_sla,2)
                        cmems_sla_var(loni,lati)=var(cmems_sla(loni,lati,:).*100);
                    end
                end
%                 cmems_msl=mean(mean_data(:), 'omitnan');
                m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                hold on;
                
                mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
                mask_model(mask_model==0)=NaN;
%                 mean_data= mean_data .* mask_model;
                cmems_sla_var=cmems_sla_var .* mask_model;
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-cmems_msl);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-0.45);
%                 m_pcolor(cut_lon_rho',cut_lat_rho',mean_data'-mean(mean_data(:),'omitnan'));
                m_pcolor(cut_lon_rho',cut_lat_rho',sqrt(cmems_sla_var)');

                shading(gca,m_pcolor_shading_method);   
                
                m_gshhs_i('color',m_gshhs_line_color)  
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
                titlename = strcat(variable, ' STD, ','cmems',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') ');  %% + glacier contribution
                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                % set colorbar 
                h = colorbar;
                colormap(jet);
                set(h,'fontsize',colorbar_fontsize);
                title(h,'(cm)','fontsize',colorbar_title_fontsize);
                caxis([0 15]);
            
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
                saveas(gcf,pngname,'tif');
                close all;
                clear lon_rho mean_data
            end
            fig_flag=0;
        end
        
% start-------------------- model pattern correlation    
        fig_flag=fig_flags{11,2};
        while (fig_flag)
            cmems_adt=ncread(filename,'cmems_adt');
            interped_ssh=ncread(filename,'interped_ssh');
            cmems_sla=ncread(filename,'cmems_sla');
            interped_sla=ncread(filename,'interped_sla');
            m_cmems_adt=mean(cmems_adt,3,'omitnan');
            m_interped_ssh=mean(interped_ssh,3,'omitnan');
            m_cmems_sla=mean(cmems_sla,3,'omitnan');
            m_interped_sla=mean(interped_sla,3,'omitnan');
            m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
            m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));

            corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))))
            corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))))
            corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
                cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)

        % end-------------------- model pattern correlation


        % start-------------------- model pattern correlation (yearly)
            pngname=strcat(outfile, '_', testname,'_',regionname, '_', nc_varname,'_','yearly_pattern_corr_',num2str(min(inputyear),'%04i'), ...
                                '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                for yearij=1:length(inputyear)
                    tempyear=inputyear(yearij);

                    yearly_cmems_adt=ncread(filename,'cmems_adt',[1 1 (yearij-1)*12+1], [inf inf 12]);
                    yearly_interped_ssh=ncread(filename,'interped_ssh',[1 1 (yearij-1)*12+1], [inf inf 12]);
                    yearly_cmems_sla=ncread(filename,'cmems_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);
                    yearly_interped_sla=ncread(filename,'interped_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);

                    m_cmems_adt=mean(yearly_cmems_adt,3,'omitnan');
                    m_interped_ssh=mean(yearly_interped_ssh,3,'omitnan');
                    m_cmems_sla=mean(yearly_cmems_sla,3,'omitnan');
                    m_interped_sla=mean(yearly_interped_sla,3,'omitnan');
        %             m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
        %             m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));

                    ssh_corr(yearij) = corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))));
                    sla_corr(yearij) = corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))));
        %             corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
        %                 cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)
                end
                plot(ssh_corr);
                hold on
                plot(sla_corr);
                hold off
                close all
            end
            fig_flag=0;
        end

% start-------------------- mean wind speed time series
        fig_flag=fig_flags{12,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_mspeed_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                
                msl=mean_speed;
                p2=polyfit(1:length(msl),msl,1);
                p2=p2(1)*12.0;

                p=polyfit(xData,msl,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl,'k');
                hold on
                mslplot2=plot(xData,msl2,'Color','r');
                xlabel('year')
                ylabel('Mean wind speed (m/s)')
    %             mean_trend=mean(trend(:), 'omitnan');
                title([regionname, ', Wind speed(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' m/s/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev_wind)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- mean shflux time series
        fig_flag=fig_flags{13,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_mshflux_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                
                msl=mean_shflux;
                p2=polyfit(1:length(msl),msl,1);
                p2=p2(1)*12.0;

                p=polyfit(xData,msl,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl,'k');
                hold on
                mslplot2=plot(xData,msl2,'Color','r');
                xlabel('year')
                ylabel('Mean shflux (W/m^2)')
    %             mean_trend=mean(trend(:), 'omitnan');
                title([regionname, ', shflux(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' W/m^2/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev_shflux)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end
        
% start-------------------- mean latent time series
        fig_flag=fig_flags{14,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_mlatent_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                
                msl=mean_latent;
                p2=polyfit(1:length(msl),msl,1);
                p2=p2(1)*12.0;

                p=polyfit(xData,msl,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl,'k');
                hold on
                mslplot2=plot(xData,msl2,'Color','r');
                xlabel('year')
                ylabel('Mean latent (W/m^2)')
    %             mean_trend=mean(trend(:), 'omitnan');
                title([regionname, ', latent(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' W/m^2/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev_shflux)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end
        
% start-------------------- mean sensible time series
        fig_flag=fig_flags{15,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_msensible_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                
                msl=mean_sensible;
                p2=polyfit(1:length(msl),msl,1);
                p2=p2(1)*12.0;

                p=polyfit(xData,msl,1);
                msl2=xData*p(1)+p(2);
                mslplot=plot(xData,msl,'k');
                hold on
                mslplot2=plot(xData,msl2,'Color','r');
                xlabel('year')
                ylabel('Mean sensible (W/m^2)')
    %             mean_trend=mean(trend(:), 'omitnan');
                title([regionname, ', sensible(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(p2(1),2)), ' W/m^2/y'])
                datetick('x','yyyy','keepticks')
                axis tight;
                ylim(meanplotlev_shflux)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- seasonal wind speed time series
        fig_flag=fig_flags{16,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\speed\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_mwind_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                mean_seasonal_data=reshape(mean_speed, [12, length(mean_speed)/12]);
                mean_seasonal_data=mean(mean_seasonal_data,2);

                mslplot=plot(1:12,mean_seasonal_data,'k');
                hold on
                xlabel('month')
                ylabel('wind (m/s)')
                title([regionname, ', wind speed (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev_wind)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end
        
% start-------------------- seasonal shflux time series
        fig_flag=fig_flags{17,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\shflux\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_mshflux_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                mean_seasonal_data=reshape(mean_shflux, [12, length(mean_shflux)/12]);
                mean_seasonal_data=mean(mean_seasonal_data,2);

                mslplot=plot(1:12,mean_seasonal_data,'k');
                hold on
                xlabel('month')
                ylabel('shflux (W/m^2)')
                title([regionname, ', shflux (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev_shflux)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end
        
% start-------------------- seasonal latent time series
        fig_flag=fig_flags{18,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\latent\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_mlatent_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                mean_seasonal_data=reshape(mean_latent, [12, length(mean_latent)/12]);
                mean_seasonal_data=mean(mean_seasonal_data,2);

                mslplot=plot(1:12,mean_seasonal_data,'k');
                hold on
                xlabel('month')
                ylabel('latent (W/m^2)')
                title([regionname, ', latent (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev_shflux)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- yearly wind speed time series
        fig_flag=fig_flags{19,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_mwind_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                mean_seasonal_data=reshape(mean_speed, [12, length(mean_speed)/12]);
                mean_seasonal_data=mean(mean_seasonal_data,1);

                mslplot=plot(inputyear,mean_seasonal_data,'k');
                hold on
                xlabel('year')
                ylabel('wind (m/s)')
                title([regionname, ', wind speed (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev_yearly_wind)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- yearly shflux time series
        fig_flag=fig_flags{20,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_mshflux_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                mean_seasonal_data=reshape(mean_shflux, [12, length(mean_shflux)/12]);
                mean_seasonal_data=mean(mean_seasonal_data,1);

                mslplot=plot(inputyear,mean_seasonal_data,'k');
                hold on
                xlabel('year')
                ylabel('shflux (W/m^2)')
                title([regionname, ', shflux (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev_yearly_shflux)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- yearly latent time series
        fig_flag=fig_flags{21,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_mlatent_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                mean_seasonal_data=reshape(mean_latent, [12, length(mean_latent)/12]);
                mean_seasonal_data=mean(mean_seasonal_data,1);

                mslplot=plot(inputyear,mean_seasonal_data,'k');
                hold on
                xlabel('year')
                ylabel('latent (W/m^2)')
                title([regionname, ', latent (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                axis tight;
                ylim(meanplotlev_yearly_latent)
                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);
                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
            end
            fig_flag=0;
        end

% start-------------------- monthly wind speed interannual time series
        fig_flag=fig_flags{22,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\speed\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_monthly_mwind_', num2str(monthij,'%02i'), '_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    mean_seasonal_data2=reshape(mean_speed, [12, length(mean_speed)/12]);
                    mean_seasonal_data=mean_seasonal_data2(monthij,:); 
                    mslplot=plot(inputyear,mean_seasonal_data,'k');
                    hold on
                    xlabel(['year, ', num2str(monthij, '%02i'), 'm'])
                    ylabel('wind (m/s)')
                    title([regionname, ', wind speed (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                    axis tight;
                    ylim(meanplotlev_monthly_wind)
                    set(mslplot,'LineWidth',2);
                    set(gca,'FontSize',20);
                    grid on
                    hold off
                    saveas(gcf,jpgname,'jpg');
                    grid off
                    close all;
                end
            end
            fig_flag=0;
        end
        
% start-------------------- monthly shflux interannual time series
        fig_flag=fig_flags{23,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\shflux\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_monthly_mshflux_', num2str(monthij,'%02i'), '_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg

                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    mean_seasonal_data2=reshape(mean_shflux, [12, length(mean_shflux)/12]);
                    mean_seasonal_data=mean_seasonal_data2(monthij,:); 
                    mslplot=plot(inputyear,mean_seasonal_data,'k');
                    hold on
                    xlabel(['year, ', num2str(monthij, '%02i'), 'm'])
                    ylabel('shflux (W/m^2)')
                    title([regionname, ', shflux (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                    axis tight;
                    ylim(meanplotlev_monthly_shflux)
                    set(mslplot,'LineWidth',2);
                    set(gca,'FontSize',20);
                    grid on
                    hold off
                    saveas(gcf,jpgname,'jpg');
                    grid off
                    close all;
                end
            end
            fig_flag=0;
        end
        
% start-------------------- climatological speed std plot
        fig_flag=fig_flags{24,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\speed\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_std_speed_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    lon_rho=lon; lat_rho=lat;
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    comb_spatial_model_speed=reshape(comb_speed, [size(comb_speed,1) size(comb_speed,2) 12 size(comb_speed,3)/12]);
                    for loni=1:size(comb_speed,1)
                        for lati=1:size(comb_speed,2)
                            std_clim_speed(loni,lati,monthij)=std(comb_spatial_model_speed(loni,lati,monthij,:));
                        end
                    end
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(std_clim_speed(:,:,monthij)'));
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                    titlename = strcat('std speed, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                        ', ', calendarname{monthij}(1:3), ') ');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
    %                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
                    caxis(stdspeedlev);

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
            fig_flag=0;
        end
        
% start-------------------- climatological speed ratio based on std plot
        fig_flag=fig_flags{25,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\speed\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_ratio_std_speed_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    lon_rho=lon; lat_rho=lat;
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    comb_spatial_model_speed=reshape(comb_speed, [size(comb_speed,1) size(comb_speed,2) 12 size(comb_speed,3)/12]);
                    for loni=1:size(comb_speed,1)
                        for lati=1:size(comb_speed,2)
                            std_clim_speed(loni,lati,monthij)=std(comb_spatial_model_speed(loni,lati,monthij,:));
                        end
                    end
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(comb_spatial_meanmodel_speed(:,:,monthij)' ./ std_clim_speed(:,:,monthij)'));
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                    titlename = strcat('std speed, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                        ', ', calendarname{monthij}(1:3), ') ');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
    %                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
%                     caxis(stdspeedlev);

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
            fig_flag=0;
        end
        
% start-------------------- climatological shflux std plot
        fig_flag=fig_flags{26,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\shflux\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_std_shflux_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    lon_rho=lon; lat_rho=lat;
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    comb_spatial_model_shflux=reshape(comb_shflux, [size(comb_shflux,1) size(comb_shflux,2) 12 size(comb_shflux,3)/12]);
                    for loni=1:size(comb_shflux,1)
                        for lati=1:size(comb_shflux,2)
                            std_clim_shflux(loni,lati,monthij)=std(comb_spatial_model_shflux(loni,lati,monthij,:));
                        end
                    end
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(std_clim_shflux(:,:,monthij)'));
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                    titlename = strcat('std shflux, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                        ', ', calendarname{monthij}(1:3), ') ');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
    %                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
                    caxis(stdshfluxlev);

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
            fig_flag=0;
        end
        
% start-------------------- climatological shflux ratio based on std plot
        fig_flag=fig_flags{27,2};
        while (fig_flag)
            climdir = [figdir,'\CLIM\shflux\'];
            if (exist(strcat(climdir) , 'dir') ~= 7)
                mkdir(strcat(climdir));
            end 
            climoutfile = strcat(climdir,regionname);
            for monthij = 1:12
                jpgname=strcat(climoutfile, '_', testname,'_',regionname, '_ratio_std_shflux_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_', num2str(monthij, '%02i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    lon_rho=lon; lat_rho=lat;
                    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                    hold on;
                    comb_spatial_model_shflux=reshape(comb_shflux, [size(comb_shflux,1) size(comb_shflux,2) 12 size(comb_shflux,3)/12]);
                    for loni=1:size(comb_shflux,1)
                        for lati=1:size(comb_shflux,2)
                            std_clim_shflux(loni,lati,monthij)=std(comb_spatial_model_shflux(loni,lati,monthij,:));
                        end
                    end
                    m_pcolor(double(lon_rho)',lat_rho',squeeze(comb_spatial_meanmodel_shflux(:,:,monthij)' ./ std_clim_shflux(:,:,monthij)'));
                    shading(gca,m_pcolor_shading_method);
                    m_gshhs_i('color',m_gshhs_line_color);
                    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
            %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                    titlename = strcat('std shflux, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                        ', ', calendarname{monthij}(1:3), ') ');  %% + glacier contribution

                    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                    % set colorbar 
                    h = colorbar;
                    colormap(byrmap);
                    set(h,'fontsize',colorbar_fontsize);
    %                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
%                     caxis(stdshfluxlev);

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
            fig_flag=0;
        end
        
% start-------------------- speed std plot
        fig_flag=fig_flags{28,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_std_speed_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                lon_rho=lon; lat_rho=lat;
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                for loni=1:size(comb_speed,1)
                    for lati=1:size(comb_speed,2)
                        std_speed(loni,lati)=std(comb_speed(loni,lati,:));
                    end
                end
                m_pcolor(double(lon_rho)',lat_rho',squeeze(std_speed(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('std speed, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ') '); 

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(byrmap);
                set(h,'fontsize',colorbar_fontsize);
%                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
%                 caxis(stdspeedlev2);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['speed std', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            close all;
            fig_flag=0;
        end

% start-------------------- shflux std plot
        fig_flag=fig_flags{29,2};
        while (fig_flag)
            jpgname=strcat(outfile, '_', testname,'_',regionname, '_std_shflux_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                lon_rho=lon; lat_rho=lat;
                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                hold on;
                for loni=1:size(comb_shflux,1)
                    for lati=1:size(comb_shflux,2)
                        std_shflux(loni,lati)=std(comb_shflux(loni,lati,:));
                    end
                end
                m_pcolor(double(lon_rho)',lat_rho',squeeze(std_shflux(:,:)'));
                shading(gca,m_pcolor_shading_method);
                m_gshhs_i('color',m_gshhs_line_color);
                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        %         titlename = strcat(regionname,', SSH trend, ','Mean=',num2str(round(mean_trend_filtered,2)),'mm/y');
                titlename = strcat('std shflux, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ') '); 

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(byrmap);
                set(h,'fontsize',colorbar_fontsize);
%                 title(h,'mm/y','fontsize',colorbar_title_fontsize);
%                 caxis(stdshfluxlev2);

                % set grid
                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp(['shflux std', ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')

                hold off
                close all;
            end
            close all;
            fig_flag=0;
        end
        
% start-------------------- mean atmvar and energies time series (corr) plot
        fig_flag=fig_flags{30,2};
        while (fig_flag)
            for var_atmi=1:length(variable_atm)
                var_atmname=variable_atm{var_atmi};
                for var_ei=1:length(variable_energy)
                    var_ename=variable_energy{var_ei};
                    for depthi=1:length(ana_stddepth)
                        tempdepth=ana_stddepth(depthi);
                        tempdepthstr=num2str(-tempdepth, '%04i');
                        efigdir=[figrawdir,'ATM\', var_atmname, '\', var_ename, '\', tempdepthstr];
                        if (exist(strcat(efigdir) , 'dir') ~= 7)
                            mkdir(strcat(efigdir));
                        end
                        eoutfile = [efigdir, '\', regionname];
                        jpgname=strcat(eoutfile, '_', testname, '_ts_', var_atmname, '_', var_ename, '_', tempdepthstr, 'm_',  ...
                            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                        if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                            clear mslplot
                            lon_rho=lon; lat_rho=lat;
                            corrdir=[filedir,'\corr'];
                            efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                            '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                            var_e=ncread(efilename, 'var_e');
                            var_e_2d=reshape(var_e, [len_lon_model*len_lat_model, size(var_e,3)]);
                            mean_var_e=mean(var_e_2d,1,'omitnan');
                            
                            axLH = gca;
                            axRH = axes('color','none');
                            eval(['mslplot{1}=plot(xData,mean_', var_atmname ', ''b'',''parent'',axLH);']);
                            mslplot{2}=plot(xData,mean_var_e, 'k','parent',axRH);
                            ylabel(axLH,var_atmname)
                            ylabel(axRH,[var_ename, '(J/m^3)'])
                            ax_pos = get(axLH,'position');
                            set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
                            set(axRH,'color','none','yaxislocation','right','xtick', xData(1:12:end), 'position', ax_pos+[0 0.02 -0.01 -0.02]);
                %             set(axRH,'color','none','yaxislocation','right');

                %             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
                            set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
                            set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
                            set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
                            xlabel(axRH, 'YearMonth');

                            title([var_atmname, ',', var_ename, ',', tempdepthstr, 'm,', num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')]);
                            datetick(axLH, 'x','yymm','keepticks');
                            datetick(axRH, 'x','yymm','keepticks');
                            set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');

%                             axis tight;
                            % ylim(meanplotlev2)
                            set(mslplot{1},'LineWidth',2);
                            set(mslplot{2},'LineWidth',2);
                            grid on

                            lgd=legend([mslplot{1} mslplot{2}], var_atmname, var_ename);
                            eval(['tempcorr=corrcoef(mean_', var_atmname ',mean_var_e);'])
                            
                            txt1=text(xData(5), min(double(mean_var_e))+diff(double([min(mean_var_e), max(mean_var_e)]))/64.0 ,['R = ', num2str(round(tempcorr(1,2),2))], 'FontSize', 20); 
                            
                %             lgd=legend('Model','TG-UST');
                            set(lgd,'FontSize',15);
                            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                            set(lgd,'Orientation','horizontal');

                            set(gcf,'PaperPosition', [0 0 36 12]) 
                            saveas(gcf,jpgname,'tif');
                            grid off
                            close all;
                            
                        end
                    end
                end
            end
            close all;
            fig_flag=0;
        end

% start-------------------- atmvar corr plot  ATM-var-e-depth
        fig_flag=fig_flags{31,2};
        while (fig_flag)
            for var_atmi=1:length(variable_atm)
                var_atmname=variable_atm{var_atmi};
                for var_ei=1:length(variable_energy)
                    var_ename=variable_energy{var_ei};
                    for depthi=1:length(ana_stddepth)
                        tempdepth=ana_stddepth(depthi);
                        tempdepthstr=num2str(-tempdepth, '%04i');
                        efigdir=[figrawdir,'ATM\', var_atmname, '\', var_ename, '\', tempdepthstr];
                        if (exist(strcat(efigdir) , 'dir') ~= 7)
                            mkdir(strcat(efigdir));
                        end
                        eoutfile = [efigdir, '\', regionname];
                        jpgname=strcat(eoutfile, '_', testname, '_corr_', var_atmname, '_', var_ename, '_', tempdepthstr, 'm_',  ...
                            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                        if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                            lon_rho=lon; lat_rho=lat;
                            corrdir=[filedir,'\corr'];
                            efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                            '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                            
                            temp_corr=ncread(efilename, 'corr_var');
                            mean_temp_corr=mean(temp_corr(:),'omitnan');
                            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                            hold on;
                            m_pcolor(double(lon_rho)',lat_rho',squeeze(temp_corr(:,:)'));
                            shading(gca,m_pcolor_shading_method);
                            m_gshhs_i('color',m_gshhs_line_color);
                            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
  
                            title(['corr,', var_atmname, ',', var_ename, ',', tempdepthstr, 'm,', ...
                                num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                                ',M=',num2str(round(mean_temp_corr,2))],'fontsize',m_pcolor_title_fontsize-5)

                            % set colorbar 
                            h = colorbar;
                            colormap(byrmap);
                            set(h,'fontsize',colorbar_fontsize);
                %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                            caxis([-1 1]);

                            % set grid
                            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                            set(gcf, 'PaperUnits', 'points');
                            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                            saveas(gcf,jpgname,'tif');

                            disp(' ')
                            disp([ fig_flags{31,1}, ' plot is created.'])
                            disp(' ')
                            disp([' File path is : ',jpgname])
                            disp(' ')

                            hold off
                            close all;
                        
                        end
                    end
                end
            end
            close all;
            fig_flag=0;
        end
            
% start-------------------- atmvar seasonal corr plot
        fig_flag=fig_flags{32,2};
        while (fig_flag)
            for var_atmi=1:length(variable_atm)
                var_atmname=variable_atm{var_atmi};
                for var_ei=1:length(variable_energy)
                    var_ename=variable_energy{var_ei};
                    for depthi=1:length(ana_stddepth)
                        tempdepth=ana_stddepth(depthi);
                        tempdepthstr=num2str(-tempdepth, '%04i');
                        efigdir=[figrawdir,'ATM\', var_atmname, '\', var_ename, '\', tempdepthstr];
                        if (exist(strcat(efigdir) , 'dir') ~= 7)
                            mkdir(strcat(efigdir));
                        end
                        eoutfile = [efigdir, '\', regionname];
                        jpgname=strcat(eoutfile, '_', testname, '_season_corr_', var_atmname, '_', var_ename, '_', tempdepthstr, 'm_',  ...
                            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                        if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                            lon_rho=lon; lat_rho=lat;
                            corrdir=[filedir,'\corr'];
                            efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                            '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                            
                            temp_corr=ncread(efilename, 'corr_spatial_mean');
                            mean_temp_corr=mean(temp_corr(:),'omitnan');
                            m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                            hold on;
                            m_pcolor(double(lon_rho)',lat_rho',squeeze(temp_corr(:,:)'));
                            shading(gca,m_pcolor_shading_method);
                            m_gshhs_i('color',m_gshhs_line_color);
                            m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
  
                            title(['season_r,', var_atmname, ',', var_ename, ',', tempdepthstr, 'm,', ...
                                num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                                ',M=',num2str(round(mean_temp_corr,2))],'fontsize',m_pcolor_title_fontsize-5)

                            % set colorbar 
                            h = colorbar;
                            colormap(byrmap);
                            set(h,'fontsize',colorbar_fontsize);
                %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                            caxis([-1 1]);

                            % set grid
                            m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                            set(gcf, 'PaperUnits', 'points');
                            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                            saveas(gcf,jpgname,'tif');

                            disp(' ')
                            disp([ fig_flags{31,1}, ' plot is created.'])
                            disp(' ')
                            disp([' File path is : ',jpgname])
                            disp(' ')

                            hold off
                            close all;
                        
                        end
                    end
                end
            end
            close all;
            fig_flag=0;
        end
        
% start-------------------- atmvar clim corr plot
        fig_flag=fig_flags{33,2};
        while (fig_flag)
            for var_atmi=1:length(variable_atm)
                var_atmname=variable_atm{var_atmi};
                for var_ei=1:length(variable_energy)
                    var_ename=variable_energy{var_ei};
                    for depthi=1:length(ana_stddepth)
                        tempdepth=ana_stddepth(depthi);
                        tempdepthstr=num2str(-tempdepth, '%04i');
                        efigdir=[figrawdir,'ATM\', var_atmname, '\', var_ename, '\', tempdepthstr];
                        if (exist(strcat(efigdir) , 'dir') ~= 7)
                            mkdir(strcat(efigdir));
                        end
                        eoutfile = [efigdir, '\', regionname];
                        for climi=1:length(inputmonth)
                            templag=inputmonth(climi);
                            jpgname=strcat(eoutfile, '_', testname, '_clim_corr_', num2str(templag, '%03i'), 'm_', var_atmname, '_', var_ename, '_', tempdepthstr, 'm_',  ...
                                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                                lon_rho=lon; lat_rho=lat;
                                corrdir=[filedir,'\corr'];
                                efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                                '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                                
                                temp_corr=ncread(efilename, 'corr_clim', [1 1 climi], [inf inf 1]);
%                                 temp_corr=clim_corr(:,:,climi);
                                mean_temp_corr=mean(temp_corr(:),'omitnan');
                                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                                hold on;
                                m_pcolor(double(lon_rho)',lat_rho',squeeze(temp_corr(:,:)'));
                                shading(gca,m_pcolor_shading_method);
                                m_gshhs_i('color',m_gshhs_line_color);
                                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                                title([calendarname{climi}(1:3),', clim r,', var_atmname, ',', var_ename, ',', tempdepthstr, 'm,', ...
                                    num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                                    ',M=',num2str(round(mean_temp_corr,2))],'fontsize',m_pcolor_title_fontsize-5)
                                calendarname{climi}(1:3)
                                % set colorbar 
                                h = colorbar;
                                colormap(byrmap);
                                set(h,'fontsize',colorbar_fontsize);
                    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                                caxis([-1 1]);

                                % set grid
                                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                                set(gcf, 'PaperUnits', 'points');
                                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                                saveas(gcf,jpgname,'tif');

                                disp(' ')
                                disp([ fig_flags{31,1}, ' plot is created.'])
                                disp(' ')
                                disp([' File path is : ',jpgname])
                                disp(' ')

                                hold off
                                close all;
                            end
                        end
                    end
                end
            end
            close all;
            fig_flag=0;
        end
        
% start-------------------- atmvar lag corr plot
        fig_flag=fig_flags{34,2};
        while (fig_flag)
            for var_atmi=1:length(variable_atm)
                var_atmname=variable_atm{var_atmi};
                for var_ei=1:length(variable_energy)
                    var_ename=variable_energy{var_ei};
                    for depthi=1:length(ana_stddepth)
                        tempdepth=ana_stddepth(depthi);
                        tempdepthstr=num2str(-tempdepth, '%04i');
                        efigdir=[figrawdir,'ATM\', var_atmname, '\', var_ename, '\', tempdepthstr];
                        if (exist(strcat(efigdir) , 'dir') ~= 7)
                            mkdir(strcat(efigdir));
                        end
                        eoutfile = [efigdir, '\', regionname];
                        for lagi=1:length(lag_month)
                            templag=lag_month(lagi);
                            jpgname=strcat(eoutfile, '_', testname, '_lag_corr_', num2str(templag, '%03i'), 'm_', var_atmname, '_', var_ename, '_', tempdepthstr, 'm_',  ...
                                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                                lon_rho=lon; lat_rho=lat;
                                corrdir=[filedir,'\corr'];
                                efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                                '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                                
                                lag_corr=ncread(efilename, 'corr_lag');
                                temp_corr=lag_corr(:,:,lagi);
                                mean_temp_corr=mean(temp_corr(:),'omitnan');
                                m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
                                hold on;
                                m_pcolor(double(lon_rho)',lat_rho',squeeze(temp_corr(:,:)'));
                                shading(gca,m_pcolor_shading_method);
                                m_gshhs_i('color',m_gshhs_line_color);
                                m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land

                                title([num2str(templag, '%03i'), 'm,','lag_r,', var_atmname, ',', var_ename, ',', tempdepthstr, 'm,', ...
                                    num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
                                    ',M=',num2str(round(mean_temp_corr,2))],'fontsize',m_pcolor_title_fontsize-5)

                                % set colorbar 
                                h = colorbar;
                                colormap(byrmap);
                                set(h,'fontsize',colorbar_fontsize);
                    %             title(h,'mm/y','fontsize',colorbar_title_fontsize);
                                caxis([-1 1]);

                                % set grid
                                m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

                                set(gcf, 'PaperUnits', 'points');
                                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                                saveas(gcf,jpgname,'tif');

                                disp(' ')
                                disp([ fig_flags{31,1}, ' plot is created.'])
                                disp(' ')
                                disp([' File path is : ',jpgname])
                                disp(' ')

                                hold off
                                close all;
                            end
                        end
                    end
                end
            end
            close all;
            fig_flag=0;
        end

% start-------------------- mean atmvar lag corr plot
        fig_flag=fig_flags{35,2};
        while (fig_flag)
            for var_atmi=1:length(variable_atm)
                var_atmname=variable_atm{var_atmi};
                for var_ei=1:length(variable_energy)
                    var_ename=variable_energy{var_ei};
                    for depthi=1:length(ana_stddepth)
                        tempdepth=ana_stddepth(depthi);
                        tempdepthstr=num2str(-tempdepth, '%04i');
                        efigdir=[figrawdir,'ATM\', var_atmname, '\', var_ename, '\', tempdepthstr];
                        if (exist(strcat(efigdir) , 'dir') ~= 7)
                            mkdir(strcat(efigdir));
                        end
                        eoutfile = [efigdir, '\', regionname];
%                         templag=lag_month(lagi);
                        jpgname=strcat(eoutfile, '_', testname, '_plot_lag_corr_', var_atmname, '_', var_ename, '_', tempdepthstr, 'm_',  ...
                            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                        if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                            lon_rho=lon; lat_rho=lat;
                            corrdir=[filedir,'\corr'];
                            efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                            '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

                            lag_corr=ncread(efilename, 'corr_lag');
                            corr_var=ncread(efilename, 'corr_var');
                            clear mean_temp_corr
                            for lagi=1:length(lag_month)
                                templag=lag_month(lagi);
                                temp_corr=lag_corr(:,:,lagi);
                                mean_temp_corr(lagi)=mean(temp_corr(:),'omitnan');
                            end
                            mean_temp_corr(end+1)=mean(corr_var(:), 'omitnan');
                            temp_lag_month=[lag_month, 0];
                            [new_lag_month, sortind]=sort(temp_lag_month);
                            mean_temp_corr2=mean_temp_corr(sortind);
                            
                            mslplot=plot(new_lag_month,mean_temp_corr2,'k');
                            hold on
                            xlabel('lag month')
                            xticks(new_lag_month)
                            ylabel('r')
                            title(['lag, r,', var_atmname, ',', var_ename, ',', tempdepthstr, 'm,', ...
                                    num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')], ...
                                    'fontsize',m_pcolor_title_fontsize-5)

%                             title([regionname, ', lag corr (',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
                            axis tight;
                            ylim([-1 1])
                            set(mslplot,'LineWidth',2);
                            set(gca,'FontSize',20);
                            grid on
                            hold off
                            saveas(gcf,jpgname,'jpg');
                            grid off
                            close all;
                        end
                    end
                end
            end
            close all;
            fig_flag=0;
        end

% start-------------------- clim mean atmvar and energies time series (corr) plot
        fig_flag=fig_flags{36,2};
        while (fig_flag)
            for var_atmi=1:length(variable_atm)
                var_atmname=variable_atm{var_atmi};
                for var_ei=1:length(variable_energy)
                    var_ename=variable_energy{var_ei};
                    for depthi=1:length(ana_stddepth)
                        tempdepth=ana_stddepth(depthi);
                        tempdepthstr=num2str(-tempdepth, '%04i');
                        efigdir=[figrawdir,'ATM\', var_atmname, '\', var_ename, '\', tempdepthstr];
                        if (exist(strcat(efigdir) , 'dir') ~= 7)
                            mkdir(strcat(efigdir));
                        end
                        eoutfile = [efigdir, '\', regionname];
                        for climi=1:length(inputmonth)
                            tempmonth=inputmonth(climi);
                            monthstr=num2str(tempmonth,'%02i');
                            jpgname=strcat(eoutfile, '_', testname, '_clim_', monthstr, 'm_ts_', var_atmname, '_', var_ename, '_', tempdepthstr, 'm_',  ...
                                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                                lon_rho=lon; lat_rho=lat;
                                corrdir=[filedir,'\corr'];
                                efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                                '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                                    
                                var_e=ncread(efilename, 'var_e');
                                var_e_3d=reshape(var_e, [len_lon_model*len_lat_model, 12, size(var_e,3)/12]);
                                
                                mean_var_e=squeeze(mean(var_e_3d(:,tempmonth,:),1,'omitnan'));
                                eval(['temp_atmvar_3d=reshape(comb_', var_atmname ', [len_lon_model*len_lat_model, 12, size(var_e,3)/12]);'])
                                mean_temp_atmvar=squeeze(mean(temp_atmvar_3d(:,tempmonth,:),1,'omitnan'));
                                
                                clear mslplot
                                axLH = gca;
                                axRH = axes('color','none');
%                                 eval(['mslplot{1}=plot(xData,mean_', var_atmname ', ''b'',''parent'',axLH);']);
                                mslplot{1}=plot(xData2,mean_temp_atmvar, 'b','parent',axLH);
                                mslplot{2}=plot(xData2,mean_var_e, 'k','parent',axRH);
                                ylabel(axLH,var_atmname)
                                ylabel(axRH,[var_ename, '(J/m^3)'])
                                ax_pos = get(axLH,'position');
                                set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
                                set(axRH,'color','none','yaxislocation','right','xtick', xData(1:12:end), 'position', ax_pos+[0 0.02 -0.01 -0.02]);
                    %             set(axRH,'color','none','yaxislocation','right');

                    %             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
                                set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
                                set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
                                set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
                                xlabel(axRH, 'YearMonth');

                                title([calendarname{climi}(1:3), ', ', var_atmname, ',', var_ename, ',', tempdepthstr, 'm,', num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')]);
                                datetick(axLH, 'x','yymm','keepticks');
                                datetick(axRH, 'x','yymm','keepticks');
                                set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');

    %                             axis tight;
                                % ylim(meanplotlev2)
                                set(mslplot{1},'LineWidth',2);
                                set(mslplot{2},'LineWidth',2);
                                grid on

                                lgd=legend([mslplot{1} mslplot{2}], var_atmname, var_ename);
                                eval(['tempcorr=corrcoef(mean_temp_atmvar,mean_var_e);'])

                                txt1=text(xData(5), min(double(mean_var_e))+diff(double([min(mean_var_e), max(mean_var_e)]))/64.0 ,['R = ', num2str(round(tempcorr(1,2),2))], 'FontSize', 20); 

                    %             lgd=legend('Model','TG-UST');
                                set(lgd,'FontSize',15);
                                set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                                set(lgd,'Orientation','horizontal');

                                set(gcf,'PaperPosition', [0 0 36 12]) 
                                saveas(gcf,jpgname,'tif');
                                grid off
                                close all;
                            end
                        end
                    end
                end
            end
            close all;
            fig_flag=0;
        end
        
% start-------------------- clim mean atmvar(wind speed, shflux) time series (corr) plot
        fig_flag=fig_flags{37,2};
        while (fig_flag)
            var_atmname1=variable_atm{1};
            var_atmname2=variable_atm{2};
            efigdir=[figrawdir,'ATM\CLIM\'];
            if (exist(strcat(efigdir) , 'dir') ~= 7)
                mkdir(strcat(efigdir));
            end
            eoutfile = [efigdir, '\', regionname];
            for climi=1:length(inputmonth)
                tempmonth=inputmonth(climi);
                monthstr=num2str(tempmonth,'%02i');
                jpgname=strcat(eoutfile, '_', testname, '_clim_', monthstr, 'm_ts_', var_atmname1, '_', var_atmname2, '_',  ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                    lon_rho=lon; lat_rho=lat;
                    corrdir=[filedir,'\corr'];
%                     efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
%                     '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
% 
%                     var_e=ncread(efilename, 'var_e');
%                     var_e_3d=reshape(var_e, [len_lon_model*len_lat_model, 12, size(var_e,3)/12]);
% 
%                     mean_var_e=squeeze(mean(var_e_3d(:,tempmonth,:),1,'omitnan'));
                    eval(['temp_atmvar_3d=reshape(comb_', var_atmname1 ', [len_lon_model*len_lat_model, 12, size(comb_', var_atmname1 ',3)/12]);'])
                    mean_temp_atmvar=squeeze(mean(temp_atmvar_3d(:,tempmonth,:),1,'omitnan'));
                    eval(['temp_atmvar_3d=reshape(comb_', var_atmname2 ', [len_lon_model*len_lat_model, 12, size(comb_', var_atmname2 ',3)/12]);'])
                    mean_temp_atmvar2=squeeze(mean(temp_atmvar_3d(:,tempmonth,:),1,'omitnan'));

                    clear mslplot
                    axLH = gca;
                    axRH = axes('color','none');
%                                 eval(['mslplot{1}=plot(xData,mean_', var_atmname ', ''b'',''parent'',axLH);']);
                    mslplot{1}=plot(xData2,mean_temp_atmvar, 'b','parent',axLH);
                    mslplot{2}=plot(xData2,mean_temp_atmvar2, 'k','parent',axRH);
                    ylabel(axLH,var_atmname1)
                    ylabel(axRH,[var_atmname2])
                    ax_pos = get(axLH,'position');
                    set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
                    set(axRH,'color','none','yaxislocation','right','xtick', xData(1:12:end), 'position', ax_pos+[0 0.02 -0.01 -0.02]);
        %             set(axRH,'color','none','yaxislocation','right');

        %             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
                    set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
                    set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
                    set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
                    xlabel(axRH, 'YearMonth');

                    title([calendarname{climi}(1:3), ', ', var_atmname1, ',', var_atmname2, ',', num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')]);
                    datetick(axLH, 'x','yymm','keepticks');
                    datetick(axRH, 'x','yymm','keepticks');
                    set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');

%                             axis tight;
                    % ylim(meanplotlev2)
                    set(mslplot{1},'LineWidth',2);
                    set(mslplot{2},'LineWidth',2);
                    grid on

                    lgd=legend([mslplot{1} mslplot{2}], var_atmname1, var_atmname2);
                    eval(['tempcorr=corrcoef(mean_temp_atmvar,mean_temp_atmvar2);'])

                    txt1=text(xData(5), min(double(mean_temp_atmvar2))+diff(double([min(mean_temp_atmvar2), max(mean_temp_atmvar2)]))/64.0 ,['R = ', num2str(round(tempcorr(1,2),2))], 'FontSize', 20); 

        %             lgd=legend('Model','TG-UST');
                    set(lgd,'FontSize',15);
                    set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                    set(lgd,'Orientation','horizontal');

                    set(gcf,'PaperPosition', [0 0 36 12]) 
                    saveas(gcf,jpgname,'tif');
                    grid off
                    close all;
                end
            end
            close all;
            fig_flag=0;
        end

% start-------------------- clim mean atmvar and energies corr profile plot
        fig_flag=fig_flags{38,2};
        while (fig_flag)
            for var_atmi=1:length(variable_atm)
                var_atmname=variable_atm{var_atmi};
                
                for climi=1:length(inputmonth)
                    tempmonth=inputmonth(climi);
                    monthstr=num2str(tempmonth,'%02i');
                    efigdir=[figrawdir,'ATM\', var_atmname, '\', 'all_depth'];
%                     if (exist(strcat(efigdir) , 'dir') ~= 7)
                        mkdir(strcat(efigdir));
%                     end
                    eoutfile = [efigdir, '\', regionname];
                    jpgname=strcat(eoutfile, '_', testname, '_clim_', monthstr, 'm_ts_', var_atmname, '_', 'all_E', '_',  ...
                                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
                    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                        for var_ei=1:length(variable_energy)
                            var_ename=variable_energy{var_ei};
                            for depthi=1:length(ana_stddepth)
                                tempdepth=ana_stddepth(depthi);
                                tempdepthstr=num2str(-tempdepth, '%04i');

                                lon_rho=lon; lat_rho=lat;
                                corrdir=[filedir,'\corr'];
                                efilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                                '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                                
%                                 disp('abc')
%                                 var_e=ncread(efilename, 'var_e');
%                                 var_e_3d=reshape(var_e, [len_lon_model*len_lat_model, 12, size(var_e,3)/12]);
%                                 
%                                 mean_var_e=squeeze(mean(var_e_3d(:,tempmonth,:),1,'omitnan'));
%                                 eval(['temp_atmvar_3d=reshape(comb_', var_atmname ', [len_lon_model*len_lat_model, 12, size(var_e,3)/12]);'])
%                                 mean_temp_atmvar=squeeze(mean(temp_atmvar_3d(:,tempmonth,:),1,'omitnan'));
                                
                                temp_corr_clim= ncread(efilename, 'corr_clim', [1 1 climi], [inf inf 1]);
                                mean_temp_corr_clim=mean(temp_corr_clim(:), 'omitnan');
                                corr_clim_depth(var_ei, depthi) = mean_temp_corr_clim;
                            end
                        end
                        clear mslplot
                        mslplot{1}=plot(corr_clim_depth(1,:), stddepth);
                        hold on
                        mslplot{2}=plot(corr_clim_depth(2,:), stddepth);
                        mslplot{3}=plot(corr_clim_depth(3,:), stddepth);
                        mslplot{4}=plot(corr_clim_depth(4,:), stddepth);
%                         axLH = gca;
%                         axRH = axes('color','none');
%                                 eval(['mslplot{1}=plot(xData,mean_', var_atmname ', ''b'',''parent'',axLH);']);
%                         mslplot{1}=plot(xData2,mean_temp_atmvar, 'b','parent',axLH);
%                         mslplot{2}=plot(xData2,mean_var_e, 'k','parent',axRH);
%                         ylabel(axLH,var_atmname)
                        ylabel('depth')
                        xlabel('r')
%                         ax_pos = get(axLH,'position');
%                         set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
%                         set(axRH,'color','none','yaxislocation','right','xtick', xData(1:12:end), 'position', ax_pos+[0 0.02 -0.01 -0.02]);
            %             set(axRH,'color','none','yaxislocation','right');

            %             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
%                         set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
%                         set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
%                         set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
%                         xlabel(axRH, 'YearMonth');

                        title([calendarname{climi}(1:3), ', ', var_atmname, ',', num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')]);
%                         datetick('x','yymm','keepticks');
%                         datetick(axRH, 'x','yymm','keepticks');
%                         set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');

%                             axis tight;
                        % ylim(meanplotlev2)
                        set(mslplot{1},'LineWidth',2);
                        set(mslplot{2},'LineWidth',2);
                        set(mslplot{3},'LineWidth',2);
                        set(mslplot{4},'LineWidth',2);
                        grid on

                        lgd=legend([mslplot{1} mslplot{2} mslplot{3} mslplot{4}],variable_energy{1}, variable_energy{2}, variable_energy{3}, variable_energy{4} );
%                         eval(['tempcorr=corrcoef(mean_temp_atmvar,mean_var_e);'])

%                         txt1=text(xData(5), min(double(mean_var_e))+diff(double([min(mean_var_e), max(mean_var_e)]))/64.0 ,['R = ', num2str(round(tempcorr(1,2),2))], 'FontSize', 20); 

            %             lgd=legend('Model','TG-UST');
                        set(lgd,'FontSize',15);
%                         set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                        set(lgd,'Position',[0.13 0.23, 0.1, 0.03]);
                        set(lgd,'Orientation','vertical');

                        set(gcf,'PaperPosition', [0 0 36 12]) 
                        saveas(gcf,jpgname,'tif');
                        grid off
                        close all;
                    end
                end
            end
            close all;
            fig_flag=0;
        end
        
        
    end
end