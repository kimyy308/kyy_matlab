close all; clear all;  clc;   
% %  get SST RMS, SST BIAS data and plot(climatological value time series)

% all_region ={'ES','SS', 'YS'}
% all_region ={'ES', 'SS', 'YS', 'ECS'}

% all_testname = {'test11', 'test12'};
% all_testname = {'test53'};
% all_testname = {'test53', 'test54', 'test55', 'test56'};
all_testname = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};

% all_testname = {'ens03'};
all_scenname = {'historical'};

% all_region ={'AKP', 'NWP', 'ES', 'SS', 'YS'};
all_region ={'AKP4'};

for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname all_scenname scenind
        scenname = all_scenname{1};
        % % % 
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
        rms_shadlev = [0 5];
        bias_shadlev = [-5 5];
        conlev  = 0:5:35;
        dl=1/20;
        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [1982:2005]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]


        varname ='temp';  %% reference variable -> temperature
        run('nwp_polygon_point.m');
        regionname=all_region{regionind}
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
            case('AKP3') %% Around Korea Peninsula
                refpolygon=akp3polygon;
            case('AKP4') %% Around Korea Peninsula
                refpolygon=akp4polygon;
            case('EKB') %% Around Korea Peninsula
                refpolygon=ekbpolygon;
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
        variable = 'SST';
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
%             filedir = strcat('J:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            filedir = strcat('D:\Data\Model\CMIP5\surface\'); % % where data files are  --> please check 181 line again
            savedir = filedir;
%             savedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            avhrrdir='Z:\내 드라이브\Data\Observation\OISST\monthly_kimyy\';
        elseif (strcmp(system_name,'GLNXA64'))
        end

        run(param_script);
        matname = [filedir,scenname, '_',testname,'_',regionname,'model_cmems_sst_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
        filename = strcat(filedir,  ...
                        scenname,'_', testname, '_NWP_CMIP5_sst_analysis_', num2str(inputyear(1),'%04i'), '_', num2str(inputyear(end),'%04i'), '.nc');
        comb_data = ncread(filename, 'raw_sst')-273.15;
        comb_spatial_mean_data=ncread(filename, 'clim_sst')-273.15;
        lon_rho=ncread(filename, 'lon_rho');
        lat_rho=ncread(filename, 'lat_rho');
        ind=1;
        
        for yearij = 1:length(inputyear)
            for monthij = 1:length(inputmonth)
                disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                tic;
                tempyear = inputyear(yearij);
                tempmonth = inputmonth(monthij);                

                % read OISST DATA
                avhrrfilename = strcat(avhrrdir,'avhrr_only_monthly_v2_', num2str(tempyear,'%04i'), '.nc');
                if (exist('avhrr_lon')==0)
                    avhrrinfo=ncinfo(avhrrfilename);
                    
                    avhrr_lon = ncread(avhrrfilename,'lon',[1],[avhrrinfo.Dimensions(2).Length]);
                    avhrr_lat = ncread(avhrrfilename,'lat',[1],[avhrrinfo.Dimensions(1).Length]);

                    avhrr_lon_west = abs(avhrr_lon - (lonlat(1)));
                    min_avhrr_lon_west=min(avhrr_lon_west);
                    avhrr_lon_east = abs(avhrr_lon - (lonlat(2)));
                    min_avhrr_lon_east=min(avhrr_lon_east);
                    avhrr_lat_south = abs(avhrr_lat - (lonlat(3)));
                    min_avhrr_lat_south=min(avhrr_lat_south);
                    avhrr_lat_north = abs(avhrr_lat - (lonlat(4)));
                    min_avhrr_lat_north=min(avhrr_lat_north);

                    avhrr_lon_min = find(avhrr_lon_west == min_avhrr_lon_west);
                    avhrr_lon_max = find(avhrr_lon_east == min_avhrr_lon_east);
                    avhrr_lat_min = find(avhrr_lat_south == min_avhrr_lat_south);
                    avhrr_lat_max = find(avhrr_lat_north == min_avhrr_lat_north);

            %         ncinfo('E:\Data\Observation\OISST\monthly\avhrr_monthly1983_11.nc');

                    avhrr_lon2 = ncread(avhrrfilename,'lon', [avhrr_lon_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1)]);
                    avhrr_lat2 = ncread(avhrrfilename,'lat', [avhrr_lat_min(1)], [avhrr_lat_max(1)-avhrr_lat_min(1)]);
                    [avhrr_lat, avhrr_lon]=meshgrid(avhrr_lat2, avhrr_lon2);
                    
                    comb_spatial_meanrms=(zeros([size(avhrr_lon),12]));
                    comb_spatial_meanbias=(zeros([size(avhrr_lon),12]));
                    comb_spatial_meanavhrr=(zeros([size(avhrr_lon),12]));
                    comb_spatial_meanmodel=(zeros([size(avhrr_lon),12]));

                    switch(regionname)
                        case('NWP') %% North western Pacific
                            mask_avhrr(1:size(avhrr_lon,1),1:size(avhrr_lon,2))=1;
                        otherwise
                            mask_avhrr = double(inpolygon(avhrr_lon,avhrr_lat,refpolygon(:,1),refpolygon(:,2)));
                            mask_avhrr(mask_avhrr==0)=NaN;
                    end
                end
                len_lon = length(avhrr_lon(:,1));
                len_lat = length(avhrr_lat(1,:));
                len_lon_model = size(comb_data,1);
                len_lat_model = size(comb_data,2);

                avhrr_data = ncread(avhrrfilename,varname,[avhrr_lon_min(1) avhrr_lat_min(1) monthij], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);
                avhrr_data=avhrr_data.*mask_avhrr;
                avhrr_err = ncread(avhrrfilename,'err',[avhrr_lon_min(1) avhrr_lat_min(1) monthij], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);
                avhrr_err=avhrr_err.*mask_avhrr;
                
%                 interped_data = griddata(double(lon), double(lat), data,double(avhrr_lon),double(avhrr_lat));   
                interped_data = griddata(double(lon_rho), double(lat_rho), squeeze(double(comb_data(:,:,ind))),double(avhrr_lon),double(avhrr_lat));   

%                 interped_data = griddata(double(lon_rho), double(lat_rho), squeeze(double(comb_data(:,:,ind))),double(cmems_lon2),double(cmems_lat2)).*mask_cmems;   

                bias = interped_data-avhrr_data;  
                mse = (interped_data-avhrr_data).^2;  
%                 meanbias = mean(mean(bias,'omitnan'),'omitnan');
%                 meanrms = mean(mean(rms,'omitnan'),'omitnan');
                meanbias = mean(bias(:),'omitnan');

                comb_interped_data(:,:,ind) = interped_data;
                comb_avhrr_data(:,:,ind) = avhrr_data;
                comb_avhrr_err(:,:,ind) = avhrr_err;
                comb_mse_data(:,:,ind) = mse;
                comb_bias_data(:,:,ind) = bias;

%                 comb_meanmse(yearij,monthij)=meanmse;
                comb_meanbias(yearij,monthij)=meanbias;
%                 comb_spatial_meanmse(:,:,monthij)=comb_spatial_meanmse(:,:,monthij)+mse/double(length(inputyear));
                comb_spatial_meanbias(:,:,monthij)=comb_spatial_meanbias(:,:,monthij)+bias/double(length(inputyear));
                comb_spatial_meanavhrr(:,:,monthij)=comb_spatial_meanavhrr(:,:,monthij)+avhrr_data/double(length(inputyear));
                comb_spatial_meanmodel(:,:,monthij)=comb_spatial_meanmodel(:,:,monthij)+interped_data/double(length(inputyear));

%                 comb_data(:,:,ind) = data;
                comb_interped_data(:,:,ind) = interped_data;
                comb_avhrr_data(:,:,ind) = avhrr_data;
    
                hold off
                close all;
                ind = ind + 1;
                toc;
            end
        end
        rms = sqrt(mean(comb_mse_data,3));
        
        abcabc= mean(comb_avhrr_err,3)';
        abcabc(abcabc>0.4)=-1;
        pcolor(abcabc); shading flat; colorbar
        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        trend(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
                trend(i,j)=p(1);
            end
        end

        avhrr_trend(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_avhrr_data(i,j,:))',1);
                avhrr_trend(i,j)=p(1);
            end
        end


        for t=1:length(inputyear)
            comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanmodel;
            comb_avhrr_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_avhrr_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanavhrr;
        end

        trend_filtered(1:len_lon,1:len_lat)=NaN;
        avhrr_trend_filtered(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
                trend_filtered(i,j)=p(1);
            end
        end

        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_avhrr_data_filtered(i,j,:))',1);
                avhrr_trend_filtered(i,j)=p(1);
            end
        end

        mean_trend=mean(mean(trend,'omitnan'),'omitnan');
        mean_trend_filtered=mean(mean(trend_filtered,'omitnan'),'omitnan');
        mean_avhrr_trend=mean(mean(avhrr_trend,'omitnan'),'omitnan');
        mean_avhrr_trend_filtered=mean(mean(avhrr_trend_filtered,'omitnan'),'omitnan');

    %     for monthij = 1:length(inputmonth)
    %         tempmonth = inputmonth(monthij);
    %         figdir=[figrawdir,'CLIM\'];
    %         outfile = strcat(figdir,regionname);
    %         if (exist(strcat(figdir) , 'dir') ~= 7)
    %             mkdir(strcat(figdir));
    %         end 
    %         % model plot
    %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    %         hold on;
    %         m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanmodel(:,:,monthij)));
    %         shading(gca,m_pcolor_shading_method);
    %         m_gshhs_i('color',m_gshhs_line_color);
    %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat('Model', ' (',char(calendarname(tempmonth)), ', clim)');
    %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    %     
    %         % set colorbar 
    %         h = colorbar;
    %         colormap(colormap_style);
    %         set(h,'fontsize',colorbar_fontsize);
    %         title(h,'^oC','fontsize',colorbar_title_fontsize);
    %         caxis(shadlev);
    %     
    %         % contour
    %         [C,h2]=m_contour(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanmodel(:,:,monthij)), conlev, m_contour_color, 'linewidth', m_contour_linewidth);
    %         clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
    %             'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
    %     
    %         % set grid
    %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    %     
    %         set(gcf, 'PaperUnits', 'points');
    %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    %     
    %         jpgname=strcat(outfile, '_', testname, '_Model_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    %         saveas(gcf,jpgname,'jpg');
    %     
    %     
    %         disp(' ')
    %         disp(['climatology mean_', num2str(tempmonth), '_', 'Model', ' plot is created.'])
    %         disp(' ')
    %         disp([' File path is : ',jpgname])
    %         disp(' ')
    %     
    %         hold off
    %         close all;
    %     
    %     
    %         % AVHRR plot
    %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    %         hold on;
    %         m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanavhrr(:,:,monthij)));
    %         shading(gca,m_pcolor_shading_method);
    %         m_gshhs_i('color',m_gshhs_line_color);
    %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat('AVHRR', ' (',char(calendarname(tempmonth)), ', clim)');
    %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    %     
    %         % set colorbar 
    %         h = colorbar;
    %         colormap(colormap_style);
    %         set(h,'fontsize',colorbar_fontsize);
    %         title(h,'^oC','fontsize',colorbar_title_fontsize);
    %         caxis(shadlev);
    %     
    %         % contour
    %         [C,h2]=m_contour(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanavhrr(:,:,monthij)), conlev, m_contour_color, 'linewidth', m_contour_linewidth);
    %         clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
    %             'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
    %     
    %         % set grid
    %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    %     
    %         set(gcf, 'PaperUnits', 'points');
    %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    %     
    %         jpgname=strcat(outfile, '_', testname, '_AVHRR_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    %         saveas(gcf,jpgname,'jpg');
    %     
    %     
    %         disp(' ')
    %         disp(['clim_', num2str(tempmonth), '_', 'AVHRR', ' plot is created.'])
    %         disp(' ')
    %         disp([' File path is : ',jpgname])
    %         disp(' ')
    %     
    %         hold off
    %         close all;
    %     
    %     
    %         % rms plot
    %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    %         hold on;
    %         m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanrms(:,:,monthij)));
    %         shading(gca,m_pcolor_shading_method);
    %         m_gshhs_i('color',m_gshhs_line_color);
    %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat('RMS', ' (',char(calendarname(tempmonth)), ', clim)','RMSE=',num2str(mean(mean(comb_spatial_meanrms(:,:,monthij),'omitnan'),'omitnan')));
    %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    %     
    %         % set colorbar 
    %         h = colorbar;
    %         colormap(colormap_style);
    %         set(h,'fontsize',colorbar_fontsize);
    %         caxis(rms_shadlev);
    %         
    %         % set grid
    %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    %     
    %         set(gcf, 'PaperUnits', 'points');
    %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    %     
    %         jpgname=strcat(outfile, '_', testname, '_rms_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    %         saveas(gcf,jpgname,'jpg');
    %     
    %     
    %         disp(' ')
    %         disp(['clim_', num2str(tempmonth), '_', 'RMS', ' plot is created.'])
    %         disp(' ')
    %         disp([' File path is : ',jpgname])
    %         disp(' ')
    %     
    %         hold off
    %         close all;
    %     
    %     
    %         % bias plot
    %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    %         hold on;
    %         m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanbias(:,:,monthij)));
    %         shading(gca,m_pcolor_shading_method);
    %         m_gshhs_i('color',m_gshhs_line_color);
    %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    %         titlename = strcat('BIAS', ' (',char(calendarname(tempmonth)), ', clim)','Mbias=',num2str(mean(mean(comb_spatial_meanbias(:,:,monthij),'omitnan'),'omitnan')));
    %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    %     
    %         % set colorbar 
    %         h = colorbar;
    %         colormap(bwrmap);
    %         set(h,'fontsize',colorbar_fontsize);
    %         caxis(bias_shadlev);
    %         
    %         % set grid
    %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
    %     
    %         set(gcf, 'PaperUnits', 'points');
    %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    %     
    %         jpgname=strcat(outfile, '_', testname, '_bias_cilm_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    %         saveas(gcf,jpgname,'jpg');
    %     
    %     
    %         disp(' ')
    %         disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
    %         disp(' ')
    %         disp([' File path is : ',jpgname])
    %         disp(' ')
    %     
    %         hold off
    %         close all;
    %     end


        figdir=[figrawdir,'CLIM\'];
        outfile = strcat(figdir,regionname);
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 

%         rmsplot=plot(mean(comb_meanrms,1),'k')
%         jpgname=strcat(outfile, '_', testname,'_',regionname,'_climrms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('month')
%         ylabel('RMS(^oC)')
%         title(['RMS, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
%         ylim([0 4])
%         set(rmsplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         saveas(gcf,jpgname,'jpg');
%         grid off

%         biasplot=plot(mean(comb_meanbias,1) ,'k')
%         jpgname=strcat(outfile, '_', testname,'_',regionname, '_climbias_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('month')
%         ylabel('bias(^o)')
%         title(['BIAS, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
%         ylim([-4 4])
%         set(biasplot,'LineWidth',2);
%         grid on
%         saveas(gcf,jpgname,'jpg');
%         grid off

        save([savedir,regionname,'_', scenname, '_sst_rms_and_bias_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

        ncid = netcdf.create(strcat(savedir, scenname, '_', testname, '_',regionname,'_rms_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');


        onedimid = netcdf.defDim(ncid,'one', 1);
        lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
        lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
        xidimid = netcdf.defDim(ncid, 'xi_rho', len_lon_model);
        etadimid = netcdf.defDim(ncid,'eta_rho',len_lat_model);
        time_dimid = netcdf.defDim(ncid, 'time', 0);
        clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'type', ['NWP 1/20 _ ', testname, 'monthly SST RMS/BIAS file']);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'title', ' monthly SST RMS/BIAS (1982-2009) ');
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'source', [' ROMS NWP 1/20 data from _ ',testname ]);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'author', 'Created by Y.Y.Kim');
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'date', date);

        timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
        netcdf.putAtt(ncid,timevarid,'long_name','time');
        netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
        netcdf.putAtt(ncid,timevarid,'calendar','gregorian');

        clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
        netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
        netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
        netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

        lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', lon_dimid);
        netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
        netcdf.putAtt(ncid,lonvarid,'units','degree_east');

        latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', lat_dimid);
        netcdf.putAtt(ncid,latvarid,'long_name','latitude');
        netcdf.putAtt(ncid,latvarid,'units','degree_north');

        lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xidimid etadimid]);
        netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
        netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

        lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xidimid etadimid]);
        netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
        netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');

        raw_sstvarid=netcdf.defVar(ncid, 'raw_sst', 'NC_FLOAT', [xidimid etadimid time_dimid]);
        netcdf.putAtt(ncid,raw_sstvarid,'long_name','raw_sst');
        netcdf.putAtt(ncid,raw_sstvarid,'units','Celsius');

        sstvarid=netcdf.defVar(ncid, 'sst', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,sstvarid,'long_name','sst');
        netcdf.putAtt(ncid,sstvarid,'units','Celsius');

        sst_filteredvarid=netcdf.defVar(ncid, 'sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,sst_filteredvarid,'long_name','sst_filtered');
        netcdf.putAtt(ncid,sst_filteredvarid,'units','Celsius');

        avhrr_sstvarid=netcdf.defVar(ncid, 'avhrr_sst', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,avhrr_sstvarid,'long_name','avhrr_sst');
        netcdf.putAtt(ncid,avhrr_sstvarid,'units','Celsius');

        avhrr_sst_filteredvarid=netcdf.defVar(ncid, 'avhrr_sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,avhrr_sst_filteredvarid,'long_name','avhrr_sst_filtered');
        netcdf.putAtt(ncid,avhrr_sst_filteredvarid,'units','Celsius');

        trendvarid=netcdf.defVar(ncid, 'trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,trendvarid,'long_name','trend');
        netcdf.putAtt(ncid,trendvarid,'units','Celsius/year');

        trend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,trend_filteredvarid,'long_name','trend_filtered');
        netcdf.putAtt(ncid,trend_filteredvarid,'units','Celsius/year');

        avhrr_trendvarid=netcdf.defVar(ncid, 'avhrr_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,avhrr_trendvarid,'long_name','avhrr_trend');
        netcdf.putAtt(ncid,avhrr_trendvarid,'units','Celsius/year');

        avhrr_trend_filteredvarid=netcdf.defVar(ncid, 'avhrr_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,avhrr_trend_filteredvarid,'long_name','avhrr_trend_filtered');
        netcdf.putAtt(ncid,avhrr_trend_filteredvarid,'units','Celsius/year');
        
        msevarid=netcdf.defVar(ncid, 'mse', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,msevarid,'long_name','mse');
        netcdf.putAtt(ncid,msevarid,'units','Celsius^2');
        
        rmsvarid=netcdf.defVar(ncid, 'rms', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,rmsvarid,'long_name','rms');
        netcdf.putAtt(ncid,rmsvarid,'units','Celsius');

        biasvarid=netcdf.defVar(ncid, 'bias', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,biasvarid,'long_name','bias');
        netcdf.putAtt(ncid,biasvarid,'units','Celsius');

        mean_trendvarid=netcdf.defVar(ncid, 'mean_trend', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trendvarid,'long_name','mean_trend');
        netcdf.putAtt(ncid,mean_trendvarid,'units','Celsius/year');

        mean_trend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'long_name','mean_trend_filtered');
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'units','Celsius/year');

        mean_avhrr_trendvarid=netcdf.defVar(ncid, 'mean_avhrr_trend', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_avhrr_trendvarid,'long_name','mean_avhrr_trend');
        netcdf.putAtt(ncid,mean_avhrr_trendvarid,'units','Celsius/year');

        mean_avhrr_trend_filteredvarid=netcdf.defVar(ncid, 'mean_avhrr_trend_filtered', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_avhrr_trend_filteredvarid,'long_name','mean_avhrr_trend_filtered');
        netcdf.putAtt(ncid,mean_avhrr_trend_filteredvarid,'units','Celsius/year');

        clim_sstvarid=netcdf.defVar(ncid, 'clim_sst', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_sstvarid,'long_name','clim_sst');
        netcdf.putAtt(ncid,clim_sstvarid,'units','Celsius');

        clim_avhrrvarid=netcdf.defVar(ncid, 'clim_avhrr_sst', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_avhrrvarid,'long_name','clim_avhrr_sst');
        netcdf.putAtt(ncid,clim_avhrrvarid,'units','Celsius');

%         clim_rmsvarid=netcdf.defVar(ncid, 'clim_rms', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
%         netcdf.putAtt(ncid,clim_rmsvarid,'long_name','clim_rms');
%         netcdf.putAtt(ncid,clim_rmsvarid,'units','Celsius');

        clim_biasvarid=netcdf.defVar(ncid, 'clim_bias', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_biasvarid,'long_name','clim_bias');
        netcdf.putAtt(ncid,clim_biasvarid,'units','Celsius');

        netcdf.endDef(ncid);

        tind=1;
        for yearij = 1:length(inputyear)
            for month=1:12 
                tempyear = inputyear(yearij);
                ftime(tind) = datenum(tempyear,month,15) - datenum(1900,12,31);
                tind=tind+1;
            end
        end
        for month=1:12 
                tempyear = inputyear(yearij);
                climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
        end

        netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
        netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
        netcdf.putVar(ncid, lonvarid, 0, len_lon, avhrr_lon(:,1));
        netcdf.putVar(ncid, latvarid, 0, len_lat, avhrr_lat(1,:));
        netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon_rho);
        netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat_rho);
        netcdf.putVar(ncid, raw_sstvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data);
        netcdf.putVar(ncid, sstvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data);
        netcdf.putVar(ncid, sst_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data_filtered);
        netcdf.putVar(ncid, avhrr_sstvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data);
        netcdf.putVar(ncid, avhrr_sst_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_avhrr_data_filtered);
        netcdf.putVar(ncid, msevarid, [0 0 0], [len_lon len_lat length(ftime)], comb_mse_data);
        netcdf.putVar(ncid, rmsvarid, [0 0], [len_lon len_lat], rms);
        netcdf.putVar(ncid, biasvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_bias_data);
        netcdf.putVar(ncid, clim_sstvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanmodel);
        netcdf.putVar(ncid, clim_avhrrvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanavhrr);
%         netcdf.putVar(ncid, clim_rmsvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanrms);
        netcdf.putVar(ncid, clim_biasvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanbias);
        netcdf.putVar(ncid, trendvarid, [0 0], [len_lon len_lat], trend);
        netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered);
        netcdf.putVar(ncid, avhrr_trendvarid, [0 0], [len_lon len_lat], avhrr_trend);
        netcdf.putVar(ncid, avhrr_trend_filteredvarid, [0 0], [len_lon len_lat], avhrr_trend_filtered);
        netcdf.putVar(ncid, mean_trendvarid, [0], [1], mean_trend);
        netcdf.putVar(ncid, mean_trend_filteredvarid, [0], [1], mean_trend_filtered);
        netcdf.putVar(ncid, mean_avhrr_trendvarid, [0], [1], mean_avhrr_trend);
        netcdf.putVar(ncid, mean_avhrr_trend_filteredvarid, [0], [1], mean_avhrr_trend_filtered);

        netcdf.close(ncid);
    end
end

% SSH_4th_mid_report3