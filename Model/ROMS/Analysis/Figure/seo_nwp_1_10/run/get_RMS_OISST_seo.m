clear all; close all; clc;
windows=1;
if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
elseif (linux==1)
    dropboxpath='/home01/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end
% % % 
% % % Read Model SST
% % % interp
% % % get RMS
% % % get BIAS

% for snu_desktop
testname='avg_ens_10km_mean'   % % need to change
% year = [2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012]; % % put year which you want to plot [year year ...]
inputyear = [2001 2002 2003 2004 2005 2006 2007 2008 2009 2010]; % % put year which you want to plot [year year ...]
% inputyear = [2010]; % % put year which you want to plot [year year ...]
% month = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
figdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\',testname,'\'); % % where figure files will be saved
outfile = strcat(figdir,'SS\SS');
filedir = strcat('E:\Data\Reanalysis\nwp_1_10_seo\', testname, '\'); % % where data files are
lonlat = [127, 129.5, 33.5, 34.5];
varname ='temp'
% % if (exist(strcat(figdir) , 'dir') ~= 7)
% %     mkdir(strcat(figdir));
% % end 
param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_south_sea.m'

run(param_script);
for yearij = 1:length(inputyear)
    for monthij = 1:length(inputmonth)
        tempyear = inputyear(yearij);
        tempmonth = inputmonth(monthij);
        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
        filename = strcat(filedir, testname, '_monthly_',num2str(tempyear,'%04i'),'_', num2str(tempmonth,'%02i'), '.nc');
        grdname = strcat(filedir,'roms_grid_final.nc');
        % read model data
        if (exist('lon')==0)
            lon = ncread(grdname,'lon_rho',[1 1],[334,1]);
            lat = ncread(grdname,'lat_rho',[1 1],[1,332]);

            lon_west = abs(lon - (lonlat(1)-1));
            min_lon_west=min(lon_west);
            lon_east = abs(lon - (lonlat(2)+1));
            min_lon_east=min(lon_east);
            lat_south = abs(lat - (lonlat(3)-1));
            min_lat_south=min(lat_south);
            lat_north = abs(lat - (lonlat(4)+1));
            min_lat_north=min(lat_north);

            lon_min = find(lon_west == min_lon_west);
            lon_max = find(lon_east == min_lon_east);
            lat_min = find(lat_south == min_lat_south);
            lat_max = find(lat_north == min_lat_north);

            lon = ncread(grdname,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
            lat = ncread(grdname,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
        end

        data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

        data = ncread(filename,varname,[lon_min(1) lat_min(1) 20 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
        data(find(data ==0))=NaN;
        data(find(data >100))=NaN;
        data(find(data <-100))=NaN;
        % read OISST DATA
        avhrrdir='E:\Data\Observation\OISST\monthly\'
        avhrrfilename = strcat(avhrrdir,'avhrr_monthly', num2str(tempyear,'%04i'), ...
                '_',num2str(tempmonth,'%02i'), '.nc');
        if (exist('avhrr_lon')==0)
            avhrr_lon = ncread(avhrrfilename,'long',[1 1],[1440,1]);
            avhrr_lat = ncread(avhrrfilename,'lat',[1 1],[1,720]);

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

            avhrr_lon = ncread(avhrrfilename,'long', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
            avhrr_lat = ncread(avhrrfilename,'lat', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
        end

        avhrr_data = ncread(avhrrfilename,varname,[avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);

        interped_data = griddata(double(lon), double(lat), double(data),double(avhrr_lon),double(avhrr_lat)); 
        
        bias = avhrr_data-interped_data;
        rms = sqrt((avhrr_data-interped_data).^2);
        meanbias = mean(mean(bias,'omitnan'),'omitnan');
        meanrms = mean(mean(rms,'omitnan'),'omitnan');
        comb_meanrms(yearij,monthij)=meanrms;
        comb_meanbias(yearij,monthij)=meanbias;
        
        
        
% % %         % rms plot
% % %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% % %         hold on;
% % %         m_pcolor(double(avhrr_lon),avhrr_lat,rms);
% % %         shading(gca,m_pcolor_shading_method);
% % %         m_gshhs_i('color',m_gshhs_line_color);
% % %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % %         titlename = strcat('RMS', ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')','RMSE=',num2str(meanrms));
% % %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % 
% % %         % set colorbar 
% % %         h = colorbar;
% % %         colormap(colormap_style);
% % %         set(h,'fontsize',colorbar_fontsize);
% % %     %     title(h,'RMS','fontsize',colorbar_title_fontsize);
% % % 
% % %     %     caxis(shadlev);
% % % 
% % %         % set grid
% % %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % % 
% % %         set(gcf, 'PaperUnits', 'points');
% % %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % 
% % %         jpgname=strcat(outfile, '_', testname, '_rms_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
% % %         saveas(gcf,jpgname,'jpg');
% % % 
% % % 
% % %         disp(' ')
% % %         disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'RMS', ' plot is created.'])
% % %         disp(' ')
% % %         disp([' File path is : ',jpgname])
% % %         disp(' ')
% % % 
% % %         hold off
% % %         close all;
% % %         
% % %         
% % %         % bias plot
% % %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% % %         hold on;
% % %         m_pcolor(double(avhrr_lon),avhrr_lat,bias);
% % %         shading(gca,m_pcolor_shading_method);
% % %         m_gshhs_i('color',m_gshhs_line_color);
% % %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % %         titlename = strcat('BIAS', ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')','RMSE=',num2str(meanrms));
% % %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % 
% % %         % set colorbar 
% % %         h = colorbar;
% % %         colormap(colormap_style);
% % %         set(h,'fontsize',colorbar_fontsize);
% % %     %     title(h,'RMS','fontsize',colorbar_title_fontsize);
% % % 
% % %     %     caxis(shadlev);
% % % 
% % %         % set grid
% % %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % % 
% % %         set(gcf, 'PaperUnits', 'points');
% % %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % 
% % %         jpgname=strcat(outfile, '_', testname, '_bias_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
% % %         saveas(gcf,jpgname,'jpg');
% % % 
% % % 
% % %         disp(' ')
% % %         disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
% % %         disp(' ')
% % %         disp([' File path is : ',jpgname])
% % %         disp(' ')
% % % 
% % %         hold off
% % %         close all;
    end
end

plot(mean(comb_meanrms,1))
jpgname=strcat(outfile, '_', testname, '_climrms', '.jpg'); %% ~_year_month.jpg
xlabel('month')
ylabel('rms(^o)')
title('rms, climatology(2001-2010)')
ylim([0 2])
saveas(gcf,jpgname,'jpg');
plot(mean(comb_meanbias,1))
jpgname=strcat(outfile, '_', testname, '_climbias', '.jpg'); %% ~_year_month.jpg
xlabel('month')
ylabel('bias(^o)')
title('bias, climatology(2001-2010)')
ylim([-2 2])
saveas(gcf,jpgname,'jpg');