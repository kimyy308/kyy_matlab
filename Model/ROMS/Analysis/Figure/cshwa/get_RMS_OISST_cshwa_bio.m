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

shadlev = [0 35];
rms_shadlev = [0 4];
bias_shadlev = [-4 4];
conlev  = 0:5:35;

% for snu_desktop
testname='cshwa_bio'   % % need to change
inputyear = [1982:1985]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

filedir = strcat('E:\Data\Model\ROMS\cshwa_bio\'); % % where data files are
varname ='temp'

% % for NWP
lon=ncread([filedir, '\ocean_avg_0001.nc'], 'lon_rho');
lat=ncread([filedir, '\ocean_avg_0001.nc'], 'lat_rho');
lonlat(1)=min(lon(:));
lonlat(2)=max(lon(:));
lonlat(3)=min(lat(:));
lonlat(4)=max(lat(:));

regionname='NWP';

% % % for EKB
% regionname='EKB';
% lonlat = [127, 129.5, 38, 40.5];


figrawdir =strcat('E:\Data\Model\ROMS\cshwa_bio\figure\'); % % where figure files will be saved

param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'

run(param_script);
load ('E:\Data\Model\ROMS\cshwa_bio\temp_surf.mat')
for yearij = 1:length(inputyear)
    for monthij = 1:length(inputmonth)
        tempyear = inputyear(yearij);
        tempmonth = inputmonth(monthij);
        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
%         filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
%                 testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
%         % read model data
%         if (exist('lon')==0)
%             modelinfo=ncinfo(filename);
%             lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
%             lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);
% 
%             lon_west = abs(lon - (lonlat(1)-1));
%             min_lon_west=min(lon_west);
%             lon_east = abs(lon - (lonlat(2)+1));
%             min_lon_east=min(lon_east);
%             lat_south = abs(lat - (lonlat(3)-1));
%             min_lat_south=min(lat_south);
%             lat_north = abs(lat - (lonlat(4)+1));
%             min_lat_north=min(lat_north);
% 
%             lon_min = find(lon_west == min_lon_west);
%             lon_max = find(lon_east == min_lon_east);
%             lat_min = find(lat_south == min_lat_south);
%             lat_max = find(lat_north == min_lat_north);
% 
%             lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
%             lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
%         end
% 
%         data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]
% 
%         data = ncread(filename,varname,[lon_min(1) lat_min(1) 40 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
        data = temp_surf(:,:,(yearij-1)*12+monthij);
        % read OISST DATA
        avhrrdir='E:\Data\Observation\OISST\monthly\';
        avhrrfilename = strcat(avhrrdir,'avhrr_monthly', num2str(tempyear,'%04i'), ...
                '_',num2str(tempmonth,'%02i'), '.nc');
        if (exist('avhrr_lon')==0)
            avhrrinfo=ncinfo(avhrrfilename);
            avhrr_lon = ncread(avhrrfilename,'long',[1 1],[avhrrinfo.Dimensions(2).Length,1]);
            avhrr_lat = ncread(avhrrfilename,'lat',[1 1],[1,avhrrinfo.Dimensions(1).Length]);

            avhrr_lon_west = abs(avhrr_lon - (lonlat(1)));
            min_avhrr_lon_west=min(avhrr_lon_west);
            avhrr_lon_east = abs(avhrr_lon - (lonlat(2)));
            min_avhrr_lon_east=min(avhrr_lon_east);
            avhrr_lat_south = abs(avhrr_lat - (lonlat(3)));
            min_avhrr_lat_south=min(avhrr_lat_south);
            avhrr_lat_north = abs(avhrr_lat - (lonlat(4)));
            min_avhrr_lat_north=min(avhrr_lat_north);

            avhrr_lon_min = find(avhrr_lon_west == min_avhrr_lon_west)-3;
            avhrr_lon_max = find(avhrr_lon_east == min_avhrr_lon_east)+3;
            avhrr_lat_min = find(avhrr_lat_south == min_avhrr_lat_south)-3;
            avhrr_lat_max = find(avhrr_lat_north == min_avhrr_lat_north)+3;

    %         ncinfo('E:\Data\Observation\OISST\monthly\avhrr_monthly1983_11.nc');

            avhrr_lon = ncread(avhrrfilename,'long', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
            avhrr_lat = ncread(avhrrfilename,'lat', [avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);
            
            comb_spatial_meanrms=(zeros([size(avhrr_lon),12]));
            comb_spatial_meanbias=(zeros([size(avhrr_lon),12]));
            comb_spatial_meanavhrr=(zeros([size(avhrr_lon),12]));
            comb_spatial_meanmodel=(zeros([size(avhrr_lon),12]));
        end

        avhrr_data = ncread(avhrrfilename,varname,[avhrr_lon_min(1) avhrr_lat_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1)]);

        interped_data = griddata(double(lon), double(lat), data,double(avhrr_lon),double(avhrr_lat)); 
        bias = avhrr_data-interped_data;
        rms = sqrt((avhrr_data-interped_data).^2);
        meanbias = mean(mean(bias,'omitnan'),'omitnan');
        meanrms = mean(mean(rms,'omitnan'),'omitnan');
        comb_meanrms(yearij,monthij)=meanrms;
        comb_meanbias(yearij,monthij)=meanbias;
        comb_spatial_meanrms(:,:,monthij)=comb_spatial_meanrms(:,:,monthij)+rms/double(length(inputyear));
        comb_spatial_meanbias(:,:,monthij)=comb_spatial_meanbias(:,:,monthij)+bias/double(length(inputyear));
        comb_spatial_meanavhrr(:,:,monthij)=comb_spatial_meanavhrr(:,:,monthij)+avhrr_data/double(length(inputyear));
        comb_spatial_meanmodel(:,:,monthij)=comb_spatial_meanmodel(:,:,monthij)+interped_data/double(length(inputyear));

        % AVHRR plot
        figdir=[figrawdir,'AVHRR\'];
        outfile = strcat(figdir,regionname);
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 

% % % %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% % % %         hold on;
% % % %         m_pcolor(double(avhrr_lon),avhrr_lat,avhrr_data);
% % % %         shading(gca,m_pcolor_shading_method);
% % % %         m_gshhs_i('color',m_gshhs_line_color);
% % % %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % %         titlename = strcat('AVHRR', ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
% % % %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % % 
% % % %         % set colorbar 
% % % %         h = colorbar;
% % % %         colormap(colormap_style);
% % % %         set(h,'fontsize',colorbar_fontsize);
% % % %         title(h,'^oC','fontsize',colorbar_title_fontsize);
% % % %         caxis(shadlev);
% % % %         
% % % % %         % contour
% % % % %         [C,h2]=m_contour(double(avhrr_lon),avhrr_lat,avhrr_data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
% % % % %         clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
% % % % %             'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
% % % %         
% % % %         % set grid
% % % %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % % % 
% % % %         set(gcf, 'PaperUnits', 'points');
% % % %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % % %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % % 
% % % %         jpgname=strcat(outfile, '_', 'cshwa_bio', '_AVHRR_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
% % % %         saveas(gcf,jpgname,'jpg');
% % % % 
% % % % 
% % % %         disp(' ')
% % % %         disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'AVHRR', ' plot is created.'])
% % % %         disp(' ')
% % % %         disp([' File path is : ',jpgname])
% % % %         disp(' ')
% % % % 
% % % %         hold off
% % % %         close all;

        
% % % %         % rms plot
% % % %         figdir=[figrawdir,'RMS\'];
% % % %         outfile = strcat(figdir,regionname);
% % % %         if (exist(strcat(figdir) , 'dir') ~= 7)
% % % %             mkdir(strcat(figdir));
% % % %         end 
% % % %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% % % %         hold on;
% % % %         m_gshhs_i('color',m_gshhs_line_color);
% % % %         m_pcolor(double(avhrr_lon),avhrr_lat,rms);
% % % % %         pcolor(double(avhrr_lon),avhrr_lat,rms);
% % % % %         shading(gca,m_pcolor_shading_method);
% % % %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % %         titlename = strcat('RMS', ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')','RMSE=',num2str(meanrms));
% % % %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % % 
% % % %         % set colorbar 
% % % %         h = colorbar;
% % % %         colormap(colormap_style);
% % % %         set(h,'fontsize',colorbar_fontsize);
% % % % 
% % % %         % set grid
% % % %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % % % 
% % % %         set(gcf, 'PaperUnits', 'points');
% % % %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % % %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % % 
% % % %         jpgname=strcat(outfile, '_', testname, '_rms_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
% % % %         saveas(gcf,jpgname,'jpg');
% % % % 
% % % % 
% % % %         disp(' ')
% % % %         disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'RMS', ' plot is created.'])
% % % %         disp(' ')
% % % %         disp([' File path is : ',jpgname])
% % % %         disp(' ')
% % % % 
% % % %         hold off
% % % %         close all;
% % % %         
% % % %         
% % % %         % bias plot
% % % %         figdir=[figrawdir,'BIAS\'];
% % % %         outfile = strcat(figdir,regionname);
% % % %         if (exist(strcat(figdir) , 'dir') ~= 7)
% % % %             mkdir(strcat(figdir));
% % % %         end 
% % % %         m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
% % % %         hold on;
% % % %         m_pcolor(double(avhrr_lon),avhrr_lat,bias);
% % % % %         shading(gca,m_pcolor_shading_method);
% % % %         m_gshhs_i('color',m_gshhs_line_color);
% % % %         m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% % % %         titlename = strcat('BIAS', ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')','Mbias=',num2str(meanbias));
% % % %         title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % % 
% % % %         % set colorbar 
% % % %         h = colorbar;
% % % %         colormap(colormap_style);
% % % %         set(h,'fontsize',colorbar_fontsize);
% % % % 
% % % %         % set grid
% % % %         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
% % % % 
% % % %         set(gcf, 'PaperUnits', 'points');
% % % %         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % % %         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % % 
% % % %         jpgname=strcat(outfile, '_', testname, '_bias_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
% % % %         saveas(gcf,jpgname,'jpg');


% % % %         disp(' ')
% % % %         disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
% % % %         disp(' ')
% % % %         disp([' File path is : ',jpgname])
% % % %         disp(' ')
% % % % 
% % % %         hold off
% % % %         close all;
    end
end


for monthij = 1:length(inputmonth)
    tempmonth = inputmonth(monthij);
    figdir=[figrawdir,'CLIM\'];
    outfile = strcat(figdir,regionname);
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    % model plot
    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanmodel(:,:,monthij)));
%     shading(gca,m_pcolor_shading_method);
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    titlename = strcat('Model', ' (',char(calendarname(tempmonth)), ', clim)');
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(colormap_style);
    set(h,'fontsize',colorbar_fontsize);
    title(h,'^oC','fontsize',colorbar_title_fontsize);
    caxis(shadlev);

    % contour
    [C,h2]=m_contour(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanmodel(:,:,monthij)), conlev, m_contour_color, 'linewidth', m_contour_linewidth);
    clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
        'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);

    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, '_', testname, '_Model_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'jpg');


    disp(' ')
    disp(['climatology mean_', num2str(tempmonth), '_', 'Model', ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;


    % AVHRR plot
    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanavhrr(:,:,monthij)));
%     shading(gca,m_pcolor_shading_method);
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    titlename = strcat('AVHRR', ' (',char(calendarname(tempmonth)), ', clim)');
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(colormap_style);
    set(h,'fontsize',colorbar_fontsize);
    title(h,'^oC','fontsize',colorbar_title_fontsize);
    caxis(shadlev);

    % contour
    [C,h2]=m_contour(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanavhrr(:,:,monthij)), conlev, m_contour_color, 'linewidth', m_contour_linewidth);
    clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
        'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);

    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, '_', testname, '_AVHRR_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'jpg');


    disp(' ')
    disp(['clim_', num2str(tempmonth), '_', 'AVHRR', ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;


    % rms plot
    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanrms(:,:,monthij)));
%     shading(gca,m_pcolor_shading_method);
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    titlename = strcat('RMS', ' (',char(calendarname(tempmonth)), ', clim)','RMSE=',num2str(mean(mean(comb_spatial_meanrms(:,:,monthij),'omitnan'),'omitnan')));
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(colormap_style);
    set(h,'fontsize',colorbar_fontsize);
    caxis(rms_shadlev);
    
    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, '_', testname, '_rms_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'jpg');


    disp(' ')
    disp(['clim_', num2str(tempmonth), '_', 'RMS', ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;


    % bias plot
    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    m_pcolor(double(avhrr_lon),avhrr_lat,squeeze(comb_spatial_meanbias(:,:,monthij)));
%     shading(gca,m_pcolor_shading_method);
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    titlename = strcat('BIAS', ' (',char(calendarname(tempmonth)), ', clim)','Mbias=',num2str(mean(mean(comb_spatial_meanbias(:,:,monthij),'omitnan'),'omitnan')));
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(bwrmap);
    set(h,'fontsize',colorbar_fontsize);
    caxis(bias_shadlev);
    
    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, '_', testname, '_bias_cilm_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'jpg');


    disp(' ')
    disp(['clim_', num2str(tempmonth), '_', 'BIAS', ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;
end


figdir=[figrawdir,'CLIM\'];
outfile = strcat(figdir,regionname);
if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 

plot(mean(comb_meanrms,1))
mrms=mean(comb_meanrms(:));
jpgname=strcat(outfile, '_', testname, '_climrms', '.jpg'); %% ~_year_month.jpg
xlabel('month')
ylabel('rms(^o)')
title(['rms, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')', 'M=', num2str(round(mrms,2))])
ylim([0 4])
grid on
saveas(gcf,jpgname,'jpg');
grid off


plot(mean(comb_meanbias,1))
mbias=mean(comb_meanbias(:));
jpgname=strcat(outfile, '_', testname, '_climbias', '.jpg'); %% ~_year_month.jpg
xlabel('month')
ylabel('bias(^o)')
title(['bias, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')', 'M=', num2str(round(mbias,2))])
ylim([-4 4])
grid on
saveas(gcf,jpgname,'jpg');
grid off

save([filedir,regionname,'_rms_and_bias_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);