close all; clear all; clc;

[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/mat_tool']));
addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
addpath(genpath([dropboxpath,'/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/MICT_pollack/paper/subroutine/']))

param_script =[dropboxpath,'/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_ES', '.m'];
variable = 'temp';
run(param_script);
run('nwp_polygon_point.m');

study_area_polygon = ...
    [127.0, 36.0;
    132, 36.0;
    132, 44.0;
    127.0, 44.0;];
% % %      North Korean EEZ
% NK_EEZ_polygon = ...
%     [132.6, 38.7;
%     133.0, 39.0;
%     133.5, 40.0;
%     130.7 42.4];
% 130.0, 43.0;
% [133.5 40.0] + [-3.5 3.0].*0.8
lonlat = [127, 144, 33, 52];
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, 'parent');
m_gshhs_i('color',m_gshhs_line_color);
m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% m_line(SK_EEZ_polygon(:,1), SK_EEZ_polygon(:,2), 'LineStyle', '--', 'Color', 'b', 'LineWidth', 1)
hold on
m_line(study_area_polygon(:,1), study_area_polygon(:,2), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
% m_line(NK_EEZ_polygon(:,1), NK_EEZ_polygon(:,2), 'LineStyle', '--')
% m_line(JP_EEZ_polygon(:,1), JP_EEZ_polygon(:,2), 'LineStyle', '--')
hold off
tifname = '/Users/kimyy/Desktop/backup/Research/Ph_D_course/2022_pollock_future/submit_figure/Figureset_20230424/Fig01_01_04.tif';
print('-dtiff','-r500',tifname); RemoveWhiteSpace([], 'file', tifname);
close all


lonlat = [127, 132, 36, 44];
spawning_ground_polygon = ...
    [127.0, 37.0;
    129.7, 37.0;
    129.7, 41.0;
    127.0, 41.0;];

etoponame='/Volumes/kyy_raid/Data/Topography/Etopo/etopo1.nc';
ncinfo(etoponame)
etopo_lon=ncread(etoponame,'lon');
etopo_lat=ncread(etoponame,'lat');
etopo_z=ncread(etoponame,'z');
[lat_rho, lon_rho] = meshgrid(etopo_lat, etopo_lon);

study_area_mask=double(inpolygon(lon_rho,lat_rho,study_area_polygon(:,1),study_area_polygon(:,2)));
study_area_mask(study_area_mask==0)=NaN;
lon_rho=lon_rho.*study_area_mask;
lat_rho=lat_rho.*study_area_mask;
etopo_z=etopo_z.*study_area_mask;



study_area_mask=double(inpolygon(lon_rho,lat_rho,study_area_polygon(:,1),study_area_polygon(:,2)));
study_area_mask(study_area_mask==0)=NaN;

size_lon = size(lon_rho, 1);
size_lat = size(lon_rho, 2);

lon_min=min(mod(find(study_area_mask(:,:)==1),size_lon));
lon_max=max(mod(find(study_area_mask(:,:)==1),size_lon));
lat_min=min(mod(find(study_area_mask(:,:)'==1),size_lat));
lat_max=max(mod(find(study_area_mask(:,:)'==1),size_lat));

loncount= lon_max(1)-lon_min(1)+1;
latcount = lat_max(1)-lat_min(1)+1;

lon_rho2=lon_rho(lon_min(1):lon_min(1)+loncount, lat_min(1):lat_min(1)+latcount);
lat_rho2=lat_rho(lon_min(1):lon_min(1)+loncount, lat_min(1):lat_min(1)+latcount);
etopo_z2=etopo_z(lon_min(1):lon_min(1)+loncount, lat_min(1):lat_min(1)+latcount);
etopo_z2(etopo_z2>0)=NaN;
hold on
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_grid('fontsize', m_grid_fontsize-7, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, 'parent');
SK_EEZ_mask = double(inpolygon(lon_rho2,lat_rho2,SK_EEZ_polygon(:,1),SK_EEZ_polygon(:,2)));
SK_EEZ_mask(SK_EEZ_mask==0)=NaN;
SK_EEZ_mask(:)=NaN;
m_pcolor(lon_rho2,lat_rho2,SK_EEZ_mask);
colormap([37/255 141/255 239/255]);
shading flat;
m_gshhs_i('color',m_gshhs_line_color);
m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
% [etopo_con, etopo_h2]=m_contour(lon_rho2, lat_rho2, etopo_z2, [0 -50 -500 -2000 -3000], 'linecolor',[0.8 0.8 0.8], 'linewidth', 0.5);

% clabel(etopo_con,etopo_h2,'FontSize',m_contour_label_fontsize,'Color',[0.8 0.8 0.8], ...
%         'labelspacing',10000000,'Rotation',m_contour_rotation);


% m_line(SK_EEZ_polygon(:,1), SK_EEZ_polygon(:,2), 'LineStyle', '--', 'Color', 'b', 'LineWidth', 1) 
% m_line(spawning_ground_polygon(:,1), spawning_ground_polygon(:,2), 'LineStyle', '--', 'Color', 'r', 'LineWidth', 1)
% m_line(NK_EEZ_polygon(:,1), NK_EEZ_polygon(:,2), 'LineStyle', '--')
% m_line(JP_EEZ_polygon(:,1), JP_EEZ_polygon(:,2), 'LineStyle', '--')
m_grid('fontsize', m_grid_fontsize-7, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, 'parent');
hold off
tifname = '/Users/kimyy/Desktop/backup/Research/Ph_D_course/2022_pollock_future/submit_figure/Figureset_20230424/Fig01_02_04.tif';
print('-dtiff','-r500',tifname); RemoveWhiteSpace([], 'file', tifname);
close all