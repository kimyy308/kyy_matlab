close all; clear all; clc;

dropboxpath='C:\users\user/Dropbox/';
addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/mat_tool']));
addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
addpath(genpath('C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\MICT_pollack\paper\subroutine\'))

param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_ES', '.m'];
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
% m_gshhs_i('color',m_gshhs_line_color);
% m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
m_gshhs_i('color',[198/256 171/256 114/256]);
m_gshhs_i('patch',[198/256 171/256 114/256]);   % gray colored land

% m_line(SK_EEZ_polygon(:,1), SK_EEZ_polygon(:,2), 'LineStyle', '--', 'Color', 'b', 'LineWidth', 1)
hold on
m_line(study_area_polygon(:,1), study_area_polygon(:,2), 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1)
% m_line(NK_EEZ_polygon(:,1), NK_EEZ_polygon(:,2), 'LineStyle', '--')
% m_line(JP_EEZ_polygon(:,1), JP_EEZ_polygon(:,2), 'LineStyle', '--')
hold off
tifname = 'D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\Fig01_01_08.tif';
print('-dtiff','-r500',tifname); RemoveWhiteSpace([], 'file', tifname);
close all


lonlat = [127, 132, 36, 44];
spawning_ground_polygon = ...
    [127.0, 37.0;
    129.7, 37.0;
    129.7, 41.0;
    127.0, 41.0;];

etoponame='D:\Data\Topography\Etopo\etopo1.nc';
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

etopo_z3=etopo_z2;
etopo_z3(etopo_z3>-50)=NaN;
etopo_z3(etopo_z3<-500)=NaN;
etopo_z3(isfinite(etopo_z3))=1;

spawning_polygon=[127, 37.5; 127, 41; 132, 41; 132, 37.5];
sp_mask=double(inpolygon(lon_rho2,lat_rho2,spawning_polygon(:,1),spawning_polygon(:,2)));
sp_mask(sp_mask==0)=NaN;
etopo_z3=etopo_z3.*sp_mask;

hold on
% ax1=axes;
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_grid('fontsize', m_grid_fontsize-7, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type);

SK_EEZ_mask = double(inpolygon(lon_rho2,lat_rho2,SK_EEZ_polygon(:,1),SK_EEZ_polygon(:,2)));
SK_EEZ_mask(SK_EEZ_mask==0)=NaN;
% SK_EEZ_mask(isfinite(etopo_z3))=1;
% SK_EEZ_mask(isfinite(etopo_z3))=SK_EEZ_mask(isfinite(etopo_z3))+etopo_z3(isfinite(etopo_z3));

% % %  83-87 (early period, period I)
% run(param_script);
% m_grid_fontsize = m_grid_fontsize -4;
% ind=1;
clear comb_data
savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
earlyyear=[1983:1992];
inputmonth = [1,2]; % % put month which you want to plot [month month ...]
testname = 'test06';
regionname = 'pollock_egg3';
for yearij = 1:length(earlyyear)
    tempyear = earlyyear(yearij);
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        data=ncread(ncname, 'egg_mask');
        lastday_m=size(data,3);
        if (exist('comb_data')==0)
            comb_data=data;
        else
            comb_data(:,:,end+1:end+lastday_m)=data;
        end
    end
end
res_comb_data=reshape(comb_data, [size(comb_data,1)*size(comb_data,2), size(comb_data,3)]);
sum_comb_data=sum(res_comb_data,1,'omitnan');
lon_RCM = ncread(ncname, 'lon_rho');
lat_RCM = ncread(ncname, 'lat_rho');
mean_RCM_early = sum(comb_data,3, 'omitnan')/size(comb_data,3);
mean_RCM_early(mean_RCM_early==0)=NaN;
mean_RCM_early(isfinite(mean_RCM_early))=1;
SK_EEZ_mask2 = double(inpolygon(lon_RCM,lat_RCM,SK_EEZ_polygon(:,1),SK_EEZ_polygon(:,2)));
SK_EEZ_mask2(SK_EEZ_mask2==0)=NaN;
SK_EEZ_mask2(isfinite(mean_RCM_early))=1;
SK_EEZ_mask2(isfinite(mean_RCM_early))=SK_EEZ_mask2(isfinite(mean_RCM_early))+mean_RCM_early(isfinite(mean_RCM_early));

ncname2='D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_7\test06_ESmodel_pollock_1992_02.nc';
h_model=ncread(ncname2, 'h');
lon_model=ncread(ncname2, 'lon_rho');
lat_model=ncread(ncname2, 'lat_rho');

m_pcolor(lon_rho2, lat_rho2, SK_EEZ_mask); 
% m_pcolor(lon_RCM, lat_RCM, SK_EEZ_mask2); 

caxis([1.2 1.3])
shading flat;
colormap([37/255 141/255 239/255; 0.5 0 0]);
% colormap([37/255 141/255 239/255]);
% shading flat;
% m_gshhs_i('color',m_gshhs_line_color);
% m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
m_gshhs_i('color',[198/256 171/256 114/256]);
m_gshhs_i('patch',[198/256 171/256 114/256]);   % gray colored land
% [etopo_con, etopo_h2]=m_contour(lon_rho2, lat_rho2, etopo_z2, [0 -50 -500 -2000 -3000], 'linecolor',[0.8 0.8 0.8], 'linewidth', 0.5);
% 
% clabel(etopo_con,etopo_h2,'FontSize',m_contour_label_fontsize,'Color',[0.8 0.8 0.8], ...
%         'labelspacing',10000000,'Rotation',m_contour_rotation);
[model_con, model_h2]=m_contour(lon_model, lat_model, h_model, [0 50 500 2000 3000], 'linecolor',[0.8 0.8 0.8], 'linewidth', 0.5);

clabel(model_con,model_h2,'FontSize',m_contour_label_fontsize,'Color',[0.8 0.8 0.8], ...
        'labelspacing',10000000,'Rotation',m_contour_rotation);

m_line(SK_EEZ_polygon(:,1), SK_EEZ_polygon(:,2), 'LineStyle', '--', 'Color', 'b', 'LineWidth', 1) 

m_grid('fontsize', m_grid_fontsize-7, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type);
hold off
tifname = 'D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\Fig01_02_08.tif';
print('-dtiff','-r500',tifname); RemoveWhiteSpace([], 'file', tifname);
close all