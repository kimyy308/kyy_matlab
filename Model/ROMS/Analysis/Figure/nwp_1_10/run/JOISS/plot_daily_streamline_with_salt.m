close all; clc; clear all;

dropboxpath='C:\users\user/Dropbox/';
addpath(genpath([dropboxpath '/source/matlab/Common/Streamcolor']));
addpath(genpath([dropboxpath '/source/matlab/Model\ROMS\Roms_tools\Preprocessing_tools']));
addpath(genpath([dropboxpath '/source/matlab/Common\Figure']));

tempyear = 2016;
tempyearstr=num2str(tempyear, '%04i');
for tempday=121:250
% tempday=210;
tempdaystr=num2str(tempday, '%04i');
fs=filesep;
rootdir=['D:\Data', fs, 'Model', fs, 'ROMS', fs, 'nwp_1_10', fs, 'test06', fs, 'DA', fs 'backup_surf_daily'];
figdir=['D:\conference & workshop_study', fs, 'JOISS', fs, 'figure'];
xrange=30:140;
yrange=150:250;

filename.fig=[figdir, fs, 'CJD_', tempyearstr, '_', tempdaystr, '.tif'];
filename.lon_rho=[rootdir, fs, 'pck_ocean_lon_rho.nc'];
filename.lat_rho=[rootdir, fs, 'pck_ocean_lat_rho.nc'];
filename.salt=[rootdir, fs, 'salt', fs, tempyearstr, fs, 'pck_ocean_avg_', tempdaystr, '_salt.nc'];
filename.u=[rootdir, fs, 'u', fs, tempyearstr, fs, 'pck_ocean_avg_', tempdaystr, '_u.nc'];
filename.v=[rootdir, fs, 'v', fs, tempyearstr, fs, 'pck_ocean_avg_', tempdaystr, '_v.nc'];

lon_rho=ncread(filename.lon_rho, 'lon_rho', [min(xrange) min(yrange)], [length(xrange) length(yrange)]);
lat_rho=ncread(filename.lat_rho, 'lat_rho', [min(xrange) min(yrange)], [length(xrange) length(yrange)]);
salt=ncread(filename.salt, 'salt', [min(xrange) min(yrange) 1 1], [length(xrange) length(yrange) 1 1]);
u=ncread(filename.u, 'u', [min(xrange) min(yrange) 1 1], [length(xrange)-1 length(yrange) 1 1]);
v=ncread(filename.v, 'v', [min(xrange) min(yrange) 1 1], [length(xrange) length(yrange)-1 1 1]);
u_rho=u2rho_2d(u')';
v_rho=v2rho_2d(v')';


% % % % minimum value must be 0
minval=28;
salt_30=salt; salt_30(salt_30<=minval)=minval; salt_30(salt_30>=34)=34; salt_30=salt_30-minval;
% pcolor(lon_rho', lat_rho', salt'); shading flat; colorbar;
% pcolor( salt'); shading flat; colorbar;

[lon_proj, lat_proj] = projfwd(defaultm('mercator'),lat_rho,lon_rho);

% hold on;
% 
% shading flat;

data_interval=3;
m_quiver_vector_size=3;
colormap(jet)
% caxis([30 35])
hlines_col=Streamcolor(lon_rho(1:data_interval:end,1:data_interval:end)', ...
                    lat_rho(1:data_interval:end,1:data_interval:end)', ...
                    u_rho(1:data_interval:end,1:data_interval:end)' * m_quiver_vector_size/1.5, ...
                    v_rho(1:data_interval:end,1:data_interval:end)' * m_quiver_vector_size/1.5, ...
                    lon_rho(1:data_interval:end,1:data_interval:end)', ...
                    lat_rho(1:data_interval:end,1:data_interval:end)', ...
                    salt_30(1:data_interval:end,1:data_interval:end)');
% plot_google_map('MapType','hybrid', 'Alpha', 1, 'Scale', 2)
plot_google_map('MapType','satellite', 'Alpha', 1, 'Scale', 2)
set(gca, 'fontsize', 12)
xlabel('Longitude')
ylabel('Latitude')
daystr=datestr(tempday+datenum(tempyear-1,12,31));
text(117,29, daystr, 'fontsize', 16, 'color', 'w')

h=colorbar;
caxis([minval, 34])
h_title=title(h,'psu','fontsize',14);
set(gcf,'PaperPositionMode','auto');

% saveas(gcf,'C:\users\user\Desktop\test.tif','tif');
saveas(gcf,filename.fig,'tif');
RemoveWhiteSpace([], 'file', filename.fig);
close all;
end


% % % % % geographical projection
% % projinv
% hlines_col=Streamcolor(lon_proj(1:data_interval:end,1:data_interval:end)', ...
%                     lat_proj(1:data_interval:end,1:data_interval:end)', ...
%                     u_rho(1:data_interval:end,1:data_interval:end)' * m_quiver_vector_size/1.5, ...
%                     v_rho(1:data_interval:end,1:data_interval:end)' * m_quiver_vector_size/1.5, ...
%                     lon_proj(1:data_interval:end,1:data_interval:end)', ...
%                     lat_proj(1:data_interval:end,1:data_interval:end)', ...
%                     salt_30(1:data_interval:end,1:data_interval:end)');
                
% land_mask=salt;
% land_mask(isfinite(land_mask))=50000;
% land_mask(isnan(land_mask))=1;
% land_mask(land_mask>1)=NaN;
% colormap([0.8 0.8 0.8]);
% pcolor(lon_proj(1:data_interval:end,1:data_interval:end)', ...
%          lat_proj(1:data_interval:end,1:data_interval:end)', ...
%          land_mask(1:data_interval:end,1:data_interval:end)');
% % % %                  pcolor(lon_proj', lat_proj', land_mask');

% xticks(lon_proj(1:10:end,end))                
% xticklabels(round(lon_rho(1:3:end,end),1))
% yticks(lat_proj(end,1:10:end))
% yticklabels(round(lat_rho(end,1:10:end),1))