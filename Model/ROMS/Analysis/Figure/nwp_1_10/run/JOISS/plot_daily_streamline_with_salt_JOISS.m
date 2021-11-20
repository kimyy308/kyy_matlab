%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% program plot_daily_streamline_with_salt_JOISS
%
% plot colored stremline with salt pcolor for animation
%
%  required external functions:
%  plot_google_map
%  RemoveWhiteSpace
%  Streamcolor
%  u2rho_2d
%  v2rho_2d
% 
%  required variables:
%  configs              configuration for data and figures
%
%  output:
%  tif figures          tif figures
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    07-Oct-2021 by Yong-Yub Kim
%  Updated    12-Oct-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear all;
fs=filesep;

% dropboxpath='C:\users\user/Dropbox/';
% addpath(genpath([dropboxpath '/source/matlab/Common/Streamcolor']));
% addpath(genpath([dropboxpath '/source/matlab/Model\ROMS\Roms_tools\Preprocessing_tools']));
% addpath(genpath([dropboxpath '/source/matlab/Common\Figure']));

% % configuration
configs.rootdir=['D:\Data', fs, 'Model', fs, 'ROMS', fs, 'nwp_1_10', fs, 'test06', fs, 'DA', fs 'backup_surf_daily'];  % rootdir of data
configs.figdir=['D:\conference & workshop_study', fs, 'JOISS', fs, 'figure'];  % directory where figures will be saved
configs.year = 2016; % years of data
configs.day = 121:250; % days of data
configs.xrange=30:140;  % xrange of data
configs.yrange=150:250;  % yrange of data
configs.minval=28;  %% minimum salinity for colored streamline
configs.maxval=34;  %% maximum salinity for colored streamline
configs.data_interval=3; % data_interval for quiver
configs.m_quiver_vector_size=3; % vector size for quiver
configs.figcmap=jet; % colormap for figures
configs.text_lon=117;
configs.text_lat=29;

for yeari=1:length(configs.year) % yearly loop
tmps.tempyear=configs.year(yeari);
tmps.tempyearstr=num2str(tmps.tempyear, '%04i');

    for dayi=1:length(configs.day)   % daily loop
    tmps.tempday=configs.day(dayi);
    tmps.tmps.temptmps.daystr=num2str(tmps.tempday, '%04i');

%     get filenames
    tmps.filename_fig=[configs.figdir, fs, 'CJD_', tmps.tempyearstr, '_', tmps.tmps.temptmps.daystr, '.tif'];
    tmps.filename_grid.lon_rho=[configs.rootdir, fs, 'pck_ocean_lon_rho.nc'];
    tmps.filename_grid.lat_rho=[configs.rootdir, fs, 'pck_ocean_lat_rho.nc'];
    tmps.filename_salt=[configs.rootdir, fs, 'salt', fs, tmps.tempyearstr, fs, 'pck_ocean_avg_', tmps.tmps.temptmps.daystr, '_salt.nc'];
    tmps.filename_u=[configs.rootdir, fs, 'u', fs, tmps.tempyearstr, fs, 'pck_ocean_avg_', tmps.tmps.temptmps.daystr, '_u.nc'];
    tmps.filename_v=[configs.rootdir, fs, 'v', fs, tmps.tempyearstr, fs, 'pck_ocean_avg_', tmps.tmps.temptmps.daystr, '_v.nc'];

%     get data
    grid.lon_rho=ncread(tmps.filename_grid.lon_rho, 'lon_rho', [min(configs.xrange) min(configs.yrange)], [length(configs.xrange) length(configs.yrange)]);
    grid.lat_rho=ncread(tmps.filename_grid.lat_rho, 'lat_rho', [min(configs.xrange) min(configs.yrange)], [length(configs.xrange) length(configs.yrange)]);
    data.salt=ncread(tmps.filename_salt, 'salt', [min(configs.xrange) min(configs.yrange) 1 1], [length(configs.xrange) length(configs.yrange) 1 1]);
    data.u=ncread(tmps.filename_u, 'u', [min(configs.xrange) min(configs.yrange) 1 1], [length(configs.xrange)-1 length(configs.yrange) 1 1]);
    data.v=ncread(tmps.filename_v, 'v', [min(configs.xrange) min(configs.yrange) 1 1], [length(configs.xrange) length(configs.yrange)-1 1 1]);
    data.u_rho=u2rho_2d(data.u')';
    data.v_rho=v2rho_2d(data.v')';

    % % % % set salinity anomaly for colored streamline
    data.salt_30=data.salt; 
    data.salt_30(data.salt_30<=configs.minval)=configs.minval; 
    data.salt_30(data.salt_30>=configs.maxval)=configs.maxval; 
    data.salt_30=data.salt_30-configs.minval;
    
    colormap(configs.figcmap);
    figs.hlines_col=Streamcolor(grid.lon_rho(1:configs.data_interval:end,1:configs.data_interval:end)', ...
                        grid.lat_rho(1:configs.data_interval:end,1:configs.data_interval:end)', ...
                        data.u_rho(1:configs.data_interval:end,1:configs.data_interval:end)' * configs.m_quiver_vector_size/1.5, ...
                        data.v_rho(1:configs.data_interval:end,1:configs.data_interval:end)' * configs.m_quiver_vector_size/1.5, ...
                        grid.lon_rho(1:configs.data_interval:end,1:configs.data_interval:end)', ...
                        grid.lat_rho(1:configs.data_interval:end,1:configs.data_interval:end)', ...
                        data.salt_30(1:configs.data_interval:end,1:configs.data_interval:end)');
    plot_google_map('MapType','satellite', 'Alpha', 1, 'Scale', 2)
    set(gca, 'fontsize', 12);
    xlabel('Longitude');
    ylabel('Latitude');
    tmps.daystr=datestr(tmps.tempday+datenum(tmps.tempyear-1,12,31));
    text(configs.text_lon, configs.text_lat, tmps.daystr, 'fontsize', 16, 'color', 'w');

    figs.h=colorbar;
    caxis([configs.minval, configs.maxval]);
    figs.h_title=title(figs.h,'psu','fontsize',14);
    set(gcf,'PaperPositionMode','auto');

    saveas(gcf,tmps.filename_fig,'tif');
    RemoveWhiteSpace([], 'file', tmps.filename_fig);
    close all;
    end  % daily loop
    
end  % yearly loop