function status=plot_ROMS_bathy(bathyfile, workdir, lonlat, clim, level_c, testname)
% % Updated 08-Jun-2018 by Yong-Yub Kim (get jet_mod colormap for windows)


% addpath(genpath('D:\MEPL\project\NWP\m_map'))
% addpath(genpath('C:\Users\KYY\Dropbox\source\matlab\Common\m_map'));

load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod

name = 'roms_grid_nwp_1_10_';
filename_suffix = '.nc';
% ex : ~workdir\roms_grid_combine2_test37.nc
filename = strcat(workdir,name,testname,filename_suffix);

%% bathymetry (1/20)
% read data
lon = ncread(filename,'lon_rho');
lat = ncread(filename,'lat_rho');
depth = ncread(filename,'h');

% plot
figure

m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_grid('fontsize',20, 'box', 'fancy');
hold on;
m_pcolor(lon,lat,depth);
shading interp;
m_gshhs_i('color','k');
m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land
% m_gshhs_h('color','k')  
% m_gshhs_h('patch',[.8 .8 .8]);   % gray colored land
titlename = 'NWP model bathymetry (1/10^o)';
title(titlename,'fontsize',25);

% set colorbar 
h = colorbar;
colormap(jet_mod);
set(h,'fontsize',20);
title(h,'depth (m)','fontsize',15);
caxis(clim)

[C,h2]=m_contour(lon,lat,depth,level_c,'k','linewidth',1);             
clabel(C,h2,'FontSize',25,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');

% make jpg file
xscale=lonlat(2)-lonlat(1);
yscale=lonlat(4)-lonlat(3);
halt = 1;
while(halt)
    if (xscale > 1000 || yscale > 1000)
        halt = 0;
    else
        xscale = xscale * 1.2; yscale = yscale * 1.2;
    end
end
xscale = 800; yscale = 920; %% temporary scale
set(gcf,'Position',[200 100 xscale yscale])
jpgname=strcat(bathyfile,'_1_10.jpg'); % ~workdir\figure\nwp_1_10_bathy.jpg
saveas(gcf,jpgname,'jpg');


disp(' ')
disp([' Making bathymetry plot (1/10) is completed.'])
disp(' ')
disp([' File path is : ',jpgname])
disp(' ')
disp(' ')



% % %% bathymetry (1/10)
% % % read data
% % filename = strcat('D:\MEPL\project\NWP\1_10_nwp_seo\','roms_grid_final',filename_suffix);
% % lon = ncread(filename,'lon_rho');
% % lat = ncread(filename,'lat_rho');
% % h = ncread(filename,'h');
% % 
% % % plot
% % figure
% % m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % m_grid('fontsize',20);
% % hold on;
% % m_pcolor(lon,lat,h);
% % shading interp;
% % m_gshhs_i('color','k')  
% % m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land
% % titlename = 'Bathymetry of the NWP model (1/10^o)';
% % title(titlename,'fontsize',25);
% % 
% % % set colorbar 
% % h = colorbar;
% % colormap(jet_mod);
% % set(h,'fontsize',20);
% % title(h,'depth (m)','fontsize',15);
% % 
% % % make jpg file
% % set(gcf,'Position',[200 100 980 920])
% % jpgname=strcat(bathyfile,'_1_10.jpg'); % ~workdir\figure\nwp_1_10_bathy.jpg
% % saveas(gcf,jpgname,'jpg');
% % 
% % disp(' ')
% % disp([' Making bathymetry plot (1/10) is completed.'])
% % disp(' ')
% % disp([' File path is : ',jpgname])
% % disp(' ')
% % disp(' ')
% % 
% % 


close all;
status=1; 
return