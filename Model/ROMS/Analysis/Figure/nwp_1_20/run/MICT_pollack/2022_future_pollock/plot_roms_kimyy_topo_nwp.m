% %  Updated 27-Apr-2018 by Yong-Yub Kim
% %  Updated 21-May-2018 by Yong-Yub Kim
% %  Updated 08-Jun-2018 by Yong-Yub Kim
% %  Updated 09-Jan-2023 by Yong-Yub Kim

% % % % % % % plot topography


figdir='/Users/kimyy/Desktop/backup/Research/Ph_D_course/2022_pollock_future/submit_figure/';
bathydir='./';
inputdir='/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/input/test2117/';
testname='test2117';
% if (windows ==1)
%     % % for windows
%     bathydir='Bathy\'; %% SNU_desktop
%     % % set colorbar parameter
%     load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod  % % set colormap (jet_modified)
% elseif (linux==1)
%     % % for linux
%     bathydir='Bathy/'; %% Linux
% end

status=plot_ROMS_bathy([figdir, bathydir, 'bathy_nwp'],inputdir, [115 164 15 52], [0 5000], [0], testname);
close all;





function status=plot_ROMS_bathy(bathyfile, workdir, lonlat, clim, level_c, testname)
% % Updated 08-Jun-2018 by Yong-Yub Kim (get jet_mod colormap for windows)


% addpath(genpath('D:\MEPL\project\NWP\m_map'))
% addpath(genpath('C:\Users\KYY\Dropbox\source\matlab\Common\m_map'));

% load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod

name = 'roms_grid_nwp_1_20_';
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
titlename = 'NWP model bathymetry (1/20^o)';
title(titlename,'fontsize',25);

% set colorbar 
h = colorbar;
% colormap(flip(cool));
colormap(cool);

set(h,'fontsize',20);
title(h,'depth (m)','fontsize',15);
caxis(clim)

[C,h2]=m_contour(lon,lat,depth,level_c,'k','linewidth',1);             
clabel(C,h2,'FontSize',25,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');

% make jpg file
xscale=lonlat(2)-lonlat(1);
yscale=lonlat(4)-lonlat(3);
halt = 1;
while(hal
    if (xscale > 1000 || yscale > 1000)
        halt = 0;
    else
        xscale = xscale * 1.2; yscale = yscale * 1.2;
    end
end
xscale = 800; yscale = 920; %% temporary scale
set(gcf,'Position',[200 100 xscale yscale])
jpgname=strcat(bathyfile,'_1_20.jpg'); % ~workdir\figure\nwp_1_20_bathy.jpg
saveas(gcf,jpgname,'jpg');
RemoveWhiteSpace([], 'file', jpgname);

disp(' ')
disp([' Making bathymetry plot (1/20) is completed.'])
disp(' ')
disp([' File path is : ',jpgname])
disp(' ')
disp(' ')


close all;
status=1;
end