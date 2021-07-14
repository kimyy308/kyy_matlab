clear all
close all

%  Updated 10-Jul-2018 by Yong-Yub Kim


dropboxpath='C:\Users\KYY\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Run']));


romstools_param
warning off
rawname=['E:\Data\Model\ROMS\nwp_1_20\input\roms_grid_nwp_1_20_not_smoothed.nc'];
system(['copy ','E:\Data\Model\ROMS\nwp_1_20\input\test42\roms_grid_nwp_1_20_test42.nc',' ',smoothname])

disp(' ')
disp([' Start smoothing the grid: ',smoothname])
disp(' ')
disp([' Title: ',ROMS_title])
disp(' ')
disp([' Resolution: 1/',num2str(1/dl),' deg'])


%
%  Smooth the topography
%
% read not smoothed grid
% mask in the not smoothed grid is wrong -> must use mask in the new grids.
nc2=netcdf(rawname,'write');
h_raw=nc2{'h'}(:);
lon_rho=nc2{'lon_rho'}(:);
lat_rho=nc2{'lat_rho'}(:);
close(nc2);

nc=netcdf(smoothname,'write');
% h=nc{'h'}(:);
h=h_raw;
maskr=nc{'mask_rho'}(:);


%% korea strait

smootharea=[127.9, 131.5, 32.3, 37.6];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% pcolor(h_raw(inds:indn,indw:inde))
% shading interp;
h(inds:indn,indw:inde)=h(inds:indn,indw:inde)+10;

%% tsugaru strait

smootharea=[139.9, 141.5, 40.7, 41.9];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(inds:indn,indw:inde)=h(inds:indn,indw:inde)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait

%% soya strait

smootharea=[141.6, 142.3, 45.3, 46.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(inds:indn,indw:inde)=h(inds:indn,indw:inde)+25; % there are too much transport in the tsugaru strait and negative transport in the soya strait



% Slope parameter (r=grad(h)/h) maximum value for topography smoothing
% - rtarget

% Number of pass of a selective filter to reduce the isolated
% seamounts on the deep ocean     -      n_filter_deep_topo

% Number of pass of a single hanning filter at the end of the
% smooting procedure to ensure that there is no 2DX noise in the 
% topography   -     n_filter_final


%% EKWC path smoothing
smootharea=[128.5, 135.0, 35.2, 38.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.8;  n_filter_deep_topo = 30;  n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

%% Kuril islands smoothing
smootharea=[150.0, 164.0, 44.8, 52];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.2; n_filter_deep_topo = 30; n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[150.0, 157.5, 44.3, 48.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.2; n_filter_deep_topo = 30; n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[150.0, 157.5, 44.3, 48.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.2; n_filter_deep_topo = 30; n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[150.4, 151.5, 46.1, 47.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.2; n_filter_deep_topo = 30; n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[153.6, 155.5, 48.5, 50.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.2; n_filter_deep_topo = 30; n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[156.9, 158.0, 51.1, 51.7];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.2; n_filter_deep_topo = 30; n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);


%% izu-ogasawara ridge
smootharea=[134.9, 135.7, 36.0, 37.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
for i=inds:indn
    for j = indw:inde
         if (h(i,j)>300) 
             h(i,j)= h(i,j)+200.0;
         end
    end
end

smootharea=[138.4, 140.0, 33.7, 34.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
for i=inds:indn
    for j = indw:inde
         h(i,j)= 2000.0;
    end
end

%% izu-ogasawara ridge (south japan) smoothing
smootharea=[138.0, 143.0, 30.5, 34.8];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.7; n_filter_deep_topo = 30; n_filter_final = 3;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

%% Luzon strait
smootharea=[119.9, 122.5, 18.5, 22.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
for i=indw:inde
    for j = inds:indn
        if (h(i,j)<2000)
            h(i,j) = h(i,j) + 50;
        end
    end
end

%% Luzon strait smoothing
smootharea=[118.5, 122.0, 18.3, 20.9];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.8; n_filter_deep_topo = 100; n_filter_final = 10;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[118.5, 122.0, 21.1, 23.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.8; n_filter_deep_topo = 100; n_filter_final = 10;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[118.5, 122.0, 17.8, 23.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.7; n_filter_deep_topo = 50; n_filter_final = 10;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[115.0, 127.5, 15.0, 23.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.5; n_filter_deep_topo = 0; n_filter_final = 10;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

%% Whole domain smoothing
smootharea=[115.0, 164.0, 15.0, 52.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.6; n_filter_deep_topo = 0; n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

%% Kuroshio path region (east japan)
smootharea=[140.0, 143.0, 35.5, 39.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(inds:indn,indw:inde)=h_raw(inds:indn,indw:inde);

smootharea=[138.7, 141.0, 26.6, 32.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(inds:indn,indw:inde)=h_raw(inds:indn,indw:inde);

smootharea=[138.0, 144.0, 26.0, 40.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.9; n_filter_deep_topo = 5; n_filter_final = 2;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

%% EKB region
smootharea=[127.0, 128.8, 37.7, 41.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(inds:indn,indw:inde)=h_raw(inds:indn,indw:inde);

smootharea=[127.0, 129, 37.5, 41.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.9;  n_filter_deep_topo = 30;  n_filter_final = 5;
h(inds:indn,indw:inde)=smoothgrid(h(inds:indn,indw:inde),maskr(inds:indn,indw:inde),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);


%
%  Write it down
%
disp(' ')
disp(' Write it down...')
nc{'h'}(:)=h;
close(nc);

% 
%  plot
% 
% pcolor(h)

% figure;

% pcolor(h(700:920,700:980))
% pcolor(h(700:800,700:850))


% lonlat = [132 143 27 35];  %%kuro, japan south
% m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
% hold on;
% m_pcolor(lon_rho,lat_rho,h);
% shading interp;
% m_gshhs_i('color','k')
% m_gshhs_i('patch',[.8 .8 .8]);  
% % set colorbar 
% ccc= colorbar;
% load jet_mod
% colormap(jet_mod);
% caxis([0 5000]);
% hold off;

% figure;
% lonlat = [117 124 16 23]; %% luzon
% m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
% hold on;
% m_pcolor(lon_rho,lat_rho,h);
% shading interp;
% % set colorbar 
% ccc= colorbar;
% load D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine\jet_mod
% colormap(jet_mod);
% caxis([0 5000]);
% hold off;


% figure;
% lonlat = [127 144 33 52]; %% East Sea
% m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
% hold on;
% m_pcolor(lon_rho,lat_rho,h);
% shading interp;
% % set colorbar 
% ccc= colorbar;
% load D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine\jet_mod
% colormap(jet_mod);
% caxis([0 5000]);
% hold off;


figure;
lonlat = [115 164 15 52]; %% nwp
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
hold on;
m_pcolor(lon_rho,lat_rho,h);
shading interp;
% set colorbar 
ccc= colorbar;
load jet_mod
colormap(jet_mod);
caxis([0 5000]);
hold off;



% figure; 
% pcolor(450:550,350:450,h(350:450,450:550))  %% izu ridge
% shading interp;
% colormap jet;
% colorbar;
% 
% figure;
% pcolor(h(1:170,1:250))  %% luzon
% shading interp
% colormap jet;
% colorbar;
% 
% figure;
% pcolor(h);
% shading interp;
% colormap jet;
% colorbar;