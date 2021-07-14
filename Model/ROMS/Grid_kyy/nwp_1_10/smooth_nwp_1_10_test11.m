clear all
close all
clc;
%  Updated 17-Nov-2019 by Yong-Yub Kim


dropboxpath='C:\Users\KYY\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Run']));

romstools_param
warning off
rawname=['E:\Data\Model\ROMS\nwp_1_10\input\test11\roms_grid_nwp_1_10_test11_not_smoothed.nc'];
system(['copy ',rawname,' ','E:\Data\Model\ROMS\nwp_1_10\input\test11\roms_grid_nwp_1_10_test11.nc'])

disp(' ')
disp([' Start smoothing the grid: ','E:\Data\Model\ROMS\nwp_1_10\input\test11\roms_grid_nwp_1_10_test11.nc'])
disp(' ')
disp([' Title: ',ROMS_title])
disp(' ')
disp([' Resolution: 1/',num2str(1/dl),' deg'])

ncid=netcdf.open('E:\Data\Model\ROMS\nwp_1_10\input\test11\roms_grid_nwp_1_10_test11.nc', 'WRITE');
sphevarid=netcdf.inqVarID(ncid,'spherical');
sphe=netcdf.getVar(ncid,sphevarid);
if sphe==1
    onedimid=netcdf.inqDimID(ncid,'one')
    netcdf.renameVar(ncid,sphevarid,'spherical2');
    netcdf.defVar(ncid,'spherical', 'NC_CHAR', onedimid);
    sphevarid=netcdf.inqVarID(ncid,'spherical')
    netcdf.putAtt(ncid,sphevarid, 'long_name', 'Grid type logical switch')
    netcdf.putAtt(ncid,sphevarid, 'option_T', 'spherical')
    netcdf.endDef(ncid)
    netcdf.putVar(ncid, sphevarid, 'T')
end
netcdf.close(ncid);
%
%  Smooth the topography
%
% read not smoothed grid
% mask in the not smoothed grid is wrong -> must use mask in the new grids.
nc2=netcdf(rawname,'write');
h_raw=nc2{'h'}(:);
lon_rho=nc2{'lon_rho'}(:);
lat_rho=nc2{'lat_rho'}(:);
maskr=nc2{'mask_rho'}(:);
masku=nc2{'mask_u'}(:);
maskv=nc2{'mask_v'}(:);
maskpsi=nc2{'mask_psi'}(:);
close(nc2);

nc=netcdf('E:\Data\Model\ROMS\nwp_1_10\input\test11\roms_grid_nwp_1_10_test11.nc','write');
h=nc{'h'}(:);
h=h_raw;
nc{'mask_rho'}(:)=maskr;
nc{'mask_u'}(:)=masku;
nc{'mask_v'}(:)=maskv;
nc{'mask_psi'}(:)=maskpsi;

%% korea strait

smootharea=[127.9, 133, 32.3, 36];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.8;  n_filter_deep_topo = 10;  n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);
% pcolor(h_raw(indw:inde,inds:indn))
% shading interp;
h(indw:inde,inds:indn)=h(indw:inde,inds:indn)+40;   %% +20 -> error2

%% tsugaru strait

smootharea=[139.9, 141.5, 40.7, 41.9];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h(indw:inde,inds:indn)+40;   %% + 20 -> error(spinup) +0 -> error2

%% soya strait

smootharea=[141.6, 142.3, 45.3, 46.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h(indw:inde,inds:indn)+40;   %% +35 -> error2



% Slope parameter (r=grad(h)/h) maximum value for topography smoothing
% - rtarget

% Number of pass of a selective filter to reduce the isolated
% seamounts on the deep ocean     -      n_filter_deep_topo

% Number of pass of a single hanning filter at the end of the
% smooting procedure to ensure that there is no 2DX noise in the 
% topography   -     n_filter_final


% %% EKWC path smoothing
% smootharea=[129.3, 135.0, 35.2, 38.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.6;  n_filter_deep_topo = 2;  n_filter_final = 4;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

%% Kuril islands smoothing
smootharea=[150.0, 162.0, 44.8, 52];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.7; n_filter_deep_topo = 30; n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[150.0, 157.5, 44.3, 48.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.7; n_filter_deep_topo = 30; n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[150.0, 157.5, 44.3, 48.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.7; n_filter_deep_topo = 30; n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[150.4, 151.5, 46.1, 47.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.7; n_filter_deep_topo = 30; n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[153.6, 155.5, 48.5, 50.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.7; n_filter_deep_topo = 30; n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[156.9, 158.0, 51.1, 51.7];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.7; n_filter_deep_topo = 30; n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);


%% Luzon strait
% smootharea=[119.9, 122.5, 18.5, 22.0];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% for i=inds:indn
%     for j = inds:indn
%         if (h(i,j)<2000)
%             h(i,j) = h(i,j) + 50;
%         end
%     end
% end

%% Luzon strait smoothing
% smootharea=[118.5, 122.0, 18.3, 20.9];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.8; n_filter_deep_topo = 50; n_filter_final = 10;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);
% 
% smootharea=[118.5, 122.0, 21.1, 23.0];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.8; n_filter_deep_topo = 50; n_filter_final = 10;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);
% 
% smootharea=[118.5, 122.0, 17.8, 23.0];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.7; n_filter_deep_topo = 50; n_filter_final = 10;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[115.0, 127.5, 15.0, 23.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.5; n_filter_deep_topo = 0; n_filter_final = 10;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

% % %  West japan
smootharea=[136.7, 139.0, 36.6, 38.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.6; n_filter_deep_topo = 10; n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[136.7, 141.0, 36.0, 39.6];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.6; n_filter_deep_topo = 10; n_filter_final = 5;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);


%% Whole domain smoothing
smootharea=[115.0, 162.0, 15.0, 52.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);

% % n_filter_deep_topo --> 4, n_filter_filnal --> 2  are default

rtarget = 0.8; n_filter_deep_topo = 2; n_filter_final = 4;  
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

% 


%% Kuroshio path region (east japan)
smootharea=[140.0, 143.0, 35.5, 39.2];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

smootharea=[138.7, 141.0, 26.6, 32.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

% smootharea=[138.0, 144.0, 26.0, 35.0];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.4; n_filter_deep_topo = 2; n_filter_final = 4;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);


% smootharea=[138.0, 144.0, 36.0, 36.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.7; n_filter_deep_topo = 20; n_filter_final = 4;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);
% 
% smootharea=[139.0, 144.0, 36.5, 37.0];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.8; n_filter_deep_topo = 15; n_filter_final = 4;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);
% 
% smootharea=[140.0, 144.0, 37.0, 37.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.9; n_filter_deep_topo = 10; n_filter_final = 4;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);
% 
% smootharea=[140.0, 144.0, 37.5, 38.0];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.9; n_filter_deep_topo = 5; n_filter_final = 4;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[138.0, 144.0, 36.0, 36.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.8; n_filter_deep_topo = 4; n_filter_final = 2;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[139.0, 144.0, 36.5, 37.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.8; n_filter_deep_topo = 4; n_filter_final = 2;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[140.0, 144.0, 37.0, 37.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.9; n_filter_deep_topo = 4; n_filter_final = 2;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[140.0, 144.0, 37.5, 38.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.9; n_filter_deep_topo = 4; n_filter_final = 2;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[140.5, 144.0, 38.0, 39.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.9; n_filter_deep_topo = 4; n_filter_final = 2;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

% smootharea=[138.0, 144.0, 26.0, 36];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.9; n_filter_deep_topo = 4; n_filter_final = 2;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[138.0, 144.0, 26.0, 35];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.9; n_filter_deep_topo = 20; n_filter_final = 2;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[140.0, 140.7, 34.7, 35.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.3; n_filter_deep_topo = 4; n_filter_final = 2;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[138.0, 144.0, 33.0, 35.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.5; n_filter_deep_topo = 4; n_filter_final = 2;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);


%% EKB region
smootharea=[127.0, 129.4, 37.2, 41.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

% smootharea=[127.0, 129.4, 37.0, 41.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.99;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[127.0, 129.0, 37.2, 41.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

% smootharea=[127.0, 129.2, 37.0, 41.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.99;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[127.0, 128.8, 37.2, 41.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

% smootharea=[127.0, 129.0, 37.4, 41.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.99;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[127.0, 128.6, 37.2, 41.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

% smootharea=[127.0, 128.8, 37.6, 41.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.99;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[127.0, 128.4, 37.2, 41.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

% smootharea=[127.0, 128.6, 37.8, 41.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.99;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[127.0, 128.2, 37.2, 41.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

% smootharea=[127.0, 128.4, 38.0, 41.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.99;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[127.0, 129.7, 37.2, 41.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

% smootharea=[127.0, 128.2, 38.2, 41.5];
% [indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
% rtarget = 0.99;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
% h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
%   rtarget, n_filter_deep_topo, n_filter_final);

smootharea=[127.0, 129.7, 36.5, 41.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
h(indw:inde,inds:indn)=h_raw(indw:inde,inds:indn);

smootharea=[127.0, 129.7, 36.5, 41.5];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.99;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);

% % % Ulleung Island smoothing
smootharea=[129.2, 135, 36.0, 40.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);
rtarget = 0.8;  n_filter_deep_topo = 1;  n_filter_final = 0.1;
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
  rtarget, n_filter_deep_topo, n_filter_final);



%% Final whole domain smoothing
smootharea=[115.0, 162.e0, 15.0, 52.0];
[indw, inde, inds, indn]=findind_Y(dl,smootharea,lon_rho,lat_rho);

% % n_filter_deep_topo --> 4, n_filter_filnal --> 2  are default

% rtarget = 0.8; n_filter_deep_topo = 5; n_filter_final = 5;  %% 0.9 --> blow up 0.8 --> blow up 0.7 --> blow up  0.6 --> blowup(67d)
rtarget = 0.95; n_filter_deep_topo = 1; n_filter_final = 0.1;  
h(indw:inde,inds:indn)=smoothgrid(h(indw:inde,inds:indn),maskr(indw:inde,inds:indn),hmin,hmax_coast,hmax, ...
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


figure;
lonlat = [127 144 33 52]; %% East Sea
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
hold on;
m_pcolor(lon_rho,lat_rho,h);
shading interp;
% set colorbar 
ccc= colorbar;
colormap(jet);
caxis([0 500]);
hold off;


% figure;
% lonlat = [115 162 15 52]; %% nwp
% m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
% hold on;
% m_pcolor(lon_rho,lat_rho,h);
% shading interp;
% % set colorbar 
% ccc= colorbar;
% load jet_mod
% colormap(jet_mod);
% caxis([0 5000]);
% hold off;



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