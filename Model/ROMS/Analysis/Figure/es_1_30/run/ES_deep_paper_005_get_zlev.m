close all; clear all; clc;

warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\user\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
end

testname='v52';
workdir = strcat('D:\Data\Model\ROMS\es_1_30\', testname, '\output\run\short_monthly\');
% workdir =[inputdir, 'input\', testname, '\'];

% load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod
% % name = 'roms_grid_combine2_';
name = 'grid_deep_eastsea_';
filename_suffix = '.nc';
% ex : ~workdir\roms_grid_combine2_test37.nc
filename = strcat(workdir,name,testname,filename_suffix);
outputname = [workdir, '\1994\short_v52_monthly_1994_01.nc'];

grd.lon_rho = ncread(filename,'lon_rho')';
grd.lat_rho = ncread(filename,'lat_rho')';
grd.mask_rho = ncread(filename,'mask_rho')';
grd.h = ncread(filename,'h')';


Vtransform = ncread(outputname, 'Vtransform');
Vstretching = ncread(outputname, 'Vstretching');
theta_s = ncread(outputname, 'theta_s');
theta_b =ncread(outputname, 'theta_b');
hc = ncread(outputname, 'hc');
N = length(ncread(outputname, 's_rho'));

% h = grd.h(495,281:381);  %   129E ~ 132E, 37N 

% zeta  
% grd.zeta = zeros(size(grd.h));
grd.zeta = ncread(outputname, 'zeta')';

temp_z_r=zlevs(Vtransform, Vstretching, grd.h,grd.zeta,theta_s,theta_b,hc,N,'r');
z_r=permute(temp_z_r, [3 2 1]);
temp_z_w=zlevs(Vtransform, Vstretching, grd.h,grd.zeta,theta_s,theta_b,hc,N+1,'w');
z_w=permute(temp_z_w, [3 2 1]);

ncoutfilename = strcat(workdir,'z_deep_eastsea_1994_01',testname,filename_suffix);

len_xi=size(z_r,1);
len_eta=size(z_r,2);
len_s_rho=size(z_r,3);
len_s_w=size(z_w,3);
ncid = netcdf.create(ncoutfilename,'NETCDF4');
xi_dimid = netcdf.defDim(ncid, 'xi', len_xi);
eta_dimid = netcdf.defDim(ncid,'eta', len_eta);
s_rho_dimid = netcdf.defDim(ncid,'s_rho', len_s_rho);
s_w_dimid = netcdf.defDim(ncid,'s_w', len_s_w);

z_rvarid=netcdf.defVar(ncid, 'z_r', 'NC_DOUBLE', [xi_dimid eta_dimid s_rho_dimid]);
netcdf.putAtt(ncid,z_rvarid,'long_name','z_r');
netcdf.putAtt(ncid,z_rvarid,'units','m');

z_wvarid=netcdf.defVar(ncid, 'z_w', 'NC_DOUBLE', [xi_dimid eta_dimid s_w_dimid]);
netcdf.putAtt(ncid,z_wvarid,'long_name','z_w');
netcdf.putAtt(ncid,z_wvarid,'units','m');

netcdf.endDef(ncid);

netcdf.putVar(ncid, z_rvarid, [0 0 0], [len_xi len_eta len_s_rho], z_r);
netcdf.putVar(ncid, z_wvarid, [0 0 0], [len_xi len_eta len_s_w], z_w);

netcdf.close(ncid);

% % %   grd.zeta = zeta;


% for theta_s = 10:10
%     for theta_b = 0:0.1:0.1
%         for hc = [5 100 250]
%             grd.z_r=squeeze(zlevs(h,grd.zeta,theta_s,theta_b,hc,N,'r'));  
%             grd.z_w=squeeze(zlevs(h,grd.zeta,theta_s,theta_b,hc,N+1,'w'));  