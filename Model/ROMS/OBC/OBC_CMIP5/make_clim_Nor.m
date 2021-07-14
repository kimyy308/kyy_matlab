clear all; close all; clc;


hisdir=['/data1/CMIP/cmip5/historical_extHDD/CMIP5/zos/historical/interp/NorESM1-M/'];
rcpdir=['/data1/CMIP/cmip5/rcp_ocn_extHDD/CMIP5/zos/rcp85/interp/NorESM1-M/'];


inputyear=1993:2019
for year=inputyear(1):inputyear(end)
    if year<=2005
        filename=[hisdir, 'zos_interp_NorESM1-M_historical_r1i1p1_', num2str(year), '.nc'];
        zos=ncread(filename, 'zos');
    else
        filename=[rcpdir, 'zos_interp_NorESM1-M_rcp85_r1i1p1_', num2str(year), '.nc'];
        zos=ncread(filename, 'zos');
    end
    if (exist('clim_zos')==0)
        clim_zos=zeros(size(zos));
    end
    clim_zos=clim_zos+zos/length(inputyear);
end

ncoutfilename=['/home/kimyy/SSH_project_source/data','/', ...
    'NorESM1-M', '_', 'climssh_', num2str(inputyear(1), '%04i'), '_', num2str(inputyear(end)), '.nc'];
% %%         make ncfile
ncid = netcdf.create(ncoutfilename,'NETCDF4');
xi_rho_dimid = netcdf.defDim(ncid, 'lon_rho', size(clim_zos,1));
eta_rho_dimid = netcdf.defDim(ncid,'lat_rho', size(clim_zos,2));
time_dimid = netcdf.defDim(ncid, 'time', 12);
clim_sshvarid=netcdf.defVar(ncid, 'clim_ssh', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid time_dimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [size(clim_zos,1) size(clim_zos,2) 12], clim_zos);
netcdf.close(ncid);


mean_zos=mean(clim_zos,3);
ncoutfilename=['/home/kimyy/SSH_project_source/data','/', ...
    'NorESM1-M', '_', 'meanssh_', num2str(inputyear(1), '%04i'), '_', num2str(inputyear(end)), '.nc'];
% %%         make ncfile
ncid = netcdf.create(ncoutfilename,'NETCDF4');
xi_rho_dimid = netcdf.defDim(ncid, 'lon_rho', size(clim_zos,1));
eta_rho_dimid = netcdf.defDim(ncid,'lat_rho', size(clim_zos,2));
time_dimid = netcdf.defDim(ncid, 'time', 1);
clim_sshvarid=netcdf.defVar(ncid, 'clim_ssh', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid time_dimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [size(clim_zos,1) size(clim_zos,2) 1], mean_zos);
netcdf.close(ncid);
