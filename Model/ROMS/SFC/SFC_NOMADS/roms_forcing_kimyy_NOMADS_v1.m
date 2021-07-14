% % % This script is based on MATLAB R2017a
% % Updated 14-Feb-2019 by Y.Y.Kim
% % NOMADS CFSR(v1) to ROMS Surface forcing file
% % NOMADS CFSR(v1) --> available for 01 Jan 1979 - 31 Mar 2011
% % 
% % dataset description
% % https://www.ncdc.noaa.gov/data-access/model-data/model-datasets/climate-forecast-system-version2-cfsv2
% % 
% % ftp link
% % ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/
% % 
% % and we make roms forcing files with that netcdf file 
% % orginated from grib2 file using kwgrib
% % spatial domain : global, 0.5 degree resolution
% % time frequency : 1 hour

clear all;close all; clc;
warning off;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
% % for windows
elseif (strcmp(system_name,'GLNXA64'))
% % for linux server
    addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
end
% start

expname='nwp_1_10';
testname='test07';

grdfile=['/data1/kimyy/Model/ROMS/',expname, ...
  '/input/',testname,'/spinup/roms_grid_',expname,'_',testname,'.nc'];
disp(['grdfile : ', grdfile])

% varname={'TMP_2maboveground',  ...  %% 2m temperature
%     'PRES_surface', ...   %% surface atmospheric pressure
%     'SPFH_2maboveground', ...  %% humidity
%     'DSWRF_surface', ...  %% downward short wave radiation
%     'UGRD_10maboveground', ...  %% U direction wind
%     'VGRD_10maboveground'};  %% V direction wind

% varname={'PRES_surface'};
varname={'UGRD_10maboveground', ...  %% U direction wind
    'VGRD_10maboveground'};  %% V direction wind

year=[1993:1993];
tinterval = 1;  %% for 1 day, do not touch it

atm_name='gdas';
expname_atm = [atm_name,'_',expname];
output_dir=['/data1/kimyy/Model/ROMS/',expname,'/forcing_matlab/NOMADS/'];

for i=1:length(year)
    for j=1:length(varname)
        status = ROMS_surface_forcing_NOMADS_v1(grdfile, year(i), varname{j}, ...
            tinterval,output_dir,expname_atm, atm_name);
    end
end