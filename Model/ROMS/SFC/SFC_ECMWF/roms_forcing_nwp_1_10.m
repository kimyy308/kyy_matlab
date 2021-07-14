% % % This script is based on MATLAB R2017a
% % 05-Oct-2018

clear all;close all; clc;
warning off;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
%     dropboxpath='C:\Users\KYY\Dropbox';
%     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_ktotalday']));
    addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

elseif (strcmp(system_name,'GLNXA64'))
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_ktotalday']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Run']));
    addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
end
% start


grdfile='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_10/input/test06/spinup/roms_grid_nwp_1_10_test06.nc';

rawdata_dir='/data1/temp/ECMWF_interim/';
varname={'airT', 'msl', 'dewt', 'ssrd', 'u10', 'v10'};
year=[1981:2015];
tinterval = 1;

testname = 'nwp_1_10';
output_dir=['/data1/kimyy/Model/ROMS/roms_nwp/',testname,'/forcing_matlab/'];
    `
% ncread('/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_2016.nc',ssrd);
for i=1:length(year)
    for j=1:length(varname)
        status = ROMS_surface_forcing_ECMWF(grdfile, year(i), varname{j}, tinterval,output_dir,testname);
    end
    status = ROMS_surface_forcing_Qair(grdfile, year(i), tinterval,output_dir,testname);
end