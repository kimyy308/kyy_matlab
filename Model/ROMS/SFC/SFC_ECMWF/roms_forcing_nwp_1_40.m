% % % This script is based on MATLAB R2017a

clear all;close all; clc;
warning off;
linux=1; windows=0;
if (windows ==1)
    % % for windows
%     dropboxpath='C:\Users\KYY\Dropbox';
%     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_ktotalday']));
elseif (linux==1)
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_ktotalday']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Run']));
end
% start


% grdfile='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_10/input/test03/spinup/2017/roms_grid_nwp_1_10_test03.nc';
grdfile='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_40/input/test01/spinup/1980/roms_grid_nwp_1_40_test01.nc';

rawdata_dir='/data1/temp/ECMWF_interim/';
varname={'airT', 'msl', 'dewt', 'ssrd', 'u10', 'v10'};
year=[2010];
tinterval = 1;

testname = 'nwp_1_40';
output_dir=['/data1/kimyy/Model/ROMS/roms_nwp/',testname,'/forcing_matlab/'];
    
% ncread('/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_2016.nc',ssrd);
for i=1:length(year)
    for j=1:length(varname)
        status = ROMS_surface_forcing_ECMWF(grdfile, year(i), varname{j}, tinterval,output_dir,testname);
    end
    status = ROMS_surface_forcing_Qair(grdfile, year(i), tinterval,output_dir,testname);
end