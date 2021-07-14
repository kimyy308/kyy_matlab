% % % This script is based on MATLAB R2017a

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


% grdfile='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_10/input/test03/spinup/2017/roms_grid_nwp_1_10_test03.nc';
grdfile='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/input/test37/2001/roms_grid_combine2_test37.nc';

% rawdata_dir='/data1/temp/ECMWF_interim/';
% varname={'tas', 'psl', 'hur', 'rsds', 'uo_wind', 'vo_wind'};
% varname={'tas', 'psl', 'hur', 'rsds', 'uo_wind', 'vo_wind'};
% varname={'vo_wind'};
varname={'tas', 'psl', 'hur', 'rsds', 'uo_wind'};


year=[1976];
tinterval = 1;

CMIP5_name='total';
% CMIP5_name='CanESM2';
% CMIP5_name='NorESM-M';
% CMIP5_name='MPI-ESM-MR';

testname = [CMIP5_name,'_nwp_1_20'];
output_dir=['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20','/forcing_matlab/'];
    
% ncread('/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_2016.nc',ssrd);
for i=1:length(year)
    for j=1:length(varname)
        status = ROMS_surface_forcing_CMIP5(grdfile, year(i), varname{j}, tinterval,output_dir,testname, CMIP5_name);
    end
%     status = ROMS_surface_forcing_Qair(grdfile, year(i), tinterval,output_dir,testname);
end