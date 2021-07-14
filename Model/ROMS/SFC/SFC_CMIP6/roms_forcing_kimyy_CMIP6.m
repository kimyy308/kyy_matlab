% % % This script is based on MATLAB R2019a
% % Updated 17-Mar-2021 by Y.Y.Kim (CMIP6 update)


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


grdfile='/home/kimyy/Model/ROMS/nwp_1_20/input/test2102/spinup/roms_grid_nwp_1_20_test2102.nc';

% rawdata_dir='/data1/temp/ECMWF_interim/';
varname={'tas', 'psl', 'huss', 'rsds', 'uas', 'vas'};
% varname={'psl', 'huss', 'rsds', 'uas', 'vas'};
% varname={'tas'};

% varname={'huss', 'rsds', 'uas', 'vas'};

% varname={'ua','va'};
% varname={'tas', 'psl', 'rsds'};

% hur -> CMCC-CM2-HR4
% huss -> 'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR'

% year=[1976];
year=[2093:2100];
scenario_name='ssp585'; % historical, ssp585
% ensemble_name='r1i1p1f1';  % CMCC-CM2-HR4, EC-Earth3-Veg, ACCESS-CM2, 
% ensemble_name='r1i1p1f2';  % CNRM-ESM2-1, CNRM-CM6-1-HR
tinterval = 1;

% CMIP6_name='ens10';
% CMIP6_names={'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
% CMIP6_names={'MPI-ESM-LR'};

% CMIP6_names = {'CMCC-CM2-HR4', 'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR'};
CMIP6_names = {'CMCC-ESM2', 'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR'};

% CMIP6_names = {'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR'};
% CMIP6_names = {'CMCC-ESM2'};

   
for i=1:length(year)
    for k=1:length(CMIP6_names)
        for j=1:length(varname)
            CMIP6_name=CMIP6_names{k};
            switch CMIP6_name
                case {'CMCC-CM2-HR4', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CMCC-ESM2'}
                    ensemble_name='r1i1p1f1';
                case {'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
                    ensemble_name='r1i1p1f2';
            end
            
            testname = [CMIP6_name,'_nwp_1_20'];
            output_dir=['/data1/RCM/CMIP6/input_hdd/nwp_1_20/input/SBC/',scenario_name,'/',CMIP6_name,'/'];
            if (exist(output_dir,'dir')~=7)
                mkdir(output_dir);
            end
            status = ROMS_surface_forcing_CMIP6(grdfile, year(i), varname{j}, tinterval,output_dir,testname, CMIP6_name, scenario_name, ensemble_name);
%             status = ROMS_surface_forcing_CMIP5_time_interp(grdfile, year(i), varname{j}, tinterval,output_dir,testname, CMIP6_name, scenario_name);
        end
    end
%     status = ROMS_surface_forcing_Qair(grdfile, year(i), tinterval,output_dir,testname);
end
