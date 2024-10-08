% % % This script is based on MATLAB R2014a
% % Updated 04-Mar-2019 by Y.Y.Kim (mepl2)
% % Updated 25-May-2019 by Y.Y.Kim (using uas, vas in IPSL models)


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


grdfile='/home/kimyy/Model/ROMS/nwp_1_20/input/test51/spinup/roms_grid_nwp_1_20_test51.nc';

% rawdata_dir='/data1/temp/ECMWF_interim/';
varname={'tas', 'psl', 'hur', 'rsds', 'ua', 'va'};
% varname={'hur','ua','va'};
% varname={'tas', 'psl', 'rsds'};

% year=[1976];
year=[2087:2100];
scenario_name='rcp45'; % historical, rcp45, rcp26, rcp85
tinterval = 1;

% CMIP5_name='ens10';
CMIP5_names={'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
% CMIP5_names={'MPI-ESM-LR'};

% CMIP5_names={'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M'};

% CMIP5_name='CanESM2';
% CMIP5_name='NorESM-M';
% CMIP5_name='MPI-ESM-MR';

    
% ncread('/data1/temp/ECMWF_interim/ssrd/ECMWF_Interim_ssrd_2016.nc',ssrd);
for i=1:length(year)
    for j=1:length(varname)
        for k=1:length(CMIP5_names)
            CMIP5_name=CMIP5_names{k};
            if(strcmp(CMIP5_name,'IPSL-CM5A-LR')==1 || strcmp(CMIP5_name,'IPSL-CM5A-MR')==1)
                if(strcmp(varname{j},'ua')==1)
                    varname{j}='uas';
                elseif(strcmp(varname{j},'va')==1)
                    varname{j}='vas';
                elseif(strcmp(varname{j},'hur')==1)
                    varname{j}='huss';
                end
            else
                if(strcmp(varname{j},'uas')==1)
                    varname{j}='ua';
                elseif(strcmp(varname{j},'vas')==1)
                    varname{j}='va';
                elseif(strcmp(varname{j},'huss')==1)
                    varname{j}='hur';
                end
            end
            testname = [CMIP5_name,'_nwp_1_20'];
            output_dir=['/home/auto/ext_hdf/nwp_1_20','/forcing_matlab/SBC/',scenario_name,'/',CMIP5_name,'/'];
            if (exist(output_dir,'dir')~=7)
                mkdir(output_dir);
            end
            status = ROMS_surface_forcing_CMIP5_ens10(grdfile, year(i), varname{j}, tinterval,output_dir,testname, CMIP5_name, scenario_name);
%             status = ROMS_surface_forcing_CMIP5_time_interp(grdfile, year(i), varname{j}, tinterval,output_dir,testname, CMIP5_name, scenario_name);
        end
    end
%     status = ROMS_surface_forcing_Qair(grdfile, year(i), tinterval,output_dir,testname);
end