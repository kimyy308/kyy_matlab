function [error_status, dropboxpath] = Func_0008_set_dropbox_path(system_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function status=set_dropbox_path(system_name);
%
% get the polygon data from nwp_polygon_point.m corresponding to regionname
%
%  input:
%  system_name             operation system name (string, result of computer)
%
%  output:
%  refpolygon             Reference polygon point (2-D array, [lon; lat])
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    20-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp.fs = filesep;

if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\User\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
%     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
    addpath(genpath([dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
    error_status=1;
elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/home/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
    error_status=2;
end


end

