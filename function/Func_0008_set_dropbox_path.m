function [dropboxpath, error_status] = Func_0008_set_dropbox_path(system_name)
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
%  Updated    07-Oct-2022 by Yong-Yub Kim  operated by hostname?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp.fs = filesep;

switch system_name
    case 'PCWIN64'
        tmp.hostname = getenv('COMPUTERNAME'); % for windows
        % tmp.hostname = getenv('HOSTNAME')
    case {'GLNXA64', 'MACI64'}
        [error_status, tmp.hostname] = system('hostname');
        tmp.hostname=tmp.hostname(1:end-1);
end


switch tmp.hostname
    case {'ROMS', 'DAMO', 'BONG'}
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/cptcmap']));
    %     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
        error_status=2;
    case 'Yong-Yubs-iMac-Pro.local'
        dropboxpath='/Volumes/kyy_raid/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/cptcmap']));
        addpath(genpath([dropboxpath '/source/matlab/Common/hatchfill2_r8']));
    %     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
        error_status=2;
    case 'Yong-Yubui-MacBookPro.local'
        dropboxpath = '/Users/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/cptcmap']));
        addpath(genpath([dropboxpath '/source/matlab/Common/hatchfill2_r8']));
    %     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));

    case {'da1', 'da2', 'da3', 'da4'}
        dropboxpath='/proj/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/cptcmap']));
        addpath(genpath([dropboxpath '/source/matlab/Common/hatchfill2_r8']));
    %     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
        error_status=2;
    otherwise % windows
        dropboxpath='C:\Users\User\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
        addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
    %     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
        addpath(genpath([dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
        addpath(genpath([dropboxpath '\source\matlab\Common\gsw_matlab_v3_06_14']));
        error_status=1;
end

%% old set
% if (strcmp(system_name,'PCWIN64'))
%     % % for windows
%     dropboxpath='C:\Users\User\Dropbox';
%     addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
%     addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
%     addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
% %     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
%     addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
%     addpath(genpath([dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
%     addpath(genpath([dropboxpath '\source\matlab\Common\gsw_matlab_v3_06_14']));
%     error_status=1;
% elseif (strcmp(system_name,'GLNXA64'))
%     dropboxpath='/home/kimyy/Dropbox';
%     addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/cptcmap']));
% %     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
%     error_status=2;
% elseif (strcmp(system_name,'MACI64'))
%     dropboxpath='/Volumes/kyy_raid/kimyy/Dropbox';
%     addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/cptcmap']));
% %     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
%     error_status=2;
% end


end

