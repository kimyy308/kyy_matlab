% % This code based on MATLAB R2016b.
% % Updated 27-Apr-2018 by Yong-Yub Kim

clc;close all;clear all;
warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
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


fid=fopen('E:\Data\Model\ROMS\nwp_1_20\test42\spinup\1989\files\modelinfo');
modelinfo=textscan(fid,'%s');
fclose(fid);

% for snu_desktop_auto
testname= modelinfo{1,1}{1,1};
year = str2num(modelinfo{1,1}{3,1}); 


% testname= 'test42';
% year = [1985 1986 1987 1988]; % % put year which you want to plot [year year ...]


month = [1]; % % put month which you want to plot [month month ...]
inputdir = strcat(['E:\Data\Model\ROMS\nwp_1_20\input\',testname,'\']);
outputdir = strcat(['E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',num2str(year,'%04i'),'\']); % % where data files are
figdir =strcat([outputdir,'\figures\']); % % where figure files will be saved

if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 


% % % % % % plot topography

% etopodir ='D:\MEPL\project\SSH\2nd_year\figure\etopo1\';
% plot_roms_kimyy_etopo

if (exist(strcat(figdir,'\Bathy\bathy_nwp_1_20.jpg'),'file') ~= 2)
    if (exist(strcat(figdir,'Bathy') , 'dir') ~= 7)
        mkdir(strcat(figdir, 'Bathy'));
    end 
    plot_roms_kimyy_topo;
end


% % % % % %  plot NWP

lonlat_NWP = [115 164 15 52]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   
    tempyear=year(i);
    if (exist(strcat(figdir,'NWP'),'dir')~=7)
        mkdir(strcat(figdir,'NWP'));
    end

    if (exist(strcat(figdir,'NWP\','TUV_200m') , 'dir') ~= 7)
        mkdir(strcat(figdir,'NWP\','TUV_200m'));
    end
    suvfile =[figdir, 'NWP\','TUV_200m\', 'NWP_200m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
    status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [115 164 15 52 40 40], tempyear, month, ...
                                    [0 30], 0:5:30, 'temp', 'fig_param_kyy_NWP');
                          
end

