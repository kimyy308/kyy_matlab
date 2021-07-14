% % This code based on MATLAB R2016b.

clc;close all;clear all;
warning off;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:/Users/KYY/Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param']));
    addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
elseif (strcmp(system_name,'GLNXA64'))    % % for linux
%     dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/seo_nwp_1_10/run']));
end




% for snu_desktop
testname='avhrr_only_monthly_v2_' 
% year = [2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012]; % % put year which you want to plot [year year ...]
% year = [1982:2016]; % % put year which you want to plot [year year ...]
year = [2018:2019]; % % put year which you want to plot [year year ...]

% month = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
month = [1 2 3]; % % put month which you want to plot [month month ...]
% figdir =strcat('D:/OneDrive - ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿????/MEPL/project/MICT_pollack/3rd year/figure/avg_ens_10km_mean/'); % % where figure files will be saved
% inputdir = strcat('E:/Data/Reanalysis/nwp_1_10_seo/', testname, '/'); % % where data files are

figdir =strcat('D:\OneDrive - ¼­¿ï´ëÇÐ±³\MEPL\project\MICT_pollack\3rd_year\figure\AVHRR\EJS\EK\hor_temp\'); % % where figure files will be saved
inputdir = strcat('E:\Data\Observation\OISST\monthly_kimyy/'); % % where data files are
% figdir =strcat('/data2/kimyy/Observation/OISST/monthly_kimyy/figures/'); % % where figure files will be saved
% inputdir = strcat('/data2/kimyy/Observation/OISST/monthly_kimyy/'); % % where data files are



if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 

% % % % % % %  plot EK
lonlat_EK = [127 134 35 43]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   
    tempyear=year(i);
%     if (exist(strcat(figdir,'EJS','/EK'),'dir')~=7)
%         mkdir(strcat(figdir,'EJS','/EK'));
%     end
%     if (exist(strcat(figdir,'EJS','/EK/hor_temp'),'dir')~=7)
%         mkdir(strcat(figdir,'EJS','/EK/hor_temp'));
%     end
%     hor_tuv_file =[figdir, 'EJS/EK/', 'hor_temp/', 'EK_hor_temp_0000m']; % % ~/SUV_EJS_testname_year_month.jpg
    hor_tuv_file =[figdir, 'EK_hor_temp_0000m']; % % ~/SUV_EJS_testname_year_month.jpg

    status=plot_OISST_monthly_hor_t_2(testname, hor_tuv_file, inputdir, lonlat_EK, tempyear, month, [0 15], [-20 2 5 100], 'fig_param_kyy_EKB_10_2');                                 
    
%     status=plot_OISST_monthly_hor_t(testname, hor_tuv_file, inputdir, [115 164 15 52], tempyear, month, [0 30], [-20 0 10 20 30 100], 'fig_param_kyy_EKB_10');                                 

%     status=plot_OISST_monthly_hor_t(testname, hor_tuv_file, inputdir, lonlat_EK, tempyear, month, [0 15], [-20 0 2 5 10 20 30 100], 'fig_param_kyy_EKB_10');                                 
% 
end
