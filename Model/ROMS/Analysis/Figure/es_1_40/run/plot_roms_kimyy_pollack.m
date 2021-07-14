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

% % % for snu_desktop_auto
% % testname= modelinfo{1,1}{1,1};
% % year = str2num(modelinfo{1,1}{3,1}); 


testname= 'test42';
year = [2006]; % % put year which you want to plot [year year ...]

month = [1,2,12,16]; % % put month which you want to plot [month month ...]
inputdir = strcat(['E:\Data\Model\ROMS\nwp_1_20\input\',testname,'\']);
% figdir =strcat([outputdir,'\figures\']); % % where figure files will be saved

figdir =strcat(['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd year\figure\nwp_1_20\']); % % where figure files will be saved


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


% % % % %  plot EJS
lonlat_EJS = [127 144 33 52]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   

    tempyear=year(i);
    outputdir = strcat(['E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',num2str(tempyear,'%04i'),'\']); % % where data files are

    if (exist(strcat(figdir,'EJS'),'dir')~=7)
        mkdir(strcat(figdir,'EJS'));
    end
    
    % %  plot SST (EJS)
    if (exist(strcat(figdir,'EJS\','SST') , 'dir') ~= 7)
        mkdir(strcat(figdir,'EJS\','SST'));
    end
    sstfile =[figdir, 'EJS\', 'SST\', 'EJS_SST']; % % ~/SST_EJS_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_data(testname, sstfile, outputdir, lonlat_EJS, tempyear, month, ...
                                    [0 30], 0:10:30, 'SST', 'fig_param_kyy_pollack');    
                                
%     plot TUV (pollack)                            
    if (exist(strcat(figdir,'EJS\','TUV\','TUV_0m') , 'dir') ~= 7)
        mkdir(strcat(figdir,'EJS\','TUV\','TUV_0m'));
    end
    suvfile =[figdir, 'EJS\','TUV\TUV_0m\', 'pollack_0m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_UV_data(testname, suvfile, outputdir, [125 135 34 42 -0 -0], tempyear, month, ...
                                    [0 30], 0:10:30, 'SST', 'fig_param_kyy_pollack');
                                
% %     plot TUV_diff (pollack)                            
%     if (exist(strcat(figdir,'EJS\','TUV\','TUV_0m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS\','TUV\','TUV_0m'));
%     end
%     suvfile =[figdir, 'EJS\','TUV\TUV_0m\', 'pollack_0m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_UV_diff(testname, suvfile, outputdir, [125 135 34 42 -0 -0], tempyear, month, ...
%                                     [-3 3], -3:1:3, 'diff', 'fig_param_kyy_pollack');           
    
                                
%     suvfile =[figdir, 'EJS\','TUV\TUV_0m\', 'pollack_0m_TUV_con2']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_UV_data_2y_con(testname, suvfile, outputdir, [125 135 34 42 -0 -0], tempyear,tempyear-2,month, ...
%                                     [0 30], 0:10:30, 'SST', 'fig_param_kyy_pollack');

    if (exist(strcat(figdir,'EJS\','EKB') , 'dir') ~= 7)
        mkdir(strcat(figdir,'EJS\','EKB'));
    end
% lonlat_EKB = [127 130 37 42]; % [lon_start lon_end lat_start lat_end]
    hor_temp_file =[figdir, 'EJS\','EKB\', 'pollack']; % % ~/SUV_EJS_testname_year_month.jpg
    status=plot_ROMS_monthly_horizontal_data_pollack(testname, hor_temp_file, outputdir, [127 130 38 41 -50 -50], tempyear, month, ...
                                    [2 5], 2:1:5, 'temp', 'fig_param_kyy_pollack');
    
end





