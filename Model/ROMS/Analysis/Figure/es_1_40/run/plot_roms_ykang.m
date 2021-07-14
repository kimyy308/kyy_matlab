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
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\es_1_40\run']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/es_1_40/run']));
end


% fid=fopen('E:\Data\Model\ROMS\es_1_40\test42\spinup\1989\files\modelinfo');
% modelinfo=textscan(fid,'%s');
% fclose(fid);

% for snu_desktop_auto
% testname= modelinfo{1,1}{1,1};
% year = str2num(modelinfo{1,1}{3,1}); 


testname= 'test01';
year = [2012]; % % put year which you want to plot [year year ...]


month = [1]; % % put month which you want to plot [month month ...]
inputdir = strcat(['E:\Data\Model\ROMS\es_1_40\input\',testname,'\']);
outputdir = strcat(['E:\Data\Model\ROMS\es_1_40\',testname,'\run\',num2str(year,'%04i'),'\']); % % where data files are
figdir =strcat([outputdir,'\figures\']); % % where figure files will be saved

if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 


% % % % % % plot topography

% etopodir ='D:\MEPL\project\SSH\2nd_year\figure\etopo1\';
% plot_roms_kimyy_etopo

% if (exist(strcat(figdir,'\Bathy\bathy_es_1_40.jpg'),'file') ~= 2)
%     if (exist(strcat(figdir,'Bathy') , 'dir') ~= 7)
%         mkdir(strcat(figdir, 'Bathy'));
%     end 
%     plot_roms_kimyy_topo;
% end


% % % % %  plot EJS
lonlat_EJS = [127 144 33 52]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   
    tempyear=year(i);
    if (exist(strcat(figdir,'EJS'),'dir')~=7)
        mkdir(strcat(figdir,'EJS'));
    end
    
% %     % %  plot SST (EJS)
% %     if (exist(strcat(figdir,'EJS\','SST') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','SST'));
% %     end
% %     sstfile =[figdir, 'EJS\', 'SST\', 'EJS_SST']; % % ~/SST_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_surface_data(testname, sstfile, outputdir, lonlat_EJS, tempyear, month, ...
% %                                     [0 35], 0:5:35, 'SST', 'fig_param_kyy_EJS');                                
% % 
% %     % %  plot SSS (EJS)
% %     if (exist(strcat(figdir,'EJS\','SSS') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','SSS'));
% %     end
% %     sssfile =[figdir, 'EJS\', 'SSS\', 'EJS_SSS']; % % ~/SSS_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_surface_data(testname, sssfile, outputdir, lonlat_EJS, tempyear, month, ...
% %                                     [30 35], 0:0:0, 'SSS', 'fig_param_kyy_EJS');
% %                                 
% %     % %  plot SSH (EJS)
% %     if (exist(strcat(figdir,'EJS\','SSH') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','SSH'));
% %     end
% %     sshfile =[figdir, 'EJS\', 'SSH\', 'EJS_SSH']; % % ~/SSH_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_surface_data(testname, sshfile, outputdir, lonlat_EJS, tempyear, month, ...
% %                                     [-0.5 0.5], -0.5:0.1:0.5, 'SSH', 'fig_param_kyy_EJS');
% % 
% %     % %  plot Surface Current (SUV)
% %     if (exist(strcat(figdir,'EJS\','SUV') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','SUV'));
% %     end
% %     suvfile =[figdir, 'EJS\', 'SUV\', 'EJS_SUV']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=vec_ROMS_monthly_surface_UV(testname, suvfile, outputdir, lonlat_EJS, tempyear, month, ...
% %                                     'UV', 'EJS', 'fig_param_kyy_EJS');

    lonlat_EKWC = [129 132 34 39]; % [lon_start lon_end lat_start lat_end]
    if (exist(strcat(figdir,'EJS\','EKWC') , 'dir') ~= 7)
        mkdir(strcat(figdir,'EJS\','EKWC'));
    end
    if (exist(strcat(figdir,'EJS\','EKWC\','SUV') , 'dir') ~= 7)
        mkdir(strcat(figdir,'EJS\','EKWC\','SUV'));
    end
% %     suvfile =[figdir, 'EJS\', 'EKWC\','SUV\', 'EJS_SUV_EKWC']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=vec_ROMS_monthly_surface_UV(testname, suvfile, outputdir, lonlat_EKWC, tempyear, month, ...
% %                                     'UV', 'EJS', 'fig_param_kyy_EKWC');
% %     if (exist(strcat(figdir,'EJS\','EKWC\','TUV_100m') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','EKWC\','TUV_100m'));
% %     end
% %     suvfile =[figdir, 'EJS\', 'EKWC\','TUV_100m\', 'EKWC_100m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_horizontal_t_uv(testname, suvfile, outputdir, [129 132 34 39 -100 -100], tempyear, month, ...
% %                                     [0 30], 0:5:30, 'temp', 'fig_param_kyy_EKWC');
% %     if (exist(strcat(figdir,'EJS\','EKWC\','TUV_200m') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','EKWC\','TUV_200m'));
% %     end
% %     suvfile =[figdir, 'EJS\', 'EKWC\','TUV_200m\', 'EKWC_200m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_horizontal_t_uv(testname, suvfile, outputdir, [129 132 34 39 -200 -200], tempyear, month, ...
% %                                     [0 30], 0:5:30, 'temp', 'fig_param_kyy_EKWC');
% %     if (exist(strcat(figdir,'EJS\','EKWC\','TUV_350m') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','EKWC\','TUV_350m'));
% %     end
% %     suvfile =[figdir, 'EJS\', 'EKWC\','TUV_350m\', 'EKWC_350m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_horizontal_t_uv(testname, suvfile, outputdir, [129 132 34 39 -350 -350], tempyear, month, ...
% %                                     [0 30], 0:5:30, 'temp', 'fig_param_kyy_EKWC');
    if (exist(strcat(figdir,'EJS\','EKWC\','TUV_500m') , 'dir') ~= 7)
        mkdir(strcat(figdir,'EJS\','EKWC\','TUV_500m'));
    end
    suvfile =[figdir, 'EJS\', 'EJS_500m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
    status=plot_ROMS_monthly_horizontal_t_uv_ykang(testname, suvfile, outputdir, [127 144 34 52 -500 -500], tempyear, month, ...
                                    [0 30], 0:5:30, 'temp', 'fig_param_kyy_EJS_deep');

%     suvfile =[figdir, 'EJS\', 'SUV\', 'EJS_SUV_temp_102_EKWC']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_UV_data(testname, suvfile, outputdir, [129.8 131 36.0767 39 -200 0], tempyear, month, ...
%                                    [6 14], 0:10:10, 'UV', 'EJS', 'fig_param_kyy_schematic');
% 
% %     if (exist(strcat(figdir,'EJS\','EKWC\','VORTICITY') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','EKWC\','VORTICITY'));
% %     end
% %     suvfile =[figdir, 'EJS\', 'EKWC\','VORTICITY\', 'EJS_SUV_relative_vorticity_EKWC']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_surface_vorticity(testname, suvfile, outputdir, lonlat_EKWC, tempyear, month, ...
% %                                    [-0.0001 0.0001], 0:10:10, 'vorticity', 'EJS', 'fig_param_kyy_EKWC');


% %      % %  plot vertical temperature
% %     if (exist(strcat(figdir,'EJS\','VERTICAL') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','VERTICAL'));
% %     end
% %     if (exist(strcat(figdir,'EJS\','VERTICAL\','LAT_FIXED') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','VERTICAL\','LAT_FIXED'));
% %     end
% %     if (exist(strcat(figdir,'EJS\','VERTICAL\','LAT_FIXED\TEMP') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','VERTICAL\','LAT_FIXED\TEMP'));
% %     end
% %     vert_temp_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\TEMP\', 'temp_1000m_lon_129_136_lat_38']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 136 38 38 -1000 0], tempyear, month, ...
% %                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
% %     vert_temp_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\TEMP\', 'temp_1000m_lon_129_133_lat_37']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 133 37 37 -1000 0], tempyear, month, ...
% %                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
% %     vert_temp_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\TEMP\', 'temp_1000m_lon_129_133_lat_36']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 133 36 36 -1000 0], tempyear, month, ...
% %                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
% %     vert_temp_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\TEMP\', 'temp_1000m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 133 35 35 -1000 0], tempyear, month, ...
% %                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
% %     vert_temp_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\TEMP\', 'temp_200m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 133 35 35 -200 0], tempyear, month, ...
% %                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
                                
% %     if (exist(strcat(figdir,'EJS\','VERTICAL\','LON_FIXED') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','VERTICAL\','LON_FIXED'));
% %     end
% %     if (exist(strcat(figdir,'EJS\','VERTICAL\','LON_FIXED','\TEMP') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','VERTICAL\','LON_FIXED','\TEMP'));
% %     end
% %     vert_temp_file =[figdir, 'EJS\', 'VERTICAL\LON_FIXED\TEMP\', 'temp_1000m_lon_129_8_lat_35_43']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129.8 129.8 35 41 -1000 0], tempyear, month, ...
% %                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
% %                                 
% %     vert_temp_file =[figdir, 'EJS\', 'VERTICAL\LON_FIXED\TEMP\', 'temp_1000m_lon_132_lat_35_43']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [132 132 35 43 -1000 0], tempyear, month, ...
% %                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
                                
    
                     

% %     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_102']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129.8 131 36.0767 36.0767 -200 0], tempyear, month, ...
% %                                     [0 25], 0:5:25, 'vert_temp', 'fig_param_kyy_EJS');  %% for 102 line
% 
% %     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_bwr_102']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129.8 131 36.0767 36.0767 -200 0], tempyear, month, ...
% %                                     [6 14], 0:10:10, 'vert_temp', 'fig_param_kyy_EJS');  %% for 102 line_bwr_map
% %         vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_bwr_36_39']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129.8 129.8 36.0767 39 -200 0], tempyear, month, ...
% %                                     [6 14], 0:10:10, 'vert_temp', 'fig_param_kyy_EJS');  %% for 102, NS_bwr_map
%                                 
% %     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
% %                                     [0 25], 0:5:25, 'vert_temp', 'fig_param_kyy_EJS');  %% 35.5N
%                                 
% %     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_WOA']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [130 130 35 39 -200 0], tempyear, month, ...
% %                                     [2 20], 0:5:20, 'vert_temp', 'fig_param_kyy_EJS');   %% for WOA


% %      % %  plot vertical salinity
% %     if (exist(strcat(figdir,'EJS\','VERTICAL\','LAT_FIXED','\SALT') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','VERTICAL\','LAT_FIXED','\SALT'));
% %     end
% %     vert_salt_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\SALT\', 'salt_1000m_lon_129_136_lat_38']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 136 38 38 -1000 0], tempyear, month, ...
% %                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
% %     vert_salt_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\SALT\', 'salt_1000m_lon_129_133_lat_37']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 133 37 37 -1000 0], tempyear, month, ...
% %                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
% %     vert_salt_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\SALT\', 'salt_1000m_lon_129_133_lat_36']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 133 36 36 -1000 0], tempyear, month, ...
% %                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
% %     vert_salt_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\SALT\', 'salt_1000m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 133 35 35 -1000 0], tempyear, month, ...
% %                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
% %     vert_salt_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\SALT\', 'salt_200m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 133 35 35 -200 0], tempyear, month, ...
% %                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
% %     
% %     
% %     if (exist(strcat(figdir,'EJS\','VERTICAL\','LON_FIXED','\SALT') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','VERTICAL\','LON_FIXED','\SALT'));
% %     end
% %     
% %     vert_salt_file =[figdir, 'EJS\', 'VERTICAL\LON_FIXED\SALT\', 'salt_1000m_lon_129_8_lat_35_43']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129.8 129.8 35 41 -1000 0], tempyear, month, ...
% %                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
% %                                 
% %     vert_salt_file =[figdir, 'EJS\', 'VERTICAL\LON_FIXED\SALT\', 'salt_1000m_lon_132_lat_35_43']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [132 132 35 43 -1000 0], tempyear, month, ...
% %                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
                                
                                
                                
% %      % %  plot vertical v
% %     if (exist(strcat(figdir,'EJS\','VERTICAL\','LAT_FIXED','\V') , 'dir') ~= 7)
% %         mkdir(strcat(figdir,'EJS\','VERTICAL\','LAT_FIXED','\V'));
% %     end
% %     vert_v_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\V\', 'v_1000m_lon_129_136_lat_38']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 136 38 38 -1000 0], tempyear, month, ...
% %                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
% %     vert_v_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\V\', 'v_1000m_lon_129_133_lat_37']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 133 37 37 -1000 0], tempyear, month, ...
% %                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
% %     vert_v_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\V\', 'v_1000m_lon_129_133_lat_36']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 133 36 36 -1000 0], tempyear, month, ...
% %                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
% %     vert_v_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\V\', 'v_1000m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 133 35 35 -1000 0], tempyear, month, ...
% %                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
% %     vert_v_file =[figdir, 'EJS\', 'VERTICAL\LAT_FIXED\V\', 'v_200m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 133 35 35 -200 0], tempyear, month, ...
% %                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%     vert_v_file =[figdir, 'EJS\', 'vert_v\', 'EJS_vert_v_102']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129.8 131 36.0767 36.0767 -200 0], tempyear, month, ...
%                                     [-1.0 1.0], -1.0:0.1:1.0, 'vert_v', 'fig_param_kyy_EJS');  %% for 102 line
                                
%     vert_v_file =[figdir, 'EJS\', 'vert_v\', 'EJS_vert_v_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
%                                     [-1.0 1.0], -1.0:0.1:1.0, 'vert_v', 'fig_param_kyy_EJS');  %% for 35.5N
                   
end


% % % % % % plot vertical section
% % % % 
% % % % NIFS KHOA, South Sea & East Sea
% % nifs_num=[207, 208, 209, 102, 103, 104, 105, 106, 107];
% % nifs_area=[ 129.12   130.6717  35.0217  34.15;    ...
% %             129.455  130.83    35.475   34.5417;  ...
% %             129.5317 131.6517  35.795   34.9817;  ...
% %             129.5883 132.455   36.0767  36.0767;  ...
% %             129.4833 132.475   36.505   36.505;   ...
% %             129.46   134.34    37.0567  37.0567;  ...
% %             129.15   132.8067  37.5533  37.5533;  ...
% %             128.875  132.8267  37.895   37.895;   ...
% %             128.6233 132.313   38.126   38.126;   ];
% %                     
% % for i=1:length(year)   
% %     tempyear=year(i);
% %     for var=1:2
% %         if (var==1)
% %             varname='vert_temp'
% %             var_lim = [0 30];
% %             var_con_lim= [0 2 5 10 15 20 25 30];
% %         elseif (var==2)
% %             varname='vert_salt'
% %             var_lim = [33 34.5];
% %             var_con_lim = [33:0.1:34.5];
% %         end
% % 
% %         for i=2:5
% %             maxdepth = -(i*100)
% %             for nifs_ind=1:length(nifs_num)
% %                 if (exist(strcat(figdir,'NIFS/',num2str(nifs_num(nifs_ind))) , 'dir') ~= 7)
% %                     mkdir(strcat(figdir,'NIFS/',num2str(nifs_num(nifs_ind))));
% %                 end
% %                 nifs_file =[figdir,'NIFS/',num2str(nifs_num(nifs_ind)),'/', ...
% %                             varname,'_',num2str(abs(maxdepth),'%04i'),'m_',num2str(nifs_num(nifs_ind)),'_line']; 
% %                 nifs_section(1:4)=nifs_area(nifs_ind,1:4);
% %                 nifs_section(5)=maxdepth;  nifs_section(6) =0;
% %                 status=plot_ROMS_monthly_vertical_data(testname, nifs_file, outputdir, nifs_section, tempyear, month, ...
% %                                                 var_lim, var_con_lim, varname, 'fig_param_kyy_EJS');    
% %             end
% %         end 
% %     end
% % end
