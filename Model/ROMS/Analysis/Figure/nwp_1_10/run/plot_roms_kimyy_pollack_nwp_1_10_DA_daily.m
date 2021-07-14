% % This code based on MATLAB R2016b.
% % Updated 15-Jun-2018 by Yong-Yub Kim

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
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run']));
%     addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run']));
end


% fid=fopen('E:\Data\Model\ROMS\nwp_1_10\test42\spinup\1989\files\modelinfo');
% modelinfo=textscan(fid,'%s');
% fclose(fid);

% % % for snu_desktop_auto
% % testname= modelinfo{1,1}{1,1};
% % year = str2num(modelinfo{1,1}{3,1}); 

testname= 'test06'
% year = [8791]; % % put year which you want to plot [year year ...]
year = [1983:2019];

month = [1:1]; % % put month which you want to plot [month month ...]
day= [1:1];
inputdir = strcat(['E:\Data\Model\ROMS\nwp_1_10\input\',testname,'\']);
% figdir =strcat([outputdir,'\figures\']); % % where figure files will be saved

figdir =strcat(['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\nwp_1_10\',testname,'\DA\']); % % where figure files will be saved


if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 


% % % % % % % % % % plot topography
% % % % 
% % % % % etopodir ='D:\MEPL\project\SSH\2nd_year\figure\etopo1\';
% % % % % plot_roms_kimyy_etopo
% % % % 
% % % % if (exist(strcat(figdir,'\Bathy\bathy_nwp_1_10.jpg'),'file') ~= 2)
% % % %     if (exist(strcat(figdir,'Bathy') , 'dir') ~= 7)
% % % %         mkdir(strcat(figdir, 'Bathy'));
% % % %     end 
% % % %     plot_roms_kimyy_topo_nwp_1_10;
% % % % end
% % % % 
% % % % 
% % % % %  plot EJS
lonlat_EJS = [127 144 33 52]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   

    tempyear=year(i);
    outputdir = strcat(['E:\Data\Model\ROMS\nwp_1_10\',testname,'\DA\',num2str(tempyear,'%04i'),'\']); % % where data files are

%     if (exist(strcat(figdir,'EJS'),'dir')~=7)
%         mkdir(strcat(figdir,'EJS'));
%     end
    
    % %  plot SST (EJS)
%     if (exist(strcat(figdir,'EJS\','SST') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS\','SST'));
%     end
%     sstfile =[figdir, 'EJS\', 'SST\', 'EJS_SST']; % % ~/SST_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_data(testname, sstfile, outputdir, [115 145 24 50], tempyear, month, ...
%                                     [0 35], 0:10:35, 'SST', 'fig_param_kyy_pollack');    
%     sstfile =[figdir, 'EJS\', 'SST\', 'EJS_SST_diff']; % % ~/SST_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_diff(testname, sstfile, outputdir, [115 145 24 50 0 0], tempyear, tempyear+1, month, ...
%                                     [-5 5], -100:100:100, 'SST', 'fig_param_kyy_pollack');    
%                                 
%     plot TUV (pollack)                            
%     if (exist(strcat(figdir,'EJS\','TUV\','TUV_0m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS\','TUV\','TUV_0m'));
%     end
%     suvfile =[figdir, 'EJS\','TUV\TUV_0m\', 'pollack_0m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_UV_data(testname, suvfile, outputdir, [125 135 34 42 -0 -0], tempyear, month, ...
%                                     [0 30], 0:10:30, 'SST', 'fig_param_kyy_pollack');
                                
% %     plot TUV_diff (pollack)                            
%     if (exist(strcat(figdir,'EJS\','TUV\','TUV_0m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS\','TUV\','TUV_0m'));
%     end
%     suvfile =[figdir, 'EJS\','TUV\TUV_0m\', 'pollack_0m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_UV_diff(testname, suvfile, outputdir, [125 135 34 42 -0 -0], tempyear, month, ...
%                                     [-3 3], -3:1:3, 'diff', 'fig_param_kyy_pollack');           
    
                                
%     suvfile =[figdir, 'EJS\','TUV\TUV_0m\', 'pollack_0m_TUV_con2']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_UV_data_2y_con(testname, suvfile, outputdir, [125 135 34 42 -0 -0], tempyear,tempyear-202,month, ...
%                                     [0 30], 0:10:30, 'SST', 'fig_param_kyy_pollack');

%     if (exist(strcat(figdir,'EJS\','EKB') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS\','EKB'));
%     end
%     lonlat_EKB = [127 130 37 42]; % [lon_start lon_end lat_start lat_end]
%     hor_temp_file =[figdir, 'EJS\','EKB\', 'pollack']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_horizontal_data_pollack(testname, hor_temp_file, outputdir, [127 130 38 41 -50 -50], tempyear, month, ...
%                                     [2 5], 2:1:5, 'temp', 'fig_param_kyy_pollack');
%     
%     if (exist(strcat(figdir,'EJS','/hor_temp'),'dir')~=7)
%         mkdir(strcat(figdir,'EJS','/hor_temp'));
%     end
%     hor_tuv_file =[figdir, 'EJS/', 'hor_temp/', 'EJS_hor_temp_0000m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, hor_tuv_file, outputdir, [127 144 33 52 0 0], tempyear,month, [-2 33], [-20 100], 'temp',  'fig_param_kyy_EJS_10');                                 
%     
%     hor_tuv_file =[figdir, 'EJS/', 'hor_temp/', 'EJS_hor_temp_0050m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, hor_tuv_file, outputdir, [127 144 33 52 -50 -50], tempyear,month, [-2 33], [-20  100], 'temp',  'fig_param_kyy_EJS_10');        
%     
%     hor_tuv_file =[figdir, 'EJS/', 'hor_temp/', 'EJS_hor_temp_0100m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, hor_tuv_file, outputdir, [127 144 33 52 -100 -100], tempyear,month, [-2 33], [-20  100], 'temp',  'fig_param_kyy_EJS_10');   
%     
%     hor_tuv_file =[figdir, 'EJS/', 'hor_temp/', 'EJS_hor_temp_0200m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, hor_tuv_file, outputdir, [127 144 33 52 -200 -200], tempyear,month, [-2 33], [-20 100], 'temp',  'fig_param_kyy_EJS_10');   
end




% % % % % % %  plot EKB
lonlat_EKB = [127 134 35 43]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   
    tempyear=year(i);
    outputdir = strcat(['E:\Data\Model\ROMS\nwp_1_10\',testname,'\DA\daily\',num2str(tempyear,'%04i'),'\']); % % where data files are
%     if (tempyear==2006)
%         tempyear2=tempyear+3;
%     else
        tempyear2=tempyear+4;
%     end
    if (exist(strcat(figdir,'EJS','/EK'),'dir')~=7)
        mkdir(strcat(figdir,'EJS','/EK'));
    end
%                        
% 
    if (exist(strcat(figdir,'EJS','/EK/hor_temp_daily'),'dir')~=7)
        mkdir(strcat(figdir,'EJS','/EK/hor_temp_daily'));
    end
   
%     hor_tuv_file =[figdir, 'EJS/EK/', 'hor_temp_daily/', 'EK_hor_uv_0000m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_uv(testname, hor_tuv_file, outputdir, [127 134 35 43 0 0], tempyear,month, [0 15], [-20 2 5 100], 'temp',  'fig_param_kyy_EKB_10');                                 

    hor_tuv_file =[figdir, 'EJS/EK/', 'hor_temp_daily/', 'EK_hor_temp_0000m']; % % ~/SUV_EJS_testname_year_month.jpg
    status=plot_ROMS_monthly_hor_t_uv_daily(testname, hor_tuv_file, outputdir, [127 134 35 43 0 0], tempyear,month, day, [0 15], [-20 2 5 10 100], 'temp',  'fig_param_kyy_EKB_10');                                 
  
    hor_tuv_file =[figdir, 'EJS/EK/', 'hor_temp_daily/', 'EK_hor_temp_0050m']; % % ~/SUV_EJS_testname_year_month.jpg
    status=plot_ROMS_monthly_hor_t_uv_daily(testname, hor_tuv_file, outputdir, [127 134 35 43 -50 -50], tempyear,month, day, [0 15], [-20 2 5 10 100], 'temp',  'fig_param_kyy_EKB_10');                                 
%  
%     hor_tuv_file =[figdir, 'EJS/EK/', 'hor_temp/', 'EK_hor_temp_diff_0000m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv_con_diff(testname, hor_tuv_file, outputdir, [127 134 35 43 0 0], tempyear, 8386, month, [0 15], [-20 2 5 100], 'temp',  'fig_param_kyy_EKB_10');                                 
%     
%     hor_tuv_file =[figdir, 'EJS/EK/', 'hor_temp/', 'EK_hor_temp_diff_0050m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv_con_diff(testname, hor_tuv_file, outputdir, [127 134 35 43 -50 -50], tempyear, 8386, month, [0 15], [-20 2 5 100], 'temp',  'fig_param_kyy_EKB_10');                                 
% 
%     hor_tuv_file =[figdir, 'EJS/EK/', 'hor_temp/', 'EK_hor_temp_diff2_0000m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv_con_diff2(testname, hor_tuv_file, outputdir, [127 134 35 43 0 0], tempyear, 8386, 0711, month, [0 15], [-20 2 5 100], 'temp',  'fig_param_kyy_EKB_10');        
%     
%     hor_tuv_file =[figdir, 'EJS/EK/', 'hor_temp/', 'EK_hor_temp_diff3_1619_shade_0000m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv_con_diff3(testname, hor_tuv_file, outputdir, [127 134 35 43 0 0], tempyear, 8386, 1619, month, [0 15], [-20 2 5 100], 'temp',  'fig_param_kyy_EKB_10_2');        
end
