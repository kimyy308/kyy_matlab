% % This code based on MATLAB R2016b.
% % Updated 21-May-2018 by Yong-Yub Kim

clc;close all;clear all;
warning off;

linux=1; windows=0;

% % for DGIST 3rd
dropboxpath='/scratch/snu01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));

fid=fopen('/scratch/snu01/kimyy/roms_nwp/nwp_1_20/output/modelinfo_fig');
modelinfo=textscan(fid,'%s');
fclose(fid);

testname= modelinfo{1,1}{1,1};
year = str2num(modelinfo{1,1}{3,1}); 
month = [1:12]; % % put month which you want to plot [month month ...]
inputdir = strcat(['/scratch/snu01/kimyy/roms_nwp/nwp_1_20/input/',testname,'/run/',num2str(year,'%04i'),'/']);
outputdir = strcat(['/scratch/snu01/kimyy/roms_nwp/nwp_1_20/output/',testname,'/run/',num2str(year,'%04i'),'/']); % % where data files are
figdir =strcat([outputdir,'figures/']); % % where figure files will be saved

if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 


% % % % % % plot topography

% etopodir ='D:/MEPL/project/SSH/2nd_year/figure/etopo1/';
% plot_roms_kimyy_etopo

if (exist(strcat(figdir,'/Bathy/bathy_nwp_1_20.jpg'),'file') ~= 2)
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
    
    % %  plot SST (NWP)
    if (exist(strcat(figdir,'NWP/','SST') , 'dir') ~= 7)
        mkdir(strcat(figdir,'NWP/','SST'));
    end
    sstfile =[figdir, 'NWP/', 'SST/', 'NWP_SST']; % % ~/SST_NWP_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_data(testname, sstfile, outputdir, lonlat_NWP, tempyear, month, ...
                                    [0 35], 0:5:35, 'SST', 'fig_param_kyy_NWP');

    % %  plot SSS (NWP)
    if (exist(strcat(figdir,'NWP/','SSS') , 'dir') ~= 7)
        mkdir(strcat(figdir,'NWP/','SSS'));
    end
    sssfile =[figdir, 'NWP/', 'SSS/', 'NWP_SSS']; % % ~/SSS_NWP_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_data(testname, sssfile, outputdir, lonlat_NWP, tempyear, month, ...
                                    [30 35], 0:0:0, 'SSS', 'fig_param_kyy_NWP');
                                
%     % %  plot SSH (NWP)
%     if (exist(strcat(figdir,'NWP/','SSH') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','SSH'));
%     end
%     sshfile =[figdir, 'NWP/', 'SSH/', 'NWP_SSH']; % % ~/SSH_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_data(testname, sshfile, outputdir, lonlat_NWP, tempyear, month, ...
%                                     [-1.5 1.5], 0:0:0, 'SSH', 'fig_param_kyy_NWP');
% 
%      % % plot Surface Current (SUV)
%     if (exist(strcat(figdir,'NWP/','SUV') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','SUV'));
%     end
%     suvfile =[figdir, 'NWP/', 'SUV/', 'NWP_SUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=vec_ROMS_monthly_surface_UV(testname, suvfile, outputdir, lonlat_NWP, tempyear, month, ...
%                                     'UV', 'NWP', 'fig_param_kyy_NWP');


%     % %  plot TUV (NWP)
%     if (exist(strcat(figdir,'NWP/','TUV_100m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_100m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_100m/', 'NWP_100m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -100 -100], tempyear, month, ...
%                                     [0 30], [1 2 3 5 10 15 20 30], 'temp', 'fig_param_kyy_NWP');
%     
%     if (exist(strcat(figdir,'NWP/','TUV_200m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_200m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_200m/', 'NWP_200m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -200 -200], tempyear, month, ...
%                                     [0 20], [1 2 3 5 10 15 20], 'temp', 'fig_param_kyy_NWP');
%     
%     if (exist(strcat(figdir,'NWP/','TUV_300m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_300m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_300m/', 'NWP_300m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -300 -300], tempyear, month, ...
%                                     [0 15], [1 2 3 5 10 15], 'temp', 'fig_param_kyy_NWP');
%     
%     if (exist(strcat(figdir,'NWP/','TUV_400m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_400m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_400m/', 'NWP_400m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -400 -400], tempyear, month, ...
%                                     [0 10], [1 2 3 5 10], 'temp', 'fig_param_kyy_NWP');
%     
% 
%     if (exist(strcat(figdir,'NWP/','TUV_500m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_500m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_500m/', 'NWP_500m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -500 -500], tempyear, month, ...
%                                     [0 5], [1 2 3 4 5 ], 'temp', 'fig_param_kyy_NWP');
% 
%     if (exist(strcat(figdir,'NWP/','TUV_700m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_700m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_700m/', 'NWP_700m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -700 -700], tempyear, month, ...
%                                     [0 3], [0:0.3:3], 'temp', 'fig_param_kyy_NWP');
% 
% 
%     if (exist(strcat(figdir,'NWP/','TUV_1000m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_1000m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_1000m/', 'NWP_1000m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -1000 -1000], tempyear, month, ...
%                                     [0 3], [0:0.3:3], 'temp', 'fig_param_kyy_NWP');
%     
%     if (exist(strcat(figdir,'NWP/','TUV_1500m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_1500m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_1500m/', 'NWP_1500m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -1500 -1500], tempyear, month, ...
%                                     [0 3], [0:0.3:3], 'temp', 'fig_param_kyy_NWP');
%     
%     if (exist(strcat(figdir,'NWP/','TUV_2000m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_2000m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_2000m/', 'NWP_2000m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -2000 -2000], tempyear, month, ...
%                                     [0 2], 0:0.2:2, 'temp', 'fig_param_kyy_NWP');
%     
% 
%     if (exist(strcat(figdir,'NWP/','TUV_2500m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_2500m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_2500m/', 'NWP_2500m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -2500 -2500], tempyear, month, ...
%                                     [0 2], 0:0.2:2, 'temp', 'fig_param_kyy_NWP');
% 
%     if (exist(strcat(figdir,'NWP/','TUV_3000m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_3000m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_3000m/', 'NWP_3000m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -3000 -3000], tempyear, month, ...
%                                     [0 2], 0:0.2:2, 'temp', 'fig_param_kyy_NWP');
% 
%                                 
%     if (exist(strcat(figdir,'NWP/','TUV_4000m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'NWP/','TUV_4000m'));
%     end
%     suvfile =[figdir, 'NWP/', 'TUV_4000m/', 'NWP_4000m_TUV']; % % ~/SUV_NWP_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -4000 -4000], tempyear, month, ...
%                                     [0 2], 0:0.2:2, 'temp', 'fig_param_kyy_NWP');
                          
end


% % % % % % %  plot YS
% 
lonlat_YS = [116 128 32 42]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   
    tempyear=year(i);
    if (exist(strcat(figdir,'YS'),'dir')~=7)
        mkdir(strcat(figdir,'YS'));
    end
    
%     % %  plot SST (YS)
%     if (exist(strcat(figdir,'YS/','SST') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'YS/','SST'));
%     end
%     sstfile =[figdir, 'YS/', 'SST/', 'YS_SST']; % % ~/SST_YS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_data(testname, sstfile, outputdir, lonlat_YS, tempyear, month, ...
%                                     [0 35], 0:5:35, 'SST', 'fig_param_kyy_YS');
                                
    % %  plot YSBCW (YS)
    if (exist(strcat(figdir,'YS/','YSBCW') , 'dir') ~= 7)
        mkdir(strcat(figdir,'YS/','YSBCW'));
    end
    sstfile =[figdir, 'YS/', 'YSBCW/', 'YSBCW']; % % ~/SST_YS_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_data(testname, sstfile, outputdir, lonlat_YS, tempyear, month, ...
                                    [0 35], 6:2:10, 'YSBCW', 'fig_param_kyy_YS');

    % %  plot SSS (YS)
    if (exist(strcat(figdir,'YS/','SSS') , 'dir') ~= 7)
        mkdir(strcat(figdir,'YS/','SSS'));
    end
    sssfile =[figdir, 'YS/', 'SSS/', 'YS_SSS']; % % ~/SSS_YS_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_data(testname, sssfile, outputdir, lonlat_YS, tempyear, month, ...
                                    [30 35], 0:0:0, 'SSS', 'fig_param_kyy_YS');
                                
%     % %  plot SSH (YS)
%     if (exist(strcat(figdir,'YS/','SSH') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'YS/','SSH'));
%     end
%     sshfile =[figdir, 'YS/', 'SSH/', 'YS_SSH']; % % ~/SSH_YS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_data(testname, sshfile, outputdir, lonlat_YS, tempyear, month, ...
%                                     [-0.5 0.5], -0.5:0.1:0.5, 'SSH', 'fig_param_kyy_YS');

%     % %  plot Surface Current (SUV)
%     if (exist(strcat(figdir,'YS/','SUV') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'YS/','SUV'));
%     end
%     suvfile =[figdir, 'YS/', 'SUV/', 'YS_SUV']; % % ~/SUV_YS_testname_year_month.jpg
%     status=vec_ROMS_monthly_surface_UV(testname, suvfile, outputdir, lonlat_YS, tempyear, month, ...
%                                     'UV', 'YS', 'fig_param_kyy_YS');
%                           
%     if (exist(strcat(figdir,'YS/','TUV_30m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'YS/','TUV_30m'));
%     end
%     suvfile =[figdir, 'YS/','TUV_30m/', 'YS_30m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [116 128 32 42 -30 -30], tempyear, month, ...
%                                     [0 30], 0:5:30, 'temp', 'fig_param_kyy_YS');
%                                 
%     if (exist(strcat(figdir,'YS/','TUV_50m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'YS/','TUV_50m'));
%     end
%     suvfile =[figdir, 'YS/','TUV_50m/', 'YS_50m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [116 128 32 42 -50 -50], tempyear, month, ...
%                                     [0 30], 0:5:30, 'temp', 'fig_param_kyy_YS');
% 
%     if (exist(strcat(figdir,'YS/','TUV_70m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'YS/','TUV_70m'));
%     end
%     suvfile =[figdir, 'YS/','TUV_70m/', 'YS_70m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [116 128 32 42 -70 -70], tempyear, month, ...
%                                     [0 30], 0:5:30, 'temp', 'fig_param_kyy_YS');
    % %  plot vertical temperature
    if (exist(strcat(figdir,'YS/','VERTICAL') , 'dir') ~= 7)
        mkdir(strcat(figdir,'YS/','VERTICAL'));
    end
    if (exist(strcat(figdir,'YS/','VERTICAL/','LON_FIXED') , 'dir') ~= 7)
        mkdir(strcat(figdir,'YS/','VERTICAL/','LON_FIXED'));
    end
    if (exist(strcat(figdir,'YS/','VERTICAL/','LAT_FIXED') , 'dir') ~= 7)
        mkdir(strcat(figdir,'YS/','VERTICAL/','LAT_FIXED'));
    end
    if (exist(strcat(figdir,'YS/','VERTICAL/','LON_FIXED','/TEMP') , 'dir') ~= 7)
        mkdir(strcat(figdir,'YS/','VERTICAL/','LON_FIXED','/TEMP'));
    end
    vert_temp_file =[figdir, 'YS/', 'VERTICAL/LON_FIXED/TEMP/', 'temp_100m_lon_124_4_lat_33_5_36_5']; 
    status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [124.4 124.4 33.5 36.5 -100 0], tempyear, month, ...
                                    [5 20], [0 7 8 9 10 15 20], 'vert_temp', 'fig_param_kyy_YS');
    if (exist(strcat(figdir,'YS/','VERTICAL/','LON_FIXED','/SALT') , 'dir') ~= 7)
        mkdir(strcat(figdir,'YS/','VERTICAL/','LON_FIXED','/SALT'));
    end
    vert_salt_file =[figdir, 'YS/', 'VERTICAL/LON_FIXED/SALT/', 'salt_100m_lon_124_4_lat_33_5_36_5']; 
    status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [124.4 124.4 33.5 36.5 -100 0], tempyear, month, ...
                                    [31 35], [32:0.5:35], 'vert_temp', 'fig_param_kyy_YS');
                                
    if (exist(strcat(figdir,'YS/','VERTICAL/','LAT_FIXED','/V') , 'dir') ~= 7)
        mkdir(strcat(figdir,'YS/','VERTICAL/','LAT_FIXED','/V'));
    end
    vert_v_file =[figdir, 'YS/', 'VERTICAL/LAT_FIXED/V/', 'v_100m_lon_119_127_lat_34']; 
    status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [119 127 34 34 -100 0], tempyear, month, ...
                                    [-0.1 0.1], [-0.1:0.02:0.1], 'vert_v', 'fig_param_kyy_YS');
    vert_v_file =[figdir, 'YS/', 'VERTICAL/LAT_FIXED/V/', 'v_100m_lon_119_127_lat_35']; 
    status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [119 127 35 35 -100 0], tempyear, month, ...
                                    [-0.1 0.1], [-0.1:0.02:0.1], 'vert_v', 'fig_param_kyy_YS');
    vert_v_file =[figdir, 'YS/', 'VERTICAL/LAT_FIXED/V/', 'v_100m_lon_119_127_lat_36']; 
    status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [119 127 36 36 -100 0], tempyear, month, ...
                                    [-0.1 0.1], [-0.1:0.02:0.1], 'vert_v', 'fig_param_kyy_YS');                            
end
% 
% % % % % %  plot EJS
% lonlat_EJS = [127 144 33 52]; % [lon_start lon_end lat_start lat_end]
% for i=1:length(year)   
%     tempyear=year(i);
%     if (exist(strcat(figdir,'EJS'),'dir')~=7)
%         mkdir(strcat(figdir,'EJS'));
%     end
%     
%     % %  plot SST (EJS)
%     if (exist(strcat(figdir,'EJS/','SST') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','SST'));
%     end
%     sstfile =[figdir, 'EJS/', 'SST/', 'EJS_SST']; % % ~/SST_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_data(testname, sstfile, outputdir, lonlat_EJS, tempyear, month, ...
%                                     [0 35], 0:5:35, 'SST', 'fig_param_kyy_EJS');                                
% 
%     % %  plot SSS (EJS)
%     if (exist(strcat(figdir,'EJS/','SSS') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','SSS'));
%     end
%     sssfile =[figdir, 'EJS/', 'SSS/', 'EJS_SSS']; % % ~/SSS_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_data(testname, sssfile, outputdir, lonlat_EJS, tempyear, month, ...
%                                     [30 35], 0:0:0, 'SSS', 'fig_param_kyy_EJS');
%                                 
%     % %  plot SSH (EJS)
%     if (exist(strcat(figdir,'EJS/','SSH') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','SSH'));
%     end
%     sshfile =[figdir, 'EJS/', 'SSH/', 'EJS_SSH']; % % ~/SSH_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_data(testname, sshfile, outputdir, lonlat_EJS, tempyear, month, ...
%                                     [-0.5 0.5], -0.5:0.1:0.5, 'SSH', 'fig_param_kyy_EJS');
%     
%     % %  plot Surface Current (SUV)
%     if (exist(strcat(figdir,'EJS/','SUV') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','SUV'));
%     end
%     suvfile =[figdir, 'EJS/', 'SUV/', 'EJS_SUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=vec_ROMS_monthly_surface_UV(testname, suvfile, outputdir, lonlat_EJS, tempyear, month, ...
%                                     'UV', 'EJS', 'fig_param_kyy_EJS');
%    
%     % %  plot TUV (EJS)
%     if (exist(strcat(figdir,'EJS/','TUV_100m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_100m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_100m/', 'EJS_100m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -100 -100], tempyear, month, ...
%                                     [0 30], [1 2 3 5 10 15 20 30], 'temp', 'fig_param_kyy_EJS');
%     
%     if (exist(strcat(figdir,'EJS/','TUV_200m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_200m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_200m/', 'EJS_200m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -200 -200], tempyear, month, ...
%                                     [0 20], [1 2 3 5 10 15 20], 'temp', 'fig_param_kyy_EJS');
%     
%     if (exist(strcat(figdir,'EJS/','TUV_300m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_300m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_300m/', 'EJS_300m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -300 -300], tempyear, month, ...
%                                     [0 15], [1 2 3 5 10 15], 'temp', 'fig_param_kyy_EJS');
%     
%     if (exist(strcat(figdir,'EJS/','TUV_400m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_400m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_400m/', 'EJS_400m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -400 -400], tempyear, month, ...
%                                     [0 10], [1 2 3 5 10], 'temp', 'fig_param_kyy_EJS');
%     
% 
%     if (exist(strcat(figdir,'EJS/','TUV_500m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_500m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_500m/', 'EJS_500m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -500 -500], tempyear, month, ...
%                                     [0 5], [1 2 3 4 5 ], 'temp', 'fig_param_kyy_EJS');
% 
%     if (exist(strcat(figdir,'EJS/','TUV_700m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_700m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_700m/', 'EJS_700m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -700 -700], tempyear, month, ...
%                                     [0 3], [0:0.3:3], 'temp', 'fig_param_kyy_EJS');
% 
% 
%     if (exist(strcat(figdir,'EJS/','TUV_1000m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_1000m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_1000m/', 'EJS_1000m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -1000 -1000], tempyear, month, ...
%                                     [0 3], [0:0.3:3], 'temp', 'fig_param_kyy_EJS');
%     
%     if (exist(strcat(figdir,'EJS/','TUV_1500m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_1500m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_1500m/', 'EJS_1500m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -1500 -1500], tempyear, month, ...
%                                     [0 3], [0:0.3:3], 'temp', 'fig_param_kyy_EJS');
%     
%     if (exist(strcat(figdir,'EJS/','TUV_2000m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_2000m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_2000m/', 'EJS_2000m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -2000 -2000], tempyear, month, ...
%                                     [0 2], 0:0.2:2, 'temp', 'fig_param_kyy_EJS');
%     
% 
%     if (exist(strcat(figdir,'EJS/','TUV_2500m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_2500m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_2500m/', 'EJS_2500m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -2500 -2500], tempyear, month, ...
%                                     [0 2], 0:0.2:2, 'temp', 'fig_param_kyy_EJS');
% 
%     if (exist(strcat(figdir,'EJS/','TUV_3000m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','TUV_3000m'));
%     end
%     suvfile =[figdir, 'EJS/', 'TUV_3000m/', 'EJS_3000m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [127 144 33 52 -3000 -3000], tempyear, month, ...
%                                     [0 2], 0:0.2:2, 'temp', 'fig_param_kyy_EJS');
% 
%     lonlat_EKWC = [129 132 34 39]; % [lon_start lon_end lat_start lat_end]
%     if (exist(strcat(figdir,'EJS/','EKWC') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','EKWC'));
%     end
%     if (exist(strcat(figdir,'EJS/','EKWC/','SUV') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','EKWC/','SUV'));
%     end
%     suvfile =[figdir, 'EJS/', 'EKWC/','SUV/', 'EJS_SUV_EKWC']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=vec_ROMS_monthly_surface_UV(testname, suvfile, outputdir, lonlat_EKWC, tempyear, month, ...
%                                     'UV', 'EJS', 'fig_param_kyy_EKWC');
%     if (exist(strcat(figdir,'EJS/','EKWC/','TUV_100m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','EKWC/','TUV_100m'));
%     end
%     suvfile =[figdir, 'EJS/', 'EKWC/','TUV_100m/', 'EKWC_100m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [129 132 34 39 -100 -100], tempyear, month, ...
%                                     [0 30], 0:5:30, 'temp', 'fig_param_kyy_EKWC');
%     if (exist(strcat(figdir,'EJS/','EKWC/','TUV_200m') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','EKWC/','TUV_200m'));
%     end
%     suvfile =[figdir, 'EJS/', 'EKWC/','TUV_200m/', 'EKWC_200m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [129 132 34 39 -200 -200], tempyear, month, ...
%                                     [0 30], 0:5:30, 'temp', 'fig_param_kyy_EKWC');
% 
% %    if (exist(strcat(figdir,'EJS/','EKWC/','TUV_350m') , 'dir') ~= 7)
% %        mkdir(strcat(figdir,'EJS/','EKWC/','TUV_350m'));
% %    end
% %    suvfile =[figdir, 'EJS/', 'EKWC/','TUV_350m/', 'EKWC_350m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
% %    status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [129 132 34 39 -350 -350], tempyear, month, ...
% %                                    [0 30], 0:5:30, 'temp', 'fig_param_kyy_EKWC');
% 
% %    if (exist(strcat(figdir,'EJS/','EKWC/','TUV_500m') , 'dir') ~= 7)
% %        mkdir(strcat(figdir,'EJS/','EKWC/','TUV_500m'));
% %    end
% %    suvfile =[figdir, 'EJS/', 'EKWC/','TUV_500m/', 'EKWC_500m_TUV']; % % ~/SUV_EJS_testname_year_month.jpg
% %    status=plot_ROMS_monthly_hor_t_uv(testname, suvfile, outputdir, [129 132 34 39 -500 -500], tempyear, month, ...
% %                                    [0 30], 0:5:30, 'temp', 'fig_param_kyy_EKWC');
% 
%     if (exist(strcat(figdir,'EJS/','EKWC/','VORTICITY') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','EKWC/','VORTICITY'));
%     end
%     suvfile =[figdir, 'EJS/', 'EKWC/','VORTICITY/', 'EJS_SUV_relative_vorticity_EKWC']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_vorticity(testname, suvfile, outputdir, lonlat_EKWC, tempyear, month, ...
%                                    [-0.0001 0.0001], 0:10:10, 'vorticity', 'EJS', 'fig_param_kyy_EKWC');
% 
% 
%      % %  plot vertical temperature
%     if (exist(strcat(figdir,'EJS/','VERTICAL') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','VERTICAL'));
%     end
%     if (exist(strcat(figdir,'EJS/','VERTICAL/','LAT_FIXED') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','VERTICAL/','LAT_FIXED'));
%     end
%     if (exist(strcat(figdir,'EJS/','VERTICAL/','LAT_FIXED/TEMP') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','VERTICAL/','LAT_FIXED/TEMP'));
%     end
%     vert_temp_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/TEMP/', 'temp_1000m_lon_129_136_lat_38']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 136 38 38 -1000 0], tempyear, month, ...
%                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
%     vert_temp_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/TEMP/', 'temp_1000m_lon_129_133_lat_37']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 133 37 37 -1000 0], tempyear, month, ...
%                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
%     vert_temp_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/TEMP/', 'temp_1000m_lon_129_133_lat_36']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 133 36 36 -1000 0], tempyear, month, ...
%                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
%     vert_temp_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/TEMP/', 'temp_1000m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 133 35 35 -1000 0], tempyear, month, ...
%                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
%     vert_temp_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/TEMP/', 'temp_200m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129 133 35 35 -200 0], tempyear, month, ...
%                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
%                                 
%     if (exist(strcat(figdir,'EJS/','VERTICAL/','LON_FIXED') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','VERTICAL/','LON_FIXED'));
%     end
%     if (exist(strcat(figdir,'EJS/','VERTICAL/','LON_FIXED','/TEMP') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','VERTICAL/','LON_FIXED','/TEMP'));
%     end
%     vert_temp_file =[figdir, 'EJS/', 'VERTICAL/LON_FIXED/TEMP/', 'temp_1000m_lon_129_8_lat_35_43']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [129.8 129.8 35 41 -1000 0], tempyear, month, ...
%                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
%                                 
%     vert_temp_file =[figdir, 'EJS/', 'VERTICAL/LON_FIXED/TEMP/', 'temp_1000m_lon_132_lat_35_43']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, outputdir, [132 132 35 43 -1000 0], tempyear, month, ...
%                                     [0 30], [0 2 5 10 15 20 25 30], 'vert_temp', 'fig_param_kyy_EJS');
% 
%      % %  plot vertical salinity
%     if (exist(strcat(figdir,'EJS/','VERTICAL/','LAT_FIXED','/SALT') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','VERTICAL/','LAT_FIXED','/SALT'));
%     end
%     vert_salt_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/SALT/', 'salt_1000m_lon_129_136_lat_38']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 136 38 38 -1000 0], tempyear, month, ...
%                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
%     vert_salt_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/SALT/', 'salt_1000m_lon_129_133_lat_37']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 133 37 37 -1000 0], tempyear, month, ...
%                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
%     vert_salt_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/SALT/', 'salt_1000m_lon_129_133_lat_36']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 133 36 36 -1000 0], tempyear, month, ...
%                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
%     vert_salt_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/SALT/', 'salt_1000m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 133 35 35 -1000 0], tempyear, month, ...
%                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
%     vert_salt_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/SALT/', 'salt_200m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129 133 35 35 -200 0], tempyear, month, ...
%                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
%     
%     
%     if (exist(strcat(figdir,'EJS/','VERTICAL/','LON_FIXED','/SALT') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','VERTICAL/','LON_FIXED','/SALT'));
%     end
%     
%     vert_salt_file =[figdir, 'EJS/', 'VERTICAL/LON_FIXED/SALT/', 'salt_1000m_lon_129_8_lat_35_43']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [129.8 129.8 35 41 -1000 0], tempyear, month, ...
%                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
%                                 
%     vert_salt_file =[figdir, 'EJS/', 'VERTICAL/LON_FIXED/SALT/', 'salt_1000m_lon_132_lat_35_43']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_salt_file, outputdir, [132 132 35 43 -1000 0], tempyear, month, ...
%                                     [33 34.5], [33:0.1:34.5], 'vert_salt', 'fig_param_kyy_EJS');
%                                 
%                                 
%                                 
%      % %  plot vertical v
%     if (exist(strcat(figdir,'EJS/','VERTICAL/','LAT_FIXED','/V') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS/','VERTICAL/','LAT_FIXED','/V'));
%     end
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_1000m_lon_129_136_lat_38']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 136 38 38 -1000 0], tempyear, month, ...
%                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_1000m_lon_129_133_lat_37']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 133 37 37 -1000 0], tempyear, month, ...
%                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_1000m_lon_129_133_lat_36']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 133 36 36 -1000 0], tempyear, month, ...
%                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_1000m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 133 35 35 -1000 0], tempyear, month, ...
%                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_200m_lon_129_133_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 133 35 35 -200 0], tempyear, month, ...
%                                     [-0.5 0.5], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
% %     vert_v_file =[figdir, 'EJS/', 'vert_v/', 'EJS_vert_v_102']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129.8 131 36.0767 36.0767 -200 0], tempyear, month, ...
% %                                     [-1.0 1.0], -1.0:0.1:1.0, 'vert_v', 'fig_param_kyy_EJS');  %% for 102 line
%                                 
% %     vert_v_file =[figdir, 'EJS/', 'vert_v/', 'EJS_vert_v_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
% %     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
% %                                     [-1.0 1.0], -1.0:0.1:1.0, 'vert_v', 'fig_param_kyy_EJS');  %% for 35.5N
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_141_lat_35']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 141 35 35 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_141_lat_36']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 141 36 36 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
% 
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_141_lat_37']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 141 37 37 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_141_lat_38']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 141 38 38 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_141_lat_39']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 141 39 39 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_141_lat_40']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 141 40 40 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_141_lat_41']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 141 41 41 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_143_lat_42']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 143 42 42 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_127_143_lat_43']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [127 141 43 43 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_135_143_lat_44']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [135 143 44 44 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_135_143_lat_45']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [135 143 45 45 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_135_143_lat_46']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [135 143 46 46 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_135_143_lat_47']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [135 143 47 47 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
%     vert_v_file =[figdir, 'EJS/', 'VERTICAL/LAT_FIXED/V/', 'v_4000m_lon_135_143_lat_48']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, outputdir, [135 143 48 48 -4000 0], tempyear, month, ...
%                                     [-0.3 0.3], -0.5:0.1:0.5, 'vert_v', 'fig_param_kyy_EJS');
%                    
% end

% % % % % plot vertical section
% % % 
% % % NIFS KHOA, South Sea & East Sea
% nifs_num=[207, 208, 209, 102, 103, 104, 105, 106, 107];
% nifs_area=[ 129.12   130.6717  35.0217  34.15;    ...
%             129.455  130.83    35.475   34.5417;  ...
%             129.5317 131.6517  35.795   34.9817;  ...
%             129.5883 132.455   36.0767  36.0767;  ...
%             129.4833 132.475   36.505   36.505;   ...
%             129.46   134.34    37.0567  37.0567;  ...
%             129.15   132.8067  37.5533  37.5533;  ...
%             128.875  132.8267  37.895   37.895;   ...
%             128.6233 132.313   38.126   38.126;   ];
%                     
% for i=1:length(year) 
%     tempyear=year(i);
%     for var=1:2
%         if (var==1)
%             varname='vert_temp'
%             var_lim = [0 30];
%             var_con_lim= [0 1 2 3 5 10 15 20 25 30];
%         elseif (var==2)
%             varname='vert_salt'
%             var_lim = [33 34.5];
%             var_con_lim = [33:0.1:34.5];
%         end
% 
%         for i=1:5
%             depthlist=[-100, -200, -350, -500, -1000]
%             maxdepth = depthlist(i);
%             for nifs_ind=1:length(nifs_num)
%                 if (exist(strcat(figdir,'NIFS/',num2str(nifs_num(nifs_ind))) , 'dir') ~= 7)
%                     mkdir(strcat(figdir,'NIFS/',num2str(nifs_num(nifs_ind))));
%                 end
%                 nifs_file =[figdir,'NIFS/',num2str(nifs_num(nifs_ind)),'/', ...
%                             varname,'_',num2str(abs(maxdepth),'%04i'),'m_',num2str(nifs_num(nifs_ind)),'_line']; 
%                 nifs_section(1:4)=nifs_area(nifs_ind,1:4);
%                 nifs_section(5)=maxdepth;  nifs_section(6) =0;
%                 status=plot_ROMS_monthly_vertical_data(testname, nifs_file, outputdir, nifs_section, tempyear, month, ...
%                                                 var_lim, var_con_lim, varname, 'fig_param_kyy_EJS');    
%             end
%         end 
%     end
% end
% 
