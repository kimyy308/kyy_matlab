% % This code based on MATLAB R2016b.
% % Updated 05-Apr-2018 by Yong-Yub Kim

clc;close all;clear all;
warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\HYCOM\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Reanalysis\MOM\kiost_nwp_oda_1_2']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/HYCOM/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Reanalysis/MOM/kiost_nwp_oda_1_2']));
end


% for snu_desktop
testname='kiost_nwp_oda_1_2'   % % need to change
year = [1985:2010]; % % put year which you want to plot [year year ...]
% year = [110]; % % put year which you want to plot [year year ...]
% month = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
month = [13]; % % put month which you want to plot [month month ...]
figdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd year\figure\',testname,'\'); % % where figure files will be saved
inputdir = strcat('E:\Data\Reanalysis\',testname,'\'); % % where data files are

if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 


% % % % % % plot topography

% etopodir ='D:\MEPL\project\SSH\2nd_year\figure\etopo1\';
% plot_hycom_kimyy_etopo

% if (exist(strcat(figdir,'\Bathy\bathy_nwp_1_20.jpg'),'file') ~= 2)
%     if (exist(strcat(figdir,'Bathy') , 'dir') ~= 7)
%         mkdir(strcat(figdir, 'Bathy'));
%     end 
%     plot_hycom_kimyy_topo;
% end

% % % % % %  plot EJS
lonlat_EJS = [127 144 33 52]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   
    tempyear=year(i);
    if (exist(strcat(figdir,'EJS'),'dir')~=7)
        mkdir(strcat(figdir,'EJS'));
    end
    if (exist(strcat(figdir,'EJS','\hor_temp'),'dir')~=7)
        mkdir(strcat(figdir,'EJS','\hor_temp'));
    end                            
    hor_tuv_file =[figdir, 'EJS\', 'hor_temp\', 'EJS_hor_temp_0100m']; % % ~/SUV_EJS_testname_year_month.jpg
    status=plot_MOM_monthly_horizontal_t_uv(testname, hor_tuv_file, inputdir, [125 135 34 42 -100 -100], tempyear, month, [0 30], [5:5:10], 'temp',  'fig_param_kyy_MICT_EKWC');
    
    hor_tuv_file =[figdir, 'EJS\', 'hor_temp\', 'EJS_hor_temp_0000m']; % % ~/SUV_EJS_testname_year_month.jpg
    status=plot_MOM_monthly_horizontal_t_uv(testname, hor_tuv_file, inputdir, [125 135 34 42 0 0], tempyear, month, [0 30], [5:5:10], 'temp',  'fig_param_kyy_MICT_EKWC');
end




% % % % % %  plot EKB
lonlat_EKB = [127 130 37 42]; % [lon_start lon_end lat_start lat_end]
for i=1:length(year)   
    tempyear=year(i);
    if (exist(strcat(figdir,'EJS','\EKB'),'dir')~=7)
        mkdir(strcat(figdir,'EJS','\EKB'));
    end
    
%     % %  plot SST (EJS)
%     if (exist(strcat(figdir,'EJS\','SST') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS\','SST'));
%     end
%     sstfile =[figdir, 'EJS\', 'SST\', 'EJS_SST']; % % ~/SST_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_surface_data(testname, sstfile, inputdir, lonlat_EJS, tempyear, month, ...
%                                     [0 35], 0:5:35, 'SST', 'fig_param_kyy_EJS');                                

    if (exist(strcat(figdir,'EJS','\EKB\hor_temp'),'dir')~=7)
        mkdir(strcat(figdir,'EJS','\EKB\hor_temp'));
    end
     % %  plot horizontal temperature
    hor_temp_file =[figdir, 'EJS\EKB\', 'hor_temp\', 'EKB_hor_temp_50m']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_MOM_monthly_horizontal_data_pollack(testname, hor_temp_file, inputdir, [127 130 38 41 -50 -50], tempyear, month, ...
%                                     [2 5], 2:1:5, 'temp', 'fig_param_kyy_EKB');
%     status=plot_MOM_monthly_horizontal_data(testname, hor_temp_file, inputdir, [127 130 38 41 -50 -50], tempyear, month, ...
%                                     [0 15], 2:1:5, 'temp', 'fig_param_kyy_EKB');
                                

     % %  plot vertical temperature
    if (exist(strcat(figdir,'EJS\','vert_temp') , 'dir') ~= 7)
        mkdir(strcat(figdir,'EJS\','vert_temp'));
    end
%     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_ykang_lat']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [132 132 35 43 -1000 0], tempyear, month, ...
%                                     [0 15], 0:1:15, 'vert_temp', 'fig_param_kyy_EJS');
%                                 
%     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_ykang_lon']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [129 136 38 38 -1000 0], tempyear, month, ...
%                                     [0 30], 0:1:35, 'vert_temp', 'fig_param_kyy_EJS');

%     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_102']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [129.8 131 36.0767 36.0767 -200 0], tempyear, month, ...
%                                     [0 25], 0:5:25, 'vert_temp', 'fig_param_kyy_EJS');  %% for 102 line

%     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_bwr_102']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [129.8 131 36.0767 36.0767 -200 0], tempyear, month, ...
%                                     [6 14], 0:10:10, 'vert_temp', 'fig_param_kyy_EJS');  %% for 102 line_bwr_map
%         vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_bwr_36_39']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [129.8 129.8 36.0767 39 -200 0], tempyear, month, ...
%                                     [6 14], 0:10:10, 'vert_temp', 'fig_param_kyy_EJS');  %% for 102, NS_bwr_map
                                
%     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
%                                     [0 25], 0:5:25, 'vert_temp', 'fig_param_kyy_EJS');  %% 35.5N
                                
%     vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_WOA']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [130 130 35 39 -200 0], tempyear, month, ...
%                                     [2 20], 0:5:20, 'vert_temp', 'fig_param_kyy_EJS');   %% for WOA
                                
%      % %  plot vertical salinity
%     if (exist(strcat(figdir,'EJS\','vert_salt') , 'dir') ~= 7)
%         mkdir(strcat(figdir,'EJS\','vert_salt'));
%     end
%     vert_temp_file =[figdir, 'EJS\', 'vert_salt\', 'EJS_vert_salt_ykang_lat']; 
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [132 132 35 43 -1000 0], tempyear, month, ...
%                                     [33 34.5], 33:0.1:34.5, 'vert_salt', 'fig_param_kyy_EJS');
%     
%     vert_temp_file =[figdir, 'EJS\', 'vert_salt\', 'EJS_vert_salt_ykang_lon']; 
%     status=plot_ROMS_monthly_vertical_data(testname, vert_temp_file, inputdir, [129 136 38 38 -1000 0], tempyear, month, ...
%                                     [33 34.5], 33:0.1:34.5, 'vert_salt', 'fig_param_kyy_EJS');


     % %  plot vertical v
    if (exist(strcat(figdir,'EJS\','vert_v') , 'dir') ~= 7)
        mkdir(strcat(figdir,'EJS\','vert_v'));
    end
    
%     vert_v_file =[figdir, 'EJS\', 'vert_v\', 'EJS_vert_v_102']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, inputdir, [129.8 131 36.0767 36.0767 -200 0], tempyear, month, ...
%                                     [-1.0 1.0], -1.0:0.1:1.0, 'vert_v', 'fig_param_kyy_EJS');  %% for 102 line
                                
%     vert_v_file =[figdir, 'EJS\', 'vert_v\', 'EJS_vert_v_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
%     status=plot_ROMS_monthly_vertical_data(testname, vert_v_file, inputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
%                                     [-1.0 1.0], -1.0:0.1:1.0, 'vert_v', 'fig_param_kyy_EJS');  %% for 35.5N
                                   
end


















% 
% % % plot SSS
% sssfile =[figdir, '\EJS\ES_1_20_SSS_'] %% East Sea
% saltlonlat = [127 144 33 52]; % [lon_start lon_end lat_start lat_end] (East Sea)
% status=plot_SSS(sssfile, workdir, saltlonlat, year, month, inputdir);
% 
% sssfile =[figdir, 'SSS_nwp_1_20_']; % [lon_start lon_end lat_start lat_end] (NWP)
% saltlonlat = [115 164 15 52]; % [lon_start lon_end lat_start lat_end] (NWP)
% status=plot_SSS(sssfile, workdir, saltlonlat, year, month, inputdir);
% 
% % % saltlonlat = [117 135 27 44]; % [lon_start lon_end lat_start lat_end]
% % close all;
% % % % % 


% % % % plot UV
% % % % % must change about Reference vector (value, text)
% % % 
% for regionflag=1:2;  %% NWP =1, ES =2
%     if (regionflag==1)
%         uvfile =[figdir, 'UV_nwp_1_20_']  %% NWP
%         uvlonlat = [115 164 15 52]; % [lon_start lon_end lat_start lat_end] (NWP)
%         status=plot_surf_UV(uvfile, workdir, uvlonlat, year, month, inputdir, regionflag);
%     elseif (regionflag==2)
%         uvfile =[figdir, '\EJS\ES_1_20_UV']
%         uvlonlat = [127 144 33 52]; % [lon_start lon_end lat_start lat_end] (East Sea)          
%         status=plot_surf_UV(uvfile, workdir, uvlonlat, year, month, inputdir, regionflag);
%     end
% end
% % close all;
% % %


% for year=1995:1995
%     uvfile =[figdir, 'model_nwp_uv_nwp_1_20_UV']
%     uvlonlat = [115 164 15 52]; % [lon_start lon_end lat_start lat_end] (NWP)
%     year = [1995 1995]; % [year_start year_end]
% %     month = [2 2] % [first month of the year_start, last month of the year_end]
%     month = [8 8] % [first month of the year_start, last month of the year_end]
%     status=plot_surf_UV(uvfile, workdir, uvlonlat, year, month, inputdir);
% end
% % 




% h_year = year(2)-1991;
% h_uvfile =[figdir, 'H_UV_kuro_1_20_', num2str(h_year),'_']
% h_uvlonlat = [132 143 27 35]; % [lon_start lon_end lat_start lat_end] (Kuro)
% status=plot_h_UV(h_uvfile, workdir, h_uvlonlat, year, month, inputdir);
% h_uvfile =[figdir, 'H_UV_nwp_1_20_', num2str(h_year),'_']
% h_uvlonlat = [115 164 15 52]; % [lon_start lon_end lat_start lat_end] (nwp)
% status=plot_h_UV(h_uvfile, workdir, h_uvlonlat, year, month, inputdir);
% h_uvfile =[figdir, 'H_UV_luzon_1_20_', num2str(h_year),'_']
% h_uvlonlat = [117 124 16 23]; % [lon_start lon_end lat_start lat_end] (luzon)
% status=plot_h_UV(h_uvfile, workdir, h_uvlonlat, year, month, inputdir);
% h_uvfile =[figdir, 'H_UV_ES_1_20_', num2str(h_year),'_']
% h_uvlonlat = [127 144 33 52]; % [lon_start lon_end lat_start lat_end] (East Sea)
% status=plot_h_UV(h_uvfile, workdir, h_uvlonlat, year, month, inputdir);


% % % % plot h_uv_seo
% h_uvfile =['D:\MEPL\project\SSH\figure\seo_avg\', 'H_UV_kuro_1_10_']
% h_uvlonlat = [132 143 27 35]; % [lon_start lon_end lat_start lat_end] (kuro)
% status=plot_h_UV_1_10(h_uvfile, 'D:\MEPL\project\SSH\data\seo_avg\', h_uvlonlat, [2001 2001], month, 'D:\MEPL\project\SSH\data\seo_avg\');
% h_uvfile =['D:\MEPL\project\SSH\figure\seo_avg\', 'H_UV_luzon_1_10_']
% h_uvlonlat = [117 124 16 23]; % [lon_start lon_end lat_start lat_end] (kuro)
% status=plot_h_UV_1_10(h_uvfile, 'D:\MEPL\project\SSH\data\seo_avg\', h_uvlonlat, [2001 2001], month, 'D:\MEPL\project\SSH\data\seo_avg\');


% % % % Yellow Sea Bottom Cold Water
% btfile =[figdir,'YSBCW_nwp_1_20']; 
% botlonlat = [115 130 30 42]; % [lon_start lon_end lat_start lat_end]  %%YSBCW
% status=plot_bottom_temp(btfile, workdir, botlonlat, year, month, inputdir);
% % botlonlat = [117 127 33 42]; % [lon_start lon_end lat_start lat_end]



% % % % % plot vertical section
% % % 
% % % NIFS KHOA, East Sea

% for var=1:1
%     if (var==1)
%         varname='temp'
%         var_lim = [2 34]
%     elseif (var==2)
%         varname='salt'
%         var_lim = [31 35]
%     end
% 
%     for i=2:5
%         maxdepth = -(i*100)
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_207_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.12 130.6717 35.0217 34.15 maxdepth 0]; %% 207 (1-13) (south sea)
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
% 
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_208_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.455 130.83 35.475 34.5417 maxdepth 0]; %% 208 (1-10)
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_209_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.5317 131.6517 35.795 34.9817 maxdepth 0]; %% 209 (00-13)
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_102_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.5883 132.455 36.0767 36.0767 maxdepth 0]; %% 102 (00-15)
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         	
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_103_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.4833 132.475 36.505 36.505 maxdepth 0]; %% 103 (00-15)
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_104_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.46 134.34 37.0567 37.0567 maxdepth 0]; %% 104 (00-21)
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_105_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.15 132.8067 37.5533 37.5533 maxdepth 0]; %% 105 (00-16)
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_106_line_',num2str(abs(maxdepth)),'_']; 
%         section = [128.875 132.8267 37.895 37.895 maxdepth 0]; %% 106 (02-16) 
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_107_line_',num2str(abs(maxdepth)),'_']; 
%         section = [128.6233 132.313	38.126 38.126 maxdepth 0]; %% 107 (01-15) 
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_EKWC_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.7 129.7 34.15 38.16 maxdepth 0]; %% 129.7 E
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%         
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_EKWC2_line_',num2str(abs(maxdepth)),'_']; 
%         section = [130.0 130.0 34.15 38.16 maxdepth 0]; %% 130 E
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
% 
%         vertfile =[figdir,'\EJS\ES_1_20_vert_',varname,'_EKWC3_line_',num2str(abs(maxdepth)),'_']; 
%         section = [129.7 129.7 34.5 36.0 maxdepth 0]; %% 129.7 E
%         status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
%     end 
% end


% % vertfile =[figdir,'model_nwp_ys_vert_temp_35N']; 
% % vertfile =[figdir,'model_nwp_ys_vert_temp_124_4E']; 
% % vertfile =[figdir,'model_nwp_ES_vert_temp_129_8E']; 
% % section = [119 127 35 35.05 -100 0];
% % section = [124.4 124.45 34 37 -100 0];
% % section = [129.8 129.85 33 38 -300 0];
% section = [129.8 129.85 35 39 -200 0];
% var = 1;
% % var_lim = [31 35];
% var_lim = [2 20];

% status=plot_vertical_var(vertfile, section, year, month, var, var_lim, inputdir);
% 
% for i=-50:-50:-200
%     horsection = [129 132 34 40 i i];
%     horfile =[figdir,'\EJS\EKWC_',num2str(abs(horsection(5))),'m_']; 
%     horvar = 1;
%     horvar_lim = [-2 33];
%     status=plot_horizontal_var(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
% end
%     horsection = [129 132 34 40 -300 -300];
%     horfile =[figdir,'\EJS\EKWC_',num2str(abs(horsection(5))),'m_']; 
%     horvar = 1;
%     horvar_lim = [-2 33];
%     status=plot_horizontal_var(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
%     
%     horsection = [129 132 34 40 -500 -500];
%     horfile =[figdir,'\EJS\EKWC_',num2str(abs(horsection(5))),'m_']; 
%     horvar = 1;
%     horvar_lim = [-2 33];
%     status=plot_horizontal_var(horfile, horsection, year, month, horvar, horvar_lim, inputdir);


% % % % horizontal temp + uv
%     horsection = [129 132 34 40 -100 -100];
%     horfile =[figdir,'\EJS\EKWC_',num2str(abs(horsection(5))),'m_temp_uv_']; 
%     horvar = 1;
%     horvar_lim = [-2 33];
%     status=plot_horizontal_t_uv(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
%     
%     horsection = [129 132 34 40 -200 -200];
%     horfile =[figdir,'\EJS\EKWC_',num2str(abs(horsection(5))),'m_temp_uv_']; 
%     horvar = 1;
%     horvar_lim = [-2 33];
%     status=plot_horizontal_t_uv(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
%     
%     horsection = [128 132.5 35.5 40 -340 -340];
%     horfile =[figdir,'\EJS\NKCC_',num2str(abs(horsection(5))),'m_temp_uv_']; 
%     horvar = 1;
%     horvar_lim = [-2 33];
%     status=plot_horizontal_t_uv(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
%     
%     horsection = [128 135 35.5 44 -235 -235];
%     horfile =[figdir,'\EJS\NKCC_',num2str(abs(horsection(5))),'m_temp_uv_']; 
%     horvar = 2;
%     horvar_lim = [33.95 34.15];
%     status=plot_horizontal_t_uv(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
%     
%     get_EKWC_info;
%     
%     
% end
