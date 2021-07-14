% % This code based on MATLAB R2016b.

clc;close all;clear all;
warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
elseif (linux==1)
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end

% % for linux

% for snu_desktop
figdir ='D:\MEPL\project\SSH\1st_year\figure\test37\';
workdir = 'D:\MEPL\project\SSH\1st_year\data\test37\';
inputdir = 'D:\MEPL\project\SSH\1st_year\data\test37\';
% outfile ='D:\MEPL\project\SSH\figure\test24\model_nwp_uv_'; 


% workdir ='/data/kimyy';
% inputdir = '/data/kimyy/';
% outfile ='/home01/kimyy/figure/ROMS/test23/yellow_temp_35_'; 


% filename example : ~workdir\monthly_spinup_0001.nc
% name             = 'monthly_spinup_';
% filename_suffix  = '.nc';
% filename         = strcat(workdir,name,num2str(1,'%04i'),filename_suffix);

% lonlat = [127 144 33 52]; % [lon_start lon_end lat_start lat_end]
lonlat = [115 164 15 52]; % [lon_start lon_end lat_start lat_end]

year = [1996 1996] % [year_start year_end]
month = [1 12] % [first month of the year_start, last month of the year_end]
fast = 0;

% % % % % % plot topography
% plot_roms_kimyy_topo;
% 

% % plot SST
sstfile =[figdir, '\EJS\ES_1_20_SST_'] %% East Sea
sstlonlat = [127 144 33 52]; % [lon_start lon_end lat_start lat_end] (East Sea)
status=plot_SST(sstfile, workdir, sstlonlat, year, month, inputdir);

sstfile =[figdir, '\SST_nwp_1_20_'] %% NWP
sstlonlat = [115 164 15 52]; % [lon_start lon_end lat_start lat_end] (NWP)
status=plot_SST(sstfile, workdir, sstlonlat, year, month, inputdir);
% close all;

% % plot SSS
sssfile =[figdir, '\EJS\ES_1_20_SSS_'] %% East Sea
saltlonlat = [127 144 33 52]; % [lon_start lon_end lat_start lat_end] (East Sea)
status=plot_SSS(sssfile, workdir, saltlonlat, year, month, inputdir);

sssfile =[figdir, 'SSS_nwp_1_20_']; % [lon_start lon_end lat_start lat_end] (NWP)
saltlonlat = [115 164 15 52]; % [lon_start lon_end lat_start lat_end] (NWP)
status=plot_SSS(sssfile, workdir, saltlonlat, year, month, inputdir);

% % saltlonlat = [117 135 27 44]; % [lon_start lon_end lat_start lat_end]
% close all;
% % % % 


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


% % % horizontal temp + uv
    horsection = [129 132 34 40 -100 -100];
    horfile =[figdir,'\EJS\EKWC_',num2str(abs(horsection(5))),'m_temp_uv_']; 
    horvar = 1;
    horvar_lim = [-2 33];
    status=plot_horizontal_t_uv(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
    
    horsection = [129 132 34 40 -200 -200];
    horfile =[figdir,'\EJS\EKWC_',num2str(abs(horsection(5))),'m_temp_uv_']; 
    horvar = 1;
    horvar_lim = [-2 33];
    status=plot_horizontal_t_uv(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
    
    horsection = [128 132.5 35.5 40 -340 -340];
    horfile =[figdir,'\EJS\NKCC_',num2str(abs(horsection(5))),'m_temp_uv_']; 
    horvar = 1;
    horvar_lim = [-2 33];
    status=plot_horizontal_t_uv(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
    
    horsection = [128 135 35.5 44 -235 -235];
    horfile =[figdir,'\EJS\NKCC_',num2str(abs(horsection(5))),'m_temp_uv_']; 
    horvar = 2;
    horvar_lim = [33.95 34.15];
    status=plot_horizontal_t_uv(horfile, horsection, year, month, horvar, horvar_lim, inputdir);
    
    get_EKWC_info;
