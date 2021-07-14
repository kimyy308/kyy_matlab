% % This code based on MATLAB R2016b.

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
    addpath(genpath([dropboxpath '\source\matlab\Reanalysis\HYCOM\Analysis\expt19_1']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/HYCOM/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Reanalysis/HYCOM/Analysis/expt19_1']));
end


%% winter surface current

testname='GLBb_0_08_expt19_1'   % % need to change
tempyear = [0110]; % % put year which you want to plot [year year ...]
tempmonth = 2;
figdir =strcat('D:\OneDrive - 서울대학교\Ph_D_course\Figure\Reanalysis\HYCOM\',testname,'\'); % % where figure files will be saved
inputdir = strcat('E:\Data\Reanalysis\HYCOM\', testname, '\'); % % where data files are                          
                            
suvfile =[figdir, 'EJS\', 'SUV\', 'EJS_SUV']; % % ~/SUV_EJS_testname_year_month.jpg
% status=vec_HYCOM_monthly_surface_UV(testname, suvfile, inputdir, lonlat_EJS, tempyear, month, ...
%                                 'UV', 'EJS', 'fig_param_HYCOM_kyy_EJS');
var = 'UV';
param_script = 'fig_param_HYCOM_kyy_EKWC';
filedir = inputdir;
lonlat = [129 133 34 41];

run(param_script);
filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                'HYCOM_', num2str(tempyear,'%04i'), num2str(tempmonth,'%02i'), '.nc');

% read data
if (exist('lon')==0)
    lon = ncread(filename,'lon');
    lat = ncread(filename,'lat');

    lon_west = abs(lon - (lonlat(1)-1));
    min_lon_west=min(lon_west);
    lon_east = abs(lon - (lonlat(2)+1));
    min_lon_east=min(lon_east);
    lat_south = abs(lat - (lonlat(3)-1));
    min_lat_south=min(lat_south);
    lat_north = abs(lat - (lonlat(4)+1));
    min_lat_north=min(lat_north);

    lon_min = find(lon_west == min_lon_west);
    lon_max = find(lon_east == min_lon_east);
    lat_min = find(lat_south == min_lat_south);
    if (lonlat(4)>47)
        lat_max=length(lat);
    else
        lat_max = find(lat_north == min_lat_north);
    end

%             lon = lon(lon_min(1):lon_max(1));
%             lat = lat(lat_min(1):lat_max(1));
    lon = ncread(filename,'longitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
    lat = ncread(filename,'latitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
end

data_info_u = ncinfo(filename, 'u');
data_info_v = ncinfo(filename, 'v');

u_rho = ncread(filename,'u',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
v_rho = ncread(filename,'v',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
u_rho(find(abs(u_rho)>10000))=NaN;
v_rho(find(abs(v_rho)>10000))=NaN;
% Reference vector value setting

u_rho(m_quiver_ref_vec_x_range2,m_quiver_ref_vec_y_range2) = m_quiver_ref_u_value2;
v_rho(m_quiver_ref_vec_x_range2,m_quiver_ref_vec_y_range2) = m_quiver_ref_v_value2; 

% plot
subplot(3,4,[1 2 5 6 9 10]);
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
hold on;
uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)*m_quiver_vector_size, ...
                v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)*m_quiver_vector_size, ...
                m_quiver_vector_color,'AutoScale','off');
m_gshhs_i('color',m_gshhs_line_color)  
m_gshhs_i('patch',m_gshhs_land_color);
m_text(m_quiver_ref_text_x_location2, m_quiver_ref_text_y_location2, m_quiver_ref_text2, 'FontSize', m_quiver_ref_text_fontsize); 

% set grid
m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'xtick',(130:1:132), 'xticklabels',[130 131 132], 'ytick', (35:1:41) );

set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

disp(' ')
disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'UV', ' plot is created.'])
disp(' ')

hold off


%% summer surface current

testname='GLBb_0_08_expt19_1'   % % need to change
tempyear = [0110]; % % put year which you want to plot [year year ...]
tempmonth = 8;
figdir =strcat('D:\OneDrive - 서울대학교\Ph_D_course\Figure\Reanalysis\HYCOM\',testname,'\'); % % where figure files will be saved
inputdir = strcat('E:\Data\Reanalysis\HYCOM\', testname, '\'); % % where data files are                          
                            
suvfile =[figdir, 'EJS\', 'SUV\', 'EJS_SUV']; % % ~/SUV_EJS_testname_year_month.jpg
% status=vec_HYCOM_monthly_surface_UV(testname, suvfile, inputdir, lonlat_EJS, tempyear, month, ...
%                                 'UV', 'EJS', 'fig_param_HYCOM_kyy_EJS');
var = 'UV';
param_script = 'fig_param_HYCOM_kyy_EKWC';
filedir = inputdir;
lonlat = [129 133 34 41];

run(param_script);
filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                'HYCOM_', num2str(tempyear,'%04i'), num2str(tempmonth,'%02i'), '.nc');

% read data
if (exist('lon')==0)
    lon = ncread(filename,'lon');
    lat = ncread(filename,'lat');

    lon_west = abs(lon - (lonlat(1)-1));
    min_lon_west=min(lon_west);
    lon_east = abs(lon - (lonlat(2)+1));
    min_lon_east=min(lon_east);
    lat_south = abs(lat - (lonlat(3)-1));
    min_lat_south=min(lat_south);
    lat_north = abs(lat - (lonlat(4)+1));
    min_lat_north=min(lat_north);

    lon_min = find(lon_west == min_lon_west);
    lon_max = find(lon_east == min_lon_east);
    lat_min = find(lat_south == min_lat_south);
    if (lonlat(4)>47)
        lat_max=length(lat);
    else
        lat_max = find(lat_north == min_lat_north);
    end

%             lon = lon(lon_min(1):lon_max(1));
%             lat = lat(lat_min(1):lat_max(1));
    lon = ncread(filename,'longitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
    lat = ncread(filename,'latitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
end

data_info_u = ncinfo(filename, 'u');
data_info_v = ncinfo(filename, 'v');

u_rho = ncread(filename,'u',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
v_rho = ncread(filename,'v',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
u_rho(find(abs(u_rho)>10000))=NaN;
v_rho(find(abs(v_rho)>10000))=NaN;
% Reference vector value setting

u_rho(m_quiver_ref_vec_x_range2,m_quiver_ref_vec_y_range2) = m_quiver_ref_u_value2;
v_rho(m_quiver_ref_vec_x_range2,m_quiver_ref_vec_y_range2) = m_quiver_ref_v_value2; 

% plot
subplot(3,4,[3 4 7 8 11 12]);
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
hold on;
uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                m_quiver_vector_color,'AutoScale','off');
m_gshhs_i('color',m_gshhs_line_color)  
m_gshhs_i('patch',m_gshhs_land_color);
m_text(m_quiver_ref_text_x_location2, m_quiver_ref_text_y_location2, m_quiver_ref_text2, 'FontSize', m_quiver_ref_text_fontsize); 


% set grid
m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'xtick',(130:1:132), 'ytick', (35:1:41),  'xticklabels',[130 131 132] , 'yaxislocation', 'right');


set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

disp(' ')
disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'UV', ' plot is created.'])
disp(' ')

hold off




fig=gcf;
fig.InvertHardcopy='off';

hold off

set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [550, 550]); %% [xsize ysize]
set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 550, 550]) 
jpgname=strcat('D:\OneDrive - 서울대학교\research\Master_course\for_publish\Figure\Fig25', '.jpg'); %% ~_year_month.jpg
saveas(gcf,'D:\OneDrive - 서울대학교\research\Master_course\for_publish\new_figure\fig_0008','jpg');