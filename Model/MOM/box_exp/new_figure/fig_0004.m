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

u_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_u_value;
v_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_v_value; 

% plot
subplot(5,4,[1 2 5 6 9 10]);
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
hold on;
uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                m_quiver_vector_size,m_quiver_vector_color);
m_gshhs_i('color',m_gshhs_line_color)  
m_gshhs_i('patch',m_gshhs_land_color);
m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 

% set grid
m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'xtick',(129:2:133), 'xticklabels',[129 131 133], 'ytick', (35:1:41) );

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

u_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_u_value;
v_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_v_value; 

% plot
subplot(5,4,[3 4 7 8 11 12]);
m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
hold on;
uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                m_quiver_vector_size,m_quiver_vector_color);
m_gshhs_i('color',m_gshhs_line_color)  
m_gshhs_i('patch',m_gshhs_land_color);
m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 


% set grid
m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type, 'xtick',(129:2:133), 'ytick', (35:1:41),  'xticklabels',[129 131 133] , 'yaxislocation', 'right');


set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

disp(' ')
disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'UV', ' plot is created.'])
disp(' ')

hold off



%% winter vertical temperature

% vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
% status=plot_HYCOM_monthly_vertical_data(testname, vert_temp_file, inputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
%                                 [0 25], 0:5:25, 'vert_temp', 'fig_param_HYCOM_kyy_EJS');  %% 35.5N
testname='GLBb_0_08_expt19_1'   % % need to change
tempyear = [0110]; % % put year which you want to plot [year year ...]
tempmonth = 2;
var = 'vert_temp';
section = [129 131 35.5 35.5 -200 0];
shadlev = [0 25];
conlev = 0:5:25;
figdir =strcat('D:\OneDrive - 서울대학교\Ph_D_course\Figure\Reanalysis\HYCOM\',testname,'\'); % % where figure files will be saved
inputdir = strcat('E:\Data\Reanalysis\HYCOM\', testname, '\'); % % where data files are 
param_script = 'fig_param_HYCOM_kyy_EKWC';
filedir = inputdir;
filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
            'HYCOM_', num2str(tempyear,'%04i'), num2str(tempmonth,'%02i'), '.nc');
run(param_script);

% if (exist('lon')==0)
    lon = ncread(filename,'lon');
    lat = ncread(filename,'lat');

    lon_west = abs(lon - (section(1)-2));
    min_lon_west=min(lon_west);
    lon_east = abs(lon - (section(2)+2));
    min_lon_east=min(lon_east);
    lat_south = abs(lat - (section(3)-2));
    min_lat_south=min(lat_south);
    lat_north = abs(lat - (section(4)+2));
    min_lat_north=min(lat_north);

    lon_min = find(lon_west == min_lon_west);
    lon_max = find(lon_east == min_lon_east);
    lat_min = find(lat_south == min_lat_south);
    if (section(4)>47)
        lat_max=length(lat);
    else
        lat_max = find(lat_north == min_lat_north);
    end

%             lon = lon(lon_min(1):lon_max(1));
%             lat = lat(lat_min(1):lat_max(1));
    cut_lon_rho = ncread(filename,'longitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1])';
    cut_lat_rho = ncread(filename,'latitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1])';
% end
depth_1d = ncread(filename,'depth');

data_info = ncinfo(filename, varname); 
data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3) 1]);
data = permute(data, [3 2 1]); %% permute data to [z y x]
cut_depth = -repmat(depth_1d,1,lat_max(1)-lat_min(1)+1,lon_max(1)-lon_min(1)+1);

dist=sqrt((cut_lon_rho-section(1)).^2+(cut_lat_rho-section(3)).^2); %% get distance from station 1
min_dist=min(min(dist));
dist2=sqrt((cut_lon_rho-section(2)).^2+(cut_lat_rho-section(4)).^2);  %% get distance from station 2
min_dist2=min(min(dist2));                
[x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
[x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]

lat1=cut_lat_rho(x1(1),y1(1));  lon1=cut_lon_rho(x1(1),y1(1));
lat2=cut_lat_rho(x2(1),y2(1));  lon2=cut_lon_rho(x2(1),y2(1));
if (lon2-lon1) >= (lat2-lat1)
    lon_line = lon1:min(gradient(cut_lon_rho(1,:))):lon2;  
    lat_line = (lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1; %% weight for latitude index (lat1 : 0.05/((lon2-lon1) * (lat2-lat1) : lat2 
    x=repmat(lon_line,length(depth_1d),1);  %% copy lon_line (length(depth_1d) times) to make matrix 
    x_label='Longitude(^oE)';
%                 domaxis=[domaxis(3) domaxis(4) domaxis(5) domaxis(6) domaxis(5) domaxis(6)];
    Temp=zeros(length(depth_1d),length(lon_line)); %% initialize temp matrix that size is same with x
else
    lat_line=[min(lat1,lat2):mean(gradient(cut_lat_rho(:,1))):max(lat1,lat2)];
    lon_line=(lat_line-lat1)*((lon2-lon1)/(lat2-lat1))+lon1;
    x=repmat(lat_line,length(depth_1d),1);
    x_label='Latitude(^oN)';
    Temp=zeros(length(depth_1d),length(lat_line)); %% initialize temp matrix that size is same with x
end

if (x1(1)==x2(1) || y1(1)==y2(1)) %% fixed lon or lat
    Temp(:,:) = squeeze(data(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
    Yi(:,:)= squeeze(cut_depth(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
else
    for k=1:1:length(depth_1d)
        lon_range=cut_lon_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
        lat_range=cut_lat_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
        data_range=squeeze(data(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));  %%data(zlevel, latmin : latmax, lonmin : lonmax)
        depth_range=squeeze(cut_depth(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)))); 
        Temp(k,:)=griddata(lon_range,lat_range,data_range,lon_line,lat_line); %%get interpolated data section
        Yi(k,:)=griddata(lon_range,lat_range,depth_range,lon_line,lat_line); %%get interpolated depth section
    end
end

data=Temp;
data(find(abs(data)>10000))=NaN;

max(section(2)-section(1), section(4)-section(3));

%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [vert_paper_size_y*max(section(2)-section(1),section(4)-section(3))/2.0, vert_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, vert_paper_size_y*max(section(2)-section(1),section(4)-section(3))/2.0, vert_paper_size_y]) 

set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [500, vert_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 500, vert_paper_size_y]) 

%         text_posi_x=(section(4)-section(3))/20+section(3);
%         text_posi_y1=(section(6)-section(5))/20+section(5);
%         text_posi_y2=2*(section(6)-section(5))/20+section(5);
%         text_posi_y3=3*(section(6)-section(5))/20+section(5);

subplot(5,4,[13 14 15 16]);
hold on
pcolor(x,Yi,data)
if (lon2-lon1) >= (lat2-lat1)        
    axis([section(1) section(2) section(5) section(6)]);
    xlim([section(1) section(2)]);
else
    axis([section(3) section(4) section(5) section(6)]);
    xlim([section(3) section(4)]);
end
shading(gca,m_pcolor_shading_method);
set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)

set(gca,'yTickLabel',{'-200m', '-100m', '0m'}, ...
    'YTick', [-200, -100, -0])
set(gca,'xTickLabel',{' ','129.5^oE','130^oE','130.5^oE', ' '}, ...
    'XTick', (129:0.5:131))
set(gca, 'color', [0.8,0.8,0.8]);
colormap(gca,jet);

hold on
[C,h]=contour(x,Yi,data,conlev,vert_contour_color,'linewidth',vert_contour_linewidth);              
if strcmp(param_script,'fig_param_HYCOM_kyy_schematic')
%         no contour label
else
    clabel(C,h,'FontSize',vert_contour_label_fontsize,'Color',vert_contour_label_color,'labelspacing',vert_contour_labelspacing,'Rotation',vert_contour_rotation,'fontweight', vert_contour_fontweight);
end

hold off




%% summer vertical temperature

% vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
% status=plot_HYCOM_monthly_vertical_data(testname, vert_temp_file, inputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
%                                 [0 25], 0:5:25, 'vert_temp', 'fig_param_HYCOM_kyy_EJS');  %% 35.5N
testname='GLBb_0_08_expt19_1'   % % need to change
tempyear = [0110]; % % put year which you want to plot [year year ...]
tempmonth = 8;
var = 'vert_temp';
section = [129 131 35.5 35.5 -200 0];
shadlev = [0 25];
conlev = 0:5:25;
figdir =strcat('D:\OneDrive - 서울대학교\Ph_D_course\Figure\Reanalysis\HYCOM\',testname,'\'); % % where figure files will be saved
inputdir = strcat('E:\Data\Reanalysis\HYCOM\', testname, '\'); % % where data files are 
param_script = 'fig_param_HYCOM_kyy_EKWC';
filedir = inputdir;
filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
            'HYCOM_', num2str(tempyear,'%04i'), num2str(tempmonth,'%02i'), '.nc');
run(param_script);

% if (exist('lon')==0)
    lon = ncread(filename,'lon');
    lat = ncread(filename,'lat');

    lon_west = abs(lon - (section(1)-2));
    min_lon_west=min(lon_west);
    lon_east = abs(lon - (section(2)+2));
    min_lon_east=min(lon_east);
    lat_south = abs(lat - (section(3)-2));
    min_lat_south=min(lat_south);
    lat_north = abs(lat - (section(4)+2));
    min_lat_north=min(lat_north);

    lon_min = find(lon_west == min_lon_west);
    lon_max = find(lon_east == min_lon_east);
    lat_min = find(lat_south == min_lat_south);
    if (section(4)>47)
        lat_max=length(lat);
    else
        lat_max = find(lat_north == min_lat_north);
    end

%             lon = lon(lon_min(1):lon_max(1));
%             lat = lat(lat_min(1):lat_max(1));
    cut_lon_rho = ncread(filename,'longitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1])';
    cut_lat_rho = ncread(filename,'latitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1])';
% end
depth_1d = ncread(filename,'depth');

data_info = ncinfo(filename, varname); 
data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3) 1]);
data = permute(data, [3 2 1]); %% permute data to [z y x]
cut_depth = -repmat(depth_1d,1,lat_max(1)-lat_min(1)+1,lon_max(1)-lon_min(1)+1);

dist=sqrt((cut_lon_rho-section(1)).^2+(cut_lat_rho-section(3)).^2); %% get distance from station 1
min_dist=min(min(dist));
dist2=sqrt((cut_lon_rho-section(2)).^2+(cut_lat_rho-section(4)).^2);  %% get distance from station 2
min_dist2=min(min(dist2));                
[x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
[x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]

lat1=cut_lat_rho(x1(1),y1(1));  lon1=cut_lon_rho(x1(1),y1(1));
lat2=cut_lat_rho(x2(1),y2(1));  lon2=cut_lon_rho(x2(1),y2(1));
if (lon2-lon1) >= (lat2-lat1)
    lon_line = lon1:min(gradient(cut_lon_rho(1,:))):lon2;  
    lat_line = (lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1; %% weight for latitude index (lat1 : 0.05/((lon2-lon1) * (lat2-lat1) : lat2 
    x=repmat(lon_line,length(depth_1d),1);  %% copy lon_line (length(depth_1d) times) to make matrix 
    x_label='Longitude(^oE)';
%                 domaxis=[domaxis(3) domaxis(4) domaxis(5) domaxis(6) domaxis(5) domaxis(6)];
    Temp=zeros(length(depth_1d),length(lon_line)); %% initialize temp matrix that size is same with x
else
    lat_line=[min(lat1,lat2):mean(gradient(cut_lat_rho(:,1))):max(lat1,lat2)];
    lon_line=(lat_line-lat1)*((lon2-lon1)/(lat2-lat1))+lon1;
    x=repmat(lat_line,length(depth_1d),1);
    x_label='Latitude(^oN)';
    Temp=zeros(length(depth_1d),length(lat_line)); %% initialize temp matrix that size is same with x
end

if (x1(1)==x2(1) || y1(1)==y2(1)) %% fixed lon or lat
    Temp(:,:) = squeeze(data(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
    Yi(:,:)= squeeze(cut_depth(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
else
    for k=1:1:length(depth_1d)
        lon_range=cut_lon_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
        lat_range=cut_lat_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
        data_range=squeeze(data(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));  %%data(zlevel, latmin : latmax, lonmin : lonmax)
        depth_range=squeeze(cut_depth(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)))); 
        Temp(k,:)=griddata(lon_range,lat_range,data_range,lon_line,lat_line); %%get interpolated data section
        Yi(k,:)=griddata(lon_range,lat_range,depth_range,lon_line,lat_line); %%get interpolated depth section
    end
end

data=Temp;
data(find(abs(data)>10000))=NaN;

max(section(2)-section(1), section(4)-section(3));

%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [vert_paper_size_y*max(section(2)-section(1),section(4)-section(3))/2.0, vert_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, vert_paper_size_y*max(section(2)-section(1),section(4)-section(3))/2.0, vert_paper_size_y]) 

set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [500, vert_paper_size_y]);
set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 500, vert_paper_size_y]) 

%         text_posi_x=(section(4)-section(3))/20+section(3);
%         text_posi_y1=(section(6)-section(5))/20+section(5);
%         text_posi_y2=2*(section(6)-section(5))/20+section(5);
%         text_posi_y3=3*(section(6)-section(5))/20+section(5);

subplot(5,4,[17 18 19 20]);
hold on
pcolor(x,Yi,data)
if (lon2-lon1) >= (lat2-lat1)        
    axis([section(1) section(2) section(5) section(6)]);
    xlim([section(1) section(2)]);
else
    axis([section(3) section(4) section(5) section(6)]);
    xlim([section(3) section(4)]);
end
shading(gca,m_pcolor_shading_method);
set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
   
set(gca,'yTickLabel',{'-200m', '-100m', '0m'}, ...
    'YTick', [-200, -100, -0])
set(gca,'xTickLabel',{'129^oE','129.5^oE','130^oE','130.5^oE', '131^oE'}, ...
    'XTick', (129:0.5:131))
set(gca, 'color', [0.8,0.8,0.8]);
colormap(gca,jet);

hold on
[C,h]=contour(x,Yi,data,conlev,vert_contour_color,'linewidth',vert_contour_linewidth);              
if strcmp(param_script,'fig_param_HYCOM_kyy_schematic')
%         no contour label
else
    clabel(C,h,'FontSize',vert_contour_label_fontsize,'Color',vert_contour_label_color,'labelspacing',vert_contour_labelspacing,'Rotation',vert_contour_rotation,'fontweight', vert_contour_fontweight);
end
h=colorbar;
caxis([0,20]);
colormap(jet);
set(h, 'Position', [.9125 .11 .0131 .2970])  % right, up, width, height
set(h,...
    'YTickLabel',{'0^oC', '5^oC','10^oC','15^oC','20^oC'},...
            'YTick',[0 5 10 15 20]);
h.Label.FontSize=7;
set(get(h,'title'),'string',' ');


fig=gcf;
fig.InvertHardcopy='off';

hold off

set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [500, 700]); %% [xsize ysize]
set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 500, 700]) 
jpgname=strcat('D:\OneDrive - 서울대학교\research\Master_course\for_publish\Figure\Fig25', '.jpg'); %% ~_year_month.jpg
saveas(gcf,'D:\OneDrive - 서울대학교\research\Master_course\for_publish\new_figure\fig_0004','jpg');