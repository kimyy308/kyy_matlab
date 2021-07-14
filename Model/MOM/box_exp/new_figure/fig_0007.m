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





%% winter vertical temperature

% vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
% status=plot_HYCOM_monthly_vertical_data(testname, vert_temp_file, inputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
%                                 [0 25], 0:5:25, 'vert_temp', 'fig_param_HYCOM_kyy_EJS');  %% 35.5N
testname='GLBb_0_08_expt19_1'   % % need to change
tempyear = [0110]; % % put year which you want to plot [year year ...]
tempmonth = 2;
var = 'vert_temp';
section = [129.875 129.875 34.875 38.875 -200 0];
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

subplot(2,1,1);
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
set(gca,'xTickLabel',{'35^oN','36^oN','37^oN','38^oN'}, ...
    'XTick', (35:1:38))
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
caxis([0,25]);
colormap(jet);
% set(h, 'Position', [.9125 .11 .0131 .2970])  % right, up, width, height
set(h, 'Position', [.9150 .122 .0131 .8050])  % right, up, width, height
set(h,...
    'YTickLabel',{'0^oC', '5^oC','10^oC','15^oC','20^oC','25^oC'},...
            'YTick',[0 5 10 15 20 25]);
h.Label.FontSize=7;
set(get(h,'title'),'string',' ');

tex2size =18;
tex2x = 38;
tex2y = -170.5;
texFontWeight= 'bold';   %%'normal or bold'
tex2=text(tex2x,tex2y,'February');
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight,'Color','w');

hold off




%% summer vertical temperature

% vert_temp_file =[figdir, 'EJS\', 'vert_temp\', 'EJS_vert_temp_35_5']; % % ~/SUV_EJS_testname_year_month.jpg
% status=plot_HYCOM_monthly_vertical_data(testname, vert_temp_file, inputdir, [129 131 35.5 35.5 -200 0], tempyear, month, ...
%                                 [0 25], 0:5:25, 'vert_temp', 'fig_param_HYCOM_kyy_EJS');  %% 35.5N
testname='GLBb_0_08_expt19_1'   % % need to change
tempyear = [0110]; % % put year which you want to plot [year year ...]
tempmonth = 8;
var = 'vert_temp';
section = [129.875 129.875 34.875 38.875 -200 0];
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

subplot(2,1,2);
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
set(gca,'xTickLabel',{'35^oN','36^oN','37^oN','38^oN'}, ...
    'XTick', (35:1:38))
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
caxis([0,25]);
colormap(jet);
% set(h, 'Position', [.9125 .11 .0131 .2970])  % right, up, width, height
set(h, 'Position', [.9150 .122 .0131 .8050])  % right, up, width, height
set(h,...
    'YTickLabel',{'0^oC', '5^oC','10^oC','15^oC','20^oC','25^oC'},...
            'YTick',[0 5 10 15 20 25]);
h.Label.FontSize=7;
set(get(h,'title'),'string',' ');


fig=gcf;
fig.InvertHardcopy='off';


tex2size =18;
tex2x = 38;
tex2y = -170.5;
texFontWeight= 'bold';   %%'normal or bold'
tex2=text(tex2x,tex2y,'August');
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight,'Color','w');

hold off



set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [500, 700]); %% [xsize ysize]
fig=gcf;
set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 500, 500]) ;
saveas(gcf,'D:\OneDrive - 서울대학교\research\Master_course\for_publish\new_figure\fig_0007','jpg');