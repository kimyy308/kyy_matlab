% %  NIFS climatology temperature W-E section in Feb and in Aug at 102 line



clear all; clc; close all;
warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\seawater_ver3_2']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Model\MOM\box_exp\Figure']));
elseif (linux==1)
    dropboxpath='/home01/kimyy/Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\seawater_ver3_2']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Model/MOM/box_exp/Figure']));
end

load 'D:\OneDrive - 서울대학교\research\Master_course\KODC\KODC_rawvar.mat'  %%for snu desktop
figdir = 'D:\MEPL\Master_course\for_publish\Figure\';

latsize=length(unique(lat))
lonsize=length(unique(lon))
stddepth=[0 10 20 30 50 75 100 125 150 200 250 300 400 500];
tempstddepth=[-5 0 10 20 30 50 75 100 125 150 200 250 300 400 500 600];

i=12;
% % % extract 102-06 station data from the all data
temp_10206=temp(find(line_st==10206));
salt_10206=salt(find(line_st==10206));
lon_10206=lon(find(line_st==10206));
lat_10206=lat(find(line_st==10206));
depth_10206=depth(find(line_st==10206));
year_10206=year(find(line_st==10206));
month_10206=month(find(line_st==10206));

j=1;
for tempyear=2001:2015
ind_year=find(year_10206==tempyear);
y_month=month_10206(ind_year);
y_temp=temp_10206(ind_year);
y_salt=salt_10206(ind_year);
y_lat=lat_10206(ind_year);
y_lon=lon_10206(ind_year);
y_depth=depth_10206(ind_year);
    for tempmonth=1:12
        ind_month=find(y_month==tempmonth);
        if (isfinite(ind_month))
            if (i<12)
                i=i+1;
            else
                i=1;
            end
%             m_time(i) = mjuliandate(tempyear,tempmonth,15);
            m_temp=y_temp(ind_month);
            m_salt=y_salt(ind_month);
            m_lat=y_lat(ind_month);
            m_lon=y_lon(ind_month);
            m_depth=y_depth(ind_month);
            if (length(m_temp)>13)
                vtemp_10206(1:14,j)=m_temp(1:14);
                vsalt_10206(1:14,j)=m_salt(1:14);
                vdepth_10206(1:14,j)=m_depth(1:14);
                vlat_10206(1:14,j)=m_lat(1:14);
                vlon_10206(1:14,j)=m_lon(1:14);
%                 vtime_10206(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10206(j)=tempyear;
                vmonth_10206(j)=tempmonth;
                j=j+1;
            elseif(length(m_temp)>0)
                m_temp(length(m_temp)+1:14)=NaN;
                m_salt(length(m_salt)+1:14)=NaN;
                m_lat(length(m_lat)+1:14)=NaN;
                m_lon(length(m_lon)+1:14)=NaN;
                m_depth(length(m_depth)+1:14)=NaN;
                
                vtemp_10206(1:14,j)=m_temp;
                vsalt_10206(1:14,j)=m_salt;
                vdepth_10206(1:14,j)=m_depth;
                vlat_10206(1:14,j)=m_lat;
                vlon_10206(1:14,j)=m_lon;
%                 vtime_10206(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10206(j)=tempyear;
                vmonth_10206(j)=tempmonth;
                j=j+1;
            else
                m_temp(1:14)=NaN;
                m_salt(1:14)=NaN;
                m_lat(1:14)=NaN;
                m_lon(1:14)=NaN;
                m_depth(1:14)=NaN;
                
                vtemp_10206(1:14,j)=m_temp;
                vsalt_10206(1:14,j)=m_salt;
                vdepth_10206(1:14,j)=m_depth;
                vlat_10206(1:14,j)=m_lat;
                vlon_10206(1:14,j)=m_lon;
%                 vtime_10206(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10206(j)=tempyear;
                vmonth_10206(j)=tempmonth;
                j=j+1;
            end
        end
    end
end
j=1;

temp_10207=temp(find(line_st==10207));
salt_10207=salt(find(line_st==10207));
lon_10207=lon(find(line_st==10207));
lat_10207=lat(find(line_st==10207));
depth_10207=depth(find(line_st==10207));
year_10207=year(find(line_st==10207));
month_10207=month(find(line_st==10207));

for tempyear=2001:2015
ind_year=find(year_10207==tempyear);
y_month=month_10207(ind_year);
y_temp=temp_10207(ind_year);
y_salt=salt_10207(ind_year);
y_lat=lat_10207(ind_year);
y_lon=lon_10207(ind_year);
y_depth=depth_10207(ind_year);
    for tempmonth=1:12
        ind_month=find(y_month==tempmonth);
        if (isfinite(ind_month))
            if (i<12)
                i=i+1;
            else
                i=1;
            end
%             m_time(i) = mjuliandate(tempyear,tempmonth,15);
            m_temp=y_temp(ind_month);
            m_salt=y_salt(ind_month);
            m_lat=y_lat(ind_month);
            m_lon=y_lon(ind_month);
            m_depth=y_depth(ind_month);
            if (length(m_temp)>13)
                vtemp_10207(1:14,j)=m_temp(1:14);
                vsalt_10207(1:14,j)=m_salt(1:14);
                vdepth_10207(1:14,j)=m_depth(1:14);
                vlat_10207(1:14,j)=m_lat(1:14);
                vlon_10207(1:14,j)=m_lon(1:14);
%                 vtime_10207(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10207(j)=tempyear;
                vmonth_10207(j)=tempmonth;
                j=j+1;
            elseif(length(m_temp)>0)
                m_temp(length(m_temp)+1:14)=NaN;
                m_salt(length(m_salt)+1:14)=NaN;
                m_lat(length(m_lat)+1:14)=NaN;
                m_lon(length(m_lon)+1:14)=NaN;
                m_depth(length(m_depth)+1:14)=NaN;
                
                vtemp_10207(1:14,j)=m_temp;
                vsalt_10207(1:14,j)=m_salt;
                vdepth_10207(1:14,j)=m_depth;
                vlat_10207(1:14,j)=m_lat;
                vlon_10207(1:14,j)=m_lon;
%                 vtime_10207(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10207(j)=tempyear;
                vmonth_10207(j)=tempmonth;
                j=j+1;
            else
                m_temp(1:14)=NaN;
                m_salt(1:14)=NaN;
                m_lat(1:14)=NaN;
                m_lon(1:14)=NaN;
                m_depth(1:14)=NaN;
                
                vtemp_10207(1:14,j)=m_temp;
                vsalt_10207(1:14,j)=m_salt;
                vdepth_10207(1:14,j)=m_depth;
                vlat_10207(1:14,j)=m_lat;
                vlon_10207(1:14,j)=m_lon;
%                 vtime_10207(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10207(j)=tempyear;
                vmonth_10207(j)=tempmonth;
                j=j+1;
            end
        end
    end
end
j=1;

temp_10208=temp(find(line_st==10208));
salt_10208=salt(find(line_st==10208));
lon_10208=lon(find(line_st==10208));
lat_10208=lat(find(line_st==10208));
depth_10208=depth(find(line_st==10208));
year_10208=year(find(line_st==10208));
month_10208=month(find(line_st==10208));

for tempyear=2001:2015
ind_year=find(year_10208==tempyear);
y_month=month_10208(ind_year);
y_temp=temp_10208(ind_year);
y_salt=salt_10208(ind_year);
y_lat=lat_10208(ind_year);
y_lon=lon_10208(ind_year);
y_depth=depth_10208(ind_year);
    for tempmonth=1:12
        ind_month=find(y_month==tempmonth);
        if (isfinite(ind_month))
            if (i<12)
                i=i+1;
            else
                i=1;
            end
%             m_time(i) = mjuliandate(tempyear,tempmonth,15);
            m_temp=y_temp(ind_month);
            m_salt=y_salt(ind_month);
            m_lat=y_lat(ind_month);
            m_lon=y_lon(ind_month);
            m_depth=y_depth(ind_month);
            if (length(m_temp)>13)
                vtemp_10208(1:14,j)=m_temp(1:14);
                vsalt_10208(1:14,j)=m_salt(1:14);
                vdepth_10208(1:14,j)=m_depth(1:14);
                vlat_10208(1:14,j)=m_lat(1:14);
                vlon_10208(1:14,j)=m_lon(1:14);
%                 vtime_10208(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10208(j)=tempyear;
                vmonth_10208(j)=tempmonth;
                j=j+1;
            elseif(length(m_temp)>0)
                m_temp(length(m_temp)+1:14)=NaN;
                m_salt(length(m_salt)+1:14)=NaN;
                m_lat(length(m_lat)+1:14)=NaN;
                m_lon(length(m_lon)+1:14)=NaN;
                m_depth(length(m_depth)+1:14)=NaN;
                
                vtemp_10208(1:14,j)=m_temp;
                vsalt_10208(1:14,j)=m_salt;
                vdepth_10208(1:14,j)=m_depth;
                vlat_10208(1:14,j)=m_lat;
                vlon_10208(1:14,j)=m_lon;
%                 vtime_10208(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10208(j)=tempyear;
                vmonth_10208(j)=tempmonth;
                j=j+1;
            else
                m_temp(1:14)=NaN;
                m_salt(1:14)=NaN;
                m_lat(1:14)=NaN;
                m_lon(1:14)=NaN;
                m_depth(1:14)=NaN;
                
                vtemp_10208(1:14,j)=m_temp;
                vsalt_10208(1:14,j)=m_salt;
                vdepth_10208(1:14,j)=m_depth;
                vlat_10208(1:14,j)=m_lat;
                vlon_10208(1:14,j)=m_lon;
%                 vtime_10208(j)= mjuliandate(tempyear,tempmonth,15);
                vyear_10208(j)=tempyear;
                vmonth_10208(j)=tempmonth;
                j=j+1;
            end
        end
    end
end
j=1;


temps(1:14,1,1:88)=vtemp_10206(1:14,1:88);
temps(1:14,2,1:88)=vtemp_10207(1:14,1:88);
temps(1:14,3,1:88)=vtemp_10208(1:14,1:88);
depths(1:14,1,1:88)=vdepth_10206(1:14,1:88);
depths(1:14,2,1:88)=vdepth_10207(1:14,1:88);
depths(1:14,3,1:88)=vdepth_10208(1:14,1:88);
lats(1:14,1,1:88)=vlat_10206(1:14,1:88);
lats(1:14,2,1:88)=vlat_10207(1:14,1:88);
lats(1:14,3,1:88)=vlat_10208(1:14,1:88);
lons(1:14,1,1:88)=vlon_10206(1:14,1:88);
lons(1:14,2,1:88)=vlon_10207(1:14,1:88);
lons(1:14,3,1:88)=vlon_10208(1:14,1:88);

temps2=temps;
temps2(:,:,89:90)=NaN;

tsz=90;
clim_temps2=reshape(temps2,14,3,6,tsz/6);
clima_temps2=squeeze(nanmean(clim_temps2,4)); 

% % % draw circle and mask grid which distance from center of circle is
% % % closer than 500m
% % % this uses haversine formula to calculate the circle distance between
% % % two points (shortest distance over the earth's surface)
R = 6373.0 *1000  %%%% earth radius(m)
dist(1:14,1)= 0;
for i=1:2
clat=deg2rad(lats(1,i,1));
clon=deg2rad(lons(1,i,1));
tlat=deg2rad(lats(1,i+1,1));
tlon=deg2rad(lons(1,i+1,1));
dlon=tlon-clon;
dlat=tlat-clat;
a=sin(dlat / 2)^2 + cos(clat) * cos(tlat) * sin(dlon / 2)^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));
dist(1:14,i+1)= R * c;
end
% cal distance
distfromstart(:,1)=dist(:,1);
dist_gvel_from_start(:,1)=(dist(:,2)-dist(:,1))/2;
for i=1:2
distfromstart(:,i+1)=distfromstart(:,i)+dist(:,i+1);
end
for i=1:1
dist_gvel_from_start(:,i+1)=distfromstart(:,i+1)+(distfromstart(:,i+2)-distfromstart(:,i+1))/2;
end


%% raw col contour 2,4,8,10
lonmaxind=3;
lonminind=1;
lon2=squeeze(lons(1,:,1));
contwidth=2;
confontsize=20;
level_c =5:5:25;
tex2size =25;
tex2x = 129.83;
tex2y = -8.5;
texFontWeight= 'bold';   %%'normal or bold'
figname = strcat(figdir,'nifs_temp_2_4_8_10.tif');


sbp1=subplot(2,1,1);
pcolor(lon2,-stddepth(1:10)/20,clima_temps2(1:10,:,1));

shading interp;
xlim([min(lon2(lonminind:lonmaxind)) max(lon2(lonminind:lonmaxind))]);
set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
    'YTick', [-8, -6, -4, -2, 0])
set(gca,'xTickLabel',{'129.8^oE','130.0^oE','130.2^oE','130.4^oE','130.6^oE','130.8^oE'}, ...
    'XTick', (129.8:0.2:130.8))
set(gca, 'color', [0.8,0.8,0.8]);
colormap(gca,jet);
set(gcf,'PaperPosition',[0 0 28 12]);
set(gca,'fontsize',20);
set(gca, 'TickDir','out');
    
caxis(gca, [0,25]);
hold on;
[C,h]=contour(lon2(lonminind:lonmaxind),-stddepth(1:10)/20,clima_temps2(1:10,:,1),level_c,'k','linewidth',contwidth);
clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',90,'fontweight','bold');
tex2=text(tex2x,tex2y,'February');
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight,'Color','w');
hold off;

sbp1=subplot(2,1,2);
pcolor(lon2,-stddepth(1:10)/20,clima_temps2(1:10,:,4));
shading interp;
xlim([min(lon2(lonminind:lonmaxind)) max(lon2(lonminind:lonmaxind))]);
set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
    'YTick', [-8, -6, -4, -2, 0])
set(gca,'xTickLabel',{'129.8^oE','130.0^oE','130.2^oE','130.4^oE','130.6^oE','130.8^oE'}, ...
    'XTick', (129.8:0.2:130.8))
set(gca, 'color', [0.8,0.8,0.8]);
colormap(gca,jet);
set(gcf,'PaperPosition',[0 0 28 12]);
set(gca,'fontsize',20);
set(gca, 'TickDir','out');
caxis(gca, [0,25]);
hold on;
[C,h]=contour(lon2(lonminind:lonmaxind),-stddepth(1:10)/20,clima_temps2(1:10,:,4),level_c,'k','linewidth',contwidth);
clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',90,'fontweight','bold');
tex2=text(tex2x,tex2y,'August');
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight,'Color','w');
hold off;

h=colorbar;
caxis([0,25]);
colormap(jet);
set(h, 'Position', [.9150 .122 .0131 .8050])  % right, up, width, height
set(h,...
    'YTickLabel',{'2^oC', '5^oC','8^oC','11^oC','14^oC','17^oC','20^oC','23^oC'},...
            'YTick',[2 5 8 11 14 17 20 23]);
h.Label.FontSize=7;
set(get(h,'title'),'string',' ');

fig=gcf;
set(gcf,'PaperPosition',[0 0 24 15]);
fig.InvertHardcopy='off';

saveas(gcf,'D:\OneDrive - 서울대학교\research\Master_course\for_publish\new_figure\fig_0002','jpg');


