clc;close all;clear all;
name = 'ocean_avg_0060';
filename_suffix = '.nc';
filename = strcat(name,filename_suffix);

%% bathymetry
lon = ncread(filename,'lon_rho');
lat = ncread(filename,'lat_rho');
h = ncread(filename,'h');


figure
m_proj('mercator','lon',[98 284],'lat',[-20 65]);
m_grid('fontsize',15)
hold on
m_pcolor(lon,lat,h)
shading interp
colormap(jet);
h = colorbar;
m_gshhs_i('color','k')  
m_gshhs_i('patch',[.8 .8 .8]);  
set(h,'fontsize',15)
title(h,'m','fontsize',15)

%% SST
lon = ncread(filename,'lon_rho');
lat = ncread(filename,'lat_rho');
temp = ncread(filename,'temp');

figure
m_proj('mercator','lon',[98 284],'lat',[-20 65]);
m_grid('fontsize',15)
hold on
m_pcolor(lon,lat,temp(:,:,30))
shading interp
colormap(jet);
[C h] = contour(lon,lat,temp(:,:,30),0:3:35,'k','showtext','on');
h = colorbar;
m_gshhs_i('color','k')  
m_gshhs_i('patch',[.8 .8 .8]);  
set(h,'fontsize',15)
caxis([0 35]);
title(h,'m','fontsize',15)


%% plot temp
close all;clear all;clc
grid_file = 'roms_grid02_2.nc';
lon_rho = ncread(grid_file,'lon_rho');
lon_u = ncread(grid_file,'lon_u');
lon_v = ncread(grid_file,'lon_v');
lat_rho = ncread(grid_file,'lat_rho');
lat_u = ncread(grid_file,'lat_u');
lat_v = ncread(grid_file,'lat_v');
load 30year_mean.mat


mean_temp_7(find(mean_temp_7 == 0)) = NaN;
mean_temp_8(find(mean_temp_8 == 0)) = NaN;
mean_salt_7(find(mean_salt_7 == 0)) = NaN;
mean_salt_8(find(mean_salt_8 == 0)) = NaN;
mean_u_7(find(mean_u_7 == 0)) = NaN;
mean_u_8(find(mean_u_8 == 0)) = NaN;
mean_v_7(find(mean_v_7 == 0)) = NaN;
mean_v_8(find(mean_v_8 == 0)) = NaN;

coast=load('D:\matlab code\fvcom 磊丰贸府\coast_sin.dat');
% geoshow(coast(:,2),coast(:,1),'displaytype','polygon','FaceColor', [0.5 0.5 0.5]);
close all
var = mean_salt_8;
varname = 'SSS';
MM = 'Aug';

hold off
pcolor(lon_rho',lat_rho',var(:,:,20)');
shading interp
c = colorbar;
ylabel(c,'degree');
set(gca,'fontsize',15)
titlename = strcat('Averaged-',varname,'(1980-2009,',MM,')');
title(titlename,'fontsize',25);

% xlim([123 132])
% ylim([31 36])
caxis([0 35]);


hold on
% [C h] = contour(lon_rho',lat_rho',var(:,:,20)','k','showtext','on');
[C h] = contour(lon_rho',lat_rho',var(:,:,20)',0:3:35,'k','showtext','on');

geoshow(coast(:,2),coast(:,1),'color','black');
set(gcf,'Position',[200 100 650 350])

% saveas(gcf,titlename,'jpg');

saveas(gcf,strcat(titlename,'_zoom'),'jpg');


%% u,v (water)
close all;clear all;clc
grid_file = 'roms_grid02_2.nc';
lon_rho = ncread(grid_file,'lon_rho');
lon_u = ncread(grid_file,'lon_u');
lon_v = ncread(grid_file,'lon_v');
lat_rho = ncread(grid_file,'lat_rho');
lat_u = ncread(grid_file,'lat_u');
lat_v = ncread(grid_file,'lat_v');
load 30year_mean.mat

length_u = length(lon_u(1,:));
length_v = length(lat_v(:,1));
a = double(mean_u_8(:,:,20));
b = double(mean_v_8(:,:,20));

u(1,1) = a(1,1);
u(334,1) = a(333,1);
for i = 1:332
    for j = 1:333
        if (j ~= 1)
            u(j,i) = (a(j-1,i) + a(j,i))/2;
        end
    end
end

v(1,1) = b(1,1);
v(1,332) = b(1,331);
for j = 1:334
    for i = 1:331
        if (i ~= 1)
            v(j,i) = (b(j,i-1) + b(j,i))/2;
        end
    end
end

close all
interval = 2;
quiver(lon_rho(1:interval:end,1:interval:end),lat_rho(1:interval:end,1:interval:end),...
    u(1:interval:end,1:interval:end),v(1:interval:end,1:interval:end),2,'k')
coast=load('D:\matlab code\fvcom 磊丰贸府\coast_sin.dat');
hold on
geoshow(coast(:,2),coast(:,1),'color','black');
% xlim([min(lon_rho(:,1)) max(lon_rho(:,1))])
% ylim([min(lat_rho(1,:)) max(lat_rho(1,:))])
xlim([124 130])
ylim([30 38])
set(gca,'fontsize',15)
titlename = strcat('u,v velocity(1980-2009)-Aug');
title(titlename,'fontsize',20);

interval = 20;
quiver(lon_rho,lat_rho,u',v')












