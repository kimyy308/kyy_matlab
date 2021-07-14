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


figname1='D:\OneDrive - 서울대학교\research\Master_course\for_publish\Figure\';
figname2='schematics_surf_beta.tif';
expname2='no_restore_04_21';
expname4='oman_restore_04_18';
windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop
%     addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
%     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\HYCOM\Grid_kyy']));
%     addpath(genpath([dropboxpath '\source\matlab\Reanalysis\HYCOM\Analysis\expt19_1']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
%     addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/HYCOM/Grid_kyy']));
%     addpath(genpath([dropboxpath '/source/matlab/Reanalysis/HYCOM/Analysis/expt19_1']));
end

bwr_map;


readname='E:\Data\Model\MOM\Box_eastsea\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname2,readname2);
u2_2=ncread(rname,'u');
v2_2=ncread(rname,'v');
temp2_2=ncread(rname,'temp');
lon2=ncread(rname,'xt_ocean');
lat2=ncread(rname,'yt_ocean');
level2=ncread(rname,'st_ocean');
ilevel2 = -level2;

'reading data for figure_schematics is completed'

readname='E:\Data\Model\MOM\Box_eastsea\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname4,readname2);
u2_4=ncread(rname,'u');
v2_4=ncread(rname,'v');
temp2_4=ncread(rname,'temp');

'reading data for figure_schematics is completed'

readname='E:\Data\Model\MOM\Box_eastsea\';
readname2='\197912\ocean_gridinfo_197912.nc';
rname=strcat(readname,expname2,readname2);
ht2=ncread(rname,'ht');
ht3=ht2;

[lon, lat] = meshgrid(lon2',lat2);

for i=1:120
    for j=1:110
        for k=1:30
            for l=1:12
                if (u2_2(i,j,k,l)< -1.0e+5)
                    u2_2(i,j,k,l)=NaN;
                    u2_4(i,j,k,l)=NaN;
                end
                if (v2_2(i,j,k,l)< -1.0e+5)
                    v2_2(i,j,k,l)=NaN;
                    v2_4(i,j,k,l)=NaN;
                end
                if (temp2_2(i,j,k,l)< -1.0e+5)
                    temp2_2(i,j,k,l)=NaN;
                    temp2_4(i,j,k,l)=NaN;
                end
            end
        end
        if (ht3(i,j)< -1.0e+5)
          ht3(i,j)=NaN;
        else
          ht3(i,j)=0;
        end
    end
end

% % for i=1:12
% %     u2_1(82:84,14:17,1:30,i)=0.2;
% %     v2_1(82:84,14:17,1:30,i)=0.00001;
% %     u2_2(82:84,14:17,1:30,i)=0.2;
% %     v2_2(82:84,14:17,1:30,i)=0.00001;
% %     u2_3(82:84,14:17,1:30,i)=0.2;
% %     v2_3(82:84,14:17,1:30,i)=0.00001;
% %     u2_4(82:84,14:17,1:30,i)=0.2;
% %     v2_4(82:84,14:17,1:30,i)=0.00001;
% % end
% 
% for i=1:12
%     u2_1(29:31,29:31,1:30,i)=0.2;
%     v2_1(29:31,29:31,1:30,i)=0.00001;
%     u2_2(29:31,29:31,1:30,i)=0.2;
%     v2_2(29:31,29:31,1:30,i)=0.00001;
%     u2_3(29:31,29:31,1:30,i)=0.2;
%     v2_3(29:31,29:31,1:30,i)=0.00001;
%     u2_4(29:31,29:31,1:30,i)=0.2;
%     v2_4(29:31,29:31,1:30,i)=0.00001;
% end
% 
% 
for i=1:120
  for j=1:110
    for k=1:30
%       u2mean_1(i,j,k,1)=mean(u2_1(i,j,k,1:12));
%       v2mean_1(i,j,k,1)=mean(v2_1(i,j,k,1:12));
      u2mean_2(i,j,k,1)=mean(u2_2(i,j,k,1:12));
      v2mean_2(i,j,k,1)=mean(v2_2(i,j,k,1:12));
%       u2mean_3(i,j,k,1)=mean(u2_3(i,j,k,1:12));
%       v2mean_3(i,j,k,1)=mean(v2_3(i,j,k,1:12));
      u2mean_4(i,j,k,1)=mean(u2_4(i,j,k,1:12));
      v2mean_4(i,j,k,1)=mean(v2_4(i,j,k,1:12));
%       temp2mean_1(i,j,k,1)=mean(temp2_1(i,j,k,1:12));
      temp2mean_2(i,j,k,1)=mean(temp2_2(i,j,k,1:12));
%       temp2mean_3(i,j,k,1)=mean(temp2_3(i,j,k,1:12));
      temp2mean_4(i,j,k,1)=mean(temp2_4(i,j,k,1:12));
    end
  end
end

% maxekwc_1=max(max(v2mean_1(1:15,21:80,1,1)))
% maxekwc_2=max(max(v2mean_2(1:15,21:80,1,1)))
% maxekwc_3=max(max(v2mean_3(1:15,21:80,1,1)))
% maxekwc_4=max(max(v2mean_4(1:15,21:80,1,1)))
% [i1, i2, i3, i4] = ind2sub(size(v2mean_1), find(v2mean_1==maxekwc_1));
% [i1, i2, i3, i4] = ind2sub(size(v2mean_2), find(v2mean_2==maxekwc_2));
% [i1, i2, i3, i4] = ind2sub(size(v2mean_3), find(v2mean_3==maxekwc_3));
% [i1, i2, i3, i4] = ind2sub(size(v2mean_4), find(v2mean_4==maxekwc_4));
% 'reading data for figure 19 is completed'

readname2='\isodepth.nc';
rname=strcat(readname,expname2,readname2);
isod5_2=ncread(rname,'isodepth');
vort5_2=ncread(rname,'sfc_rel_vort');
f5_2=ncread(rname,'f');
pot_vort5_2=ncread(rname,'pot_vort');
lon5_2=ncread(rname,'XT_OCEAN');
lat5_2=ncread(rname,'YT_OCEAN');
vort5mean_2=mean(vort5_2(:,:,109:120),3);
isod5mean_2=mean(isod5_2(:,:,109:120),3);
pot_vort5mean_2=mean(pot_vort5_2(:,:,109:120),3);

rname=strcat(readname,expname4,readname2);

isod5_4=ncread(rname,'isodepth');
vort5_4=ncread(rname,'sfc_rel_vort');
f5_4=ncread(rname,'f');
pot_vort5_4=ncread(rname,'pot_vort');
lon5=ncread(rname,'XT_OCEAN');
lat5=ncread(rname,'YT_OCEAN');
vort5mean_4=mean(vort5_4(:,:,109:120),3);
isod5mean_4=mean(isod5_4(:,:,109:120),3);
pot_vort5mean_4=mean(pot_vort5_4(:,:,109:120),3);

'reading data for figure_schematics is completed'



[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
pcolindex=1;
lon_min_ind = 4;  %% 128.3 E
lon_max_ind = 23; %% 133 E
lat_min_ind = 25; %% 36.5 N
lat_max_ind = 90; %% 43 N
skipind = 3;
vecwidth = 1.2;
xlength=lon_max_ind-lon_min_ind+1;
ylength=lat_max_ind-lat_min_ind+1;
texx = 130.6;
texy = 36.7;
texsize =15;
tex2size =25;
tex2x = 132.0;
tex2y = 37.0;
texFontWeight= 'bold';   %%'normal or bold'
colordepthindex=1;
quivdepthindex=1;


pcol=pcolor(lon5(lon_min_ind:lon_max_ind),lat5(lat_min_ind:lat_max_ind),(isod5mean_2(lon_min_ind:lon_max_ind,lat_min_ind:lat_max_ind,1))');
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
colormap(gca,jet);
set(gca,'fontsize',13);

axis equal;
axis tight;
caxis(gca, [0,200]);
quiver(lon(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind),lat(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind), ...
  squeeze(u2mean_2(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  squeeze(v2mean_2(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  'k','AutoScale','off', 'LineWidth',vecwidth);

xlim([min(lon2(lon_min_ind:lon_max_ind)) max(lon2(lon_min_ind:lon_max_ind))]);
ylim([min(lat2(lat_min_ind:lat_max_ind)) max(lat2(lat_min_ind:lat_max_ind))]);
xdivvar=4.0;
ydivvar=6.0;
lontick_min_ind=10; %%128.5;
lontick_max_ind=50; %%132.5;
lattick_min_ind=30;  %%37.0;
lattick_max_ind=90;  %%43.0;
xtickvar=lon2(lontick_min_ind):(lon2(lontick_max_ind)-lon2(lontick_min_ind))/xdivvar:lon2(lontick_max_ind);
ytickvar=lat2(lattick_min_ind):(lat2(lattick_max_ind)-lat2(lattick_min_ind))/ydivvar:lat2(lattick_max_ind);
set(gca,...
    'XTickLabel',{'129^oE','','131^oE','','133^oE'},...
    'XTick',xtickvar);
set(gca,...
    'YTickLabel',{'37^oN','38^oN','39^oN','40^oN','41^oN', '42^oN', '43^oN'},...
        'YTick',ytickvar);
set(gca, 'TickDir','out');
set(gca,'fontsize',20);
set(gca, 'color', [0.8,0.8,0.8]);

% tex2=text(tex2x,tex2y,'(a)');
% tex=text(texx,texy,'0.2m/s');
% set(tex,'fontsize',texsize, 'FontWeight',texFontWeight);
% set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight);
% fig=gcf;

% set(gcf, 'renderer', 'zbuffer')
% set(pcol,'ZData',0.1+zeros(size(v2mean_2(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))))  % Move the surface plot to Z = -1

hold off;




%% vertical temp
figure; 

latlimind=90;
[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
[lat3,ilevel3] = meshgrid(ilevel2(1:20),lat2(1:latlimind));
lonind = 17; %% 129E
zlimind=20; %% ~ 200m
skipind=5;
tex2x=34.55;
tex2y=-8.5;
vecwidth=1.2;

pcolindex=1;

pcolor(lat2(lat_min_ind:lat_max_ind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_2(4,lat_min_ind:lat_max_ind,1:zlimind),lat_max_ind-lat_min_ind+1,zlimind))');
    shading interp;
%     xlim([min(lat2(1:latlimind)) max(lat2(1:latlimind))]);
%     set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
%         'YTick', [-8, -6, -4, -2, 0])
%     set(gca,'xTickLabel',{'35^oN', '36^oN', '37^oN', '38^oN', '39^oN', '40^oN', '41^oN', '42^oN'}, ...
%         'XTick', [35, 36, 37, 38, 39, 40, 41, 42])
    set(gca, 'color', [0.8,0.8,0.8]);
    colormap(gca,bwrmap);
    set(gcf,'PaperPosition',[0 0 30 12]);
    set(gca, 'TickDir','out');
    set(gca,'fontsize',20);
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);

hold on;
% quiver(ilevel3(1:skipind:latlimind,1:zlimind),lat3(1:skipind:latlimind,1:zlimind)/zlimind, ...
%     (reshape(v2mean_2(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*3, ...
%     (reshape(w2mean_2(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*2, ...
%     '-k','filled','AutoScale','off', 'LineWidth',vecwidth);tex2=text(tex2x,tex2y,'0.2m/s');
% set(tex2,'fontsize',20);
hold off;










% % vert_we section
figure;

lonmaxind=23;
lonminind=4;
[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
[lon3,ilevel3] = meshgrid(ilevel2(1:20),lon2(lonminind:lonmaxind));
latind = 25; %% 20 : 36N
zlimind=20; %% ~ 200m
skipind=5;
vecwidth=1.2;
xlength=lonmaxind-lonminind+1;
level_c =-0.5:0.05:0.5;
pcolindex=1;
contwidth=2;
confontsize=20;

pcolor(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_2(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))');
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);

% hold on;
% [C,h]=contour(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(v2mean_2(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))',level_c,'k','linewidth',contwidth);
% clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
% hold off;

%     hold on;
    shading interp;
%     xlim([min(lon2(lonminind:lonmaxind)) max(lon2(lonminind:lonmaxind))]);
%     set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
%         'YTick', [-160, -120, -80, -40, 0])
%     set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
%         'YTick', [-8, -6, -4, -2, 0])
%     set(gca,'xTickLabel',{'128.3^oE', '128.8^oE', '129.3^oE', '129.8^oE', '130.3^oE'}, ...
%     set(gca,'xTickLabel',{'128.5^oE','129.0^oE','129.5^oE','130.0^oE'}, ...
%         'XTick', (128.5:0.5:130.0))
    set(gca, 'color', [0.8,0.8,0.8]);
%     load jet_mod;
%     colormap(gca,jet_mod);
    colormap(gca,bwrmap);
    set(gcf,'PaperPosition',[0 0 28 12]);
%     axis equal;
    set(gca,'fontsize',20);
    set(gca, 'TickDir','out');