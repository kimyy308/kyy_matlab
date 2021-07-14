clear all; clc; close all;
warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
%     addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
%     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\MOM\box_exp\Figure']));
elseif (linux==1)
    dropboxpath='/home01/kimyy/Dropbox';
%     addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/MOM/box_exp/Figure']));
end


% % % 
% % %  Figure 2. Meridional vertical mean temperature at 129°E
% % % (WOA 2013 (2005-2012 mean), upper panel ? summer, lower panel ? winter)

rname ='D:\need_to_presentation\WOA2013/T_1_7.nc';
temp2=ncread(rname,'T_MOD');
lon2=ncread(rname,'LON1240_1240');
lat2=ncread(rname,'LAT500_516');
level2=ncread(rname,'DEPTH1_25');
ilevel2 = -level2;

'reading data for figure 2 is completed'

% 
% Figure 2.1. July (upper panel)
%

filename=('D:\MEPL\Master_course\for_publish\Figure\Fig2_1.tif');
gifname=('Temperature at 130°E, JJA');
% pcolor(lat2(1:17),ilevel2(1:25),(reshape((temp2(1,1:17,1:25,6)+temp2(1,1:17,1:25,7)+temp2(1,1:17,1:25,8))/3,17,25))');
pcolor(lat2(1:17),ilevel2(1:25),squeeze((temp2(1,1:17,1:25,6)+temp2(1,1:17,1:25,7)+temp2(1,1:17,1:25,8))/3)');
hold on;
shading interp;
c=colorbar;
c.Label.String= 'Temp (°C)';
c.Label.FontSize= 12;
colormap jet;
% title(gifname,'fontsize',15);
xlabel('Latitude (°N)','fontsize',15);
ylabel('Depth (M)','fontsize',15);
zindex = 10;
caxis([0.2,20]);
c=colorbar;
c.Label.String= 'Temp (°C)';
c.Label.FontSize= 12;
% con=contour(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))','k', 'ShowText','on','LevelStep',5, 'LineWidth',2);
[C,h]=contour(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))','k','LevelStep',5, 'LineWidth',2);
t=clabel(C,h,'FontSize',15, 'LabelSpacing', 300)
% contourc(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))',4);
set(gcf,'PaperPosition',[0 0 30 18]);  %왼쪽 여백, 위쪽 여백, 가로길이, 세로 길이
set(gca,'fontsize',20);
set(gca, 'color', [0.8,0.8,0.8]);
fig=gcf;
fig.InvertHardcopy='off';
saveas(gcf,filename,'tiff');
hold off
'making figure 2(upper) is completed'

% 
% Figure 2.2. January (lower panel)
%

filename=('D:\MEPL\Master_course\for_publish\Figure\Fig2_2.tif');
gifname=('Temperature at 130°E, DJF');
% pcolor(lat2(1:17),ilevel2(1:25),(reshape((temp2(1,1:17,1:25,12)+temp2(1,1:17,1:25,1)+temp2(1,1:17,1:25,2))/3,17,25))');
pcolor(lat2(1:17),ilevel2(1:25),squeeze((temp2(1,1:17,1:25,12)+temp2(1,1:17,1:25,1)+temp2(1,1:17,1:25,2))/3)');
hold on;
shading interp;
c=colorbar;
c.Label.String= 'Temp (°C)';
c.Label.FontSize= 12;
colormap jet;
% title(gifname,'fontsize',15);
xlabel('Latitude (°N)','fontsize',15);
ylabel('Depth (M)','fontsize',15);
zindex = 10;

caxis([0.2,20]);
c=colorbar;
c.Label.String= 'Temp (°C)';
c.Label.FontSize= 12;
% con=contour(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))','k', 'ShowText','on','LevelStep',5, 'LineWidth',2);
[C,h]=contour(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,1),17,25))','k','LevelStep',5, 'LineWidth',2);
t=clabel(C,h,'FontSize',15, 'LabelSpacing', 300)
% contourc(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))',4);
set(gcf,'PaperPosition',[0 0 30 18]);  %왼쪽 여백, 위쪽 여백, 가로길이, 세로 길이
set(gca,'fontsize',20);
set(gca, 'color', [0.8,0.8,0.8]);
fig=gcf;
fig.InvertHardcopy='off';
saveas(gcf,filename,'tiff');
hold off
'making figure 2(lower) is completed'

% % % 
% % % Figure 3. Horizontal depth distribution in the study area the East/Japan Sea.
% % %

clear all; clc; close all;

expname='oman_restore_04_18';
readname='D:\need_to_presentation\';
readname2='\197912\ocean_gridinfo_197912.nc';
rname=strcat(readname,expname,readname2);
ht2=ncread(rname,'ht');
lon2=ncread(rname,'xt_ocean');
lat2=ncread(rname,'yt_ocean');

for i=1:120
    for j=1:110
        if (ht2(i,j)< -1.0e+5)
          ht2(i,j)=NaN;
        end
    end
end

'reading data for figure 3 is completed'

 readname='D:\MEPL\Master_course\for_publish\Figure';
 readname2='\Fig3.tif';
 filename=strcat(readname,readname2);
    pcolor(lon2,lat2,(ht2'));
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
    caxis([0,2000]);
    xlim([128 140]);
    ylim([34 45]);
    c=colorbar;
    c.Label.String= 'Depth (M)';
    c.Label.FontSize= 12;
    colormap jet;
    xlabel('Longitude (^oE)','fontsize',27);
    ylabel('Latitude (^oN)','fontsize',27);
    set(gcf,'PaperPosition',[0 0 20 15]);

    
% %  %     right edge
% %     set(gca,'fontsize',15); 
% %     ar=annotation(gcf,'rectangle',...
% %     [0.455 0.8275 0.00353378378378377 0.075],...
% %     'FaceColor','red');
% %  %     left edge
% %     set(gca,'fontsize',15); 
% %     ar=annotation(gcf,'rectangle',...
% %     [0.345 0.8275 0.00353378378378377 0.075],...
% %     'FaceColor','red');
% %  %    lower edge
% %     set(gca,'fontsize',15); 
% %     ar=annotation(gcf,'rectangle',...
% %     [0.345 0.8275 0.11 0.00353378378378377],...
% %     'FaceColor','red');
% %  %    upper edge
% %     set(gca,'fontsize',15); 
% %     ar=annotation(gcf,'rectangle',...
% %     [0.345 0.9 0.11 0.00353378378378377],...
% %     'FaceColor','red');
annotation('rectangle',[.3542 .83 .11 .07],'Color','white','LineWidth',3);

    tex=text(134.315,44.2,'Cooling Area');
    set(tex,'fontsize',15);
%     set(tex,'Color',[0.5 0.5 0.5]);
  hold off;
    fig=gcf;
    fig.InvertHardcopy='off';
    set(gca, 'color', [0.8,0.8,0.8]);
    saveas(gcf,filename,'tiff');
    
'making figure 3 is completed'




% % % 
% % % Figure 4. Horizontal mean Kinetic Energy & Temperature
% % %

clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
rdata_fig4;

%
% figure4_1, mean kinetic energy on each depth
%

 readname='D:\MEPL\Master_course\for_publish\Figure\';
 readname2='Fig4_1.tif';
 filename=strcat(readname,readname2);
 
 startDate = datenum('01-01-0001');
 endDate = datenum('12-31-0010');
 xData = linspace(startDate,endDate,120);
 plot(xData,kine(1,1:120),'k','LineWidth',1.5);
 datetick('x','yy')
 
 hold on
 plot(xData,kine(5,1:120),'b');
 plot(xData,kine(10,1:120),'y');
 plot(xData,kine(15,1:120),'r');
 plot(xData,kine(20,1:120),'g');
 ylim([0 7])
 legend('5m','45m','95m','145m','195m','Location','southoutside','Orientation','horizontal')
 hold off   
 xlabel('Years','fontsize',27);
 ylabel('K. E.','fontsize',27);
 set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
 set(gca,'fontsize',15);
    fig=gcf;
    fig.InvertHardcopy='off';
 saveas(gcf,filename,'tiff');
'making figure 4_1 is completed'

%
% figure 4_2, mean temperature on each depth
%

 readname='D:\MEPL\Master_course\for_publish\Figure\';
 readname2='Fig4_2.tif';
 filename=strcat(readname,readname2);
 
 startDate = datenum('01-01-0001');
 endDate = datenum('12-31-0010');
 xData = linspace(startDate,endDate,120);
 plot(xData,tempe(1,1:120),'k','LineWidth',1.5);
 datetick('x','yy')
 
 hold on
 plot(xData,tempe(5,1:120),'b');
 plot(xData,tempe(10,1:120),'y');
 plot(xData,tempe(15,1:120),'r');
 plot(xData,tempe(20,1:120),'g');
 ylim([0 21])
 legend('5m','45m','95m','145m','195m','Location','southoutside','Orientation','horizontal')
 hold off 
 xlabel('Years','fontsize',27);
 ylabel('Temp (^oC)','fontsize',27);
 set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
 set(gca,'fontsize',15);
    fig=gcf;
    fig.InvertHardcopy='off';
 saveas(gcf,filename,'tiff');
'making figure 4_2 is completed'




% % % 
% % % Figure 5. Surface Vector (Exp. 1)
% % %
clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
rdata_fig5;
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig5.tif';
filename=strcat(readname,readname2);
vector_common1;
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean(1:3:120,1:3:110,1,1))'*2,squeeze(v2mean(1:3:120,1:3:110,1,1))'*2,'k','AutoScale','off');
hold off;
vector_common2;
'making figure 5 is completed'



% % % 
% % % Figure 6. Surface Vector (Exp. 2) & vertical section of A & B
% % %
clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
expname='no_restore_04_21';

rdata_fig6_8_10;

%
% figure 6_1, Surface Vector (Exp. 2)
%
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig6_1.tif';
filename=strcat(readname,readname2);
vector_common1;
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean(1:3:120,1:3:110,1,1))'*2,squeeze(v2mean(1:3:120,1:3:110,1,1))'*2,'k','AutoScale','off');
ar=annotation(gcf,'rectangle',...
[0.465 0.3815 0.00353378378378377 0.07],...
'FaceColor','red');
tex=text(132.815,37.1,'A');
set(tex,'fontsize',20);
br=annotation(gcf,'rectangle',...
[0.225 0.55 0.07010055865921788 0.00353378378378377 ],...    
 'FaceColor','red');
tex=text(127.8,40.05,'B');
set(tex,'fontsize',20);
vector_common2;
hold off;
delete(ar);
delete(br); 
'making figure 6_1 is completed'

%
% figure 6_2, Vertical Velocity section of A & B (Exp. 2)
%

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig6_2.tif';
filename=strcat(readname,readname2);
fig6_2;

% % % 
% % % Figure 7. 100m Temperature (Exp. 3) & vector
% % %

clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
rdata_fig7;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig7.tif';
filename=strcat(readname,readname2);

pcolor(lon2,lat2,(temp2mean(:,:,10,1)'));
temp_pcolor_common1;
qv=  quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
squeeze(u2mean(1:3:120,1:3:110,1,1))'*2,squeeze(v2mean(1:3:120,1:3:110,1,1))'*2,'k','AutoScale','off');    
set(qv,'Color',[0 0 0]);  
hold off;
tex=text(135.9,35.2,'0.2m/s');
set(tex,'fontsize',10);
    fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
saveas(gcf,filename,'tiff');



% % % 
% % % Figure 8. Surface Vector (Exp. 3) & vertical section of A & B
% % %
clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
expname='restore_f_plane_07_01';
rdata_fig6_8_10;

%
% figure 8_1, Surface Vector (Exp. 3)
%
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig8_1.tif';
filename=strcat(readname,readname2);
vector_common1;
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean(1:3:120,1:3:110,1,1))'*2,squeeze(v2mean(1:3:120,1:3:110,1,1))'*2,'k','AutoScale','off');
ar=annotation(gcf,'rectangle',...
[0.465 0.3815 0.00353378378378377 0.07],...
'FaceColor','red');
tex=text(132.815,37.1,'A');
set(tex,'fontsize',20);
br=annotation(gcf,'rectangle',...
[0.225 0.55 0.07010055865921788 0.00353378378378377 ],...    
 'FaceColor','red');
tex=text(127.8,40.05,'B');
set(tex,'fontsize',20);
pcolindex=1
vector_common2;
hold off;
delete(ar);
delete(br); 
'making figure 8_1 is completed'

%
% figure 8_2, Vertical Velocity section of A & B (Exp. 3)
%

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig8_2.tif';
filename=strcat(readname,readname2);
fig6_2;   %%%same with 8_2



% % % 
% % % Figure 9. Isothermal surface (19.5^o) (Exp. 3)
% % %

clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
rdata_fig9;

view(-24,8);
hold on
surf(lon2(4:23),lat2(1:110),isod1(4:23,1:110)')
shading interp;
surf(lon2(4:23),lat2(1:110),isod5(4:23,1:110)','FaceColor',[0.8 0.8 0.8])
contour3(lon2(4:23),lat2(1:110),isod1(4:23,1:110)', 10, 'color', 'k', 'LineWidth',1)
surf_common_3d;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig9_1.tif';
filename=strcat(readname,readname2);
saveas(gcf,filename,'tiff');
hold off

view(-154,21);
hold on
surf(lon2(4:23),lat2(1:110),isod1(4:23,1:110)')
shading interp;
surf(lon2(4:23),lat2(1:110),isod5(4:23,1:110)','FaceColor',[0.8 0.8 0.8])
contour3(lon2(4:23),lat2(1:110),isod1(4:23,1:110)', 10, 'color', 'k', 'LineWidth',1)
surf_common_3d;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig9_2.tif';
filename=strcat(readname,readname2);
saveas(gcf,filename,'tiff');
hold off


view(-24,8);
hold on
surf(lon2(4:23),lat2(1:110),isod3(4:23,1:110)')
shading interp;
surf(lon2(4:23),lat2(1:110),isod6(4:23,1:110)','FaceColor',[0.8 0.8 0.8])
contour3(lon2(4:23),lat2(1:110),isod3(4:23,1:110)', 10, 'color', 'k', 'LineWidth',1)
surf_common_3d;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig9_1_cf.tif';
filename=strcat(readname,readname2);
saveas(gcf,filename,'tiff');
hold off

view(-154,21);
hold on
surf(lon2(4:23),lat2(1:110),isod3(4:23,1:110)')
shading interp;
surf(lon2(4:23),lat2(1:110),isod6(4:23,1:110)','FaceColor',[0.8 0.8 0.8])
contour3(lon2(4:23),lat2(1:110),isod3(4:23,1:110)', 10, 'color', 'k', 'LineWidth',1)
surf_common_3d;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig9_2_cf.tif';
filename=strcat(readname,readname2);
saveas(gcf,filename,'tiff');
hold off



% % % 
% % % Figure 10. Surface Vector (Exp. 4) 
% % %

clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
expname='oman_restore_04_18';
rdata_fig6_8_10;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig10.tif';
filename=strcat(readname,readname2);
vector_common1;
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean(1:3:120,1:3:110,1,1))'*2,squeeze(v2mean(1:3:120,1:3:110,1,1))'*2,'k','AutoScale','off');
vector_common2;
hold off;
delete(ar);
delete(br); 
'making figure 10 is completed'



% % % 
% % % Figure 11. Potential Vorticity (Exp. 2) 
% % %
clear all; clc; close all;
expname='no_restore_04_21';
rdata_fig11;
rdata_fig6_8_10;
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig11.tif';
filename=strcat(readname,readname2);
pot_vort_pcolor_common;


% % % 
% % % Figure 12. Upper Layer Thickness (Exp. 2) 
% % %
clear all; clc; close all;
expname='no_restore_04_21';
rdata_fig11;
rdata_fig6_8_10;
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig12.tif';
filename=strcat(readname,readname2);
h_pcolor_common;


% % % 
% % % Figure 13. Surface Relative Vorticity (Exp. 2) 
% % %
clear all; clc; close all;
expname='no_restore_04_21';
rdata_fig11;
rdata_fig6_8_10;
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig13.tif';
filename=strcat(readname,readname2);
vort_pcolor_common;
hold off;


% % % 
% % % Figure 14. Potential Vorticity (Exp. 3) 
% % %
clear all; clc; close all;
expname='restore_f_plane_07_01';
rdata_fig11;
rdata_fig6_8_10;
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig14.tif';
filename=strcat(readname,readname2);
pot_vort_pcolor_common;


% % % 
% % % Figure 15. Upper Layer Thickness (Exp. 3) 
% % %
clear all; clc; close all;
expname='restore_f_plane_07_01';
rdata_fig11;
rdata_fig6_8_10;
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig15.tif';
filename=strcat(readname,readname2);
h_pcolor_common;


% % % 
% % % Figure 16. Surface Relative Vorticity (Exp. 3) 
% % %
clear all; clc; close all;
expname='restore_f_plane_07_01';
rdata_fig11;
rdata_fig6_8_10;
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig16.tif';
filename=strcat(readname,readname2);
vort_pcolor_common;
hold off;



% % % Figure 19. Temp + Vec 1,2,3,4  (surf, 100m, 200m)
clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
expname1='no_restore_f_plane_06_29';
expname2='no_restore_04_21';
expname3='restore_f_plane_07_01';
expname4='oman_restore_04_18';
rdata_fig19;

%
% figure 8_1, Surface Vector (Exp. 3)
%
readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19.tif';
colordepthindex=1;
quivdepthindex=1;
filename=strcat(readname,readname2);
fig19_1;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19_2.tif';
colordepthindex=10;
quivdepthindex=10;
filename=strcat(readname,readname2);
fig19_1;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19_3.tif';
colordepthindex=20;
quivdepthindex=20;
filename=strcat(readname,readname2);
fig19_1;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19_4.tif';
colordepthindex=10;
quivdepthindex=1;
filename=strcat(readname,readname2);
fig19_1;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19_5.tif';
filename=strcat(readname,readname2);
fig19_2;
'making figure 19 is completed'

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19_6.tif';
filename=strcat(readname,readname2);
fig19_3;
'making figure 19 is completed'

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19_7.tif';
filename=strcat(readname,readname2);
fig19_7;
'making figure 19_7 is completed'

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19_7_2.tif';
filename=strcat(readname,readname2);
fig19_7_2;
'making figure 19_7_2 is completed'



readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig19_8.tif';
filename=strcat(readname,readname2);
fig19_8;
'making figure 19_8 is completed'

figname1='D:\MEPL\Master_course\for_publish\Figure\';
figname2='Fig19_9.tif';
expname1='no_restore_f_plane_06_29';
expname2='no_restore_04_21';
expname3='restore_f_plane_07_01';
expname4='oman_restore_04_18';
rdata_fig19;
rdata_fig19_9;
colordepthindex=10;
quivdepthindex=10;
figname=strcat(figname1,figname2);
fig19_9;
'making figure 19_9 is completed'


figname1='D:\MEPL\Master_course\for_publish\Figure\';
figname2='Fig20.tif';
expname1='no_restore_f_plane_06_29';
expname2='no_restore_04_21';
expname3='restore_f_plane_07_01';
expname4='oman_restore_04_18';
rdata_fig19;
rdata_fig19_9;
colordepthindex=10;
quivdepthindex=10;
figname=strcat(figname1,figname2);
fig20;
'making figure 20 is completed'




expname1='no_restore_f_plane_06_29';
expname2='no_restore_04_21';
expname3='restore_f_plane_07_01';
expname4='oman_restore_04_18';
rdata_fig21;

figname1='D:\MEPL\Master_course\for_publish\Figure\';
% figname2='Fig21.tif';
figname2='ns_section_x129_7.tif';
figname=strcat(figname1,figname2);
fig21;
'making figure 21 is completed'




expname1='no_restore_f_plane_06_29';
expname2='no_restore_04_21';
expname3='restore_f_plane_07_01';
expname4='oman_restore_04_18';
rdata_fig21;

figname1='D:\MEPL\Master_course\for_publish\Figure\';
% figname2='Fig22.tif';
figname2='we_setction_y36.tif';
figname=strcat(figname1,figname2);
fig22;
'making figure 22 is completed'




expname1='no_restore_f_plane_06_29';
expname2='no_restore_04_21';
expname3='restore_f_plane_07_01';
expname4='oman_restore_04_18';
rdata_fig21;

figname1='D:\MEPL\Master_course\for_publish\Figure\';
figname2='fig22_2.tif';
figname=strcat(figname1,figname2);
fig22_2;
'making figure 22_2 is completed'







figname1='D:\MEPL\Master_course\for_publish\Figure\';
figname2='Fig23.tif';
expname1='no_restore_f_plane_06_29';
expname2='no_restore_04_21';
expname3='restore_f_plane_07_01';
expname4='oman_restore_04_18';
rdata_fig23;
rdata_fig19_9;
colordepthindex=1;
quivdepthindex=1;
figname=strcat(figname1,figname2);
fig23;
'making figure 23 is completed'


figname1='D:\MEPL\Master_course\for_publish\Figure\';
figname2='Fig24.tif';
expname1='no_restore_f_plane_06_29';
expname2='no_restore_04_21';
expname3='restore_f_plane_07_01';
expname4='oman_restore_04_18';
rdata_fig23;
rdata_fig19_9;
colordepthindex=1;
quivdepthindex=1;
figname=strcat(figname1,figname2);
fig24;
'making figure 24 is completed'








rdata_fig_woa;

figname1='D:\MEPL\Master_course\for_publish\Figure\';
figname_summer2='woa_summer.tif';
figname_winter2='woa_winter.tif';
figname_woa2='woa_summer_winter.tif';
figname_summer=strcat(figname1,figname_summer2);
figname_winter=strcat(figname1,figname_winter2);
figname_woa=strcat(figname1,figname_woa2);
fig_woa;
'making figure_woa is completed'



rdata_fig_woa;

figname1='D:\MEPL\Master_course\for_publish\Figure\';
figname_summer2='woa_summer.tif';
figname_winter2='woa_winter.tif';
figname_woa2='woa_summer_winter_fnr.tif';
figname_summer=strcat(figname1,figname_summer2);
figname_winter=strcat(figname1,figname_winter2);
figname_woa=strcat(figname1,figname_woa2);
fig_woa_fnr;
'making figure_woa_fnr is completed'




%
% figure 8_2, Vertical Velocity section of A & B (Exp. 3)
%

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig8_2.tif';
filename=strcat(readname,readname2);
fig6_2;   %%%same with 8_2








% % % Figure 17_2. Vort + Vec 1,2,3,4 (surf)

clear all; clc; close all;
addpath(genpath('D:\MEPL\Master_course\for_publish'))
rdata_fig7;

readname='D:\MEPL\Master_course\for_publish\Figure\';
readname2='Fig7.tif';
filename=strcat(readname,readname2);

pcolor(lon2,lat2,(temp2mean(:,:,10,1)'));
temp_pcolor_common1;
qv=  quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
squeeze(u2mean(1:3:120,1:3:110,1,1))'*2,squeeze(v2mean(1:3:120,1:3:110,1,1))'*2,'k','AutoScale','off');    
set(qv,'Color',[0 0 0]);  
hold off;
tex=text(135.9,35.2,'0.2m/s');
set(tex,'fontsize',10);
    fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
saveas(gcf,filename,'tiff');