clc;clear all;close all

%openfile

% filename='eastsea.nc';
% filename='..\Eastsea\grid\grid_etopo1_Eastsea.nc';
filename='F:\Nesting\make_grid\nwp2_eastsea.nc';
% filename='f:\Nesting\input\roms_grid_combine2.nc';

nc=netcdf(filename,'nowrite');
h=nc{'h'}(:);
maskr=nc{'mask_rho'}(:);
Lonr=nc{'lon_rho'}(:);
Latr=nc{'lat_rho'}(:);

figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.9]);
set(gca,'Position',[0.13 0.06 0.8 0.82]);
hold on;
  themask=ones(size(maskr));
  themask(maskr==0)=NaN; 
  domaxis=[min(min(Lonr)) max(max(Lonr)) min(min(Latr)) max(max(Latr))];
%   domaxis=[112  165 20 52 ];
%         domaxis=[116 128 33 42]; %%- È²ÇØ
  colaxis=[min(min(h)) max(max(h))];
% set(gca, 'box','on','linewidth',1.5,'layer','top')
% xlabel('Longitude(¢ªE)');ylabel('Latitude(¢ªN)')
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 20])
 m_proj('mercator',...
            'lon',[domaxis(1) domaxis(2)],...
            'lat',[domaxis(3) domaxis(4)]);
        
%  [c1,h1]=m_contour(Lonr,Latr,h,[20 40 60 80 ],'color',[.6 .6 .6]);
%  [c2,h2]=m_contour(Lonr,Latr,h,[100 110 120 140 160 180],'k');
%  [c3,h3]=m_contour(Lonr,Latr,h,[200 220 240 260 280 ],'k');
%  [c4,h4]=contour(Lonr,Latr,h,[1000 1500 2000 3000 4000],'k');
% [c5,h5]=m_contour(Lonr,Latr,h,[100,1000 3000],'k');
% back_x=[min(min(Lonr)) max(max(Lonr)) max(max(Lonr)) min(min(Lonr))];
% back_y=[min(min(Latr)) min(min(Latr)) max(max(Latr)) max(max(Latr))];
% patch(back_x,back_y,[.9 .9 .9])
m_pcolor(Lonr,Latr,h.*themask);shading flat;
caxis(colaxis);bar_h=colorbar('fontsize',16,'fontweight','bold');
set(bar_h,'YTick',[0:2000:8000]);
A = cptcmap('GMT_ocean');
colormap(flipud(A))
% A = cptcmap('GMT_globe');
colormap(flipud(A(1:end/2,:)))
caxis([0 7000])
% caxis([0 5000])
% color=cbrewer('div','RdYlBu',50);
% color=flipud(color);
% colormap(color)
title(bar_h,'Depth(m)','FontSize',16,'fontweight','bold');



        m_gshhs_l('patch',[.95 .95 .8]);
%         m_grid('linestyle','none','box','fancy','tickdir','out','fontsize',14);
                m_grid('linestyle','none','tickdir','out','fontsize',16);
% m_text(126.7,40.5,'KOREA','fontsize',7,'rotation',-70,'fontweight','bold','FontName','Palatino Linotype');
% m_text(131.3,40.2,'East/','fontsize',7.,'fontweight','bold','FontAngle','italic','FontName','Palatino Linotype');
% m_text(130.5,39.1,'Japan','fontsize',7.,'fontweight','bold','FontAngle','italic','FontName','Palatino Linotype');
% m_text(131.3,38,'Sea','fontsize',7.,'fontweight','bold','FontAngle','italic','FontName','Palatino Linotype');
% m_text(121,36,'Yellow','fontsize',7.,'fontweight','bold','FontAngle','italic','FontName','Palatino Linotype');
% m_text(122,35,'Sea','fontsize',7.,'fontweight','bold','FontAngle','italic','FontName','Palatino Linotype');
% % m_text(121.5,26,'East China Sea','fontsize',9,'rotation',40,'fontweight','bold','FontAngle','italic','FontName','Palatino Linotype');
% m_text(113.,36.5,'CHINA','fontsize',7,'fontweight','bold','FontName','Palatino Linotype');
% m_text(133,35,'JAPAN','fontsize',7,'fontweight','bold','FontName','Palatino Linotype');
Y1=34;
Y2=44;
X1=132;
X2=132;    
box_gridx=[X1 X2 X2 X1 X1];box_gridy=[Y2 Y2 Y1 Y1 Y2];
m_plot(box_gridx,box_gridy,'--','color',[1 0 0],'linewidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
ratio_x=(domaxis(4)-domaxis(3))/100;
ratio_y=(domaxis(2)-domaxis(1))/100;
% m_text(domaxis(1)+ratio_x,domaxis(4)+(4*ratio_y),'(b)','color','k','FontSize',15,'fontweight','bold')
saveas(gcf,'eastsea_nest2_grid','tif')
% print('-dpng','roms_grid_ADD_10_ep')