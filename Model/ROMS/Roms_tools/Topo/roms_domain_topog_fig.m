clear all; clc; close all;
rname='D:\MEPL\project\NWP\Roms_tools\Topo\etopo2.nc';
ht2=ncread(rname,'topo');
lon2=ncread(rname,'lon');
lat2=ncread(rname,'lat');
lons=find(lon2==115);
lone=find(lon2==162);
lats=find(lat2==15);
late=find(lat2==52);
for i=lons:lone
    for j=lats:late
        if (ht2(i,j)>=0)
          ht2(i,j)=NaN;
        end
    end
end
filename='D:\MEPL\project\NWP\Roms_tools\Topo\topog.tif';
 pcolor(lon2(lons:lone),lat2(lats:late),(ht2(lons:lone,lats:late)'));
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
    caxis([min(min(ht2(lons:lone,lats:late))) 0]);
    c=colorbar;
    c.Label.String= 'Depth (M)';
    c.Label.FontSize= 12;
    colormap jet;
    xlabel('Longitude (^oE)','fontsize',27);
    ylabel('Latitude (^oN)','fontsize',27);
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',15);
  hold off;
    ar=annotation(gcf,'rectangle',...
    [0.316595 0.50848648 0.00353378378378377 0.06],'FaceColor','k');
    ar=annotation(gcf,'rectangle',...
    [0.316595 0.50848648 0.043 0.00353378378378377],'FaceColor','k');
    ar=annotation(gcf,'rectangle',...
    [0.316595 0.56848648 0.043 0.00353378378378377],'FaceColor','k');
    ar=annotation(gcf,'rectangle',...
    [0.358995 0.50848648 0.00353378378378377 0.063],'FaceColor','k');
    fig=gcf;
    fig.InvertHardcopy='off';
    set(gca, 'color', [0.8,0.8,0.8]);
%     grid on
%     set (gca, 'layer','top');
    saveas(gcf,filename,'tiff');