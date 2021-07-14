 
rname='D:\MEPL\project\NWP\Roms_tools\Topo\etopo2.nc';
ht2=ncread(rname,'topo');
lon2=ncread(rname,'lon');
lat2=ncread(rname,'lat');
lons=find(lon2==128);
lone=find(lon2==131);
lats=find(lat2==33);
late=find(lat2==36);
for i=lons:lone
    for j=lats:late
        if (ht2(i,j)>=0)
          ht2(i,j)=NaN;
        end
    end
end
filename='D:\MEPL\project\NWP\Roms_tools\Topo\ks_topog.tif';
 fff=pcolor(lon2(lons:lone),lat2(lats:late),(ht2(lons:lone,lats:late)'));
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
%     caxis([min(min(ht2(lons:lone,lats:late))) 0]);
%     c=colorbar;
%     c.Label.String= 'Depth (M)';
%     c.Label.FontSize= 12;
%     colormap jet;
%     xlabel('Longitude (^oE)','fontsize',27);
%     ylabel('Latitude (^oN)','fontsize',27);
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',15);
%     grid on
%     grid minor

% axis square
% grid(ax,
    ax=gca;
  hold off;
    fig=gcf;
    fig.InvertHardcopy='off';
    set(gca, 'color', [0.8,0.8,0.8]);
%      set(gca, 'XTickLabel', 128:0.5:131,'Xtick',128:0.05:131)
%     set (gca, 'XTickLabel',128:0.5:131);
    set(gca,'XTickLabel',[])
    set (gca, 'XTick',128:0.05:131);
%     set (gca, 'XTickLabelMode','manual');
%     set (gca, 'XTickLabelRotation', 0);
    set (gca, 'XGrid','on');
%     set(gca, 'YTickLabel', 33:0.5:36,'Ytick',33:0.05:36)
    set (gca, 'YTick',33:0.05:36);
%     set (gca, 'YTickLabel',33:0.5:36);
    set(gca,'YTickLabel',[])
%     set (gca, 'YTickLabelMode','manual');
    set (gca, 'YGrid','on');
    set (gca, 'GridColor','k');
    set (gca, 'GridAlpha',1);
    set (gca, 'layer','top');
    set (gca, 'Linewidth',1);
    saveas(gcf,filename,'tiff');