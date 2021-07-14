clear all

expname='oman_restore_04_18';

readname='D:\need_to_presentation\';
readname2='\197912\ocean_gridinfo_197912.nc';
rname=strcat(readname,expname,readname2);
ht2=ncread(rname,'ht');
lon2=ncread(rname,'xt_ocean');
lat2=ncread(rname,'yt_ocean');

[lat,lon] = meshgrid(lat2,lon2);
for i=1:120
    for j=1:110
        if (ht2(i,j)< -1.0e+5)
          ht2(i,j)=NaN;
        end
    end
end

'reading data is completed'

%%
%% figure 6, topography field
%%
 readname='D:\need_to_presentation\';
 readname2='\paper\topog.tif';
 filename=strcat(readname,expname,readname2);
 gifname = sprintf('100m Vector & Temperature, 10th year averaged value');
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
%     axis equal;
    set(gca,'fontsize',15);
  hold off;
    fig=gcf;
    fig.InvertHardcopy='off';
    set(gca, 'color', [0.8,0.8,0.8]);
    saveas(gcf,filename,'tiff');
'making topog figure is completed'
    
    