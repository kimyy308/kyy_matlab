
lonind = 1; %% 130E
latminind=1; %% 34.8750N
latmaxind=17; %% 38.8750N
zminind=1; %% 0m
zmaxind=25; %% 200m
% [lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
% [lon3,ilevel3] = meshgrid(ilevel2(1:20),lat2(latminind:latmaxind));
% skipind=5;
% vecwidth=1.2;
ylength=latmaxind-latminind+1;
zlength=zmaxind-zminind+1;
level_c =5:5:20;
pcolindex=1;
contwidth=2;
confontsize=20;
temp_summer=(temp2_woa(lonind,latminind:latmaxind,zminind:zmaxind,6) + ...
    temp2_woa(lonind,latminind:latmaxind,zminind:zmaxind,7) + ...
    temp2_woa(lonind,latminind:latmaxind,zminind:zmaxind,8)) /3;
temp_winter=(temp2_woa(lonind,latminind:latmaxind,zminind:zmaxind,12) + ...
    temp2_woa(lonind,latminind:latmaxind,zminind:zmaxind,1) + ...
    temp2_woa(lonind,latminind:latmaxind,zminind:zmaxind,2)) /3;

sbp1=subplot(2,1,1);
pcolor(lat2_woa(latminind:latmaxind),ilevel2_woa(zminind:zmaxind),reshape(temp_summer,ylength,zlength)');
fig_woa_pcol_common;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);

[C,h]=contour(lat2_woa(latminind:latmaxind),ilevel2_woa(1:zmaxind),reshape(temp_summer,ylength,zlength)',level_c,'k','linewidth',contwidth);
clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',90,'fontweight','bold');
hold off;
 


sbp1=subplot(2,1,2);
pcolor(lat2_woa(latminind:latmaxind),ilevel2_woa(zminind:zmaxind),reshape(temp_winter,ylength,zlength)');
fig_woa_pcol_common;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);

[C,h]=contour(lat2_woa(latminind:latmaxind),ilevel2_woa(1:zmaxind),reshape(temp_winter,ylength,zlength)',level_c,'k','linewidth',contwidth);
clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',90,'fontweight','bold');
hold off;

h=colorbar;
caxis([2,20]);
colormap(jet);
set(h, 'Position', [.9150 .11 .0131 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'2^oC', '5^oC','8^oC','11^oC','14^oC','17^oC','20^oC'},...
            'YTick',[2 5 8 11 14 17 20]);
h.Label.FontSize=7;
set(get(h,'title'),'string',' ');


fig=gcf;
set(gcf,'PaperPosition',[0 0 24 24]);
fig.InvertHardcopy='off';
figname=figname_woa;
saveas(gcf,figname,'tiff');


hold off;



  