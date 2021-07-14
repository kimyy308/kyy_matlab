
lonmaxind=24;
lonminind=3;
[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
[lon3,ilevel3] = meshgrid(ilevel2(1:20),lon2(lonminind:lonmaxind));
latind = 20; %% 20 : 36N
zlimind=20; %% ~ 200m
skipind=5;
vecwidth=1.2;
xlength=lonmaxind-lonminind+1;
level_c =-0.5:0.05:0.5;
pcolindex=1;
contwidth=2;
confontsize=20;

sbp1=subplot(4,1,1);
pcolor(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_1(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))');
fig22_pcol_common;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);

[C,h]=contour(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(v2mean_1(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))',level_c,'k','linewidth',contwidth);
clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',90,'fontweight','bold');
hold off;
 



sbp1=subplot(4,1,2);

pcolor(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_2(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))');
fig22_pcol_common;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);

hold on;
[C,h]=contour(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(v2mean_2(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))',level_c,'k','linewidth',contwidth);
clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
hold off;


sbp1=subplot(4,1,3);
pcolor(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_3(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))');
fig22_pcol_common;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);
hold on;
[C,h]=contour(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(v2mean_3(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))',level_c,'k','linewidth',contwidth);
clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
hold off;
 

sbp1=subplot(4,1,4);
pcolor(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_4(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))');
fig22_pcol_common;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
hold on;
[C1,h1]=contour(lon2(lonminind:lonmaxind),ilevel2(1:zlimind)/zlimind,(reshape(v2mean_4(lonminind:lonmaxind,latind,1:zlimind),xlength,zlimind))',level_c,'k','linewidth',contwidth);
clabel(C1,h1,'FontSize',confontsize,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
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
set(gcf,'PaperPosition',[0 0 24 30]);
fig.InvertHardcopy='off';
saveas(gcf,figname,'tiff');


hold off;



  