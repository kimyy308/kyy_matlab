
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

sbp1=subplot(4,1,1);
pcolor(lat2(1:latlimind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_1(lonind,1:latlimind,1:zlimind),latlimind,zlimind))');
temp_pcolor_common4;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);
quiver(ilevel3(1:skipind:latlimind,1:zlimind),lat3(1:skipind:latlimind,1:zlimind)/zlimind, ...
    (reshape(v2mean_1(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*3, ...
    (reshape(w2mean_1(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*2, ...
    '-k','filled','AutoScale','off', 'LineWidth',vecwidth);
tex2=text(tex2x,tex2y,'0.2m/s');
set(tex2,'fontsize',20);

% fig=gcf;
hold off;
 



sbp1=subplot(4,1,2);

pcolor(lat2(1:latlimind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_2(lonind,1:latlimind,1:zlimind),latlimind,zlimind))');
temp_pcolor_common4;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);

hold on;
quiver(ilevel3(1:skipind:latlimind,1:zlimind),lat3(1:skipind:latlimind,1:zlimind)/zlimind, ...
    (reshape(v2mean_2(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*3, ...
    (reshape(w2mean_2(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*2, ...
    '-k','filled','AutoScale','off', 'LineWidth',vecwidth);tex2=text(tex2x,tex2y,'0.2m/s');
set(tex2,'fontsize',20);
hold off;


sbp1=subplot(4,1,3);
pcolor(lat2(1:latlimind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_3(lonind,1:latlimind,1:zlimind),latlimind,zlimind))');
temp_pcolor_common4;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [2,20]);
hold on;
quiver(ilevel3(1:skipind:latlimind,1:zlimind),lat3(1:skipind:latlimind,1:zlimind)/zlimind, ...
    (reshape(v2mean_3(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*3, ...
    (reshape(w2mean_3(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*2, ...
    '-k','filled','AutoScale','off', 'LineWidth',vecwidth);tex2=text(tex2x,tex2y,'0.2m/s');
set(tex2,'fontsize',20);
hold off;
 

sbp1=subplot(4,1,4);
pcolor(lat2(1:latlimind),ilevel2(1:zlimind)/zlimind,(reshape(temp2mean_4(lonind,1:latlimind,1:zlimind),latlimind,zlimind))');
temp_pcolor_common4;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
hold on;
quiver(ilevel3(1:skipind:latlimind,1:zlimind),lat3(1:skipind:latlimind,1:zlimind)/zlimind, ...
    (reshape(v2mean_4(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*3, ...
    (reshape(w2mean_4(lonind,1:skipind:latlimind,1:zlimind),latlimind/skipind,zlimind))*2, ...
    '-k','filled','AutoScale','off', 'LineWidth',vecwidth);

h=colorbar;
caxis([2,20]);
colormap(jet);
set(h, 'Position', [.9150 .11 .0131 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'2^oC', '5^oC','8^oC','11^oC','14^oC','17^oC','20^oC'},...
            'YTick',[2 5 8 11 14 17 20]);
h.Label.FontSize=7;
set(get(h,'title'),'string',' ');

tex2=text(tex2x,tex2y,'0.2m/s');
set(tex2,'fontsize',20);

fig=gcf;
set(gcf,'PaperPosition',[0 0 25 30]);
fig.InvertHardcopy='off';
saveas(gcf,figname,'tiff');


hold off;



  