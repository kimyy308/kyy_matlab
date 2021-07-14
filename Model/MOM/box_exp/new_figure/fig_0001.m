% %  World ocean atlas 2013 climatology temperature N-S section in Feb and in Aug

rname ='E:\Data\Observation\WOA\2013\WOA2013\T_1_7.nc';
temp2_woa=ncread(rname,'T_MOD');
lon2_woa=ncread(rname,'LON1240_1240');
lat2_woa=ncread(rname,'LAT500_516');
level2_woa=ncread(rname,'DEPTH1_25');
ilevel2_woa = -level2_woa/20.0;

rname ='E:\Data\Observation\WOA\2013\WOA2013\T_1_7_FNR.nc';
temp2_woa_fnr=ncread(rname,'T_MOD_FNR');

for i=1:12
    temp2_woa_fnr2(:,:,i)=griddata(lat2_woa,ilevel2_woa,squeeze(temp2_woa_fnr(1,:,:,i))',(34.875:0.025:38.875)',(0:-0.125:-10));  %% [81 161 12]
end
lat2_woa_fnr=(34.875:0.025:38.875)';
ilevel2_woa_fnr=(0:-0.125:-10);

yfilled(1:10)=-6.25;
for i=11:31
    yfilled(i)=-6.25-(i-10)*0.125;
end
yfilled(32)=-10;
yfilled(33)=-10;
yfilled(34)=-6.25;
xfilled(1:31)=lat2_woa_fnr(1:31);
xfilled(32)=lat2_woa_fnr(31);
xfilled(33)=lat2_woa_fnr(1);
xfilled(34)=lat2_woa_fnr(1);
'reading data for figure_woa is completed'


figname1='D:\OneDrive - 서울대학교\research\Master_course\for_publish\new_figure\';
figname_woa2='fig_0001';
figname_woa=strcat(figname1,figname_woa2);

lonind = 1; %% 130E
latminind=1; %% 34.8750N
latmaxind=161; %% 38.8750N
zminind=1; %% 0m
zmaxind=81; %% 200m
ylength=latmaxind-latminind+1;
zlength=zmaxind-zminind+1;
level_c =5:5:20;
pcolindex=1;
contwidth=2;
confontsize=20;
tex2size =25;
tex2x = 38;
tex2y = -8.5;
texFontWeight= 'bold';   %%'normal or bold'

temp_summer=(temp2_woa_fnr2( zminind:zmaxind,latminind:latmaxind,8));
temp_winter=(temp2_woa_fnr2(zminind:zmaxind, latminind:latmaxind,2));



sbp1=subplot(2,1,2);
pcolor(lat2_woa_fnr(latminind:latmaxind),ilevel2_woa_fnr(zminind:zmaxind),temp_summer);
    shading interp;
    xlim([min(lat2_woa_fnr(latminind:latmaxind)) max(lat2_woa_fnr(latminind:latmaxind))]);
    set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
        'YTick', [-8, -6, -4, -2, 0])
    set(gca,'xTickLabel',{'35^oN','36^oN','37^oN','38^oN'}, ...
        'XTick', (35:1:39))
    set(gca, 'color', [0.8,0.8,0.8]);
    colormap(gca,jet);
    set(gcf,'PaperPosition',[0 0 16 12]);
    set(gca,'fontsize',20);

hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [0,25]);

[C,h]=contour(lat2_woa_fnr(latminind:latmaxind),ilevel2_woa_fnr(1:zmaxind),temp_summer,level_c,'k','linewidth',contwidth);

clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',300,'Rotation',90,'fontweight','bold');
fill(xfilled,yfilled,[0.8 0.8 0.8]);
tex2=text(tex2x,tex2y,'August');
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight,'Color','w');
hold off;
 


sbp1=subplot(2,1,1);
pcolor(lat2_woa_fnr(latminind:latmaxind),ilevel2_woa_fnr(zminind:zmaxind),temp_winter);
    shading interp;
    xlim([min(lat2_woa_fnr(latminind:latmaxind)) max(lat2_woa_fnr(latminind:latmaxind))]);
    set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
        'YTick', [-8, -6, -4, -2, 0])
    set(gca,'xTickLabel',{'35^oN','36^oN','37^oN','38^oN'}, ...
        'XTick', (35:1:39))
    set(gca, 'color', [0.8,0.8,0.8]);
    colormap(gca,jet);
    set(gcf,'PaperPosition',[0 0 16 12]);
    set(gca,'fontsize',20);
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [0,25]);

[C,h]=contour(lat2_woa_fnr(latminind:latmaxind),ilevel2_woa_fnr(1:zmaxind),temp_winter,level_c,'k','linewidth',contwidth);
clabel(C,h,'FontSize',confontsize,'Color','k','labelspacing',300,'Rotation',90,'fontweight','bold');
fill(xfilled,yfilled,[0.8 0.8 0.8]);
tex2=text(tex2x,tex2y,'February');
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight,'Color','w');
hold off;

h=colorbar;
caxis([0,25]);
colormap(jet);
set(h, 'Position', [.9150 .11 .0131 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'0^oC', '5^oC','10^oC','15^oC','20^oC','25^oC'},...
            'YTick',[0 5 10 15 20 25]);
h.Label.FontSize=7;
set(get(h,'title'),'string',' ');


fig=gcf;
set(gcf,'PaperPosition',[0 0 24 24]);
fig.InvertHardcopy='off';
figname=figname_woa;
saveas(gcf,figname,'jpg');

hold off;

'making figure_woa_fnr is completed'
