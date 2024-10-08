[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
pcolindex=1;
lon_min_ind = 3;  %% 128.3 E
lon_max_ind = 50; %% 133 E
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
% % % 
% % % Exp1, temp(surf) + vec(surf)
% % % 
sbp1=subplot(2,2,1);

% vector_common1; % for white background color
% caxis([-10 10]); % for white background color

% pcolor(lon5,lat5,(isod5mean_1(:,:,1))');
pcolor(lon5(lon_min_ind:lon_max_ind),lat5(lat_min_ind:lat_max_ind),(isod5mean_1(lon_min_ind:lon_max_ind,lat_min_ind:lat_max_ind,1))');
fig23_temp_pcolor_common;
caxis(gca, [0,200]);
hold on;
% quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
%   squeeze(u2mean_1(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_1(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
%   'k','AutoScale','off');
quiver(lon(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind),lat(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind), ...
  squeeze(u2mean_1(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  squeeze(v2mean_1(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  'k','AutoScale','off', 'LineWidth',vecwidth);

fig23_vector_common;
tex2=text(tex2x,tex2y,'(a)');
tex=text(texx,texy,'0.2m/s');
set(tex,'fontsize',texsize, 'FontWeight',texFontWeight);
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight);
% fig=gcf;
hold off;
 

sbp1=subplot(2,2,2);
% vector_common1;  % for white background color
% caxis([-10 10]); % for white background color
pcolor(lon5(lon_min_ind:lon_max_ind),lat5(lat_min_ind:lat_max_ind),(isod5mean_2(lon_min_ind:lon_max_ind,lat_min_ind:lat_max_ind,1))');
fig23_temp_pcolor_common;
hold on;
caxis(gca, [0,200]);
quiver(lon(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind),lat(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind), ...
  squeeze(u2mean_2(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  squeeze(v2mean_2(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  'k','AutoScale','off', 'LineWidth',vecwidth);

fig23_vector_common;
tex2=text(tex2x,tex2y,'(b)');
tex=text(texx,texy,'0.2m/s');
set(tex,'fontsize',texsize, 'FontWeight',texFontWeight);
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight);
hold off;


sbp1=subplot(2,2,3);
pcolor(lon5(lon_min_ind:lon_max_ind),lat5(lat_min_ind:lat_max_ind),(isod5mean_3(lon_min_ind:lon_max_ind,lat_min_ind:lat_max_ind,1))');
fig23_temp_pcolor_common;
hold on;
caxis(gca, [0,200]);
quiver(lon(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind),lat(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind), ...
  squeeze(u2mean_3(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  squeeze(v2mean_3(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  'k','AutoScale','off', 'LineWidth',vecwidth);

fig23_vector_common;
tex2=text(tex2x,tex2y,'(c)');
tex=text(texx,texy,'0.2m/s');
set(tex,'fontsize',texsize, 'FontWeight',texFontWeight);
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight);
hold off;
 

sbp1=subplot(2,2,4);
pcolor(lon5(lon_min_ind:lon_max_ind),lat5(lat_min_ind:lat_max_ind),(isod5mean_4(lon_min_ind:lon_max_ind,lat_min_ind:lat_max_ind,1))');
fig23_temp_pcolor_common;
hold on;
caxis(gca, [0,200]);
quiver(lon(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind),lat(lat_min_ind:skipind:lat_max_ind,lon_min_ind:skipind:lon_max_ind), ...
  squeeze(u2mean_4(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  squeeze(v2mean_4(lon_min_ind:skipind:lon_max_ind,lat_min_ind:skipind:lat_max_ind,quivdepthindex,1))'*2, ...
  'k','AutoScale','off', 'LineWidth',vecwidth);


h=colorbar;
set(h, 'Position', [.9150 .11 .0181 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'20m','40m','60m','80m','100m','120m','140m','160m','180m','200m'},...
        'YTick',[20 40 60 80 100 120 140 160 180 200]);
h.Label.FontSize=7;
set(get(h,'title'),'string','Upper layer');
tex2=text(tex2x,tex2y,'(d)');
tex=text(texx,texy,'0.2m/s');
set(tex,'fontsize',texsize, 'FontWeight',texFontWeight);
set(tex2,'fontsize',tex2size, 'FontWeight',texFontWeight);
fig23_vector_common;
hold off;

fig=gcf;
set(gcf,'PaperPosition',[0 0 27 30]);
fig.InvertHardcopy='off';
saveas(gcf,figname,'tiff');



  