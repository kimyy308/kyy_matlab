[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
pcolindex=1;
% % % 
% % % Exp1, temp(surf) + vec(surf)
% % % 
sbp1=subplot(2,2,1);
hold on;
% vector_common1; % for white background color
% caxis([-10 10]); % for white background color

pcolor(lon5,lat5,(isod5mean_1(:,:,1))');
temp_pcolor_common2;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [0,210]);

quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_1(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_1(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
vector_common2;
tex2=text(137.5,36.0,'(a)');
set(tex2,'fontsize',20);
fig=gcf;
hold off;
 
sbp1=subplot(2,2,2);
% vector_common1;  % for white background color
% caxis([-10 10]); % for white background color

pcolor(lon5,lat5,(isod5mean_2(:,:,1))');
temp_pcolor_common2;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [0,210]);

hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_2(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_2(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
vector_common2;
tex2=text(137.5,36.0,'(b)');
set(tex2,'fontsize',20);
hold off;


sbp1=subplot(2,2,3);
pcolor(lon5,lat5,(isod5mean_3(:,:,1))');
temp_pcolor_common2;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [0,210]);
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_3(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_3(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
vector_common2;
tex2=text(137.5,36.0,'(c)');
set(tex2,'fontsize',20);
hold off;
 

sbp1=subplot(2,2,4);
pcolor(lon5,lat5,(isod5mean_4(:,:,1))');
temp_pcolor_common2;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [0,210]);
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_4(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_4(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
h=colorbar;
set(h, 'Position', [.9150 .11 .0181 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'20m','40m','60m','80m','100m','120m','140m','160m','180m','200m'},...
        'YTick',[20 40 60 80 100 120 140 160 180 200]);
h.Label.FontSize=7;
set(get(h,'title'),'string','upper layer');
tex2=text(137.5,36.0,'(d)');
set(tex2,'fontsize',20);
vector_common2;
hold off;



  