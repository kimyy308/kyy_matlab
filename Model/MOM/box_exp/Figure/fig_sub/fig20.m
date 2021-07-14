[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
pcolindex=1;
% % % 
% % % Exp1, temp(surf) + vec(surf)
% % % 
sbp1=subplot(2,2,1);
hold on;
% vector_common1; % for white background color
% caxis([-10 10]); % for white background color

pcolor(lon5,lat5,(vort5mean_1(:,:,1))');
temp_pcolor_common3;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [-1e-5,1e-5]);

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

pcolor(lon5,lat5,(vort5mean_2(:,:,1))');
temp_pcolor_common3;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [-1e-5,1e-5]);

hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_2(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_2(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
vector_common2;
tex2=text(137.5,36.0,'(b)');
set(tex2,'fontsize',20);
hold off;


sbp1=subplot(2,2,3);
pcolor(lon5,lat5,(vort5mean_3(:,:,1))');
temp_pcolor_common3;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [-1e-5,1e-5]);
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_3(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_3(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
vector_common2;
tex2=text(137.5,36.0,'(c)');
set(tex2,'fontsize',20);
hold off;
 

sbp1=subplot(2,2,4);
pcolor(lon5,lat5,(vort5mean_4(:,:,1))');
temp_pcolor_common3;
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
caxis(gca, [-1e-5,1e-5]);
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_4(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_4(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
h=colorbar;
set(h, 'Position', [.8950 .11 .0181 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'<-1.0，10^-^5','-0.8，10^-^5','-0.6，10^-^5','-0.4，10^-^5','-0.2，10^-^5','0', ...
                '0.2，10^-^5','0.4，10^-^5','0.6，10^-^5','0.8，10^-^5', '>1.0，10^-^5'}, ...
                'YTick',[-1.0e-5 -0.8e-5 -0.6e-5 -0.4e-5 -0.2e-5 0 0.2e-5 0.4e-5 0.6e-5 0.8e-5 1.0e-5]);
h.Label.FontSize=7;
set(get(h,'title'),'string','Relative vorticity (s^-^1)');
% readname='D:\MEPL\Master_course\for_publish\Figure\';
% readname2='Fig19_9.tif';
% filename=strcat(readname,readname2);
tex2=text(137.5,36.0,'(d)');
set(tex2,'fontsize',20);
vector_common2;
hold off;



  