clear all
% rname ='D:\need_to_presentation\WOA2013/woa2013_temp.nc';
rname ='D:\need_to_presentation\WOA2013/T_1_7.nc';
temp2=ncread(rname,'T_MOD');
lon2=ncread(rname,'LON1240_1240');
lat2=ncread(rname,'LAT500_516');
level2=ncread(rname,'DEPTH1_25');
ilevel2 = -level2;


% 
% July
%

filename=('D:\need_to_presentation\WOA2013/july.png');
gifname=('Temperature at 130°E, JJA');
pcolor(lat2(1:17),ilevel2(1:25),(reshape((temp2(1,1:17,1:25,6)+temp2(1,1:17,1:25,7)+temp2(1,1:17,1:25,8))/3,17,25))');
hold on;
shading interp;
c=colorbar;
c.Label.String= 'Temp (°C)';
c.Label.FontSize= 12;
colormap jet;
% title(gifname,'fontsize',15);
xlabel('Latitude (°N)','fontsize',15);
ylabel('Depth (M)','fontsize',15);
zindex = 10;
caxis([0.2,20]);
c=colorbar;
c.Label.String= 'Temp (°C)';
c.Label.FontSize= 12;
% con=contour(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))','k', 'ShowText','on','LevelStep',5, 'LineWidth',2);
[C,h]=contour(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))','k','LevelStep',5, 'LineWidth',2);
t=clabel(C,h,'FontSize',15, 'LabelSpacing', 300)
% contourc(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))',4);
set(gcf,'PaperPosition',[0 0 30 18]);  %왼쪽 여백, 위쪽 여백, 가로길이, 세로 길이
set(gca,'fontsize',20);
saveas(gcf,filename,'png');
hold off


% 
% January
%

filename=('D:\need_to_presentation\WOA2013/jan.png');
gifname=('Temperature at 130°E, DJF');
pcolor(lat2(1:17),ilevel2(1:25),(reshape((temp2(1,1:17,1:25,12)+temp2(1,1:17,1:25,1)+temp2(1,1:17,1:25,2))/3,17,25))');
hold on;
shading interp;
c=colorbar;
c.Label.String= 'Temp (°C)';
c.Label.FontSize= 12;
colormap jet;
% title(gifname,'fontsize',15);
xlabel('Latitude (°N)','fontsize',15);
ylabel('Depth (M)','fontsize',15);
zindex = 10;

caxis([0.2,20]);
c=colorbar;
c.Label.String= 'Temp (°C)';
c.Label.FontSize= 12;
% con=contour(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))','k', 'ShowText','on','LevelStep',5, 'LineWidth',2);
[C,h]=contour(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,1),17,25))','k','LevelStep',5, 'LineWidth',2);
t=clabel(C,h,'FontSize',15, 'LabelSpacing', 300)
% contourc(lat2(1:17),ilevel2(1:25),(reshape(temp2(1,1:17,1:25,7),17,25))',4);
set(gcf,'PaperPosition',[0 0 30 18]);  %왼쪽 여백, 위쪽 여백, 가로길이, 세로 길이
set(gca,'fontsize',20);
saveas(gcf,filename,'png');
hold off