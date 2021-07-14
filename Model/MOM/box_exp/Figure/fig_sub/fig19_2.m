[lon44,ilevel44] = meshgrid(ilevel2(1:20),lon2(1:30));
pcolindex=1;
% % % 
% % % Exp1, temp(surf) + vec(surf)
% % % 


sbp1=subplot(2,2,1);
pcolor(ilevel44,lon44,reshape(v2mean_1(1:30,60,1:20,1),30,20));
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(ilevel44,lon44,reshape(v2mean_1(1:30,60,1:20,1),30,20),'ShowText','on');
    h.LineColor='black';
vector_common3;
hold off;

transport_1 = 0;
lat_beg = 40; lat_end = 40;
for lon_beg = 128.5 : 0.1 : 130.9
    lon_end = lon_beg + 0.1; 
    for depth_beg = 0 : 10 : 190
        if (v2mean_1(int8(5+10*(lon_beg - 128.5)), 60, (1+depth_beg/10), 1)>=0)
            lon_lat_distance;
            transport_1 = transport_1 + v2mean_1(int8(5+10*(lon_beg - 128.5)), 60, (1+depth_beg/10), 1) * dist * 10;
        end
    end
end
transport_1

sbp2=subplot(2,2,2);
pcolor(ilevel44,lon44,reshape(v2mean_2(1:30,60,1:20,1),30,20));
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(ilevel44,lon44,reshape(v2mean_2(1:30,60,1:20,1),30,20),'ShowText','on');
    h.LineColor='black';
vector_common3;
hold off;

transport_2 = 0;
lat_beg = 40; lat_end = 40;
for lon_beg = 128.5 : 0.1 : 130.9
    lon_end = lon_beg + 0.1; 
    for depth_beg = 0 : 10 : 190
        if (v2mean_2(int8(5+10*(lon_beg - 128.5)), 60, (1+depth_beg/10), 1)>=0)
            lon_lat_distance;
            transport_2 = transport_2 + v2mean_2(int8(5+10*(lon_beg - 128.5)), 60, (1+depth_beg/10), 1) * dist * 10;
        end
    end
end
transport_2

% % % 38N (40)
sbp3=subplot(2,2,3);
pcolor(ilevel44,lon44,reshape(v2mean_3(1:30,40,1:20,1),30,20));
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(ilevel44,lon44,reshape(v2mean_4(1:30,40,1:20,1),30,20),'ShowText','on');
    h.LineColor='black';
vector_common3;
hold off;

transport_3 = 0;
lat_beg = 38; lat_end = 38;
for lon_beg = 128.5 : 0.1 : 130.9
    lon_end = lon_beg + 0.1; 
    for depth_beg = 0 : 10 : 190
        if (v2mean_3(int8(5+10*(lon_beg - 128.5)), 40, (1+depth_beg/10), 1)>=0)
            lon_lat_distance;
            transport_3 = transport_3 + v2mean_3(int8(5+10*(lon_beg - 128.5)), 40, (1+depth_beg/10), 1) * dist * 10;
        end
    end
end
transport_3


% % % % % % NKCC
% % transport_3 = 0;
% % lat_beg = 38; lat_end = 38;
% % for lon_beg = 128.5 : 0.1 : 129.0
% %     lon_end = lon_beg + 0.1; 
% %     for depth_beg = 0 : 10 : 190
% %         if (v2mean_3(int8(5+10*(lon_beg - 128.5)), 40, (1+depth_beg/10), 1) < 0.)
% %             lon_lat_distance;
% %             transport_3 = transport_3 + v2mean_3(int8(5+10*(lon_beg - 128.5)), 40, (1+depth_beg/10), 1) * dist * 10
% %         end
% %     end
% % end
% % transport_3





sbp4=subplot(2,2,4);
pcolor(ilevel44,lon44,reshape(v2mean_4(1:30,60,1:20,1),30,20));
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(ilevel44,lon44,reshape(v2mean_4(1:30,60,1:20,1),30,20),'ShowText','on');
    h.LineColor='black';
h=colorbar;
caxis([-0.35,0.35]);
colormap(bwrmap);
set(h, 'Position', [.9150 .11 .0131 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'-0.2m/s', '0m/s', '0.2m/s'},...
            'YTick',[-0.2 0 0.2]);
h.Label.FontSize=7;
% set(h,'Visible','off');
vector_common3;
hold off;


transport_4 = 0;
lat_beg = 40; lat_end = 40;
for lon_beg = 128.5 : 0.1 : 130.9
    lon_end = lon_beg + 0.1; 
    for depth_beg = 0 : 10 : 190
        if (v2mean_4(int8(5+10*(lon_beg - 128.5)), 60, (1+depth_beg/10), 1)>=0)
            lon_lat_distance;
            transport_4 = transport_4 + v2mean_4(int8(5+10*(lon_beg - 128.5)), 40, (1+depth_beg/10), 1) * dist * 10;
        end
    end
end
transport_4
    
    
% % % % % % NKCC
% % transport_4 = 0;
% % lat_beg = 38; lat_end = 38;
% % for lon_beg = 128.5 : 0.1 : 129.0
% %     lon_end = lon_beg + 0.1; 
% %     for depth_beg = 0 : 10 : 190
% %         if (v2mean_4(int8(5+10*(lon_beg - 128.5)), 60, (1+depth_beg/10), 1) < 0.)
% %             lon_lat_distance;
% %             transport_4 = transport_4 + v2mean_4(int8(5+10*(lon_beg - 128.5)), 60, (1+depth_beg/10), 1) * dist * 10;
% %         end
% %     end
% % end
% % transport_4
    
    
    
  