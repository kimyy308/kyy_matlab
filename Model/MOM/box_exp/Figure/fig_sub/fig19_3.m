[lon44,ilevel44] = meshgrid(ilevel2(1:20),lon2(1:30));
pcolindex=1;
% % % 
% % % Exp1, temp(surf) + vec(surf)
% % % 


sbp1=subplot(2,2,1);
pcolor(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean_1(50,31:50,1:20,1),20,20))');
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean_1(50,31:50,1:20,1),20,20))','ShowText','on');
    h.LineColor='black';
vector_common4;
hold off;

transport_11 = 0;
lon_beg = 133; lon_end = 133;
for lat_beg = 37.0 : 0.1 : 38.9
    lat_end = lat_beg + 0.1; 
    for depth_beg = 0 : 10 : 190
        if (u2mean_1(50, int8(31+10*(lat_beg - 37.0)), (1+depth_beg/10), 1)>=0)
            lon_lat_distance;
            transport_11 = transport_11 + u2mean_1(50, int8(31+10*(lat_beg - 37.0)), (1+depth_beg/10), 1) * dist * 10;
        end
    end
end
transport_11

sbp2=subplot(2,2,2);
pcolor(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean_2(50,31:50,1:20,1),20,20))');
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean_2(50,31:50,1:20,1),20,20))','ShowText','on');
    h.LineColor='black';
vector_common4;
hold off;

transport_22 = 0;
lon_beg = 133; lon_end = 133;
for lat_beg = 37.0 : 0.1 : 38.9
    lat_end = lat_beg + 0.1; 
    for depth_beg = 0 : 10 : 190
        if (u2mean_2(50, int8(31+10*(lat_beg - 37.0)), (1+depth_beg/10), 1)>=0)
            lon_lat_distance;
            transport_22 = transport_22 + u2mean_2(50, int8(31+10*(lat_beg - 37.0)), (1+depth_beg/10), 1) * dist * 10;
        end
    end
end
transport_22

% % % 38N (40)
sbp3=subplot(2,2,3);
pcolor(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean_3(50,31:50,1:20,1),20,20))');
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean_3(50,31:50,1:20,1),20,20))','ShowText','on');
    h.LineColor='black';
vector_common4;
hold off;

transport_33 = 0;
lon_beg = 133; lon_end = 133;
for lat_beg = 37.0 : 0.1 : 38.9
    lat_end = lat_beg + 0.1; 
    for depth_beg = 0 : 10 : 190
        if (u2mean_3(50, int8(31+10*(lat_beg - 37.0)), (1+depth_beg/10), 1)>=0)
            lon_lat_distance;
            transport_33 = transport_33 + u2mean_3(50, int8(31+10*(lat_beg - 37.0)), (1+depth_beg/10), 1) * dist * 10;
        end
    end
end
transport_33


sbp4=subplot(2,2,4);
pcolor(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean_4(50,31:50,1:20,1),20,20))');
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean_4(50,31:50,1:20,1),20,20))','ShowText','on');
    h.LineColor='black';
h=colorbar;
caxis([-0.35,0.35]);
colormap(bwrmap);
set(h, 'Position', [.9150 .11 .0131 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'-0.2m/s', '0m/s', '0.2m/s'},...
        'YTick',[-0.2 0 0.2]);
h.Label.FontSize=5;
% set(h,'Visible','off');
vector_common4;
hold off;


transport_44 = 0;
lon_beg = 133; lon_end = 133;
for lat_beg = 37.0 : 0.1 : 38.9
    lat_end = lat_beg + 0.1; 
    for depth_beg = 0 : 10 : 190
        if (u2mean_4(50, int8(31+10*(lat_beg - 37.0)), (1+depth_beg/10), 1)>=0)
            lon_lat_distance;
            transport_44 = transport_44 + u2mean_4(50, int8(31+10*(lat_beg - 37.0)), (1+depth_beg/10), 1) * dist * 10;
        end
    end
end
transport_44
    
    
    
    
    
    
  