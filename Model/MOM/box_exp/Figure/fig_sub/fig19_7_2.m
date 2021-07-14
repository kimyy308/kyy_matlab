addpath(genpath('D:\MEPL\Master_course\for_publish\seawater_ver3_2'));

[lon44,ilevel44] = meshgrid(ilevel2(1:20),lon2(1:15));
pcolindex=1;
% % % 
% % % Exp1, temp(surf) + vec(surf)
% % % 

% % % 40N (j=60)
sbp1=subplot(2,2,1);
pcolor(ilevel44,lon44,reshape(temp2mean_1(1:15,60,1:20,1),15,20));
shading interp;
caxis([2,20]);
hold on;
[C,h]=contour(ilevel44,lon44,reshape(v2mean_1(1:15,60,1:20,1),15,20),'ShowText','on'); %% 
    h.LineColor='black';
vector_common5;
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
pcolor(ilevel44,lon44,reshape(temp2mean_2(1:15,60,1:20,1),15,20));
shading interp;
caxis([2,20]);
hold on;
[C,h]=contour(ilevel44,lon44,reshape(v2mean_2(1:15,60,1:20,1),15,20),'ShowText','on');
h.LineColor='black';
[C2,h2]=contour(ilevel44,lon44,reshape(v2mean_2(1:15,60,1:20,1),15,20),[-0.005 -0.005],'ShowText','on','LineStyle','--');
vector_common5;
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
pcolor(ilevel44,lon44,reshape(temp2mean_3(1:15,40,1:20,1),15,20));
shading interp;
caxis([2,20]);
hold on;
[C,h]=contour(ilevel44,lon44,reshape(v2mean_3(1:15,40,1:20,1),15,20),'ShowText','on');
h.LineColor='black';
[C2,h2]=contour(ilevel44,lon44,reshape(v2mean_3(1:15,40,1:20,1),15,20),[-0.005 -0.01 -0.015 -0.02],'ShowText','on','LineStyle','--');
vector_common5;
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
pcolor(ilevel44,lon44,reshape(temp2mean_4(1:15,60,1:20,1),15,20));
shading interp;
caxis([2,20]);
hold on;
[C,h]=contour(ilevel44,lon44,reshape(v2mean_4(1:15,60,1:20,1),15,20),'ShowText','on');
    h.LineColor='black';
[C2,h2]=contour(ilevel44,lon44,reshape(v2mean_4(1:15,60,1:20,1),15,20),[-0.005 -0.005],'ShowText','on','LineStyle','--');
h=colorbar;
caxis([2,20]);
colormap(jet);
set(h, 'Position', [.9150 .11 .0131 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'2^oC', '5^oC','8^oC','11^oC','14^oC','17^oC','20^oC'},...
            'YTick',[2 5 8 11 14 17 20]);
h.Label.FontSize=7;
% set(h,'Visible','off');
vector_common5;
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
    
    
% calculate geostrophic velocity
% squeeze(temp2mean_1(1:15,60,1:20,1));
depth=5:10:195;
for i=1:15
    depths(i,:)=depth;
end
salt2mean_1(1:15,1:20)=35.;




ga_exp1(:,:)=sw_gpan(salt2mean_1',squeeze(temp2mean_1(1:15,60,1:20,1))',depths');
ga_exp2(:,:)=sw_gpan(salt2mean_1',squeeze(temp2mean_2(1:15,60,1:20,1))',depths');
ga_exp3(:,:)=sw_gpan(salt2mean_1',squeeze(temp2mean_3(1:15,40,1:20,1))',depths');
ga_exp4(:,:)=sw_gpan(salt2mean_1',squeeze(temp2mean_4(1:15,60,1:20,1))',depths');
velp_exp1(:,:) = sw_gvel(ga_exp1(:,:),squeeze(lat(60,1:15)),squeeze(lon(60,1:15)));
velp_exp2(:,:) = sw_gvel(ga_exp2(:,:),squeeze(lat(60,1:15)),squeeze(lon(60,1:15)));
velp_exp3(:,:) = sw_gvel(ga_exp3(:,:),squeeze(lat(40,1:15)),squeeze(lon(40,1:15)));
velp_exp4(:,:) = sw_gvel(ga_exp4(:,:),squeeze(lat(60,1:15)),squeeze(lon(60,1:15)));



% geostrophic velocity of exp.1
for vv=1:length(velp_exp1(1,:))
    mnid=max(find(~isnan(velp_exp1(:,vv))));
    if isempty(mnid)
        idm(vv)=1;
    else
        idm(vv)=mnid;
    end
    gvel_exp1(:,vv)=velp_exp1(:,vv)-velp_exp1(idm(vv),vv);
end
% geostrophic velocity of exp.2
for vv=1:length(velp_exp2(1,:))
    mnid=max(find(~isnan(velp_exp2(:,vv))));
    if isempty(mnid)
        idm(vv)=1;
    else
        idm(vv)=mnid;
    end
    gvel_exp2(:,vv)=velp_exp2(:,vv)-velp_exp2(idm(vv),vv);
end
% geostrophic velocity of exp.3
for vv=1:length(velp_exp3(1,:))
    mnid=max(find(~isnan(velp_exp3(:,vv))));
    if isempty(mnid)
        idm(vv)=1;
    else
        idm(vv)=mnid;
    end
    gvel_exp3(:,vv)=velp_exp3(:,vv)-velp_exp3(idm(vv),vv);
end
% geostrophic velocity of exp.4
for vv=1:length(velp_exp4(1,:))
    mnid=max(find(~isnan(velp_exp4(:,vv))));
    if isempty(mnid)
        idm(vv)=1;
    else
        idm(vv)=mnid;
    end
    gvel_exp4(:,vv)=velp_exp4(:,vv)-velp_exp4(idm(vv),vv);
end


% % Example
% % (depth, station, time)
% % i=1;
% % for t=1:88
% %     ga(:,:,t)=sw_gpan(sals(:,:,t),temps(:,:,t),depths(:,:,t));
% %     velp(:,:,t) = sw_gvel(ga(:,:,t),squeeze(lats(1,:,t)),squeeze(lons(1,:,t)));
% % end
% % for vv=1:length(velp(1,:))
% %     mnid=max(find(~isnan(velp(:,vv))));
% %     if isempty(mnid)
% %         idm(vv)=1;
% %     else
% %         idm(vv)=mnid;
% %     end
% %     gvel(:,vv)=velp(:,vv)-velp(idm(vv),vv);
% % end