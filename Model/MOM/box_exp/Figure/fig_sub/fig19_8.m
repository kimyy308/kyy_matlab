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




addpath(genpath('D:\MEPL\Master_course\for_publish\seawater_ver3_2'));



for i=1:14
    lon3(i)=(lon2(i)+lon2(i+1))/2;
end

[lon44,ilevel44] = meshgrid(ilevel2(1:20),lon3(1:14));
pcolindex=1;

% % % 
% % % Exp1, geostrophic vel
% % % 

% % % 40N (j=60)
sbp1=subplot(2,2,1);
pcolor(ilevel44,lon44,gvel_exp1');
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(ilevel44,lon44,gvel_exp1','ShowText','on'); %% 
    h.LineColor='black';
vector_common5;
hold off;

sbp1=subplot(2,2,2);
pcolor(ilevel44,lon44,gvel_exp2');
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(ilevel44,lon44,gvel_exp2','ShowText','on'); %% 
h.LineColor='black';
[C2,h2]=contour(ilevel44,lon44,gvel_exp2',[-0.005 -0.005],'ShowText','on','LineStyle','--');
vector_common5;
hold off;

sbp1=subplot(2,2,3);
pcolor(ilevel44,lon44,gvel_exp3');
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(ilevel44,lon44,gvel_exp3','ShowText','on'); %% 
h.LineColor='black';
[C2,h2]=contour(ilevel44,lon44,gvel_exp3',[-0.005 -0.005],'ShowText','on','LineStyle','--');
vector_common5;
hold off;

sbp1=subplot(2,2,4);
pcolor(ilevel44,lon44,gvel_exp4');
shading interp;
caxis([-0.35,0.35]);
hold on;
[C,h]=contour(ilevel44,lon44,gvel_exp4','ShowText','on'); %% 
h.LineColor='black';
[C2,h2]=contour(ilevel44,lon44,gvel_exp4',[-0.005 -0.005],'ShowText','on','LineStyle','--');
h=colorbar;
caxis([-0.35,0.35]);
colormap(bwrmap);
set(h, 'Position', [.9150 .11 .0131 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'-0.3m/s', '0m/s', '0.3m/s'},...
            'YTick',[-0.3 0 0.3]);
h.Label.FontSize=7;
vector_common5;
hold off;


max_exp1=max(max(v2mean_1(1:15,60,1:20,1)));
max_exp2=max(max(v2mean_2(1:15,60,1:20,1)));
max_exp3=max(max(v2mean_3(1:15,40,1:20,1)));
max_exp4=max(max(v2mean_4(1:15,60,1:20,1)));

min_exp1=min(min(v2mean_1(1:15,60,1:20,1)));
min_exp2=min(min(v2mean_2(1:15,60,1:20,1)));
min_exp3=min(min(v2mean_3(1:15,40,1:20,1)));
min_exp4=min(min(v2mean_4(1:15,60,1:20,1)));

max_exp1_gvel=max(max(gvel_exp1));
max_exp2_gvel=max(max(gvel_exp2));
max_exp3_gvel=max(max(gvel_exp3));
max_exp4_gvel=max(max(gvel_exp4));

min_exp1_gvel=min(min(gvel_exp1));
min_exp2_gvel=min(min(gvel_exp2));
min_exp3_gvel=min(min(gvel_exp3));
min_exp4_gvel=min(min(gvel_exp4));