clear all; close all;
% clc; 
% close all;
expname='oman_restore_04_18';

readname='D:\need_to_presentation\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname,readname2);
temp2=ncread(rname,'temp');
u2=ncread(rname,'u');
v2=ncread(rname,'v');
lon2=ncread(rname,'xt_ocean');
lat2=ncread(rname,'yt_ocean');
level2=ncread(rname,'st_ocean');
ilevel2 = -level2;
[lat,lon,lev] = meshgrid(lat2,lon2,ilevel2');

readname='D:\need_to_presentation\';
readname2='\197912\ocean_diag_197912.nc';
rname=strcat(readname,expname,readname2);
w2=ncread(rname,'wt');


readname='D:\need_to_presentation\';
readname2='\isodepth_3dd.nc';
rname=strcat(readname,expname,readname2);
isod1=ncread(rname,'isodepth');
isod2=ncread(rname,'isodepth2');
isod3=ncread(rname,'isodepth3');
isod4=ncread(rname,'isodepth4');
isod11=(0-isod1);
isod22=(0-isod2);
isod33=(0-isod3);
isod44=(0-isod4);
% vort5=ncread(rname,'sfc_rel_vort');
% f5=ncread(rname,'f');
% pot_vort5=ncread(rname,'pot_vort');
lon5=ncread(rname,'XT_OCEAN');
lat5=ncread(rname,'YT_OCEAN');
% axis equal;
% aa = reshape(isod5(:,:,1),[120 110]);
aa1=mean(isod11(:,:,109:120),3); % NaN, 19.5
aa2=mean(isod22(:,:,109:120),3);
aa3=mean(isod33(:,:,109:120),3); % NaN, 10.0
aa4=mean(isod44(:,:,109:120),3);
% bb=mean(isod1(:,:,36),3);

uu2=mean(u2,4);
vv2=mean(v2,4);
ww2=mean(w2,4);

for i=1:25
    for j=1:110
        for k=1:20         
            if (uu2(i,j,k)< -1.0e+5)
                uu2(i,j,k)=NaN;
                vv2(i,j,k)=NaN;
                ww2(i,j,k)=NaN;
            end
        end
    end
end

for i=1:120
    for j=1:110
        if(isnan(aa1(i,j))==0)
            aa5(i,j)=-200.0;
        else
            aa5(i,j)=NaN;
        end
        if(isnan(aa3(i,j))==0)
            aa6(i,j)=-200.0;
        else
            aa6(i,j)=NaN;
        end
    end
end

for i=4:23
    for j=1:10
        if(isnan(aa5(i,j))==0)
            aa5(i,j)=-100;
            aa6(i,j)=-100;
        end
    end
end
for i=4:23
    for j=11:35
        if(isnan(aa5(i,j))==0)
            aa5(i,j)= -( 200 - 100*((35-j)/24.0) );
        end
        if(isnan(aa6(i,j))==0)
            aa6(i,j)= -( 200 - 100*((35-j)/24.0) );
        end
    end
end

[lat55,lon55,ilevel55] = meshgrid(lat5(1:110),lon5(4:23),ilevel2(1:20));

% i=115;
% while i<120
%     j=61;
%     while j<70
%         isod(i,j,120)=isod(114,j,120);
%         aa(i,j)=aa(114,j);
%         j=j+1
%     end
%     i=i+1
% end

% az= 37.5
% view(az, el);

% % % % % % % % % % 19.5 isosurface figure


% surf(lon5(4:23),lat5(1:110),aa2(4:23,1:110)', 'FaceColor','white')
% view([-0.2,-0.5,-0.1]);
view(-24,8);
% view(-154,21);
hold on
surf(lon5(4:23),lat5(1:110),aa1(4:23,1:110)')
shading interp;
% surf(lon5(4:23),lat5(1:110),aa3(4:23,1:110)')
% surf(lon5(4:23),lat5(1:110),aa5(4:23,1:110)', 'FaceColor','red')
% surf(lon5(4:23),lat5(1:110),aa5(4:23,1:110)', aa1(4:23,1:110)')
surf(lon5(4:23),lat5(1:110),aa5(4:23,1:110)','FaceColor',[0.8 0.8 0.8])
contour3(lon5(4:23),lat5(1:110),aa1(4:23,1:110)', 10, 'color', 'k', 'LineWidth',1)
zlim([-201 0]);
grid on
hold off
colormap jet;
colormap(flipud(colormap))
set(gca,'PlotBoxAspectRatio',[2.5 10 2])

readname='D:\need_to_presentation\';
readname2='\paper\upperlayer_surface.tif';
% readname2='\paper\upperlayer_surface2.tif';
filename=strcat(readname,expname,readname2);
set(gcf,'PaperPosition',[0 0 22 15]);
set(gca,'fontsize',15);
xlabel('Longitude (^oE)','fontsize',12);
ylabel('Latitude (^oN)','fontsize',12);
zlabel('Depth (M)','fontsize',12);
fig=gcf;
fig.InvertHardcopy='off';
c=colorbar;
c.Label.String= 'Depth (m)';
c.Label.FontSize= 12;
saveas(gcf,filename,'tiff');
% axis equal



% % % % % % % % % % % % 10.0 isosurface figure

figure;
% surf(lon5(4:23),lat5(1:110),aa4(4:23,1:110)', 'FaceColor','white')
view(-24,8);
% view(-154,21);
hold on
surf(lon5(4:23),lat5(1:110),aa3(4:23,1:110)')
shading interp;
% surf(lon5(4:23),lat5(1:110),aa6(4:23,1:110)', 'FaceColor','blue')
% surf(lon5(4:23),lat5(1:110),aa6(4:23,1:110)', aa3(4:23,1:110)')
surf(lon5(4:23),lat5(1:110),aa6(4:23,1:110)','FaceColor',[0.8 0.8 0.8])
contour3(lon5(4:23),lat5(1:110),aa3(4:23,1:110)', 10, 'color', 'k', 'LineWidth',1)
zlim([-201 0]);
grid on
% quiver3(lon55(1:4:20,10:4:110,1:9:20),lat55(1:4:20,10:4:110,1:9:20),ilevel55(1:4:20,10:4:110,1:9:20), ...
%     uu2(4:4:23,10:4:110,1:9:20)/2,vv2(4:4:23,10:4:110,1:9:20)/2,ww2(4:4:23,10:4:110,1:9:20)/2*100000000, 'color', 'black');

% quiver3(lon55(1:2:20,10:2:110,1:9:20),lat55(1:2:20,10:2:110,1:9:20),ilevel55(1:2:20,10:2:110,1:9:20), ...
%     uu2(4:2:23,10:2:110,1:9:20)/2,vv2(4:2:23,10:2:110,1:9:20)/2,ww2(4:2:23,10:2:110,1:9:20)/2*100000000, 'color', 'black');
% cccc=ww2/2*100000000;
% quiver3(lon55(1:2:20,10:2:110,1:2:20),lat55(1:2:20,10:2:110,1:2:20),ilevel55(1:2:20,10:2:110,1:2:20), ...
%     uu2(4:2:23,10:2:110,1:2:20)/2,vv2(4:2:23,10:2:110,1:2:20)/2,cccc(4:2:23,10:2:110,1:2:20), 'color', 'black');
hold off 
colormap jet;
colormap(flipud(colormap))
xlabel('Longitude (^oE)','fontsize',15);
ylabel('Latitude (^oN)','fontsize',15);
zlabel('Depth (M)','fontsize',15);
set(gca,'PlotBoxAspectRatio',[2.5 10 2])

readname='D:\need_to_presentation\';
readname2='\paper\lowerlayer_surface.tif';
% readname2='\paper\lowerlayer_surface2.tif';
filename=strcat(readname,expname,readname2);
set(gcf,'PaperPosition',[0 0 22 15]);
set(gca,'fontsize',15);
xlabel('Longitude (^oE)','fontsize',12);
ylabel('Latitude (^oN)','fontsize',12);
zlabel('Depth (M)','fontsize',12);
fig=gcf;
fig.InvertHardcopy='off';
c=colorbar;
c.Label.String= 'Depth (m)';
c.Label.FontSize= 12;
saveas(gcf,filename,'tiff');



%%% 19.5 + 10 + slope
% % figure;
% % view(-32,13);
% % hold on
% % surf(lon5(4:23),lat5(1:110),aa1(4:23,1:110)','EdgeColor','none')
% % surf(lon5(4:23),lat5(1:110),aa3(4:23,1:110)','EdgeColor','none')
% % % surf(lon5(4:23),lat5(1:110),aa6(4:23,1:110)', aa3(4:23,1:110)')
% % surf(lon5(4:23),lat5(1:110),aa6(4:23,1:110)', aa3(4:23,1:110)','FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
% % 
% % shading flat;
% % 
% % zlim([-201 0]);
% % grid on
% % hold off 
% % colormap jet;
% % colormap(flipud(colormap))
% % xlabel('Longitude (^oE)','fontsize',15);
% % ylabel('Latitude (^oN)','fontsize',15);
% % zlabel('Depth (M)','fontsize',15);
% % set(gca,'PlotBoxAspectRatio',[2.5 10 2])





% set(gca, 'color', [0.8,0.8,0.8]);
% surf(lon5,lat5,isod(:,:,120)')

% patch(lon5,lat5,isod(:,:,120)','red')
% colormap jet;
% colormap(flipud(colormap))  % reversed jet;
% caxis([0 -200]);
% colormap gray;
% surfc(lon5,lat5,isod(:,:,120)')
% mesh(lon5,lat5,isod(:,:,120)')
% meshz(lon5,lat5,isod(:,:,120)')
% hidden off;
% surf(lon5,lat5,isod(:,:,120)','FaceColor','red','EdgeColor','none')
% surf(lon5,lat5,isod(:,:,120)','EdgeColor','interp','EdgeColor','none','FaceLighting','phong')
% shading flat;
% shading interp;
% light;
% lighting phong;
% shading interp;
% camlight down; lighting phong;

% isosurface(temp2(:,:,:,1),15)

% [faces,verts,colors] = isosurface(lat2,lon2,ilevel2,temp2(:,:,:,1),10.0); 
% patch('Vertices', verts, 'Faces', faces, ... 
%     'FaceVertexCData', colors, ... 
%     'FaceColor','interp', ... 
%     'edgecolor', 'interp');
% axis vis3d;

