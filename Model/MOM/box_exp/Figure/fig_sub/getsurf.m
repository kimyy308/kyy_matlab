clear all; 
% clc; 
% close all;
expname='oman_restore_04_18';

readname='D:\need_to_presentation\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname,readname2);
temp2=ncread(rname,'temp');
lon2=ncread(rname,'xt_ocean');
lat2=ncread(rname,'yt_ocean');
level2=ncread(rname,'st_ocean');
ilevel2 = -level2;
[lat,lon,lev] = meshgrid(lat2,lon2,ilevel2');

readname='D:\need_to_presentation\';
readname2='\isodepth3.nc';
rname=strcat(readname,expname,readname2);
isod1=ncread(rname,'isodepth');
isod11=(0-isod1);


readname='D:\need_to_presentation\';
readname2='\isodepth4.nc';
rname=strcat(readname,expname,readname2);
isod5=ncread(rname,'isodepth');
isod=(0-isod5);
vort5=ncread(rname,'sfc_rel_vort');
f5=ncread(rname,'f');
% pot_vort5=ncread(rname,'pot_vort');
lon5=ncread(rname,'XT_OCEAN');
lat5=ncread(rname,'YT_OCEAN');

% aa = reshape(isod5(:,:,1),[120 110]);
aa=mean(isod(:,:,109:120),3)*0.01;
aa11=mean(isod11(:,:,109:120),3)*0.01;
bb=mean(isod(:,:,36),3);
i=115;
while i<120
    j=61;
    while j<70
        isod(i,j,120)=isod(114,j,120);
        isod11(i,j,120)=isod11(114,j,120);
        aa(i,j)=aa(114,j);
        aa11(i,j)=aa11(114,j);
        j=j+1
    end
    i=i+1
end

set(gca, 'color', [0.8,0.8,0.8]);
% surf(lon5,lat5,isod(:,:,120)')
% surf(lon5,lat5,aa')
surf(lon5(1:30),lat5(1:110),aa11(1:30,1:110)', 'FaceColor','red')
hold on
surf(lon5(1:30),lat5(1:110),aa(1:30,1:110)')
hold off
axis equal;
% patch(lon5,lat5,isod(:,:,120)','red')
colormap jet;
colormap(flipud(colormap))  % reversed jet;
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

xlabel('Longitude (^oE)','fontsize',15);
ylabel('Latitude (^oN)','fontsize',15);
zlabel('Depth (M)','fontsize',15);
