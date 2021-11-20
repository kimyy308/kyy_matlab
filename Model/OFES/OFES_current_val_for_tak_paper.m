close all; clear all; clc;

% 1. make climate
% 2. 0-30m(1~15th), 30~bot(16-28th)
% 3. vec

ufilename='D:\Research\Ph_D_course\2021_YS_biogeochemical\OFES_uo_YS_2006_2019_mon.nc';
vfilename='D:\Research\Ph_D_course\2021_YS_biogeochemical\OFES_vo_YS_2006_2019_mon.nc';

udata_info=ncinfo(ufilename);
vdata_info=ncinfo(vfilename);

uo=ncread(ufilename, 'uvel')/100.0;
vo=ncread(vfilename, 'vvel')/100.0;
depth=ncread(ufilename, 'LEV');
longitude=ncread(ufilename, 'longitude');
latitude=ncread(ufilename, 'latitude');


for i=2:length(depth)-1
    layer_thickness(i)=(depth(i+1)-depth(i-1))/2;
end
layer_thickness(1)=depth(1)*2;
layer_thickness(length(depth))=layer_thickness(length(depth)-1);

uo_5d=reshape(uo, [size(uo,1), size(uo,2), size(uo,3), 12, size(uo,4)/12]);
clim_uo=mean(uo_5d,5);
vo_5d=reshape(vo, [size(uo,1), size(uo,2), size(uo,3), 12, size(uo,4)/12]);
clim_vo=mean(vo_5d,5);
% clim_uo(:,:,end,1)


% % depth weighted mean
surf_clim_uo=zeros([size(uo,1), size(uo,2), 12]);
surf_clim_vo=zeros([size(uo,1), size(uo,2), 12]);
for x=1:size(uo,1)
    for y=1:size(uo,2)
        valid_depthnum=sum(isfinite(squeeze(clim_uo(x,y,1:5,1))));
        if valid_depthnum>0
            for i=1:valid_depthnum
                surf_clim_uo(x,y,:)=squeeze(surf_clim_uo(x,y,:))+squeeze(clim_uo(x,y,i,:).*layer_thickness(i)./sum(layer_thickness(1:valid_depthnum)));
                surf_clim_vo(x,y,:)=squeeze(surf_clim_vo(x,y,:))+squeeze(clim_vo(x,y,i,:).*layer_thickness(i)./sum(layer_thickness(1:valid_depthnum)));
            end
        else
            surf_clim_uo(x,y,:)=NaN(1,12);
            surf_clim_vo(x,y,:)=NaN(1,12);
        end
    end
end

% % depth weighted mean
subsurf_clim_uo=zeros([size(uo,1), size(uo,2), 12]);
subsurf_clim_vo=zeros([size(uo,1), size(uo,2), 12]);
for x=1:size(uo,1)
    for y=1:size(uo,2)
        valid_depthnum=sum(isfinite(squeeze(clim_uo(x,y,6:end,1))));
        if valid_depthnum>0
            for i=6:6+valid_depthnum-1
                subsurf_clim_uo(x,y,:)=squeeze(subsurf_clim_uo(x,y,:))+squeeze(clim_uo(x,y,i,:).*layer_thickness(i)./sum(layer_thickness(6:6+valid_depthnum-1)));
                subsurf_clim_vo(x,y,:)=squeeze(subsurf_clim_vo(x,y,:))+squeeze(clim_vo(x,y,i,:).*layer_thickness(i)./sum(layer_thickness(6:6+valid_depthnum-1)));
            end
        else
            subsurf_clim_uo(x,y,:)=NaN(1,12);
            subsurf_clim_vo(x,y,:)=NaN(1,12);
        end
    end
end


[lat2, lon2]=meshgrid(latitude, longitude);

save('D:\Research\Ph_D_course\2021_YS_biogeochemical\OFES_2006_2019_ys_current.mat', 'clim_uo', 'clim_vo', 'depth', 'layer_thickness', 'surf_clim_uo', 'surf_clim_vo', 'subsurf_clim_uo', 'subsurf_clim_vo', 'longitude', 'latitude', 'lon2', 'lat2')

surf_clim_uo(10,10,:)=0.1;
surf_clim_vo(10,10,:)=0.0000001;



quivsize=7;
quiver(lon2, lat2, mean(surf_clim_uo(:,:,3:5),3).*quivsize, mean(surf_clim_vo(:,:,3:5),3).*quivsize, 'AutoScale', 'off')
quiver(lon2, lat2, mean(subsurf_clim_uo(:,:,3:5),3).*quivsize, mean(subsurf_clim_vo(:,:,3:5),3).*quivsize, 'AutoScale', 'off')

quiver(lon2, lat2, mean(surf_clim_uo(:,:,6:8),3).*quivsize, mean(surf_clim_vo(:,:,6:8),3).*quivsize, 'AutoScale', 'off')
quiver(lon2, lat2, mean(subsurf_clim_uo(:,:,6:8),3).*quivsize, mean(subsurf_clim_vo(:,:,6:8),3).*quivsize, 'AutoScale', 'off')

quiver(lon2, lat2, mean(surf_clim_uo(:,:,9:11),3).*quivsize, mean(surf_clim_vo(:,:,9:11),3).*quivsize, 'AutoScale', 'off')



%               uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
%                                 cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
%                                 u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
%                                 v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
%                                 'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);