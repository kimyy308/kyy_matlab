clear all

expname='oman_restore_04_18';

readname='D:\need_to_presentation\';
readname2='\197912\ocean_snap_197912.nc';
rname=strcat(readname,expname,readname2);
lon2=ncread(rname,'xt_ocean');
lat2=ncread(rname,'yt_ocean');

pi = 3.14159265358979323846;
radius = 6371.0e3;
radian = 180./pi;
omega  = 7.292e-5;
cp_ocean = 3989.24495292815;

readname='D:\need_to_presentation\';
readname2='\197912\ocean_gridinfo_197912.nc';
rname=strcat(readname,expname,readname2);
ht2=ncread(rname,'ht');

deg2m=radius/radian;
for i=1:120
    for j=1:110
        sin1(i,j)  = sin(lat2(j)*pi/180.0);
        cos1(i,j)  = cos(lat2(j)*pi/180.0);
        beta(i,j) = 2.0*omega*cos(lat2(j)*pi/180.0)/radius;
        f2(i,j)    = 2.0*omega*sin1(i,j) + beta(i,j);
    end
end

% [lat,lon] = meshgrid(lat2,lon2);
for i=1:120
    for j=1:110
        if (ht2(i,j)< -1.0e+5)
          ht2(i,j)=NaN;
          f2(i,j)=NaN;
        end
    end
end




%%
%% figure, coriolis frequency field
%%
 readname='D:\need_to_presentation\';
 readname2='\paper\cor_freq.tif';
filename=strcat(readname,expname,readname2);
    pcolor(lon2,lat2,(f2'));
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
%     caxis([0,2000]);
    c=colorbar;
%     c.Label.String= 'sec^-^1';
    c.Label.FontSize= 12;
    colormap jet;
    xlabel('Longitude (^oE)','fontsize',27);
    ylabel('Latitude (^oN)','fontsize',27);
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',15);
  hold off;
    fig=gcf;
    fig.InvertHardcopy='off';
    set(gca, 'color', [0.8,0.8,0.8]);
    saveas(gcf,filename,'tiff');
'making coriolis frequency figure is completed'


%%
%% figure , topography field
%%
 readname='D:\need_to_presentation\';
 readname2='\paper\topog.tif';
 filename=strcat(readname,expname,readname2);
    pcolor(lon2,lat2,(ht2'));
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
    caxis([0,2000]);
    c=colorbar;
    c.Label.String= 'Depth (M)';
    c.Label.FontSize= 12;
    colormap jet;
    xlabel('Longitude (^oE)','fontsize',27);
    ylabel('Latitude (^oN)','fontsize',27);
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',15);
  hold off;
    fig=gcf;
    fig.InvertHardcopy='off';
    set(gca, 'color', [0.8,0.8,0.8]);
    saveas(gcf,filename,'tiff');
'making topog figure is completed'