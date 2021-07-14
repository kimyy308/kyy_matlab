close all; clear all; clc; 
    %%
    %% blue-white-red colormap
    %%
    i=1:20;
      bwrmap(i,1)= 0.;
      bwrmap(i,2)= 0.2:(0.3/19.):0.5;
      bwrmap(i,3)= 0.4:(0.6/19.):1.;
    
    i=21:49;
      bwrmap(i,1)= 0.:(1./28.):1.;
      bwrmap(i,2)= 0.5:(0.5/28.):1.;
      bwrmap(i,3)= 1.;   
    
    i=49:51;
      bwrmap(i,1)= 1.;
      bwrmap(i,2)= 1.;
      bwrmap(i,3)= 1.;
    
    i=51:56;
      bwrmap(i,1)= 1.;
      bwrmap(i,2)= 1.:(-0.1/5.):0.9;
      bwrmap(i,3)= 1.:(-0.1/5.):0.9;
    
    i=56:70;
      bwrmap(i,1)= 1.;
      bwrmap(i,2)= 0.9:(-0.45/14.):0.45;
      bwrmap(i,3)= 0.9:(-0.45/14.):0.45;
     
    i=70:80;
      bwrmap(i,1)= 1.;
      bwrmap(i,2)= 0.45:(-0.45/10.):0.;
      bwrmap(i,3)= 0.45:(-0.45/10.):0.;
     
    i=80:100;
      bwrmap(i,1)= 1.:(-0.6/20.):0.4;
      bwrmap(i,2)= 0.;
      bwrmap(i,3)= 0.;


filename = strcat('E:\Data\Model\CMIP5\NorESM1-M\nor_diff_thetao_rcp45_his_20y.nc');
info=ncinfo(filename);
lon=ncread(filename,'LON');
lat=ncread(filename,'LAT');
diff=ncread(filename,'DIFF_THETAO');
diff(diff<-100000)=NaN;
pcolor(lon(140:180,248:310)',lat(140:180,248:310)',diff(140:180,248:310)')
meanval=mean(mean(diff(140:180,248:310)','omitnan'),'omitnan')
shading interp;
xlabel('Lon(^oE)');
ylabel('Lat(^oN)');
colormap(bwrmap);
h=colorbar;
caxis([-4 4])
title(h,'^oC');
set(gca,'fontsize', 20);