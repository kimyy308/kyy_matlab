clc;clear all;close all

% year=1994;
% 
% 
% for mon=1:1:12;
%     filename=['SODA_Y',num2str(year),'M',num2str(mon),'.cdf'];
%     nc=netcdf(filename);
%     lat_rho=nc{'latT'}(:);lon_rho=nc{'lonT'}(:);
%     lat_u=nc{'latU'}(:);lon_u=nc{'lonU'}(:);
%     lat_v=nc{'latV'}(:);lon_v=nc{'lonV'}(:);
%     
%     
%     
%     
% end

year=1994;mon=1;
bryname=['SODA_Y',num2str(year),'M',num2str(mon),'.cdf'];
grdname=['D:\add_ini_bry_grd\1993\roms_grid_ADD_03.nc'];
title=['making roms boundary file by SODA 2.1.6'];

romstools_param

create_bryfile(bryname,grdname,title,obc,...
                        theta_s,theta_b,hc,N,...
                        time,cycle,clobber);