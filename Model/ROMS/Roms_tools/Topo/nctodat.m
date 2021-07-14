clc;clear all;close all

nc=netcdf('etopo1.nc');

lat=nc{'lat'}(:);
lon=nc{'lon'}(:);
tt=size(lat);
fid=fopen('etopo1.dat','w');
tt2=size(lon);
for i=1:1:tt2(1);
    z=nc{'z'}(:,i);
    lon2=repmat(lon(i),tt(1),1);
    fprintf(fid,'%9.4f  %9.4f  %5d \r\n',[lat,lon2,z]');
end

fclose(fid);