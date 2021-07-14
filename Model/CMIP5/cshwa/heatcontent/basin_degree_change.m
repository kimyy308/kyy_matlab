%1 degree -> 2 degree
%load data
load('basin.mat');

%lon, lat setting
lon1 = [0.5:1:360];
lat1 = [-89.5:90];
[lon1,lat1] = meshgrid(lon1,lat1);

lon2 = [0.5:2:360];
lat2 = [-89.5:2:90];
[lon2,lat2] = meshgrid(lon2,lat2);

%regrid data
basin2 = interp2(lon1,lat1,basin,lon2,lat2,'nearest');
basin = basin2;
save('basin_low.mat','basin');