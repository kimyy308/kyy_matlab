clc; clear all; close all;

minyear=1974;
maxyear=1977;
list = dir('*uo_interped_nor_*')

uo_sd_1 = load(list(2).name,'uo_sd');
uo_sd_2 = load(list(3).name,'uo_sd');
uo_sd_3 = load(list(4).name,'uo_sd');
uo_sd_4 = load(list(5).name,'uo_sd');
uo_sd_5 = load(list(6).name,'uo_sd');
uo_sd_6 = load(list(7).name,'uo_sd');
uo_sd_7 = load(list(8).name,'uo_sd');
uo_sd_8 = load(list(9).name,'uo_sd');
uo_sd_1= uo_sd_1.uo_sd;
uo_sd_2= uo_sd_2.uo_sd;
uo_sd_3= uo_sd_3.uo_sd;
uo_sd_4= uo_sd_4.uo_sd;
uo_sd_5= uo_sd_5.uo_sd;
uo_sd_6= uo_sd_6.uo_sd;
uo_sd_7= uo_sd_7.uo_sd;
uo_sd_8= uo_sd_8.uo_sd;
load('standard_grid_ocean.mat');

[size11 size12 size13 size14]= size(uo_sd_1);

%%% concatenate array along specified dimension
uo_sd = cat(4, uo_sd_1, uo_sd_2, uo_sd_3, uo_sd_4,uo_sd_5,uo_sd_6,uo_sd_7,uo_sd_8); %combine matrix into one

t_sd=size(uo_sd,4)
size(uo_sd)

timetime=minyear:1/365:maxyear+1-1/365;

tstep=0;
% for j = minyear:maxyear
for j = minyear:maxyear
tstep =tstep +1;
fstep=(tstep-1)*12+1:1:tstep*12;
ftime=timetime(fstep); % slice every year in time
uo_temp=uo_sd(:,:,:,fstep); % slice every year in uo

ncid = netcdf.create(strcat('../forcing_Nor/NorESM-M_ocean_uo_',num2str(j),'_v2.nc'),'CLOBBER');
lon_dimid = netcdf.defDim(ncid,'lon',size11);
lat_dimid = netcdf.defDim(ncid, 'lat', size12);
dep_dimid = netcdf.defDim(ncid, 'depth', size13);
time_dimid = netcdf.defDim(ncid, 'time',12);  

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' NorESM1-M interpolated data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', 'sea_water_x_velocity');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'uource', 'CMIP5 data');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by S.H Chae');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
netcdf.putAtt(ncid,lonvarid,'long_name','lon');
netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
netcdf.putAtt(ncid,latvarid,'standard_name','lat');
netcdf.putAtt(ncid,latvarid,'long_name','lat');
netcdf.putAtt(ncid,latvarid,'units','degrees_north');

depvarid = netcdf.defVar(ncid, 'depth', 'NC_float', dep_dimid);
netcdf.putAtt(ncid,depvarid,'standard_name','ocean depth');
netcdf.putAtt(ncid,depvarid,'long_name','depth_below_ocean_surface');
netcdf.putAtt(ncid,depvarid,'units','meter');

timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
netcdf.putAtt(ncid,timevarid,'standard_name','time');
netcdf.putAtt(ncid,timevarid,'long_name','time');
netcdf.putAtt(ncid,timevarid,'units','months from 1974-jan');

uovarid = netcdf.defVar(ncid, 'uo', 'NC_float', [lon_dimid lat_dimid dep_dimid time_dimid]);
netcdf.putAtt(ncid,uovarid,'standard_name','uo');
netcdf.putAtt(ncid,uovarid,'long_name','sea_water_x_velocity');
netcdf.putAtt(ncid,uovarid,'units','m/s');

netcdf.endDef(ncid);

netcdf.putVar(ncid, lonvarid, 0, length(squeeze(lon_ec_f(:,1,1))), squeeze(lon_ec_f(:,1,1)));
netcdf.putVar(ncid, latvarid, 0, length(squeeze(lat_ec_f(1,:,1))), squeeze(lat_ec_f(1,:,1)));
netcdf.putVar(ncid, depvarid, 0, length(squeeze(stan_dep(1,1,:))), squeeze(stan_dep(1,1,:)));
netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
netcdf.putVar(ncid, uovarid, [0 0 0 0], [size11 size12 size13 12], uo_temp);

netcdf.close(ncid);
end