clear all; clc; close all;
ystart = [1950, 1960, 1970, 1980, 1990, 2000];
dstart = [0103, 0103, 0103, 0103, 0103, 0103];
yend   = [1959, 1969, 1979, 1989, 1999, 2009];
dend   = [1227, 1227, 1227, 1227, 1227, 0627];

% time=21909;
% abc=datestr(datenum(1900,1,1)+time,'yyyymmdd')
% [info, struc]=read_nc_file_struct('E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\CCAR_recon_sea_level_19500103_19591227_v1.nc')

% % read separated file

tind_start=1;
for i=1:length(ystart)
    [info, struc]=read_nc_file_struct ...
                (['E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\CCAR_recon_sea_level_', ...
                    num2str(ystart(i),'%04i'),num2str(dstart(i),'%04i'),'_', ...
                    num2str(yend(i),'%04i'),num2str(dend(i),'%04i'),'_v1.nc']);
    if(exist('lon','var')==0)
        lon=struc.lon;
        lat=struc.lat;
    end
    tind_end=tind_start+length(struc.time)-1;
    ssha(:,:,tind_start:tind_end)=struc.ssha;
    time(tind_start:tind_end)=struc.time;
    tind_start=tind_end+1;
end

% ssha=single(ssha);
time_str=datestr(datenum(1900,1,1)+time,'yyyymmdd');
year=str2num(time_str(:,1:4));
month=str2num(time_str(:,5:6));
mind_3=1;
for i=min(year):max(year)
    yind_found=find(year==i);
    clear year_2 month_2
    year_2(1:length(yind_found))=year(yind_found);
    month_2(1:length(yind_found))=month(yind_found);
    ssha_2(:,:,1:length(yind_found))=ssha(:,:,yind_found);
%     ssha_3(1:info.Dimensions(1).Length,1:info.Dimensions(2).Length,1:max(month_2))=0;
%     time_3(1:max(month_2))=0;
    for j=1:max(month_2)
        mind_found=find(month_2==j);
        ssha_3(:,:,mind_3+j-1)=mean(ssha_2(:,:,mind_found),3);
        time_3(mind_3+j-1)=datenum(i,j,15) - datenum(1900,1,1);
    end
    mind_3=mind_3+max(month_2);
end

trendtime_from_80 = 1980 + 1/12 : 1/12 : max(year)+1-6/12;  %% from 361 ~ 714
trend_80(1:info.Dimensions(1).Length,1:info.Dimensions(2).Length)=NaN;
trend2_80(1:info.Dimensions(1).Length,1:info.Dimensions(2).Length)=NaN;
for i=1:info.Dimensions(1).Length
    for j=1:info.Dimensions(2).Length
        p=polyfit(trendtime_from_80,squeeze(ssha_3(i,j,361:714))',1);
        trend_80(i,j)=p(1);
        trend2_80(i,j)=p(2);
    end
end

msl = mean(squeeze(mean(ssha_3,1,'omitnan')),1,'omitnan');



ncid = netcdf.create(strcat(['E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\CCAR_recon_sea_level_','merged','.nc']),'CLOBBER');

lat_dimid = netcdf.defDim(ncid,'lat',info.Dimensions(2).Length);
lon_dimid = netcdf.defDim(ncid, 'lon', info.Dimensions(1).Length);
time_dimid = netcdf.defDim(ncid, 'time', 0);

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' NOAA CSEOF Reconstructed Sea Level merged file ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', ' Sea Surface Height and Mean Sea Level (1950.1-2009.6), sla trends (1980.1-2009.6) ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'source', ' NOAA CSEOF Reconstructed Sea Level weekly data ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.Y.Kim');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', datestr(now));


lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
for i=1:length(info.Variables(1).Attributes)
    netcdf.putAtt(ncid,lonvarid, ...
        info.Variables(1).Attributes(i).Name,info.Variables(1).Attributes(i).Value);
end

latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
for i=1:length(info.Variables(2).Attributes)
    netcdf.putAtt(ncid,latvarid, ...
        info.Variables(2).Attributes(i).Name,info.Variables(2).Attributes(i).Value);
end

timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
% netcdf.putAtt(ncid,timevarid,'long_name','time');
% netcdf.putAtt(ncid,timevarid,'units','days since 1900-01-01 00:00:00');
% netcdf.putAtt(ncid,timevarid,'calendar','gregorian');
for i=1:4
    netcdf.putAtt(ncid,timevarid, ...
        info.Variables(3).Attributes(i).Name,info.Variables(3).Attributes(i).Value);
end

sshavarid = netcdf.defVar(ncid, 'ssha', 'NC_float', [lon_dimid lat_dimid time_dimid]);
for i=1:4
    netcdf.putAtt(ncid,sshavarid, ...
        info.Variables(4).Attributes(i).Name,info.Variables(4).Attributes(i).Value);
end

trendvarid = netcdf.defVar(ncid, 'trend', 'NC_float', [lon_dimid lat_dimid]);
netcdf.putAtt(ncid,trendvarid,'long_name','sea level rise trend from Jan 1980 to Jun 2009 (cm/year)');
for i=3:4
    netcdf.putAtt(ncid,trendvarid, ...
        info.Variables(4).Attributes(i).Name,info.Variables(4).Attributes(i).Value);
end

mslvarid = netcdf.defVar(ncid, 'msl', 'NC_float', time_dimid);
netcdf.putAtt(ncid,mslvarid,'long_name','global mean sea level');
for i=3:4
    netcdf.putAtt(ncid,mslvarid, ...
        info.Variables(4).Attributes(i).Name,info.Variables(4).Attributes(i).Value);
end

netcdf.endDef(ncid);

netcdf.putVar(ncid, timevarid, 0, length(time_3), time_3);
netcdf.putVar(ncid, lonvarid, 0, length(lon), lon);
netcdf.putVar(ncid, latvarid, 0, length(lat), lat);
netcdf.putVar(ncid, mslvarid, 0, length(msl), msl);
netcdf.putVar(ncid, sshavarid, [0 0 0], [info.Dimensions(1).Length,info.Dimensions(2).Length length(time_3)], ssha_3);
netcdf.putVar(ncid, trendvarid, [0 0], [info.Dimensions(1).Length,info.Dimensions(2).Length], trend_80);

netcdf.close(ncid);
