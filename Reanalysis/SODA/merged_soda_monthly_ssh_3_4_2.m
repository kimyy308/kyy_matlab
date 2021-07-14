clc; clear all; close all;
% % Updated 26-Mar-2018 by Y.Y.Kim
% % Updated 27-Mar-2018 by Y.Y.Kim

sodadir='/data2/kimyy/Reanalysis/SODA/SODA_3_4_2/';
sodaprefix='soda3.4.2_mn_ocean_reg_';
avisodir='/data2/kimyy/Observation/AVISO/monthly/';
avisoprefix='dt_global_allsat_msla_h_y';

minyear=1993;
maxyear=2012;

i=1;
lon_aviso = double(ncread(strcat(avisodir,avisoprefix,num2str(minyear,'%04i'),'_m',num2str(12,'%02i'),'.nc'),'lon'));
lat_aviso = double(ncread(strcat(avisodir,avisoprefix,num2str(minyear,'%04i'),'_m',num2str(12,'%02i'),'.nc'),'lat'));
lon = ncread(strcat(sodadir,sodaprefix,num2str(2015,'%4i'),'.nc'),'xt_ocean');
lat = ncread(strcat(sodadir,sodaprefix,num2str(2015,'%4i'),'.nc'),'yt_ocean');
for year=minyear:maxyear
    ssh(:,:,:,year-minyear+1)=ncread(strcat(sodadir,sodaprefix,num2str(year,'%4i'),'.nc'),'ssh');
    for month=1:12
        ftime(i) = datenum(year,month,15) - datenum(1900,12,31);
        i=i+1;
        tempssha=ncread(strcat(avisodir,avisoprefix,num2str(year,'%04i'),'_m',num2str(month,'%02i'),'.nc'),'sla');
        ssha_aviso(:,:,month,year-minyear+1)=griddata(lon_aviso,lat_aviso,tempssha',lon,lat'); %% [x y v'(y,x) xv yv]
    end
    year
end
ssha_aviso=permute(ssha_aviso,[2 1 3 4]);
ssh(find(ssh<-10000))=NaN; %% unit : m
for year=minyear:maxyear
    for month=1:12
        ssha_soda(:,:,month,year-minyear+1)=ssh(:,:,month,year-minyear+1)-mean(mean(ssh,4),3);
    end
end

% % time_data_info = ncinfo(strcat(sodadir,sodaprefix,num2str(year,'%4i'),'.nc'), 'time');
lon_data_info = ncinfo(strcat(sodadir,sodaprefix,num2str(year,'%4i'),'.nc'), 'xt_ocean');
lat_data_info = ncinfo(strcat(sodadir,sodaprefix,num2str(year,'%4i'),'.nc'), 'yt_ocean');
ssh_data_info = ncinfo(strcat(sodadir,sodaprefix,num2str(year,'%4i'),'.nc'), 'ssh');



trendtime=minyear:1/12:maxyear+1-1/12;
reshaped_ssha_soda=reshape(ssha_soda,[ssh_data_info.Size(1),ssh_data_info.Size(2),(maxyear-minyear+1)*12]);
trend(1:ssh_data_info.Size(1),1:ssh_data_info.Size(2))=NaN;
trend2(1:ssh_data_info.Size(1),1:ssh_data_info.Size(2))=NaN;
for i=1:ssh_data_info.Size(1)
    for j=1:ssh_data_info.Size(2)
        p=polyfit(trendtime,squeeze(reshaped_ssha_soda(i,j,:))',1);
        trend(i,j)=p(1);
        trend2(i,j)=p(2);
    end
end

for month=1:12 
    climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
end



msl = reshape(squeeze(mean(squeeze(mean(ssh,1,'omitnan')),1,'omitnan')),[(maxyear-minyear+1)*12 1]); %% unit : m
% reshaped_ssh = reshape(ssh,ssh_data_info.Size(1),ssh_data_info.Size(2),12,36); %% [lon lat month year]
% clim_ssh = mean(reshaped_ssh,4,'omitnan');
clim_ssh = mean(ssh,4,'omitnan');
% clear reshaped_ssh
diff_ssha = ssha_soda-ssha_aviso;
ssh(isnan(ssh))=ssh_data_info.Attributes(1,4).Value;
diff_ssha(isnan(diff_ssha))=ssh_data_info.Attributes(1,4).Value;
clim_ssh(isnan(clim_ssh))=ssh_data_info.Attributes(1,4).Value;

for year=minyear:maxyear
    rm_seasonal_ssh(:,:,:,year-minyear+1) = squeeze(ssh(:,:,:,year-minyear+1)) - clim_ssh(:,:,:);
end


ncid = netcdf.create(strcat(sodadir,sodaprefix,'merged','.nc'),'CLOBBER');

lat_dimid = netcdf.defDim(ncid,'latitude',ssh_data_info.Size(2));
lon_dimid = netcdf.defDim(ncid, 'longitude', ssh_data_info.Size(1));
time_dimid = netcdf.defDim(ncid, 'time', 0);
clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);


netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' SODA 3.4.2 global monthly SSH merged file ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', ' Sea Surface Height, Mean Sea Level, Climatology SSH (1993-2012) ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'source', ' SODA 3.4.2 monthly average data with MOM5 driven by ERA-Interim forcing ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.Y.Kim');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
netcdf.putAtt(ncid,timevarid,'long_name','time');
netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
netcdf.putAtt(ncid,timevarid,'calendar','gregorian');

clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

lonvarid = netcdf.defVar(ncid, 'longitude', 'NC_float', lon_dimid);
for i=1:length(lon_data_info.Attributes)
    netcdf.putAtt(ncid,lonvarid, ...
        lon_data_info.Attributes(1,i).Name,lon_data_info.Attributes(1,i).Value);
end

latvarid = netcdf.defVar(ncid, 'latitude', 'NC_float', lat_dimid);
for i=1:length(lat_data_info.Attributes)
    netcdf.putAtt(ncid,latvarid, ...
        lat_data_info.Attributes(1,i).Name,lat_data_info.Attributes(1,i).Value);
end

sshvarid = netcdf.defVar(ncid, 'ssh', 'NC_float', [lon_dimid lat_dimid time_dimid]);
for i=1:length(ssh_data_info.Attributes)
    netcdf.putAtt(ncid,sshvarid, ...
        ssh_data_info.Attributes(1,i).Name,ssh_data_info.Attributes(1,i).Value);
end

diff_sshavarid = netcdf.defVar(ncid, 'diff_ssha', 'NC_float', [lon_dimid lat_dimid time_dimid]);
netcdf.putAtt(ncid,diff_sshavarid,'long_name','SODA - AVISO SSHA');
for i=2:length(ssh_data_info.Attributes)
    netcdf.putAtt(ncid,diff_sshavarid, ...
        ssh_data_info.Attributes(1,i).Name,ssh_data_info.Attributes(1,i).Value);
end

ssha_avisovarid = netcdf.defVar(ncid, 'aviso_ssha', 'NC_float', [lon_dimid lat_dimid time_dimid]);
netcdf.putAtt(ncid,ssha_avisovarid,'long_name','AVISO SSHA');
for i=2:length(ssh_data_info.Attributes)
    netcdf.putAtt(ncid,ssha_avisovarid, ...
        ssh_data_info.Attributes(1,i).Name,ssh_data_info.Attributes(1,i).Value);
end

ssha_sodavarid = netcdf.defVar(ncid, 'soda_ssha', 'NC_float', [lon_dimid lat_dimid time_dimid]);
netcdf.putAtt(ncid,ssha_sodavarid,'long_name','SODA SSHA');
for i=2:length(ssh_data_info.Attributes)
    netcdf.putAtt(ncid,ssha_sodavarid, ...
        ssh_data_info.Attributes(1,i).Name,ssh_data_info.Attributes(1,i).Value);
end

clim_sshvarid = netcdf.defVar(ncid, 'clim_ssh', 'NC_float', [lon_dimid lat_dimid clim_time_dimid]);
netcdf.putAtt(ncid,clim_sshvarid,'long_name','climatological mean Sea Surface Height');
for i=2:length(ssh_data_info.Attributes)-1
    netcdf.putAtt(ncid,clim_sshvarid, ...
        ssh_data_info.Attributes(1,i).Name,ssh_data_info.Attributes(1,i).Value);
end

trendvarid = netcdf.defVar(ncid, 'trend', 'NC_float', [lon_dimid lat_dimid]);
netcdf.putAtt(ncid,trendvarid,'long_name','sea level rise trend(m/year)');
for i=2:length(ssh_data_info.Attributes)-1
    netcdf.putAtt(ncid,trendvarid, ...
        ssh_data_info.Attributes(1,i).Name,ssh_data_info.Attributes(1,i).Value);
end

rm_seasonal_sshvarid = netcdf.defVar(ncid, 'rm_seasonal_ssh', 'NC_float', [lon_dimid lat_dimid time_dimid]);
netcdf.putAtt(ncid,rm_seasonal_sshvarid,'long_name','SSH - clim_SSH');
for i=2:length(ssh_data_info.Attributes)-1
    netcdf.putAtt(ncid,rm_seasonal_sshvarid, ...
        ssh_data_info.Attributes(1,i).Name,ssh_data_info.Attributes(1,i).Value);
end

mslvarid = netcdf.defVar(ncid, 'msl', 'NC_float', time_dimid);
netcdf.putAtt(ncid,mslvarid,'long_name','global mean sea level');
for i=2:4
    netcdf.putAtt(ncid,mslvarid, ...
        ssh_data_info.Attributes(1,i).Name,ssh_data_info.Attributes(1,i).Value);
end

netcdf.endDef(ncid);

netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
netcdf.putVar(ncid, lonvarid, 0, length(lon), lon);
netcdf.putVar(ncid, latvarid, 0, length(lat), lat);
netcdf.putVar(ncid, mslvarid, 0, length(msl), msl);
netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [ssh_data_info.Size(1) ssh_data_info.Size(2) length(climtime)], clim_ssh);
netcdf.putVar(ncid, sshvarid, [0 0 0], [ssh_data_info.Size(1) ssh_data_info.Size(2) length(ftime)], ssh);
netcdf.putVar(ncid, ssha_avisovarid, [0 0 0], [ssh_data_info.Size(1) ssh_data_info.Size(2) length(ftime)], ssha_aviso);
netcdf.putVar(ncid, ssha_sodavarid, [0 0 0], [ssh_data_info.Size(1) ssh_data_info.Size(2) length(ftime)], ssha_soda);
netcdf.putVar(ncid, trendvarid, [0 0], [ssh_data_info.Size(1) ssh_data_info.Size(2)], trend);
netcdf.putVar(ncid, diff_sshavarid, [0 0 0], [ssh_data_info.Size(1) ssh_data_info.Size(2) length(ftime)], diff_ssha);
netcdf.putVar(ncid, rm_seasonal_sshvarid, [0 0 0], [ssh_data_info.Size(1) ssh_data_info.Size(2) length(ftime)], rm_seasonal_ssh);

netcdf.close(ncid);