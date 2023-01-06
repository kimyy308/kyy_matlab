clear all; clc; close all;

%% set path
tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

rootdir='/Volumes/kyy_raid/kimyy/Observation/ERSST/www.ncei.noaa.gov/pub/data/cmb/ersst/v5/netcdf';

load('a.mat')
% ex) ersst.v5.196001.nc

for year=1960:2020
    yearstr=num2str(year);
    for month=1:12
        monthstr=num2str(month, '%02i');
        filename=[rootdir, filesep, 'ersst.v5.',yearstr, monthstr,'.nc'];
%         ncinfo(filename)
        if year==1960 && month==1
            lon=ncread(filename, 'lon');
            lat=ncread(filename, 'lat');
        end
        ti=(year-1960)*12+month;
        sst(:,:,ti)=ncread(filename, 'sst');
    end
end

config_obs.staname = 'PDO';
config_obs.sta_lon = [105, 360-98]; % degrees west to 0~360 degrees east
config_obs.sta_lat = [20 70]; % degrees north

[lon_ind_w, lon_ind_e, lat_ind_s, lat_ind_n] = ...
    Func_0012_findind_Y(0.1, [config_obs.sta_lon, config_obs.sta_lat], ...
    lon, ...
    lat, 'CESM2'); % find valid lon, lat index near station

pdo_sst=sst(lon_ind_w:lon_ind_e,lat_ind_s:lat_ind_n,:);


[lat2, lon2]=meshgrid(lat,lon);

y=m_lldist([0 0], [0,2]);

for i=1:89
    x(i)=m_lldist([0,2], [-88+(i-1)*2, -88+i*2]);
end

area=repmat(x,180,1)*y;

area_pdo=area(lon_ind_w:lon_ind_e,lat_ind_s:lat_ind_n);
lon2_pdo=lon2(lon_ind_w:lon_ind_e,lat_ind_s:lat_ind_n);
lat2_pdo=lat2(lon_ind_w:lon_ind_e,lat_ind_s:lat_ind_n);

for i=1:732
    gmsst(i)=f_gm_var(sst(:,:,i),area);
    gmpdosst(i)=f_gm_var(pdo_sst(:,:,i),area_pdo);
end




% plot(gmsst)

re_gmsst=reshape(gmsst,[12,61]);
clim_gmsst=mean(re_gmsst,2);
re_gmpdosst=reshape(gmpdosst,[12,61]);
clim_gmpdosst=mean(re_gmpdosst,2);

for y=1:61
    for m=1:12
        ti=(y-1)*12+m;
        gmssta(ti)=gmsst(ti)-clim_gmsst(m);
        gmpdossta(ti)=gmpdosst(ti)-clim_gmpdosst(m);
    end
end
plot(gmssta);
hold on
plot(gmpdossta);
hold off


for i=1:732
    det_pdo_sst(:,:,i)=pdo_sst(:,:,i)-gmssta(i);
end

re_pdo_sst=reshape(det_pdo_sst, [80, 26, 12, 61]);
clim_pdo_sst=mean(re_pdo_sst, 4);
for y=1:61
    for m=1:12
        ti=(y-1)*12+m;
        pdo_ssta(:,:,ti)=det_pdo_sst(:,:,ti)-clim_pdo_sst(:,:,m);
    end
end

[lv, pc, var_exp] = Func_0024_EOF_3d(pdo_ssta,4);
ERSST.lv=lv;
ERSST.pc=pc;
ERSST.var_exp=var_exp;
ERSST.lon2_pdo=lon2_pdo;
ERSST.lat2_pdo=lat2_pdo;

plot(pc(:,1))
pcolor(-lv(:,:,1)'); shading flat; colormap jet

save('/Volumes/kyy_raid/Data/Observation/ERSST/EOF.mat', 'ERSST');


[r,p]=corrcoef(-a(:,1), pc(:,1))

bbb=lv(:,:,1);
bbb=bbb(isfinite(bbb));
b=b(isfinite(b));
[r,p]=corrcoef(bbb,b)

function gmsst = f_gm_var(var_2d, area)
    mask_var=NaN(size(var_2d));
    mask_var(isfinite(var_2d))=1;
    area=area.*mask_var;
    var_2d=squeeze(var_2d);
    var_2d_sum=var_2d.*area;
    var_2d_sum=sum(var_2d_sum(:), 'omitnan');
    gmsst=var_2d_sum ./ sum(area(:), 'omitnan');
end

