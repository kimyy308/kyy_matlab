clear all; clc; close all;

rootdir = '/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/cut_ES_daily/test2117';
testname='test2117';
[dropboxpath, error_status] = Func_0008_set_dropbox_path(computer);

    refpolygon = ...
       [127.0, 42.0;
        133.0, 42.0
        133.0, 35.0;
        128.0, 35.0];
lonlat(1)=min(refpolygon(:,1));
lonlat(2)=max(refpolygon(:,1));
lonlat(3)=min(refpolygon(:,2));
lonlat(4)=max(refpolygon(:,2));

for years = 1985:1987
    tempyear=years;
    yearstr=num2str(years, '%04i');
    for days=1:90
        tempday=days;
        daystr=num2str(days, '%04i');
        fname=[rootdir, '/', yearstr, '/', 'pck_ES_', testname, '_daily_', yearstr, '_', daystr, '.nc'];
        if years==1985 & days ==1
            lon_rho=ncread(fname, 'lon_rho');
            lat_rho=ncread(fname, 'lat_rho');
%             [refpolygon, section, error_status]=Func_0007_get_polygon_data_from_regionname('EKWC2');
            [indw, inde, inds, indn]=Func_0012_findind_Y(1/20,lonlat,lon_rho,lat_rho, 1);
            comb_v=NaN(inde-indw+1, indn-inds+1, 2014-1985+1, 90);
            comb_uwind=NaN(inde-indw+1, indn-inds+1, 2014-1985+1, 90);
            comb_vwind=NaN(inde-indw+1, indn-inds+1, 2014-1985+1, 90);            
        end
        comb_zeta(:,:,years-1984, days)=ncread(fname, 'zeta', [indw, inds, 1], [inde-indw+1, indn-inds+1, 1]);
        comb_v(:,:,years-1984, days)=ncread(fname, 'v', [indw, inds, 40, 1], [inde-indw+1, indn-inds+1, 1, 1]);
        comb_uwind(:,:,years-1984, days)=ncread(fname, 'Uwind', [indw, inds, 1], [inde-indw+1, indn-inds+1, 1]);
        comb_vwind(:,:,years-1984, days)=ncread(fname, 'Vwind', [indw, inds, 1], [inde-indw+1, indn-inds+1, 1]);
    end
end

[byrmap, error_status] = Func_0009_get_colormaps('byr', dropboxpath);
for loni=1:size(comb_v, 1)
    for lati=1:size(comb_v, 2)
        if isfinite(comb_v(loni,lati,1,1))
            temp_corr=corrcoef(comb_v(loni,lati,1,1:90), comb_uwind(loni,lati,1,1:90));
            corr_v_uwin(loni,lati)=temp_corr(1,2);
        end
    end
end
pcolor(corr_v_uwin'); shading flat; colormap(byrmap); colorbar; caxis([-0.8 0.8]);


for loni=1:size(comb_v, 1)
    for lati=1:size(comb_v, 2)
        if isfinite(comb_v(loni,lati,1,1))
            temp_corr=corrcoef(comb_zeta(loni,lati,1,1:90), comb_uwind(loni,lati,1,1:90));
            corr_zeta_uwin(loni,lati)=temp_corr(1,2);
        end
    end
end
pcolor(corr_zeta_uwin'); shading flat; colormap(byrmap); colorbar; caxis([-0.7 0.7]);

pcolor(mean(comb_v(:,:,1,:),4)'); shading flat; colormap(byrmap); colorbar; caxis([-0.7 0.7]);
