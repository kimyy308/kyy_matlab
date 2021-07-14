function status=ROMS_surface_forcing_NOMADS_total(grdfile, year, varname, tinterval, outfiledir, testname, atm_name)
% % This function is based on MATLAB 2017a
% % Updated 14-Feb-2019 by Y.Y.Kim

    root_dir = ['/data1/stlee/ext_hdd/NOMADS/',num2str(year,'%04i'),'/nc_file/'];

%     outfiledir = '/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/forcing_matlab/';
%     testname = 'nwp_1_20';

% % Example of filename
% % 
% % tmp2m.gdas.199301.nc
% % pressfc.gdas.199301.nc
% % q2m.gdas.199301.nc
% % dswsfc.gdas.199301.nc
% % wnd10m.gdas.199301.nc

% % from varname, define filename, units,
% % varname in the file, varname in the ROMS forcing file,
for month = 1:12
    switch varname
        case 'TMP_2maboveground'
            varname_file = 'tmp2m';
            filename = strcat(varname_file, '.', atm_name, '.', ...
                num2str(year,'%04i'), num2str(month,'%02i'), '.nc');
            atmfile  = strcat(root_dir, filename);
            disp([atm_name, ' file is the ',atmfile]);
            varname_input = varname;
            varname_output = 'Tair';
            long_name = [atm_name, ' 2 meter Temperature'];
            units = 'Celsius';           
        case 'SPFH_2maboveground'
            varname_file = 'q2m';
            filename = strcat(varname_file, '.', atm_name, '.', ...
                num2str(year,'%04i'), num2str(month,'%02i'), '.nc');
            atmfile  = strcat(root_dir, filename);
            disp([atm_name, ' file is the ',atmfile]);
            varname_input = varname;
            varname_output = 'Qair';
            long_name = [atm_name, ' Relative humidity'];
            units = 'Percentage';   
        case 'PRES_surface'
            varname_file = 'pressfc';
            filename = strcat(varname_file, '.', atm_name, '.', ...
                num2str(year,'%04i'), num2str(month,'%02i'), '.nc');
            atmfile  = strcat(root_dir, filename);
            disp([atm_name, ' file is the ',atmfile]);
            varname_input = varname;
            varname_output = 'Pair';
            long_name = [atm_name, ' Mean Sea Level pressure'];
            units = 'mbar';   
        case 'DSWRF_surface'
            varname_file = 'dswsfc';
            filename = strcat(varname_file, '.', atm_name, '.', ...
                num2str(year,'%04i'), num2str(month,'%02i'), '.nc');
            atmfile  = strcat(root_dir, filename);
            disp([atm_name, ' file is the ',atmfile]);
            varname_input = varname;
            varname_output = 'swrad';
            long_name = [atm_name, ' Short Wave RAdiation Downwards'];
            units = 'W/m^2';
        case 'UGRD_10maboveground'
            varname_file = 'Uwnd10m';
            filename = strcat(varname_file, '.', atm_name, '.', ...
                num2str(year,'%04i'), num2str(month,'%02i'), '.nc');
            atmfile  = strcat(root_dir, filename);
            disp([atm_name, ' file is the ',atmfile]);
            varname_input = varname;
            varname_output = 'Uwind';
            long_name = [atm_name, ' 10 Meter zonal wind velocity'];
            units = 'm/s';
        case 'VGRD_10maboveground'
            varname_file = 'Vwnd10m';
            filename = strcat(varname_file, '.', atm_name, '.', ...
                num2str(year,'%04i'), num2str(month,'%02i'), '.nc');
            atmfile  = strcat(root_dir, filename);
            disp([atm_name, ' file is the ',atmfile]);
            varname_input = varname;
            varname_output = 'Vwind';
            long_name = [atm_name, ' 10 Meter meridional wind velocity'];
            units = 'm/s';
        otherwise
% %             If unknown varname is used, return
            '?'
            return
    end
    
    varname_time = 'time';
    outfile = strcat(outfiledir,testname,'_',num2str(year,'%04i'),'_',varname_output,'.nc');
    disp(['ROMS forcing file is the ',outfile])
    
    lon_raw   = ncread(atmfile,'longitude');
    lat_raw   = ncread(atmfile,'latitude');
    totalday  = yeardays(year); 
    
    lon_rho = ncread(grdfile,'lon_rho');
    lat_rho = ncread(grdfile,'lat_rho');
    data_info = ncinfo(grdfile, 'lon_rho'); 
    
    lon_grd_max = lon_rho(data_info.Size(1),data_info.Size(2));
    lat_grd_max = lat_rho(data_info.Size(1),data_info.Size(2));
    lon_grd_min = lon_rho(1,1);
    lat_grd_min = lat_rho(1,1);

    time=[0.5:1:totalday-0.5];

    lon_west = abs(lon_raw - (lon_grd_min-1));
    min_lon_west=min(lon_west);
    lon_east = abs(lon_raw - (lon_grd_max+1));
    min_lon_east=min(lon_east);
    lat_south = abs(lat_raw - (lat_grd_min-1));
    min_lat_south=min(lat_south);
    lat_north = abs(lat_raw - (lat_grd_max+1));
    min_lat_north=min(lat_north);

    lon_min = find(lon_west == min_lon_west);
    lon_max = find(lon_east == min_lon_east);
    lat_min = find(lat_south == min_lat_south);
    lat_max = find(lat_north == min_lat_north);
    
    if (year==2011 && month ==7)  %% because grid is differenct since Aug, 2011 in pressure
        lon_NOMADS_old = lon_NOMADS;
        lat_NOMADS_old = lat_NOMADS;
    elseif(year==2011 && month >7)
        lon_min=lon_min-1;
        lom_max=lon_max+1;
        lat_min=lat_min-1;
        lat_max=lat_max+1;
    end
    lon_NOMADS = ncread(atmfile,'longitude', [lon_min(1)-1], [lon_max(1)-lon_min(1)+2]);
    lat_NOMADS = ncread(atmfile,'latitude', [lat_min(1)-1], [lat_max(1)-lat_min(1)+2]); 
    
    NOMADS_data_info = ncinfo(atmfile, varname_input); 
    raw_data = ncread(atmfile, varname_input, [lon_min(1)-1 lat_min(1)-1 1], ...
        [lon_max(1)-lon_min(1)+2 lat_max(1)-lat_min(1)+2 NOMADS_data_info.Size(3)]);  %% [x y t], get data except last day's 24:00 data.
%     pcolor(squeeze(raw_data(:,:,1))')
    if (year>=1994) %% reshape raw_data to easy daily average
        monthday = (NOMADS_data_info.Size(3))/8;
        raw_data2 = reshape(raw_data, ...
            [size(raw_data,1), size(raw_data,2), 8, monthday]);
    else
        monthday = (NOMADS_data_info.Size(3))/24;
        raw_data2 = reshape(raw_data, ...
            [size(raw_data,1), size(raw_data,2), 24, monthday]);
    end
    raw_data = squeeze(mean(raw_data2,3));
    
    switch varname
        case 'TMP_2maboveground'
            raw_data = raw_data - 273.15;
        case 'PRES_surface'
            raw_data = raw_data ./ 100;
        case 'SPFH_2maboveground'
            sphum = raw_data;
            t2mfilename = strcat(root_dir,'tmp2m', '.', atm_name, '.', ...
                num2str(year,'%04i'), num2str(month,'%02i'), '.nc');
            t2m = ncread(t2mfilename, 'TMP_2maboveground', [lon_min(1)-1 lat_min(1)-1 1], ...
                [lon_max(1)-lon_min(1)+2 lat_max(1)-lat_min(1)+2 NOMADS_data_info.Size(3)]) - 273.15;  %% [x y t], K -> Celcius degree
            spfilename = strcat(root_dir,'pressfc', '.', atm_name, '.', ...
                num2str(year,'%04i'), num2str(month,'%02i'), '.nc');
            sp = ncread(spfilename, 'PRES_surface', [lon_min(1)-1 lat_min(1)-1 1], ...
                [lon_max(1)-lon_min(1)+2 lat_max(1)-lat_min(1)+2 NOMADS_data_info.Size(3)]) / 100.0;  %% [x y t], Pa -> mbar
            if (year>=1994)
                monthday = (NOMADS_data_info.Size(3))/8;
                t2m2 = reshape(t2m, ...
                    [size(t2m,1), size(t2m,2), 8, monthday]);
                sp2 = reshape(sp, ...
                    [size(sp,1), size(sp,2), 8, monthday]);
            else
                monthday = (NOMADS_data_info.Size(3))/24;
                t2m2 = reshape(t2m, ...
                    [size(t2m,1), size(t2m,2), 24, monthday]);
                sp2 = reshape(sp, ...
                    [size(sp,1), size(sp,2), 24, monthday]);
            end
            t2m = squeeze(mean(t2m2,3));
            sp = squeeze(mean(sp2,3));
            
            
            for humi=1:size(sphum,1)
                for humj=1:size(sphum,2)
                    for humt=1:size(sphum,3)
                        raw_data(humi,humj,humt)=qair2rh(sphum(humi,humj,humt), t2m(humi,humj,humt), sp(humi,humj,humt)) * 100; %% convert 1 to 100 percentage
                    end
                end
            end
    end
    
    size_data=size(raw_data);
    if (month ==1)
        data_daily = ones(size_data(1),size_data(2),totalday)*NaN;
        data_daily=raw_data;
        clear raw_data;
    else
        if (size(data_daily,2)==size(raw_data,2))
            data_daily(:,:,end+1:end+1+monthday-1)=raw_data;
            clear raw_data;
        else
            for interpi=1:monthday
                raw_data_interpolated(:,:,interpi) = griddata(double(lon_NOMADS),double(lat_NOMADS),double(raw_data(:,:,interpi)'),double(lon_NOMADS_old),double(lat_NOMADS_old)')';
            end
            data_daily(:,:,end+1:end+1+monthday-1)=raw_data_interpolated;
            clear raw_data;
            clear raw_data_interpolated
        end
    end
end   

% %  downward shortwave radiations in the ECMWF Interim data have unit of J/m^2.
% %     switch varname
% %         case 'ssrd'
% % % unit of swrad raw data is J/m^2; unit of swrad forcing data is W/m^2;
% % % these values are accumulated from the beginning of the forecast.
% % % (beginning of the forecast : 00:00, 12:00, 24:00, ...)
% % % (several peak value from the beginning : 12:00, 24:00, 36:00 ...)
% % % 1st day accumulated value = ssrd(12:00) + ssrd(24:00)
% % % convert Joules(for a day) to Watts --> divide by 86400s

% %  But dsrsfc in NOMADS is already W/m^2, 
% % So we didn't have to convert this. Just daily mean
    
ncid = netcdf.create(outfile,'CLOBBER');

eta_rho_dimid = netcdf.defDim(ncid,'eta_rho',data_info.Size(2));
xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', data_info.Size(1));
eta_u_dimid = netcdf.defDim(ncid,'eta_u',data_info.Size(2));
xi_u_dimid = netcdf.defDim(ncid, 'xi_u', data_info.Size(1)-1);
eta_v_dimid = netcdf.defDim(ncid,'eta_v',data_info.Size(2)-1);
xi_v_dimid = netcdf.defDim(ncid, 'xi_v', data_info.Size(1));
eta_psi_dimid = netcdf.defDim(ncid,'eta_psi',data_info.Size(2)-1);
xi_psi_dimid = netcdf.defDim(ncid, 'xi_psi', data_info.Size(1)-1);
time_dimid = netcdf.defDim(ncid, varname_time, 0);

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'type', ' ROMS Surface forcing file ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'title', ' Bulk Formular Forcing file ');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'source', [' NOMADS ', atm_name]);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'author', 'Created by Y.Y.Kim');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
    'date', date);

timevarid=netcdf.defVar(ncid, varname_time, 'NC_DOUBLE', time_dimid);
netcdf.putAtt(ncid,timevarid,'long_name','cyclic day');
netcdf.putAtt(ncid,timevarid,'units','DAYS');
netcdf.putAtt(ncid,timevarid,'cycle_length', totalday);

dvarid=netcdf.defVar(ncid,varname_output, 'NC_FLOAT', [xi_rho_dimid eta_rho_dimid time_dimid]);  %% [x y t]
netcdf.putAtt(ncid,dvarid,'long_name',long_name);
netcdf.putAtt(ncid,dvarid,'units',units);
netcdf.putAtt(ncid,dvarid,'time',varname_time);

netcdf.endDef(ncid);

if (year==2011)  %% because grid is differenct since Aug, 2011 in pressure
    for i=1:totalday
        netcdf.putVar(ncid, timevarid, i-1, 1, time(i));
        Zi=squeeze(data_daily(:,:,i)); %%[x y]
        Z=griddata(double(lon_NOMADS_old),double(lat_NOMADS_old),double(Zi'),lon_rho(:,1),lat_rho(1,:)); %% [x y v'(y,x) xv yv]

        netcdf.putVar(ncid, dvarid, [0 0 i-1], [data_info.Size(1) data_info.Size(2) 1], Z'); %%[x y t], Z'(x,y)
    end
else
    for i=1:totalday
        netcdf.putVar(ncid, timevarid, i-1, 1, time(i));
        Zi=squeeze(data_daily(:,:,i)); %%[x y]
        Z=griddata(double(lon_NOMADS),double(lat_NOMADS),double(Zi'),lon_rho(:,1),lat_rho(1,:)); %% [x y v'(y,x) xv yv]

        netcdf.putVar(ncid, dvarid, [0 0 i-1], [data_info.Size(1) data_info.Size(2) 1], Z'); %%[x y t], Z'(x,y)
    end
end

netcdf.close(ncid);
status=1;
end

% function rh=qair2rh(qair, temp, press = 1013.25)
function rh=qair2rh(qair, temp, press)

% % https://github.com/PecanProject/pecan/blob/master/modules/data.atmosphere/R/metutils.R
% % from this R code

% % ##' Convert specific humidity to relative humidity
% % ##'
% % ##' converting specific humidity into relative humidity
% % ##' NCEP surface flux data does not have RH
% % ##' from Bolton 1980 The computation of Equivalent Potential Temperature 
% % ##' \url{http://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html}
% % ##' @title qair2rh
% % ##' @param qair specific humidity, dimensionless (e.g. kg/kg) ratio of water mass / total air mass
% % ##' @param temp degrees C
% % ##' @param press pressure in mb
% % ##' @return rh relative humidity, ratio of actual water mixing ratio to saturation mixing ratio
% % ##' @export
% % ##' @author David LeBauer  


% %     es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
% %     e <- qair * press / (0.378 * qair + 0.622)
% %     rh <- e / es
% %     rh[rh > 1] <- 1
% %     rh[rh < 0] <- 0
% %     return(rh)
    
    es = 6.112 * exp((17.67 * temp)/(temp + 243.5));
    e = qair * press / (0.378 * qair + 0.622);
    rh = e/es;
    if (rh>=1)
        rh = 1;
    elseif (rh<=0)
        rh = 0;
    end
end