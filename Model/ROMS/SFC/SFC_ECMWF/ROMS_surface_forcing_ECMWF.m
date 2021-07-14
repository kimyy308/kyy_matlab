function status=ROMS_surface_forcing_ECMWF(grdfile, year, varname, tinterval, outfiledir, testname)
% % This function is based on MATLAB 2017a
    root_dir = '/data1/temp/ECMWF_interim/';
    filename = strcat('ECMWF_Interim_',varname,'_',num2str(year,'%04i'),'.nc');
    atmfile  = strcat(root_dir,varname,'/',filename);
%     outfiledir = '/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/forcing_matlab/';
%     testname = 'nwp_1_20';
    
    disp(['ECMWF file is the ',atmfile])
    
    switch varname
        case 'airT'
            varname_input = 't2m';
            varname_output = 'Tair';
            long_name = 'ECMWF 2 meter Temperature';
            units = 'Celsius';           
        case 'dewt'
            varname_input = 'd2m';
            varname_output = 'Dair';
            long_name = 'ECMWF 2 meter Dewpoint temperature';
            units = 'Celsius';   
        case 'msl'
            varname_input = varname;
            varname_output = 'Pair';
            long_name = 'ECMWF Mean Sea Level pressure';
            units = 'mbar';   
        case 'ssrd'
            varname_input = varname;
            varname_output = 'swrad';
            long_name = 'ECMWF Short Wave RAdiation Downwards';
            units = 'W/m^2';
        case 'u10'
            varname_input = varname;
            varname_output = 'Uwind';
            long_name = 'ECMWF 10 Meter zonal wind velocity';
            units = 'm/s';
        case 'v10'
            varname_input = varname;
            varname_output = 'Vwind';
            long_name = 'ECMWF 10 Meter meridional wind velocity';
            units = 'm/s';
        otherwise
            '?'
            return
    end
    varname_time = strcat(varname_output,'_time');
    outfile = strcat(outfiledir,testname,'_',num2str(year,'%04i'),'_',varname_output,'.nc');
    disp(['forcing file is the ',outfile])
    
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

    lon_ECMWF = ncread(atmfile,'longitude', [lon_min(1)-1], [lon_max(1)-lon_min(1)+2]);
    lat_ECMWF = ncread(atmfile,'latitude', [lat_max(1)-1], [lat_min(1)-lat_max(1)+2]); %% lat -> descending series -> max to min
    
    ECMWF_data_info = ncinfo(atmfile, varname_input); 
    raw_data = ncread(atmfile, varname_input, [lon_min(1)-1 lat_max(1)-1 1], [lon_max(1)-lon_min(1)+2 lat_min(1)-lat_max(1)+2 ECMWF_data_info.Size(3)]);  %% [x y t]
    
    switch varname
        case 'airT'
            raw_data = raw_data - 273.15;
        case 'dewt'
            raw_data = raw_data - 273.15;
        case 'msl'
            raw_data = raw_data ./ 100;
    end
    
    size_data=size(raw_data);
    data_daily = ones(size_data(1),size_data(2),totalday)*NaN;
    
    switch varname
        case 'ssrd'
% unit of swrad raw data is J/m^2; unit of swrad forcing data is W/m^2;
% these values are accumulated from the beginning of the forecast.
% (beginning of the forecast : 00:00, 12:00, 24:00, ...)
% (several peak value from the beginning : 12:00, 24:00, 36:00 ...)
% 1st day accumulated value = ssrd(12:00) + ssrd(24:00)
% convert Joules(for a day) to Watts --> divide by 86400s
            for i=1:1:totalday;
                data_daily(:,:,i)=(raw_data(:,:,4+(i-1)*8)+raw_data(:,:,8+(i-1)*8))/(24*60*60);        
            end
        otherwise
            for i=1:1:totalday;
                if (year >= 2017)
                    aa=sum(raw_data(:,:,1+8*(i-1):8*i),3);  %% sum 1~4, 5~8, ... (6-hour interval data) 
                    data_daily(:,:,i)=(aa)/8;
                else
                    aa=sum(raw_data(:,:,1+4*(i-1):4*i),3);  %% sum 1~4, 5~8, ... (6-hour interval data) 
                    data_daily(:,:,i)=(aa)/4;
                end
            end
    end  
    clear raw_data;
   
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
        'source', ' ECMWF ERA Interim ');
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

    for i=1:totalday
        netcdf.putVar(ncid, timevarid, i-1, 1, time(i));
        Zi=squeeze(data_daily(:,:,i)); %%[x y]
        Z=griddata(double(lon_ECMWF),double(lat_ECMWF),Zi',lon_rho(:,1),lat_rho(1,:)); %% [x y v'(y,x) xv yv]
        
        netcdf.putVar(ncid, dvarid, [0 0 i-1], [data_info.Size(1) data_info.Size(2) 1], Z'); %%[x y t], Z'(x,y)
    end

    netcdf.close(ncid);
    status=1;
end
