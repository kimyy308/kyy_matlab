function status=ROMS_surface_forcing_CMIP5_time_interp(grdfile, year, varname, tinterval, outfiledir, testname, CMIP5_name, scenario_name)
% % This function is based on MATLAB 2017a
% % Updated 17-Feb-2019 by Y.Y.Kim

    root_dir = ['/home/auto/ext_hdd/CMIP5/',varname,'/',scenario_name,'/interp/',CMIP5_name,'/'];
    filename = strcat(varname,'_interp_',CMIP5_name,'_',scenario_name,'_r1i1p1_',num2str(year,'%04i'),'.nc');
    atmfile  = strcat(root_dir,filename);
%     outfiledir = '/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/forcing_matlab/';
%     testname = 'nwp_1_20';
    
% % NorESM-M_hur_1976.nc
% % NorESM-M_psl_1976.nc
% % NorESM-M_rsds_1976.nc
% % NorESM-M_salinity_1976.nc
% % NorESM-M_tas_1976.nc
% % NorESM-M_temp_1976.nc
% % NorESM-M_uo_1976.nc
% % NorESM-M_uo_wind_1976.nc
% % NorESM-M_vo_1976.nc
% % NorESM-M_vo_wind_1976.nc
% % NorESM-M_zos_1976.nc

    disp([CMIP5_name, ' file is the ',atmfile])
    
    switch varname
        case 'tas'
            varname_input = varname;
            varname_output = 'Tair';
            long_name = [CMIP5_name, ' 2 meter Temperature'];
            units = 'Celsius';           
        case 'hur'
            varname_input = varname;
            varname_output = 'Qair';
            long_name = [CMIP5_name, ' Relative humidity'];
            units = 'Percentage';   
        case 'psl'
            varname_input = varname;
            varname_output = 'Pair';
            long_name = [CMIP5_name, ' Mean Sea Level pressure'];
            units = 'mbar';   
        case 'rsds'
            varname_input = varname;
            varname_output = 'swrad';
            long_name = [CMIP5_name, ' Short Wave RAdiation Downwards'];
            units = 'W/m^2';
        case 'ua'
            varname_input = varname;
            varname_output = 'Uwind';
            long_name = [CMIP5_name, ' 10 Meter zonal wind velocity'];
            units = 'm/s';
        case 'va'
            varname_input = varname;
            varname_output = 'Vwind';
            long_name = [CMIP5_name, ' 10 Meter meridional wind velocity'];
            units = 'm/s';
        otherwise
            '?'
            return
    end
    varname_time = strcat(varname_output,'_time');
    outfile = strcat(outfiledir,testname,'_',num2str(year,'%04i'),'_',varname_output,'_tinterp.nc');
    disp(['forcing file is the ',outfile])
    
%     lon_raw   = ncread(atmfile,'lon');
%     lat_raw   = ncread(atmfile,'lat');
    lon_raw   = ncread(atmfile,'lon');
    lat_raw   = ncread(atmfile,'lat');
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

%     lon_CMIP5 = ncread(atmfile,'lon', [lon_min(1)-1], [lon_max(1)-lon_min(1)+2]);
%     lat_CMIP5 = ncread(atmfile,'lat', [lat_min(1)-1], [lat_max(1)-lat_min(1)+2]); 
    lon_CMIP5 = double(ncread(atmfile,'lon', [lon_min(1)-1], [lon_max(1)-lon_min(1)+2]));
    lat_CMIP5 = double(ncread(atmfile,'lat', [lat_min(1)-1], [lat_max(1)-lat_min(1)+2])); 
    
    CMIP5_data_info = ncinfo(atmfile, varname_input); 
    raw_data = ncread(atmfile, varname_input, [lon_min(1)-1 lat_min(1)-1 1], [lon_max(1)-lon_min(1)+2 lat_max(1)-lat_min(1)+2 CMIP5_data_info.Size(3)]);  %% [x y t]
    
    switch varname
        case 'tas'
            raw_data = raw_data - 273.15;
        case 'psl'
            raw_data = raw_data ./ 100;
    end
    
    disp('reading raw file is completed')
    size_data=size(raw_data);
    
    [mesh_lon_CMIP5, mesh_lat_CMIP5]=meshgrid(lat_CMIP5 ,lon_CMIP5);
    summask=zeros(size(raw_data(:,:,1)));
    for ttttt=1:CMIP5_data_info.Size(3)
%     for ttttt=1:1
        if(isnan(sum(sum(raw_data(:,:,ttttt))))>0)
%             disp('abc')
            ismask=isnan(raw_data(:,:,ttttt));
            summask=summask+double(ismask);
            temp_raw_data=raw_data(:,:,ttttt);
            temp_raw_data(ismask)=griddata(mesh_lon_CMIP5(~ismask), ...
                mesh_lat_CMIP5(~ismask), double(temp_raw_data(~ismask)), ...
                mesh_lon_CMIP5(ismask), mesh_lat_CMIP5(ismask),'nearest');
            raw_data(:,:,ttttt)=temp_raw_data;
        end
    end
    disp('filling empty value is completed')

    data_daily = ones(size_data(1),size_data(2),totalday)*NaN;
    
% %     switch varname
% %         case 'ssrd'
% % % unit of swrad raw data is J/m^2; unit of swrad forcing data is W/m^2;
% % % these values are accumulated from the beginning of the forecast.
% % % (beginning of the forecast : 00:00, 12:00, 24:00, ...)
% % % (several peak value from the beginning : 12:00, 24:00, 36:00 ...)
% % % 1st day accumulated value = ssrd(12:00) + ssrd(24:00)
% % % convert Joules(for a day) to Watts --> divide by 86400s
% %             for i=1:1:totalday;
% %                 data_daily(:,:,i)=(raw_data(:,:,4+(i-1)*8)+raw_data(:,:,8+(i-1)*8))/(24*60*60);        
% %             end
% %         otherwise
% %             for i=1:1:totalday;
% %                 if (year >= 2017)
% %                     aa=sum(raw_data(:,:,1+8*(i-1):8*i),3);  %% sum 1~4, 5~8, ... (6-hour interval data) 
% %                     data_daily(:,:,i)=(aa)/8;
% %                 else
% %                     aa=sum(raw_data(:,:,1+4*(i-1):4*i),3);  %% sum 1~4, 5~8, ... (6-hour interval data) 
% %                     data_daily(:,:,i)=(aa)/4;
% %                 end
% %             end
% %     end  
    data_daily=raw_data;
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
        'source', [' CMIP5 ', CMIP5_name]);
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
    
    missvalue=NaN;
    if isnan(missvalue)
      ismask=isnan(data_daily(:,:,1));
    else
      ismask=(data_daily(:,:,1)==missvalue);
    end
    
    [lat_mesh, lon_mesh]=meshgrid(lat_CMIP5, lon_CMIP5);
    for i=1:totalday
        netcdf.putVar(ncid, timevarid, i-1, 1, time(i));
        Zi=squeeze(data_daily(:,:,i)); %%[x y]
        Zi(ismask)=griddata(lon_mesh(~ismask),lat_mesh(~ismask),double(Zi(~ismask)), ...
                         lon_mesh(ismask),lat_mesh(ismask),'nearest');
        Z=griddata(lon_CMIP5,lat_CMIP5,double(Zi'),lon_rho(:,1),lat_rho(1,:)); %% [x y v'(y,x) xv yv]
        
        netcdf.putVar(ncid, dvarid, [0 0 i-1], [data_info.Size(1) data_info.Size(2) 1], Z'); %%[x y t], Z'(x,y)
    end

    netcdf.close(ncid);
    status=1;
end
