%-----------------------------------------------------------
%  Create empty river data NetCDF file
%
%   Called from latte_rivers.m
%-----------------------------------------------------------

disp('  ')
disp(['The RIVER netcdf file will be ' Fname])
disp('  ')

    
% ncid = netcdf.create(Fname, 'NETCDF4');
    
ncid = netcdf.create(Fname, '64BIT_OFFSET');

% dimensions
xi_u_dimid = netcdf.defDim(ncid, 'xi_u', size(grd.lon_u,2));
xi_v_dimid = netcdf.defDim(ncid, 'xi_v', size(grd.lon_v,2));
xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', size(grd.lon_rho,2));
eta_u_dimid = netcdf.defDim(ncid, 'eta_u', size(grd.lon_u,1));
eta_v_dimid = netcdf.defDim(ncid, 'eta_v', size(grd.lon_v,1));
eta_rho_dimid = netcdf.defDim(ncid, 'eta_rho', size(grd.lon_rho,1));
s_rho_dimid = netcdf.defDim(ncid, 's_rho', N);
river_dimid = netcdf.defDim(ncid, 'river', r);
time_dimid = netcdf.defDim(ncid, 'time', length(River.time));

% the variables
river_varid = netcdf.defVar(ncid, 'river', 'float', river_dimid);
netcdf.putAtt(ncid, river_varid, 'long_name', 'river identification number');
netcdf.putAtt(ncid, river_varid, 'units', 'non-dimensional');
netcdf.putAtt(ncid, river_varid, 'field', 'river, scalar');

river_name_varid = netcdf.defVar(ncid, 'river_name', 'float', river_dimid);
netcdf.putAtt(ncid, river_name_varid, 'long_name', 'river name');

river_Xposition_varid = netcdf.defVar(ncid, 'river_Xposition', 'float', river_dimid);
netcdf.putAtt(ncid, river_Xposition_varid, 'long_name', 'river XI-position at RHO points');
netcdf.putAtt(ncid, river_Xposition_varid, 'units', 'non-dimensional');
netcdf.putAtt(ncid, river_Xposition_varid, 'valid_min', 1);
netcdf.putAtt(ncid, river_Xposition_varid, 'valid_max', size(grd.lon_u,2));
netcdf.putAtt(ncid, river_Xposition_varid, 'field', 'river_Xposition, scalar');

river_Eposition_varid = netcdf.defVar(ncid, 'river_Eposition', 'float', river_dimid);
netcdf.putAtt(ncid, river_Eposition_varid, 'long_name', 'river ETA-position at RHO points');
netcdf.putAtt(ncid, river_Eposition_varid, 'units', 'non-dimensional');
netcdf.putAtt(ncid, river_Eposition_varid, 'valid_min', 1);
netcdf.putAtt(ncid, river_Eposition_varid, 'valid_max', size(grd.lon_v,2));
netcdf.putAtt(ncid, river_Eposition_varid, 'field', 'river_Eposition, scalar');

river_direction_varid = netcdf.defVar(ncid, 'river_direction', 'float', river_dimid);
netcdf.putAtt(ncid, river_direction_varid, 'long_name', 'river runoff direction');
netcdf.putAtt(ncid, river_direction_varid, 'units', 'non-dimensional');
netcdf.putAtt(ncid, river_direction_varid, 'field', 'river_direction, scalar');

river_flag_varid = netcdf.defVar(ncid, 'river_flag', 'float', river_dimid);
netcdf.putAtt(ncid, river_flag_varid, 'long_name', 'river ETA-position at RHO points');
netcdf.putAtt(ncid, river_flag_varid, 'option_0', 'all tracers are off');
netcdf.putAtt(ncid, river_flag_varid, 'option_1', 'only temperature is on');
netcdf.putAtt(ncid, river_flag_varid, 'option_2', 'only salinity is on');
netcdf.putAtt(ncid, river_flag_varid, 'option_3', 'both temperature and salinity are on');
netcdf.putAtt(ncid, river_flag_varid, 'units', 'non-dimensional');
netcdf.putAtt(ncid, river_flag_varid, 'field', 'river_flag, scalar');

river_Vshape_varid = netcdf.defVar(ncid, 'river_Vshape', 'float', [river_dimid, s_rho_dimid]);
netcdf.putAtt(ncid, river_Vshape_varid, 'long_name', 'river runoff mass tranport vertical profile');
netcdf.putAtt(ncid, river_Vshape_varid, 'units', 'non-dimensional');
netcdf.putAtt(ncid, river_Vshape_varid, 'field', 'river_Vshape, scalar');

river_time_varid = netcdf.defVar(ncid, 'river_time', 'float', time_dimid);
netcdf.putAtt(ncid, river_time_varid, 'long_name', 'river day of year');
netcdf.putAtt(ncid, river_time_varid, 'units', River.time_units);
netcdf.putAtt(ncid, river_time_varid, 'cycle_length',cycle);
netcdf.putAtt(ncid, river_time_varid, 'field', 'river_time, scalar, series');

river_transport_varid = netcdf.defVar(ncid, 'river_transport', 'float', [river_dimid, time_dimid]);
netcdf.putAtt(ncid, river_transport_varid, 'long_name', 'river runoff volume transport');
netcdf.putAtt(ncid, river_transport_varid, 'units', 'meter^3 / sec');
netcdf.putAtt(ncid, river_transport_varid, 'field', 'river_transport, scalar, series');
netcdf.putAtt(ncid, river_transport_varid, 'time', 'river_time');

river_temp_varid = netcdf.defVar(ncid, 'river_temp', 'float', [river_dimid, s_rho_dimid, time_dimid]);
netcdf.putAtt(ncid, river_temp_varid, 'long_name', 'river runoff potential temperature');
netcdf.putAtt(ncid, river_temp_varid, 'units', 'Celsius');
netcdf.putAtt(ncid, river_temp_varid, 'field', 'river_temp, scalar, series');
netcdf.putAtt(ncid, river_temp_varid, 'time', 'river_time');

river_salt_varid = netcdf.defVar(ncid, 'river_salt', 'float', [river_dimid, s_rho_dimid, time_dimid]);
netcdf.putAtt(ncid, river_salt_varid, 'long_name', 'river runoff salinity');
netcdf.putAtt(ncid, river_salt_varid, 'units', '  ');
netcdf.putAtt(ncid, river_salt_varid, 'field', 'river_salt, scalar, series');
netcdf.putAtt(ncid, river_salt_varid, 'time', 'river_time');

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'type', 'ROMS River Forcing file ');
netcdf.putAtt(ncid,varid,'date',datestr(date, 'yyyymmdd'));
netcdf.putAtt(ncid,varid,'grid_file','grd_file');
netcdf.putAtt(ncid,varid,'source',Rname);
netcdf.putAtt(ncid,varid,'details',detailstr);
netcdf.putAtt(ncid,varid,'history', 'Created by Yong-Yub Kim');

%
% Leave define mode
%
netcdf.endDef(ncid);
netcdf.close(ncid);

