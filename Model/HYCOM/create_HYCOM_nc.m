function create_HYCOM_nc(fname, size_3d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Create empty HYCOM NetCDF file for monthly data
%
%       create_HYCOM_nc(fname, [Depth, Lat, Lon])
%       Fname: Output file name
%       Example) create_HYCOM_nc(myncfile.nc, [20, 488, 386])
%
%       J. Jung
%       edited by Y. Y. Kim (Jan 16, 2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = size_3d(1);
n = size_3d(2);
m = size_3d(3);

% Generate NetCDF file
ncid = netcdf.create(fname, 'clobber');

% Global Attributes
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'type','HYCOM monthly mean');
netcdf.putAtt(ncid,varid,'author','Created by Yong-Yub Kim');
netcdf.putAtt(ncid,varid,'date',datestr(date, 'yyyymmdd'));

% Dimensions
time_dimID = netcdf.defDim(ncid, 'time', netcdf.getConstant('NC_UNLIMITED'));
depth_dimID = netcdf.defDim(ncid,'depth', s);
lat_dimID = netcdf.defDim(ncid,'lat', n);
lon_dimID = netcdf.defDim(ncid,'lon', m);

% Attributes associated with the variable
depth_ID = netcdf.defVar(ncid, 'depth', 'double', depth_dimID);
netcdf.putAtt(ncid, depth_ID, 'standard_name', 'depth');
netcdf.putAtt(ncid, depth_ID, 'units', 'm');
netcdf.putAtt(ncid, depth_ID, 'positive', 'down');
netcdf.putAtt(ncid, depth_ID, 'axis', 'Z');

time_ID = netcdf.defVar(ncid, 'time', 'double', time_dimID);
netcdf.putAtt(ncid, time_ID, 'long_name', 'time');
netcdf.putAtt(ncid, time_ID, 'units', 'days since 1900-12-31 00:00:00');
netcdf.putAtt(ncid, time_ID, 'calendar', 'gregorian');
netcdf.putAtt(ncid, time_ID, 'field', 'time, scalar, series');

lon_ID = netcdf.defVar(ncid, 'lon', 'double', lon_dimID);
netcdf.putAtt(ncid, lon_ID, 'standard_name', 'lon');
netcdf.putAtt(ncid, lon_ID, 'units', 'degrees_east');

lat_ID = netcdf.defVar(ncid, 'lat', 'double', lat_dimID);
netcdf.putAtt(ncid, lat_ID, 'standard_name', 'lat');
netcdf.putAtt(ncid, lat_ID, 'units', 'degrees_north');

lon_ID = netcdf.defVar(ncid, 'longitude', 'double', [lon_dimID lat_dimID]);
netcdf.putAtt(ncid, lon_ID, 'standard_name', 'longitude');
netcdf.putAtt(ncid, lon_ID, 'units', 'degrees_east');

lat_ID = netcdf.defVar(ncid, 'latitude', 'double', [lon_dimID lat_dimID]);
netcdf.putAtt(ncid, lat_ID, 'standard_name', 'latitude');
netcdf.putAtt(ncid, lat_ID, 'units', 'degrees_north');

% Sea surface height
var_ID = netcdf.defVar(ncid, 'ssh', 'double', [lon_dimID lat_dimID time_dimID]);
netcdf.putAtt(ncid, var_ID, 'long_name', 'Sea Surface Height');
netcdf.putAtt(ncid, var_ID, 'units', 'm');
netcdf.putAtt(ncid, var_ID, 'missing_value', double(1.267650600228229e+30));

% Temperature
var_ID = netcdf.defVar(ncid, 'temp', 'double', [lon_dimID lat_dimID depth_dimID time_dimID]);
netcdf.putAtt(ncid, var_ID, 'long_name', 'Temperature');
netcdf.putAtt(ncid, var_ID, 'units', 'Celsius');
netcdf.putAtt(ncid, var_ID, 'missing_value', double(1.267650600228229e+30));

% Salinity
var_ID = netcdf.defVar(ncid, 'salt', 'double', [lon_dimID lat_dimID depth_dimID time_dimID]);
netcdf.putAtt(ncid, var_ID, 'long_name', 'Salinity');
netcdf.putAtt(ncid, var_ID, 'units', 'psu');
netcdf.putAtt(ncid, var_ID, 'missing_value', double(1.267650600228229e+30));

% U
var_ID = netcdf.defVar(ncid, 'u', 'double', [lon_dimID lat_dimID depth_dimID time_dimID]);
netcdf.putAtt(ncid, var_ID, 'long_name', 'Zonal Velocity');
netcdf.putAtt(ncid, var_ID, 'units', 'm/s');
netcdf.putAtt(ncid, var_ID, 'missing_value', double(1.267650600228229e+30));

% V
var_ID = netcdf.defVar(ncid, 'v', 'double', [lon_dimID lat_dimID depth_dimID time_dimID]);
netcdf.putAtt(ncid, var_ID, 'long_name', 'Meridional Velocity');
netcdf.putAtt(ncid, var_ID, 'units', 'm/s');
netcdf.putAtt(ncid, var_ID, 'missing_value', double(1.267650600228229e+30));

netcdf.endDef(ncid);

netcdf.close(ncid);