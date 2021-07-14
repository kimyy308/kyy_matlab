function create_frc_merge(fname,grdname,varname,cycle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Create an empty netcdf frc file 
%       x: total number of rho points in x direction
%       y: total number of rho points in y direction
%       varname: name of field variable
%       fname: name of the ecmwf file
%       var: variable of ecmwf file
%
%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nc=netcdf(grdname);
L=length(nc('xi_psi'));
M=length(nc('eta_psi'));
lon = nc{'lon_rho'}(:);
lat = nc{'lat_rho'}(:);
result=close(nc);
Lp=L+1;
Mp=M+1;

nw = netcdf(fname, 'clobber');
% result = redef(nw);

disp(['file is ',fname])
ocean_time='ocean_time';

%
%  Create dimensions
%
nw('eta_rho') = Mp;
nw('xi_rho') = Lp;
nw(ocean_time) = 'unlimited';

%
%  Create variables and attributes
%
nw.type = ' ROMS forcing file ';
nw.title = ' Bulk Formular Forcing file ';
nw.source = ' ROMS 4 forcing file';
nw.author = 'Created by C-S, Kim ';
nw.date = date;

nw{ocean_time} = ncdouble(ocean_time);
nw{ocean_time}.long_name = ncchar('Julian days from ref_year');
nw{ocean_time}.units = ncchar('Julian days');
nw{ocean_time}.cycle_length = cycle;

nw{'lon_rho'} = ncdouble('eta_rho','xi_rho');
nw{'lon_rho'}.long_name = ncchar('x location of RHO-points');
nw{'lon_rho'}.units = ncchar('degree');
nw{'lon_rho'}(:) = lon;

nw{'lat_rho'} = ncdouble('eta_rho','xi_rho');
nw{'lat_rho'}.long_name = ncchar('y location of RHO-points');
nw{'lat_rho'}.units = ncchar('degree');
nw{'lat_rho'}(:) = lat;

switch varname
case 'Uwind'
 nw{'Uwind'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'Uwind'}.long_name = ncchar('10 meter wind of x-direction');
 nw{'Uwind'}.units = ncchar('m/s');
 nw{'Uwind'}.time = ncchar('ocean_time') ;
close(nw);
 
case 'Vwind'
 nw{'Vwind'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'Vwind'}.long_name = ncchar('10 meter wind of y-direction');
 nw{'Vwind'}.units = ncchar('m/s');
 nw{'Vwind'}.time = ncchar('ocean_time') ;
close(nw);
 
case 'Tair'
 nw{'Tair'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'Tair'}.long_name = ncchar('2 meter Temperature');
 nw{'Tair'}.units = ncchar('Celsius');
 nw{'Tair'}.time = ncchar('ocean_time') ;
% nw{'Tair'}(:)= var;
close(nw);

case 'Pair'
 nw{'Pair'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'Pair'}.long_name = ncchar('sea level air pressure');
 nw{'Pair'}.units = ncchar('millibar');
 nw{'Pair'}.time = ncchar('ocean_time') ;
% nw{'Pair'}(:)= var ;
close(nw); 

case 'Qair'
 nw{'Qair'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'Qair'}.long_name = ncchar('relative humidity');
 nw{'Qair'}.units = ncchar(' % ');
 nw{'Qair'}.time = ncchar('ocean_time') ;
% nw{'Qair'}(:) = var ;
close(nw);

case 'swrad'
 nw{'swrad'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'swrad'}.long_name = ncchar('solar shortwave radiation flux');
 nw{'swrad'}.units = ncchar('watt meter-2');
 nw{'swrad'}.time = ncchar('ocean_time') ;
% nw{'swrad'}(:) = var ;
 close(nw);
 
case 'lwrad'
 nw{'lswrad'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'lwrad'}.long_name = ncchar('net longwave radiation flux');
 nw{'lwrad'}.units = ncchar('watt meter-2');
 nw{'lwrad'}.time = ncchar('ocean_time') ;
% nw{'swrad'}(:) = var ;
 close(nw);
 
 case 'rain'
 nw{'rain'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'rain'}.long_name = ncchar('Precipitation Rate');
 nw{'rain'}.units = ncchar('kilogram meter-2 second-1');
 nw{'rain'}.time = ncchar('ocean_time') ;
%  nw{'rain'}(:) = var ;
 close(nw);
 
 case 'cloud'
 nw{'cloud'} = ncfloat('ocean_time', 'eta_rho', 'xi_rho');
 nw{'cloud'}.long_name = ncchar('Total Cloud Cover');
 nw{'cloud'}.units = ncchar('nondimensional (0-1)');
 nw{'cloud'}.time = ncchar('ocean_time') ;
%  nw{'cloud'}(:) = var ;
 close(nw);
end

