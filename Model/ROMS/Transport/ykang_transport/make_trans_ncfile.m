function make_trans_ncfile(var, sname, fname)

[t,n,m] = size(var);

nw = netcdf.create(sname, 'clobber');
%
%  Create dimensions
%
time_dimID = netcdf.defDim(nw,'time',t) ; 
xi_rho_dimID = netcdf.defDim(nw,'xi_rho',m) ; 
level_dimID = netcdf.defDim(nw,'level',n);
wlevel_dimID = netcdf.defDim(nw,'wlevel',n+1);
one_dimID = netcdf.defDim(nw,'one',1);
%
%  Create variables and attributes
%
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(nw,varid,'type','ROMS file');
netcdf.putAtt(nw,varid,'title','calculate variables at calculated transport section');
netcdf.putAtt(nw,varid,'source','roms data');

% switch varname
var_ID = netcdf.defVar(nw, 'lon', 'double', [xi_rho_dimID, one_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','longitude at transport section');
netcdf.putAtt(nw, var_ID, 'units','degree of east');

var_ID = netcdf.defVar(nw, 'mask', 'double', [xi_rho_dimID, one_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','mask at transport section');

var_ID = netcdf.defVar(nw, 'h', 'double', [xi_rho_dimID, one_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','depth at transport section');
netcdf.putAtt(nw, var_ID, 'units','Celsius');

var_ID = netcdf.defVar(nw, 'temp', 'double', [xi_rho_dimID, level_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','potential temperature at transport section');
netcdf.putAtt(nw, var_ID, 'units','Celsius');

var_ID = netcdf.defVar(nw, 'salt', 'double', [xi_rho_dimID, level_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','salinity at transport section');

var_ID = netcdf.defVar(nw, 'rho', 'double', [xi_rho_dimID, level_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','density at transport section');
netcdf.putAtt(nw, var_ID, 'units','kilogram meter-3');

var_ID = netcdf.defVar(nw, 'zeta', 'double', [xi_rho_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','sea surface height at transport section');
netcdf.putAtt(nw, var_ID, 'units','meter');

var_ID = netcdf.defVar(nw, 'u', 'double', [xi_rho_dimID, level_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','u-momentum component at transport section');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');

var_ID = netcdf.defVar(nw, 'v', 'double', [xi_rho_dimID, level_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','v-momentum component at transport section');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');

var_ID = netcdf.defVar(nw, 'w', 'double', [xi_rho_dimID, wlevel_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','vertical momentum component at transport section');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');

var_ID = netcdf.defVar(nw, 'Hz', 'double', [xi_rho_dimID, level_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','vertical depth at transport section');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');

var_ID = netcdf.defVar(nw, 'transport', 'double', [xi_rho_dimID, level_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','transport at section');
netcdf.putAtt(nw, var_ID, 'units','Sv');

var_ID = netcdf.defVar(nw, 'Uwind', 'double', [xi_rho_dimID one_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','u-momentum wind component at transport section');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');

var_ID = netcdf.defVar(nw, 'Vwind', 'double', [xi_rho_dimID one_dimID, time_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','v-momentum wind component at transport section');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');

netcdf.endDef(nw);
netcdf.close(nw);

nw = netcdf(sname,'w');
nw{'lon'}(:)=zeros(1,m);
nw{'mask'}(:)=zeros(1,m);
nw{'h'}(:)=zeros(1,m);
nw{'temp'}(:)= zeros(t,n,m); 
nw{'salt'}(:)= zeros(t,n,m); 
nw{'rho'}(:)= zeros(t,n,m);
nw{'zeta'}(:)=zeros(t,m);
nw{'u'}(:)= zeros(t,n,m); 
nw{'v'}(:)= zeros(t,n,m); 
nw{'w'}(:)= zeros(t,n+1,m);
nw{'Hz'}(:)=zeros(t,n,m);
nw{'transport'}(:)=zeros(t,n,m);
nw{'Uwind'}(:)=zeros(t,m);
nw{'Vwind'}(:)=zeros(t,m);
close(nw)
