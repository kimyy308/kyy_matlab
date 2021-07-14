function make_stddepth_ncfile(var,stddepth, fname, gname);

nc = netcdf(gname);
temp = nc{var}(:);
[t,n,m] = size(temp);
close(nc);

t = length(stddepth);
nw = netcdf.create(fname, 'NETCDF4');
%
%  Create dimensions
%
disp(['xi_rho is ',num2str(m)])
disp(['eta_rho is ',num2str(n)])
eta_rho_dimID = netcdf.defDim(nw,'eta_rho', n); %middle
xi_rho_dimID = netcdf.defDim(nw,'xi_rho',m) ; 
level_dimID = netcdf.defDim(nw,'level',t);
one_dimID = netcdf.defDim(nw,'one',1);
%
%  Create variables and attributes
%
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(nw,varid,'type','ROMS file');
netcdf.putAtt(nw,varid,'title','averaging monthly to yearly');
netcdf.putAtt(nw,varid,'source','roms data');

% switch varname
var_ID = netcdf.defVar(nw, 'temp', 'NC_FLOAT', [xi_rho_dimID, eta_rho_dimID, level_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','time-averaged potential temperature');
netcdf.putAtt(nw, var_ID, 'units','Celsius');
netcdf.defVarChunking(nw,var_ID,'CHUNKED',[m/10, n/10, t/10]);
netcdf.defVarDeflate(nw,var_ID,true,true,1);

var_ID = netcdf.defVar(nw, 'salt', 'NC_FLOAT', [xi_rho_dimID, eta_rho_dimID, level_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','time-averaged salinity');
netcdf.defVarChunking(nw,var_ID,'CHUNKED',[m/10, n/10, t/10]);
netcdf.defVarDeflate(nw,var_ID,true,true,1);

var_ID = netcdf.defVar(nw, 'u', 'NC_FLOAT', [xi_rho_dimID, eta_rho_dimID, level_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','time-averaged u-momentum component');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');
netcdf.defVarChunking(nw,var_ID,'CHUNKED',[m/10, n/10, t/10]);
netcdf.defVarDeflate(nw,var_ID,true,true,1);

var_ID = netcdf.defVar(nw, 'v', 'NC_FLOAT', [xi_rho_dimID, eta_rho_dimID, level_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','time-averaged v-momentum component');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');
netcdf.defVarChunking(nw,var_ID,'CHUNKED',[m/10, n/10, t/10]);
netcdf.defVarDeflate(nw,var_ID,true,true,1);

var_ID = netcdf.defVar(nw, 'w', 'NC_FLOAT', [xi_rho_dimID, eta_rho_dimID, level_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','time-averaged vertical momentum component');
netcdf.putAtt(nw, var_ID, 'units','meter second-1');
netcdf.defVarChunking(nw,var_ID,'CHUNKED',[m/10, n/10, t/10]);
netcdf.defVarDeflate(nw,var_ID,true,true,1);

var_ID = netcdf.defVar(nw, 'zeta', 'NC_FLOAT', [xi_rho_dimID, eta_rho_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','time-averaged zeta');
netcdf.putAtt(nw, var_ID, 'units','meter');
netcdf.defVarChunking(nw,var_ID,'CHUNKED',[m/10, n/10]);
netcdf.defVarDeflate(nw,var_ID,true,true,1);

var_ID = netcdf.defVar(nw, 'depth', 'NC_FLOAT', [one_dimID level_dimID]);
netcdf.putAtt(nw, var_ID, 'long_name','depth');
netcdf.putAtt(nw, var_ID, 'units','meter');

netcdf.endDef(nw);
netcdf.close(nw);

nw = netcdf(fname,'w');
nw{'temp'}(:)= zeros(t,n,m); nw{'salt'}(:)= zeros(t,n,m); nw{'zeta'}(:)= zeros(n,m);
nw{'u'}(:)= zeros(t,n,m); nw{'v'}(:)= zeros(t,n,m); nw{'w'}(:)= zeros(t,n,m);
nw{'depth'}(:)= zeros(t,1);
close(nw)
