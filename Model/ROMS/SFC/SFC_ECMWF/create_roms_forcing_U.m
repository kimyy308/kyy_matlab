function create_roms_forcing_U(fname,Uwind,time,yy,cycle)%,tframe,cycle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Create an empty netcdf frc file 
%       x: total number of rho points in x direction
%       y: total number of rho points in y direction
%       varname: name of field variable
%       fname: name of the ecmwf file
%       var: mean file
%
%                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nw = netcdf(fname, 'clobber');
result = redef(nw);

disp(['file is ',fname])
[t,n,m]=size(Uwind);
var_time=['Uwind_time'];
%
%  Create dimensions
%
disp(['xi_rho is ',num2str(m)])
disp(['eta_rho is ',num2str(n)])
nw('eta_rho') = n; 
nw('xi_rho') = m; 
nw('eta_u') = n;  
nw('xi_u') = m-1;  
nw('eta_v') = n-1;
nw('xi_v') = m;   
nw('eta_psi') = n-1;
nw('xi_psi') = m-1; 
nw(var_time) = t;
%
%  Create variables and attributes
%

nw.type = ' ROMS forcing file ';
nw.title = ' Bulk Formular Forcing file ';
nw.source = 'ECMWF-interim ';
nw.author = 'Created by T';
nw.date = date;

nw{var_time} = ncdouble(var_time);
nw{var_time}.long_name = ncchar('Year Month (1: January, ... 12: December)');
nw{var_time}.long_name = 'Year Month (1: January, ... 12: December)';
nw{var_time}.units = ncchar('DAYS');
nw{var_time}.units = 'DAYS';
nw{var_time}.cycle_length = ncdouble(cycle);
nw{var_time}.cycle_length = cycle;
nw{var_time}(:)=time';

 nw{'Uwind'} = ncfloat(var_time, 'eta_rho', 'xi_rho');
 nw{'Uwind'}.long_name = ncchar('ECMWF Ten meter U');
 nw{'Uwind'}.units = ncchar('m/s');
 nw{'Uwind'}.time = ncchar(var_time);
for i=1:1:yy
 nw{'Uwind'}(i,:,:)= Uwind(i,:,:);
end
 close(nw);
