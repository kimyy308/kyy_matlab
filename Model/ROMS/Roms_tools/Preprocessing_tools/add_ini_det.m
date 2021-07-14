function add_ini_det(inifile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function add_ini_det(inifile);
%

%  input:
%    detritus = 0

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp('Add_ini_det: creating variable and attribute')
%
% open the initial file  
% 
nc=netcdf(inifile,'write');
time= nc{'scrum_time'}(:);
tlen=length(time);
redef(nc);
nc{'detritus'} = ncdouble('time','s_rho','eta_rho','xi_rho') ;
nc{'detritus'}.long_name = ncchar('detritus');
nc{'detritus'}.long_name = 'detritus';
nc{'detritus'}.units = ncchar('mMol N m-3');
nc{'detritus'}.units = 'mMol N m-3';
nc{'detritus'}.fields = ncchar('detritus, scalar, series');
nc{'detritus'}.fields = 'detritus, scalar, series';
endef(nc);
%
% loop on time
%
nc{'detritus'}(:)=0;

close(nc);
return
