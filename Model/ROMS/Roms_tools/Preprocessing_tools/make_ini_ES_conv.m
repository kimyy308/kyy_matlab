%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS initial file from Levitus Data
%
%  Extrapole and interpole temperature and salinity from a
%  Climatology to get initial conditions for
%  ROMS (initial netcdf files) .
%  Get the velocities and sea surface elevation via a 
%  geostrophic computation.
%
%  Data input format (netcdf):
%     temperature(T, Z, Y, X)
%     T : time [Months]
%     Z : Depth [m]
%     Y : Latitude [degree north]
%     X : Longitude [degree east]
%
%  Data source : IRI/LDEO Climate Data Library (World Ocean Atlas 1998)
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA98/
%
%  P. Marchesiello & P. Penven - IRD 2005
%
%  Version of 21-Sep-2005
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
romstools_param

%
%  Title 
%

title ='ROMS initial file for East Sea Convection experiment';

%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%
disp(' ')
disp([' Making initial file: ',ininame])
disp(' ')
disp([' Title: ',title])
%
% Initial file
%
% create_inifile(ininame,grdname,title,...
%                theta_s,theta_b,hc,N,...
%                tini,'clobber');

create_inifile_Y(ininame,grdname,title,...
               theta_s,theta_b,hc,N,...
               tini,'clobber');
           
           
       
           
%
% Horizontal and vertical interp/extrapolations 
%
disp(' ')
disp(' Interpolations / extrapolations')
disp(' ')
disp(' Temperature...')

nc = netcdf(ininame, 'write');
nc{'salt'}(:) =  34; 
close(nc);



nc=netcdf(grdname);
h=nc{'h'}(:);  
close(nc);

P=-1e-4*1025*9.81*zlevs(Vtransform,Vstretching,h,0.*h,theta_s,theta_b,hc,N,'r');
nc = netcdf(ininame);
temp_insitu=nc{'temp'}(:);  
S=nc{'salt'}(:);  
close(nc);

nc = netcdf(ininame, 'write');
nc{'temp'}(1,:,:,:)=theta(S,temp_insitu,P);
close(nc);


