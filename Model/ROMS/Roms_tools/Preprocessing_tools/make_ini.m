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
if (WOA_switch==0)
    title='ROMS initial file from WOA 1998';
else
    title='ROMS initial file from WOA 2013';
end


% for making initial file  --> in romstools_param
% % WOA_switch = 0; %% 0 -> WOA1998, 1 -> WOA2013


%  Data climatologies file names:
%
%    temp_month_data : monthly temperature climatology
%    temp_ann_data   : annual temperature climatology
%    salt_month_data : monthly salinity climatology
%    salt_ann_data   : annual salinity climatology
%

if (WOA_switch==0)
    temp_month_data  = [woa_dir,'temp_month.cdf'];  %% for WOA 1998
    salt_month_data  = [woa_dir,'salt_month.cdf']; %% for WOA 1998
else
    temp_month_data  = [woa_dir,'woa2013_temp.nc'];  %% for WOA 2013
    salt_month_data  = [woa_dir,'woa2013_salt.nc']; %% for WOA 2013
end

temp_ann_data    = [woa_dir,'temp_ann.cdf'];
salt_ann_data    = [woa_dir,'salt_ann.cdf'];
insitu2pot       = 1;   %1: transform in-situ temperature to potential temperature

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

if (WOA_switch==0)
    ext_tracers_ini(ininame,grdname,temp_month_data,temp_ann_data,...
                'temperature','temp','r',tini);  %% for WOA 1998
else
    ext_tracers_ini(ininame,grdname,temp_month_data,temp_ann_data,...
                'T_AN','temp','r',tini);  %% for WOA 2013
end
disp(' ')
disp(' Salinity...')
if (WOA_switch==0)
    ext_tracers_ini(ininame,grdname,salt_month_data,salt_ann_data,...
                 'salinity','salt','r',tini); %% for WOA 1998
else
     ext_tracers_ini(ininame,grdname,salt_month_data,salt_ann_data,...
         'S_AN','salt','r',tini); %% for WOA 1998
end


%
% Geostrophy
%
 disp(' ')
 disp(' Compute geostrophic currents')

geost_currents(ininame,grdname,oaname,frcname,zref,obc,15) %% 15 means January

% % disp('geostrophic current is not generated (commented)');
% % disp('geostrophic current is not generated (commented)');
% % disp('geostrophic current is not generated (commented)');
% % disp('geostrophic current is not generated (commented)');

%
% Initial file
%
if (insitu2pot)
  disp(' ')
  disp(' Compute potential temperature from in-situ...')
  getpot(Vtransform,Vstretching,ininame,grdname)
end
%
% Make a few plots
%
disp(' ')
disp(' Make a few plots...')
test_clim(ininame,grdname,'temp',1,coastfileplot)
figure
test_clim(ininame,grdname,'salt',1,coastfileplot)
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
