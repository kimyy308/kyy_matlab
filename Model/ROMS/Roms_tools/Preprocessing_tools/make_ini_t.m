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
%  Title 
%
title='Climatology';
%
% Common parameters
%
 romstools_param

%  Data climatologies file names:
%
%    temp_month_data : monthly temperature climatology
%    temp_ann_data   : annual temperature climatology
%    salt_month_data : monthly salinity climatology
%    salt_ann_data   : annual salinity climatology
%
% temp_month_data  = [woa_dir,'temp_month.cdf'];
% temp_ann_data    = [woa_dir,'temp_ann.cdf'];
temp_month_data  = [woa_dir,'temperature_monthly_1deg.nc'];
temp_ann_data    = [woa_dir,'temperature_annual_1deg.nc'];
insitu2pot       = 1;   %1: transform in-situ temperature to potential temperature
% salt_month_data  = [woa_dir,'salt_month.cdf'];
% salt_ann_data    = [woa_dir,'salt_ann.cdf'];
salt_month_data  = [woa_dir,'salinity_monthly_1deg.nc'];
salt_ann_data    = [woa_dir,'salinity_annual_1deg.nc'];
% %
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%

% ininame=['roms_add4_05_2ep_ini.nc'];
% grdname=['D:\add4_ini_bry_grd\grid\roms_grid4_ADD_05_2_ep.nc'];
% scoord = [5.0 0.4 5 40]; 
% theta_s=scoord(1);
% theta_b=scoord(2);
% hc=scoord(3);
% N=scoord(4);
tini=0;
disp(' ')
disp([' Making initial file: ',ininame])
disp(' ')
disp([' Title: ',title])
%
% Initial file
%

create_inifile(ininame,grdname,title,...
               theta_s,theta_b,hc,N,...
               tini,'clobber');
%
% Horizontal and vertical interp/extrapolations 
%
disp(' ')
disp(' Interpolations / extrapolations')
disp(' ')
disp(' Temperature...')
% ext_tracers_ini(ininame,grdname,temp_month_data,temp_ann_data,...
%             'temperature','temp','r',tini);
ext_tracers_ini_2009(ininame,grdname,temp_month_data,temp_ann_data,...
            't_an','temp','r',tini);
disp(' ')
disp(' Salinity...')
% ext_tracers_ini(ininame,grdname,salt_month_data,salt_ann_data,...
%              'salinity','salt','r',tini);
 ext_tracers_ini_2009(ininame,grdname,salt_month_data,salt_ann_data,...
    's_an','salt','r',tini);
%
% Geostrophy
%
 disp(' ')
 disp(' Compute geostrophic currents')
 frcname=[];
 oaname=['D:\Roms_tools\WOA1998\temp_ann.cdf'];
 geost_currents(ininame,grdname,oaname,frcname,zref,obc,0)

% Initial file
%
if (insitu2pot)
  disp(' ')
  disp(' Compute potential temperature from in-situ...')
  getpot(ininame,grdname)
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
