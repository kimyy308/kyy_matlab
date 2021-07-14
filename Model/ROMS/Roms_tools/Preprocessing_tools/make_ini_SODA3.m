%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS initial file from SODA 3.4.2 data
%
%  Extrapole and interpole temperature, salinity, u, v and ssh from a
%  Climatology to get initial conditions for
%  ROMS (initial netcdf files) .
%
%  Data input format (netcdf):
%     temperature(T, Z, Y, X)
%     T : time [Months]
%     Z : st_ocean(depth) [m]
%     Y : yt_ocean [degree north]
%     X : xt_ocean [degree east]
%
%     Data Source : http://dsrs.atmos.umd.edu/DATA/soda3.4.2/REGRIDED/ocean/soda3.4.2_5dy_ocean_reg_1980_01_03.nc
%     Updated 23-Apr-2018 by Yong-Yub Kim
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
    title='ROMS initial file from SODA 3.4.2. 5 day averaged result';



%  Data climatologies file names:
%
%    temp_daily_data : 5-day averaged SODA3 temperature
%    salt_daily_data : 5-day averaged SODA3 salinity
%

soda3_dir=['E:\Data\Reanalysis\SODA_3_4_2\daily\'];
soda3_year=1980;
soda3_month=1;
soda3_day=3;

soda3_daily_data  = [soda3_dir,'soda3.4.2_5dy_ocean_reg_', num2str(soda3_year,'%04i'),'_',num2str(soda3_month,'%02i'), '_', num2str(soda3_day,'%02i'),'.nc']; 

% SODA 3.4.2 data already has potential temperature. 
insitu2pot       = 0;   %1: transform in-situ temperature to potential temperature.

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

    ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
                'temp','temp','r',tini);  %% for SODA 3.4.2 daily data

disp(' ')
disp(' Salinity...')

    ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
                 'salt','salt','r',tini); %% for SODA 3.4.2 daily data

disp(' ')
disp(' Sea Surface Height...')

    ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
                'ssh','zeta','r',tini);  %% for SODA 3.4.2 daily data
            
disp(' ')
disp(' U...')
    ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
                'u','u','u',tini);  %% for SODA 3.4.2 daily data
            
disp(' ')
disp(' V...')
    ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
                'v','v','v',tini);  %% for SODA 3.4.2 daily data
            
% Initial file
%
if (insitu2pot)  %% it doesn't run for SODA3
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
