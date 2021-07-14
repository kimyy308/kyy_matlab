%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS initial file from CMIP5 interped data
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
%     Data Source : CMIP5 interped data from B-G. Kim
%     Updated 19-May-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%

% romstools_param

model_names={'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};

for modeli = 1:length(model_names)

grdname = 'D:\Data\Model\ROMS\nwp_1_20\input\test53\roms_grid_nwp_1_20_test53.nc';
model_name = model_names{modeli};
ininame = ['D:\Data\Model\ROMS\nwp_1_20\input\roms_nwp_ini_', model_name, '.nc'];
theta_s = 10;
theta_b = 1;
Vtransform =2;
Vstretching = 4;
hc= 250;
N= 40;
tini =1;

%
%  Title 
%
    title='ROMS initial file from CMIP5 interped data';



%  Data climatologies file names:
%
%    temp_daily_data : 5-day averaged SODA3 temperature
%    salt_daily_data : 5-day averaged SODA3 salinity
%

CMIP5_dir = 'D:\Data\Model\CMIP5';
CMIP5_year=1976;
CMIP5_month=1;
% soda3_day=3;

% soda3_daily_data  = [soda3_dir,'soda3.4.2_5dy_ocean_reg_', num2str(soda3_year,'%04i'),'_',num2str(soda3_month,'%02i'), '_', num2str(soda3_day,'%02i'),'.nc']; 



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
               Vtransform, Vstretching, theta_s,theta_b,hc,N,...
               tini,'clobber');
% %
% % Horizontal and vertical interp/extrapolations 
% %
variables={'thetao', 'so', 'uo', 'vo', 'zos'};
roms_variables = {'temp', 'salt', 'u', 'v', 'zeta'};
grid_locations = {'r', 'r', 'u', 'v', 'r'};
for namei=1:length(variables)
    variable = variables{namei};
    roms_variable = roms_variables{namei};
    grid_location = grid_locations{namei};
    
    CMIP5name.(variable) = [CMIP5_dir, '\', variable, '\historical\interp\', variable, '_interp_', model_name, '_historical_r1i1p1_', num2str(CMIP5_year), '.nc'];
    disp(' ')
    disp(' Interpolations / extrapolations')
    disp(' ')
    disp([roms_variable, ' ...'])
    ext_datas_ini_CMIP5_interp(ininame,grdname,CMIP5name.(variable), ...
                    variable,roms_variable,grid_location,tini, Vtransform, Vstretching);  %% for SODA 3.4.2 daily data
end


% disp(' ')
% disp(' Interpolations / extrapolations')
% disp(' ')
% disp(' Temperature...')
% 
%     ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
%                 'temp','temp','r',tini);  %% for SODA 3.4.2 daily data
% 
% disp(' ')
% disp(' Salinity...')
% 
%     ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
%                  'salt','salt','r',tini); %% for SODA 3.4.2 daily data
% 
% disp(' ')
% disp(' Sea Surface Height...')
% 
%     ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
%                 'ssh','zeta','r',tini);  %% for SODA 3.4.2 daily data
%             
% disp(' ')
% disp(' U...')
%     ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
%                 'u','u','u',tini);  %% for SODA 3.4.2 daily data
%             
% disp(' ')
% disp(' V...')
%     ext_datas_ini_SODA3(ininame,grdname,soda3_daily_data, ...
%                 'v','v','v',tini);  %% for SODA 3.4.2 daily data
            
% Initial file
%
if (insitu2pot)  %% it doesn't need to run for CMIP5 interped
  disp(' ')
  disp(' Compute potential temperature from in-situ...')
  getpot(Vtransform,Vstretching,ininame,grdname)
end
%
% Make a few plots
%
disp(' ')
disp(' Make a few plots...')
% test_clim(ininame,grdname,'temp',1,coastfileplot)
% figure
% test_clim(ininame,grdname,'salt',1,coastfileplot)
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end