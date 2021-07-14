clear all;close all;clc

% romstools_param

ROMS_title  = 'NorthPacific Model_';
ROMS_config = 'NorthPacific_';

% bry_prefix='D:\SODA2_2_4\roms_bry_SODA-';
bry_prefix='D:\data\roms_input\np\SODA\roms_bry_SODA-';
rmdepth     = 0;         % Number of bottom levels to remove 
                         %(This is usefull when there is no valid data at this level
                         %i.e if the depth in the domain is shallower than
                         %the OGCM depth)

% %Overlap parameters : before (_a) and after (_p) the months.
itolap_a=0;           %Overlap parameters

% Objective analysis decorrelation scale [m]
% (if Roa=0: simple extrapolation method; crude but much less costly)
%
%Roa=300e3;
Roa=0;

interp_method = 'linear';           % Interpolation method: 'linear' or 'cubic'

obc = [1 1 1 1]; % open boundaries (1=open , [S E N W])

theta_s = 5;
theta_b = 0.4;
hc      =5.;
N       =30;

Ymin=2007;
Ymax=2007;
Mmin=1;
Mmax=12;

makeplot=1;
makebry=1;

nc_suffix='.nc';
%
% Get the model grid
%
grdname=['D:\data\roms_input\np\grid\NP_grid_30layer.nc'];
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
angle=nc{'angle'}(:);
h=nc{'h'}(:);
close(nc)

lonmin=min(min(lon));
lonmax=max(max(lon));
latmin=min(min(lat));
latmax=max(max(lat));
%
%------------------------------------------------------------------------------------
%
% Get the OGCM grid 
% grid_name=['.\reanalysis\SODA_2.2.4_',num2str(Ymin),num2char(Mmin,2),'.cdf'];
grid_name=['.\soda3.3.1_mn_ocean_reg_2007.nc'];
nc=netcdf(grid_name);
LAT=nc{'LAT'}(:);  % from 2009 ~
LON=nc{'LON'}(:);  % from 2009 ~
% LAT=nc{'lon'}(:);  % before 2009 ~
% LON=nc{'lat'}(:);  % before 2009 ~
lonT=LON;
latT=LAT;
lonU=LON;
latU=LAT;
lonV=LON;
latV=LAT;
Z=-nc{'DEPTH'}(:);  % from 2009 ~
% Z=-nc{'depth'}(:);  % before 2009 ~
NZ=length(Z);
NZ=NZ-rmdepth;
Z=Z(1:NZ);

close(nc)

  %
  % Loop on the years and the months
  %
  for Y=Ymin:Ymax
    if Y==Ymin 
      mo_min=Mmin;
    else
      mo_min=1;
    end
    if Y==Ymax
      mo_max=Mmax;
    else
      mo_max=12;
    end
    for M=mo_min:mo_max
      disp(' ')
      disp(['Processing  year ',num2str(Y),...
	    ' - month ',num2str(M)])
      disp(' ')
%       OGCM_dir='D:\SODA2_2_4\reanalysis\';
%       OGCM_prefix = 'SODA_2.2.4_';
      OGCM_dir='D:\data\roms_input\np\SODA\';
      OGCM_prefix = 'SODA3.3.1_mn_ocean_reg_';
      %
      Mm=M-1;Ym=Y;
      if Mm==0
        Mm=12;
        Ym=Y-1;
      end
      Mp=M+1;Yp=Y;
      if Mp==13
        Mp=1;
        Yp=Y+1;
      end
      %
      % Add 2 times step in the ROMS files: 1 at the beginning and 1 at the end 
      %
      bndy_times=[15:30:365];
%       nc=netcdf(['.\reanalysis\SODA_2.2.4_',num2str(Y),num2char(M,2),'.cdf']);
      nc=netcdf(['.\soda3.3.1_mn_ocean_reg_2007.nc']);
      OGCM_time=bndy_times(M);
      ntimes=length(OGCM_time);
      if ntimes==1
        dt=30; % monthly files (SODA..)
      else
        dt=max(gradient(OGCM_time));
      end

      roms_time=0*(1:ntimes+0);
      %
      % Create and open the ROMS files
      %
    bryname=[bry_prefix,'Y',num2str(Y),...
		 'M',num2str(M),nc_suffix];
     
    create_bryfile(bryname,grdname,ROMS_title,[1 1 1 1],...
		       theta_s,theta_b,hc,N,...
		       roms_time,0,'clobber');
	nc_bry=netcdf(bryname,'write');
    warning off

	nc_clm=[];

      %
      % Perform the interpolations for the current month
      %
      disp(' Current month :')
      disp('================')
      for tndx_OGCM=1:ntimes
	disp([' Time step : ',num2str(tndx_OGCM),' of ',num2str(ntimes),' :'])
	interp_OGCM_224(OGCM_dir,OGCM_prefix,Y,M,Roa,interp_method,...
		    lonU,latU,lonV,latV,lonT,latT,Z,tndx_OGCM,...
		    nc_clm,nc_bry,lon,lat,angle,h,tndx_OGCM+itolap_a,obc)
      end
      %
      % Close the ROMS files
      %
      if ~isempty(nc_clm)
	close(nc_clm);
      end
      if ~isempty(nc_bry)
	close(nc_bry);
      end
      %
    end
  end

%
% Spin-up: (reproduce the first year 'SPIN_Long' times)
% just copy the files for the first year and change the time
%

%---------------------------------------------------------------
% Make a few plots
%---------------------------------------------------------------
if makeplot==1
  disp(' ')
  disp(' Make a few plots...')
  if makebry==1
    bryname=[bry_prefix,'Y',num2str(Y),'M',num2str(M),nc_suffix];
    figure
    test_bry(bryname,grdname,'temp',1,obc)
    figure
    test_bry(bryname,grdname,'salt',1,obc)
    figure
    test_bry(bryname,grdname,'u',1,obc)
    figure
    test_bry(bryname,grdname,'v',1,obc)
  end
end