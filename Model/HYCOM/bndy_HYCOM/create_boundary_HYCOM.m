clear; clc; close all

% addpath(genpath('./'))
addpath(genpath('/home/kimyy/Dropbox/source/matlab/Common/netcdf_old/'))

ROMS_title  = 'NorthPacific Model';
ROMS_config = 'NorthPacific';

OGCM_path = '/data1/kimyy/Model/HYCOM/hycom/GLBb/GLBb0.08_expt19.1_monthly/data/';
OGCM_prefix = 'HYCOM_';
Ymin = 0
Ymax = 0;
Mmin = 1;
Mmax = 1;

bry_prefix='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/forcing_matlab/OBC/HYCOM/HYCOM_';
bry_suffix = '.nc';

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

obc = [1 1 1 1]; % open boundaries (1 = open , [S E N W])

makeplot = 1;
%
% Get the model grid
%
% g = grd('NWP_1_20');
g = grd('ES_1_40');
grdname = g.grd_file;
theta_s = g.theta_s;
theta_b = g.theta_b;
hc      = g.hc;
N       = g.N;
Vtransform =g.Vtransform;
Vstretching =g.Vstretching;

nc = netcdf(grdname);
lon = nc{'lon_rho'}(:);
lat = nc{'lat_rho'}(:);
angle = nc{'angle'}(:);
h = nc{'h'}(:);
close(nc)

lonmin = min(min(lon));
lonmax = max(max(lon));
latmin = min(min(lat));
latmax = max(max(lat));
%
%------------------------------------------------------------------------------------
%
% Get the OGCM grid
grid_name = [OGCM_path, 'HYCOM_', num2str(Ymin,'%04i'), num2str(Mmin,'%02i'), '.nc'];
nc = netcdf(grid_name);
% % original code
LAT = nc{'latitude'}(:); 
LON = nc{'longitude'}(:);
LAT = LAT(:,1);
LON = LON(1,:)';




lonT = LON; latT = LAT;
lonU = LON; latU = LAT;
lonV = LON; latV = LAT;
Z = -nc{'depth'}(:);
NZ = length(Z);
NZ = NZ - rmdepth;
Z = Z(1:NZ);
close(nc)
%
% Loop on the years and the months
%
for Y = Ymin:Ymax
    if Y == Ymin
        mo_min = Mmin;
    else
        mo_min = 1;
    end
    if Y == Ymax
        mo_max = Mmax;
    else
        mo_max = 12;
    end
    for M = mo_min:mo_max
        disp(' ')
        disp(['Processing  year ',num2str(Y),...
            ' - month ',num2str(M)])
        disp(' ')
        %
        Mm = M-1; Ym = Y;
        if Mm == 0
            Mm = 12;
            Ym = Y-1;
        end
        Mp = M+1; Yp = Y;
        if Mp == 13
            Mp = 1;
            Yp = Y+1;
        end
        %
        % Add 2 times step in the ROMS files: 1 at the beginning and 1 at the end
        %
        bndy_times = [15:30:365];
        nc = netcdf([OGCM_path, OGCM_prefix, num2str(Y,'%04i'), num2str(M,'%02i'), '.nc']);
        OGCM_time = bndy_times(M);
        ntimes = length(OGCM_time);
        if ntimes == 1
            dt = 30; % monthly files (SODA..)
        else
            dt = max(gradient(OGCM_time));
        end
        
        roms_time = 0*(1:ntimes+0);
        %
        % Create and open the ROMS files
        %
        bryname = [bry_prefix, 'Y', num2str(Y,'%04i'), ...
            'M', num2str(M,'%02i'), bry_suffix];
        
        create_bryfile_J(bryname, grdname, ROMS_title, [1 1 1 1],...
            theta_s, theta_b, hc, N,...
            roms_time, 0, 'clobber',Vtransform,Vstretching);
        nc_bry = netcdf(bryname, 'write');
        
        nc_clm = [];
        
        %
        % Perform the interpolations for the current month
        %
        disp(' Current month :')
        disp('================')
        for tndx_OGCM = 1:ntimes
            disp([' Time step : ',num2str(tndx_OGCM),' of ',num2str(ntimes),' :'])
            interp_OGCM_HYCOM(OGCM_path, OGCM_prefix, Y, M, Roa, interp_method, ...
                lonU, latU, lonV, latV, lonT, latT, Z, tndx_OGCM, ...
                nc_clm, nc_bry, lon, lat, angle, h, tndx_OGCM + itolap_a, obc,Vtransform,Vstretching)
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
if makeplot == 1
    disp(' ')
    disp(' Make a few plots...')
    bryname=[bry_prefix, 'Y', num2str(Y,'%04i'), 'M', num2str(M,'%02i'), bry_suffix];
    figure
    test_bry(bryname,grdname,'temp',1,obc,Vtransform,Vstretching)
    figure
    test_bry(bryname,grdname,'salt',1,obc,Vtransform,Vstretching)
    figure
    test_bry(bryname,grdname,'u',1,obc,Vtransform,Vstretching)
    figure
    test_bry(bryname,grdname,'v',1,obc,Vtransform,Vstretching)
end
