clear; clc; close all
% % Updated 28-Mar-2018 by Y.Y.Kim
% % Updated 29-Mar-2018 by Y.Y.Kim
warning off;
linux=1; windows=0;
if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop\
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
    OGCM_path = 'E:\Data\Reanalysis\SODA_3_4_2\';
    bry_prefix='E:\Data\Model\ROMS\nwp_1_20\input\test38\SODA_';
elseif (linux==1)
    % % for linux
%     dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
    OGCM_path = '/data2/kimyy/Reanalysis/SODA/SODA_3_4_2/monthly/';
    bry_prefix='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_40/forcing_matlab/OBC/SODA_3_4_2/SODA_';
end

% % check please
Vtransform      = 1;
Vstretching     = 1;



%
% Get the model grid
%
testname = 'test01';
g = grd('NWP_1_40_linux',testname);
grdname = g.grd_file;
theta_s = g.theta_s;
theta_b = g.theta_b;
hc      = g.hc;
N       = g.N;




ROMS_title  = 'NorthPacific Model';
ROMS_config = 'NorthPacific';

OGCM_prefix = 'soda3.4.2_mn_ocean_reg_';
bry_suffix = ['_',testname,'.nc'];

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

Ymin = 2010;
Ymax = 2010;
Mmin = 1;
Mmax = 12;

makeplot = 1;


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
grid_name = [OGCM_path, OGCM_prefix, num2str(Ymin), '.nc'];
nc = netcdf(grid_name);

% % SODA 3.4.1
% LAT = nc{'latitude'}(:);
% LON = nc{'longitude'}(:);
% Z = -nc{'depth'}(:);

% % SODA 3.4.2
LAT = nc{'yt_ocean'}(:);
LON = nc{'xt_ocean'}(:);
Z = -nc{'st_ocean'}(:);

lonT = LON; latT = LAT;
lonU = LON; latU = LAT;
lonV = LON; latV = LAT;

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
        nc = netcdf([OGCM_path, OGCM_prefix, num2str(Y), '.nc']);
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
        bryname = [bry_prefix, 'Y', num2str(Y), ...
            'M', num2str(M), bry_suffix];
        
        create_bryfile_J(bryname, grdname, ROMS_title, [1 1 1 1],...
            theta_s, theta_b, hc, N,...
            roms_time, 0, 'clobber');
        nc_bry = netcdf(bryname, 'write');
        
        nc_clm = [];
        
        %
        % Perform the interpolations for the current month
        %
        disp(' Current month :')
        disp('================')
        for tndx_OGCM = 1:ntimes
            disp([' Time step : ',num2str(tndx_OGCM),' of ',num2str(ntimes),' :'])
            interp_OGCM_SODA3(OGCM_path, OGCM_prefix, Y, M, Roa, interp_method, ...
                lonU, latU, lonV, latV, lonT, latT, Z, M, ...
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
    bryname=[bry_prefix, 'Y', num2str(Y), 'M', num2str(M), bry_suffix];
    figure
    test_bry(bryname,grdname,'temp',1,obc,Vtransform,Vstretching)
    figure
    test_bry(bryname,grdname,'salt',1,obc,Vtransform,Vstretching)
    figure
    test_bry(bryname,grdname,'u',1,obc,Vtransform,Vstretching)
    figure
    test_bry(bryname,grdname,'v',1,obc,Vtransform,Vstretching)
end

% % % clear all;close all;clc;
% % brytype='SODA';
% % for yy=Ymin:Ymax
% %     ydays=yeardays(yy);
% % 
% %     if ydays == 366
% %         time_month=[31 29 31 30 31 30 31 31 30 31 30 31];
% %     else
% %         time_month=[31 28 31 30 31 30 31 31 30 31 30 31];
% %     end
% % 
% % 
% %     for i_ot=1:1:12
% %         if i_ot ==1
% %             time(i_ot,1)=(time_month(i_ot)/2);
% %         else
% %             time(i_ot,1)=(sum(time_month(1:i_ot-1))+time_month(i_ot)/2);
% %         end
% %     end
% %     time=[15:30:345];
% %     bryname=['roms_NWP_bry2_',brytype,'-Y',num2str(yy),'_test38.nc'];
% % %     grdname='D:\add2_ini_bry_grd\grid\roms_grid2_ADD_08_2_ep.nc';
% % %     grdname=['D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_edite
% % %     d_nwp_kyy.nc'];
% %     create_bryfile_t(bryname,brytype,grdname,yy,ydays,time,Vtransform);
% % end



ROMS_compile_bndy_SODA_nwp_1_40;