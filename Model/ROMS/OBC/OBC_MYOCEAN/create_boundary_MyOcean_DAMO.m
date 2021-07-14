clear; clc; close all
% % 
% % make lateral boundary file from MyOcean Reanalysis monthly data (t,s,u,v,zeta)
% % Updated 15-Feb-2019 by Y.Y.Kim

warning off;
system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop\
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
    OGCM_path = '???';
    bry_prefix='???';
elseif (strcmp(system_name,'GLNXA64'))
    % % for linux
%     dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
end

% % for nwp_1_10, EAST
expname='nwp_1_10';
testname='test06';
domainname='North Western Pacific';

% % Model name (MyOcean Reanalysis)
model_name='MyOcean_rys';

% % specify and make output path directory
OGCM_path = '/data1/stlee/ext_hdd/MyOcean/';
bry_prefix=['/data1/kimyy/Model/ROMS/',expname,'/forcing_matlab/OBC/'];

bry_prefix=[bry_prefix,model_name,'/',testname,'/'];
if (exist(bry_prefix,'dir')~=7)
    mkdir(bry_prefix);
end
bry_prefix=[bry_prefix,model_name,'_'];
bry_suffix = ['_',testname,'.nc'];

OGCM_prefix = 'mercatorglorys12v1_gl12_mean_';


% % check please
Vtransform      = 2;
Vstretching     = 4;

%
% Get the model grid
%
% testname = 'test06';
g = grd('nwp_1_10_linux',testname);
grdname = g.grd_file;
theta_s = g.theta_s;
theta_b = g.theta_b;
hc      = g.hc;
N       = g.N;

ROMS_title  = [expname, ' Model OBC file from ', model_name];

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

Ymin = 2019;
Ymax = 2019;
Mmin = 1;
Mmax = 4;

if (Ymin>=2018)
  OGCM_prefix = 'mercatorpsy4v3r1_gl12_mean_';
end

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
tmp_OGCM_path=OGCM_path;
OGCM_path=[tmp_OGCM_path,num2str(Ymin,'%04i'),'/'];

grid_name = [OGCM_path, ...
    OGCM_prefix, num2str(Ymin,'%04i'), '01.nc'];  % mercatorglorys12v1_gl12_mean_199301.nc
nc = netcdf(grid_name);

LAT = nc{'latitude'}(:); %% 1 dimension, -90 ~ 90 (2041)
LON = nc{'longitude'}(:);  %% 1 dimension -180 ~ 180 (4320)
Z = -nc{'depth'}(:);  %% 1 dimension, 0.494025 ~ 5727.917 (50)
Z = Z(1:end-1);  %% bottom data -> NaN, so last index of Z -> trash

lonT = LON; latT = LAT;
lonU = LON; latU = LAT;
lonV = LON; latV = LAT;
lonZ = LON; latZ = LAT;

NZ = length(Z); 
NZ = NZ - rmdepth;
Z = Z(1:NZ);
close(nc)

for Y = Ymin:Ymax
    OGCM_path=[tmp_OGCM_path,num2str(Y,'%04i'),'/'];
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
        bndy_times = [15:30:365]; % 15 45 ... 345
        nc = netcdf([OGCM_path, ...
            OGCM_prefix, num2str(Y,'%04i'), num2str(M,'%02i'), '.nc']);
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
            interp_OGCM_MyOcean(OGCM_path, OGCM_prefix, Y, M, Roa, interp_method, ...
                lonU, latU, lonV, latV, lonT, latT, lonZ, latZ, Z, 1, ...
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


% ROMS_compile_bndy_CMIP5;