clear; clc; close all
% % 
% % make lateral boundary file from CMIP5 model monthly data(yearly packed) (t,s,u,v,zeta)
% % Updated 04-Mar-2019 by Y.Y.Kim (mepl2)

warning off;
system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop\
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run']));
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
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
    
end

% testname='test53';

%testnames ={'test61', 'test62', 'test63', 'test64'};
 testnames ={'test16'};

% testnames ={'test66'};

% % CMIP5 Model name
% model_name='ens10';
model_names = {'NorESM1-M'};
% model_names = {'NorESM1-M', 'MPI-ESM-LR'};

% model_names = {'IPSL-CM5A-MR'};
% scenario_name='historical'; % historical, rcp45, rcp26, rcp85
scenario_name='rcp85'; % historical, rcp45, rcp26, rcp85

switch scenario_name
    case 'historical'
        tmp_OGCM_path = '/data1/CMIP/cmip5/historical_extHDD/CMIP5/';
    case 'rcp26'
        tmp_OGCM_path = '/data1/CMIP/cmip5/rcp_ocn_extHDD/CMIP5/';
    case 'rcp45'
        tmp_OGCM_path = '/data1/CMIP/cmip5/rcp_ocn_extHDD/CMIP5/';
    case 'rcp85'
        tmp_OGCM_path = '/data1/CMIP/cmip5/rcp_ocn_extHDD/CMIP5/';
end

tmp_bry_prefix='/data1/kimyy/Model/ROMS/nwp_1_10/forcing_matlab/OBC/';

for model_names_index=1:length(model_names)

    testname = testnames{model_names_index};
    model_name = model_names{model_names_index};
   
    bry_prefix=[tmp_bry_prefix,scenario_name,'/',model_name,'/'];
    if (exist(bry_prefix,'dir')~=7)
            mkdir(bry_prefix);
    end
    bry_prefix=[bry_prefix,model_name,'_'];
    OGCM_prefix = model_name;
    bry_suffix = ['_',testname,'.nc'];


    % % check please
    Vtransform      = 2;
    Vstretching     = 4;

    %
    % Get the model grid
    %
%     testname = 'test51';
    g = grd('NWP_1_10_linux',testname);
    grdname = g.grd_file;
    theta_s = g.theta_s;
    theta_b = g.theta_b;
    hc      = g.hc;
    N       = g.N;

    ROMS_title  = ['NWP 1/10 Model OBC file from ', model_name];
    ROMS_config = 'NorthWesternPacific';

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
    Roa=0;  % It should be 0. 0 is used in other functions

    interp_method = 'linear';           % Interpolation method: 'linear' or 'cubic'

    obc = [1 1 1 1]; % open boundaries (1 = open , [S E N W])
    

    Ymin = 2006;
    Ymax = 2018;
    Mmin = 1;
    Mmax = 12;

    makeplot = 0;


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
%     tmp_OGCM_path=OGCM_path;
%     OGCM_path=[tmp_OGCM_path,'thetao/historical/interp/ensemble/'];
    OGCM_path=[tmp_OGCM_path,'thetao/',scenario_name,'/interp/',model_name,'/'];

%     grid_name = [OGCM_path, 'thetao_interp_ensemble_historical_r1i1p1_', num2str(Ymin,'%04i'), '.nc'];
    grid_name = [OGCM_path, 'thetao_interp_',model_name,'_',scenario_name,'_r1i1p1_', num2str(Ymin,'%04i'), '.nc'];

    nc = netcdf(grid_name);

    LAT = nc{'lat'}(:);
    LON = nc{'lon'}(:);
    Z = -nc{'depth'}(:);

    lonT = LON; latT = LAT;
    lonU = LON; latU = LAT;
    lonV = LON; latV = LAT;
    lonZ = LON; latZ = LAT;

    NZ = length(Z);
    NZ = NZ - rmdepth;
    Z = Z(1:NZ);
    close(nc)

    % grid_name = [OGCM_path, OGCM_prefix, '_zos_', num2str(Ymin), '_v2.nc'];
    % nc = netcdf(grid_name);
    % latZ = nc{'lat'}(:);
    % lonZ = nc{'lon'}(:);
    % close(nc)

    for Y = Ymin:Ymax
    %     OGCM_path=[tmp_OGCM_path];
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
            nc = netcdf([OGCM_path, 'thetao_interp_',model_name,'_',scenario_name,'_r1i1p1_',num2str(Y,'%04i'), '.nc']);
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
                switch scenario_name
                    case 'historical'
                         interp_OGCM_CMIP5_climssh(tmp_OGCM_path, OGCM_prefix, Y, M, interp_method, ...
                    lonU, latU, lonV, latV, lonT, latT, lonZ, latZ, Z, M, ...
                    nc_bry, lon, lat, angle, h, tndx_OGCM + itolap_a, obc,Vtransform,Vstretching,model_name,scenario_name)

                    case 'rcp26'
                         interp_OGCM_CMIP5_climssh(tmp_OGCM_path, OGCM_prefix, Y, M, interp_method, ...
                    lonU, latU, lonV, latV, lonT, latT, lonZ, latZ, Z, M, ...
                    nc_bry, lon, lat, angle, h, tndx_OGCM + itolap_a, obc,Vtransform,Vstretching,model_name,scenario_name)

                    case 'rcp85'
                        interp_OGCM_CMIP5_climssh(tmp_OGCM_path, OGCM_prefix, Y, M, interp_method, ...
                    lonU, latU, lonV, latV, lonT, latT, lonZ, latZ, Z, M, ...
                    nc_bry, lon, lat, angle, h, tndx_OGCM + itolap_a, obc,Vtransform,Vstretching,model_name,scenario_name)

                end
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

end
% ROMS_compile_bndy_CMIP5;
