% clear; clc; close all
% % 
% % make lateral boundary file from MyOcean Reanalysis monthly data (t,s,u,v,zeta)
% Updated 14-Oct-2020 by Y.Y.Kim
clear all; close all; clc;

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
    addpath(genpath([dropboxpath '/util/m_map']));
    addpath(genpath([dropboxpath '/util/Figure']));
    addpath(genpath([dropboxpath '/util/netcdf_old']));
    addpath(genpath([dropboxpath '/util/roms_grid']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
    addpath(genpath([dropboxpath '/Roms_tools/Preprocessing_tools']));
    addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
end

% % for nwp_1_10, EAST
expname='nwp_1_10';
testname='test19';
domainname='North Western Pacific';
% % Model name (MyOcean Reanalysis)
% model_name='auto01_OBC';


% % specify and make input, output path directory
% runtime=datevec(datenum(date)-2); % year, month, day(2 days ago)

inputyear=1994:2018;
inputmonth=1:12;

for yearind=1:length(inputyear)
    year = inputyear(yearind); 
    for monthind=1:length(inputmonth)
        month = inputmonth(monthind);
        yearstr = num2str(year,'%04i');
        monthstr = num2str(month,'%02i');
        OGCM_path = '/home/kimyy/Model/CMEMS/monthly/';

        bry_path=['/data1/kimyy/Model/ROMS/nwp_1_10/forcing_matlab/OBC/MyOcean/'];
        % tempday=[num2str(year,'%04i'),num2str(month,'%02i'),num2str(day,'%02i')];

        if (exist(bry_path,'dir')~=7)
            mkdir(bry_path);
        end
        bry_prefix=[bry_path,'MyOcean'];
        bry_suffix = ['.nc'];
        ROMS.bryfilename = [bry_prefix, '_Y',yearstr, 'M',num2str(month), '_', testname, bry_suffix]; % auto01_bndy_20190608.nc

        OGCM_prefix = 'myocean_';



        ROMS.conserv = 0; % same barotropic velocities as the OGCM,
        % % if switch is on, barotropic velocity is obtained by 
        % % divided volume transport by maximum depth --> can cause error

        % % check please
        ROMS.Vtransform      = 2;
        ROMS.Vstretching     = 4;

        %
        % Get the model grid
        %
        % testname = 'test19';
        g = grd('NWP_1_10_linux',testname);
        ROMS.grdname = g.grd_file;
        ROMS.theta_s = g.theta_s;
        ROMS.theta_b = g.theta_b;
        ROMS.hc      = g.hc;
        ROMS.N       = g.N;

        ROMS.title  = ['Automated Model OBC file from MyOcean analysis forecast file'];

        rmdepth     = 0;         % Number of bottom levels to remove
        %(This is usefull when there is no valid data at this level
        %i.e if the depth in the domain is shallower than
        %the OGCM depth)

        % Objective analysis decorrelation scale [m]
        % (if Roa=0: simple extrapolation method; crude but much less costly)
        %
        %Roa=300e3;
        Roa=0;

        interp_method = 'linear';           % Interpolation method: 'linear' or 'cubic'

        ROMS.obc = [1 1 1 1]; % open boundaries (1 = open , [S E N W])

        makeplot = 0;

        nc = netcdf(ROMS.grdname);
        ROMS.lon.rho.all = nc{'lon_rho'}(:);
        ROMS.lon.u.all = nc{'lon_u'}(:);
        ROMS.lon.v.all = nc{'lon_v'}(:);
        ROMS.lat.rho.all = nc{'lat_rho'}(:);
        ROMS.lat.u.all = nc{'lat_u'}(:);
        ROMS.lat.v.all = nc{'lat_v'}(:);
        ROMS.angle = nc{'angle'}(:);
        ROMS.h = nc{'h'}(:);
        close(nc)
        ROMS.hu = rho2u_2d(ROMS.h);
        ROMS.hv = rho2v_2d(ROMS.h);

        ROMS.data(1).direction = 'east';
        ROMS.data(2).direction = 'west';
        ROMS.data(3).direction = 'south';
        ROMS.data(4).direction = 'north';


        %
        %------------------------------------------------------------------------------------
        %
        % Get the OGCM grid

        tmp_OGCM_path=OGCM_path;
        OGCM_path_day=[tmp_OGCM_path,yearstr,'/'];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ddddddd
        % % % need to be different by boundary 
        for diri=1:4
            myocean_grid_name{diri} = [OGCM_path_day, ...
                OGCM_prefix,ROMS.data(diri).direction,'_', yearstr,monthstr '.nc'];  % myocean_direction_20190608.nc
            OGCM.data(diri).lon = ncread(myocean_grid_name{diri}, 'longitude');
            OGCM.data(diri).lat = ncread(myocean_grid_name{diri}, 'latitude');
            OGCM.depth = ncread(myocean_grid_name{diri}, 'depth');
            OGCM.NZ = length(OGCM.depth);
        end
        OGCM.data(1).validgrid = OGCM.data(1).lat;
        OGCM.data(2).validgrid = OGCM.data(2).lat;
        OGCM.data(3).validgrid = OGCM.data(3).lon;
        OGCM.data(4).validgrid = OGCM.data(4).lon;

        % totalday  = 4; %% 2d ago, yesterday, today, tomorrow
        % firstdofy=datenum(['01-Jan-',num2str(year)]); % dd-mmm-yyyy
        % presenttime=datenum([num2str(month),'-',num2str(day),'-',num2str(year)]); % mm-dd-yyyy
        % runtime = presenttime - firstdofy;
        bndy_times = 0; 

        % ROMS.time = 0*(1:totalday+0); %% all zero for initializing
        ROMS.time = bndy_times; %% all zero for initializing

        %
        % Create and open the ROMS files
        %

        create_bryfile_auto(ROMS.bryfilename, ROMS.grdname, ROMS.title, ROMS.obc,...
            ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N,...
            ROMS.time, 0, 'clobber');
        ROMS.brync = netcdf(ROMS.bryfilename, 'write');

        %
        % Perform the interpolations for the current month
        %
        disp(' Current month :')
        disp('================')
        % for tndx_OGCM = 1:totalday
        %     disp([' Time step : ',num2str(tndx_OGCM),' of ',num2str(totalday),' :'])
            interp_OGCM_bndy_MyOcean_sep(OGCM_path_day, Roa, interp_method, ...
                OGCM, year, month, ROMS)
        % end

        %
        % Close the ROMS files
        %
        if ~isempty(ROMS.brync)
            close(ROMS.brync);
        end
        %

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
            figure
            test_bry(ROMS.bryfilename,ROMS.grdname,'temp',1,ROMS.obc,ROMS.Vtransform,ROMS.Vstretching)
            figure
            test_bry(ROMS.bryfilename,ROMS.grdname,'salt',1,ROMS.obc,ROMS.Vtransform,ROMS.Vstretching)
            figure
            test_bry(ROMS.bryfilename,ROMS.grdname,'u',1,ROMS.obc,ROMS.Vtransform,ROMS.Vstretching)
            figure
            test_bry(ROMS.bryfilename,ROMS.grdname,'v',1,ROMS.obc,ROMS.Vtransform,ROMS.Vstretching)
        end
    end
end

% ROMS_compile_bndy_CMIP5;