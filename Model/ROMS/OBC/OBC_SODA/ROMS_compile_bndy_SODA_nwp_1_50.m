clear all; close all; clc;
warning off;
run /home/kimyy/Dropbox/source/matlab/Model/ROMS/Roms_tools/Run/romstools_param;
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
 end


year = 2010;
mm = 1:12;
casename = 'NWP_1_50';
testname = 'test01';
addpath(genpath('/home/kimyy/Dropbox/source/matlab/Common/netcdf_old/'))
quot = '''';

for yind=1:length(year)
    yyyy=year(yind);
    ROMS_title = 'Northwestern Pacific Model';
    if (windows==1)
        bryname = ['E:\Data\Model\ROMS\nwp_1_20\input\test38\','roms_bndy_',casename,'_SODA_3_4_2_',num2str(yyyy), '_',testname,'.nc'];
    elseif (linux==1)
        bryname = ['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_50/forcing_matlab/OBC/SODA_3_4_2/','roms_bndy_',casename,'_SODA_3_4_2_',num2str(yyyy), '_',testname,'.nc'];
    end

    g = grd('NWP_1_50_linux',testname);
    grdname = g.grd_file;
    theta_s = g.theta_s;
    theta_b = g.theta_b;
    hc = g.hc;
    N = g.N;

    roms_time = [15:30:365];
    % if leapyear(yyyy)
    %     cycle_length = 366;
    % else
        cycle_length = 365.25;
    % end

    varis_3d = {'temp', 'salt', 'u', 'v'};
    varis_2d = {'zeta', 'ubar', 'vbar'};
    direction_list = {'east', 'west', 'south', 'north'};
    obc = [1 1 1 1];

    create_bryfile_J(bryname, grdname, ROMS_title, ...
        obc, theta_s, theta_b, hc, N, roms_time, cycle_length, 'clobber');

    for mi = mm(1):mm(end);
        if (windows==1)
            filename = ['E:\Data\Model\ROMS\nwp_1_20\input\test38\','SODA_Y', num2str(yyyy,4), 'M', num2str(mi), '_',testname,'.nc'];
        elseif (linux==1)
            filename = ['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_50/forcing_matlab/OBC/SODA_3_4_2/','SODA_Y', num2str(yyyy,4), 'M', num2str(mi), '_',testname,'.nc'];
        end
            nc = netcdf(filename);

        for i_dir = 1:length(direction_list);
            direction = direction_list{i_dir};

            for i_3d = 1:length(varis_3d);
                var_3d = [varis_3d{i_3d}, '_', direction];
                eval([var_3d, '_all(mi, :, :) = nc{', quot, var_3d, quot, '}(:);'])
            end
            for i_2d = 1:length(varis_2d);
                var_2d = [varis_2d{i_2d}, '_', direction];
                eval([var_2d, '_all(mi, :) = nc{', quot, var_2d, quot, '}(:);'])
            end
        end
        close(nc);
    end

    % %
    % % Compute S coordinates (when vtransform = 1 vstretching =1)
    % %
    % cff1=1./sinh(theta_s);
    % cff2=0.5/tanh(0.5*theta_s);
    % sc_r=((1:N)-N-0.5)/N;
    % Cs_r=(1.-theta_b)*cff1*sinh(theta_s*sc_r)...
    %     +theta_b*(cff2*tanh(theta_s*(sc_r+0.5))-0.5);
    % sc_w=((0:N)-N)/N;
    % Cs_w=(1.-theta_b)*cff1*sinh(theta_s*sc_w)...
    %     +theta_b*(cff2*tanh(theta_s*(sc_w+0.5))-0.5);


    kgrid=0;
    [temp_sc_r,temp_Cs_r]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);
    kgrid=1;
    [temp_sc_w,temp_Cs_w]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);
    sc_r=temp_sc_r';
    Cs_r=temp_Cs_r';
    sc_w=temp_sc_w';
    Cs_w=temp_Cs_w';


    % Write variables
    nc = netcdf(bryname, 'write');
    nc{'spherical'}(:)= 'T';
    nc{'Vtransform'}(:)=Vtransform;
    nc{'Vstretching'}(:)=Vstretching;
    nc{'tstart'}(:) = min([min(roms_time) min(roms_time) min(roms_time)]);
    nc{'tend'}(:) = max([max(roms_time) max(roms_time) max(roms_time)]);
    nc{'theta_s'}(:) = theta_s;
    nc{'theta_b'}(:) = theta_b;
    nc{'Tcline'}(:) = hc;
    nc{'hc'}(:) = hc;
    nc{'s_rho'}(:) = sc_r;
    nc{'s_w'}(:) = sc_w;
    nc{'Cs_r'}(:) = Cs_r;
    nc{'Cs_w'}(:) = Cs_w;
    nc{'tclm_time'}(:) = roms_time;
    nc{'temp_time'}(:) = roms_time;
    nc{'sclm_time'}(:) = roms_time;
    nc{'salt_time'}(:) = roms_time;
    nc{'uclm_time'}(:) = roms_time;
    nc{'vclm_time'}(:) = roms_time;
    nc{'v2d_time'}(:) =  roms_time;
    nc{'v3d_time'}(:) =  roms_time;
    nc{'ssh_time'}(:) =  roms_time;
    nc{'zeta_time'}(:) = roms_time;
    nc{'bry_time'}(:) =  roms_time;
    % South
    nc{'u_south'}(:) = u_south_all;
    nc{'v_south'}(:) = v_south_all;
    nc{'ubar_south'}(:) = ubar_south_all;
    nc{'vbar_south'}(:) = vbar_south_all;
    nc{'zeta_south'}(:) = zeta_south_all;
    nc{'temp_south'}(:) = temp_south_all;
    nc{'salt_south'}(:) = salt_south_all;
    % East
    nc{'u_east'}(:) = u_east_all;
    nc{'v_east'}(:) = v_east_all;
    nc{'ubar_east'}(:) = ubar_east_all;
    nc{'vbar_east'}(:) = vbar_east_all;
    nc{'zeta_east'}(:) = zeta_east_all;
    nc{'temp_east'}(:) = temp_east_all;
    nc{'salt_east'}(:) = salt_east_all;
    % North
    nc{'u_north'}(:) = u_north_all;
    nc{'v_north'}(:) = v_north_all;
    nc{'ubar_north'}(:) = ubar_north_all;
    nc{'vbar_north'}(:) = vbar_north_all;
    nc{'zeta_north'}(:) = zeta_north_all;
    nc{'temp_north'}(:) = temp_north_all;
    nc{'salt_north'}(:) = salt_north_all;
    % West
    nc{'u_west'}(:) = u_west_all;
    nc{'v_west'}(:) = v_west_all;
    nc{'ubar_west'}(:) = ubar_west_all;
    nc{'vbar_west'}(:) = vbar_west_all;
    nc{'zeta_west'}(:) = zeta_west_all;
    nc{'temp_west'}(:) = temp_west_all;
    nc{'salt_west'}(:) = salt_west_all;
    close(nc)
end