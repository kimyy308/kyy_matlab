close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
% all_region ={'YS', 'SS'}

% all_var ={'temp', 'salt', 'u', 'v', 'w', 'zeta', 'swrad', 'shflux', 'Uwind', 'Vwind'};
% all_var ={'temp', 'salt', 'u', 'v', 'zeta'}

% for varind=1:length(all_var)
    clearvars '*' -except varind all_var
    
    % % % 
    % % % Read Model variable at surface, bottom, 
    % % % boundary(north, south, east, west)

    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        dropboxpath='C:\Users\KYY\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
        addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
        addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end

%     shadlev = [0 35];
%     rms_shadlev = [0 4];
%     bias_shadlev = [-4 4];
%     conlev  = 0:5:35;
    dl=1/20;
    % for snu_desktop
    testname='test52'   % % need to change
    inputyear = [1980 : 2008]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

%     varname =all_var{varind};  %% reference variable -> temperature
    
    run('nwp_polygon_point.m');
%     regionname=all_region{regionind}
    regionname='NWP'
    switch(regionname)
        case('NWP') %% North western Pacific
            lonlat = [115, 164, 15, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('ES') %% East Sea
            refpolygon=espolygon;
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        otherwise
            ('?')
    end
    lonlat(1)=min(refpolygon(:,1));
    lonlat(2)=max(refpolygon(:,1));
    lonlat(3)=min(refpolygon(:,2));
    lonlat(4)=max(refpolygon(:,2));
    
%     switch(varname)
%         case('zeta') %% North western Pacific
%             ndim=3;
%             data_units = 'm';
%         case('temp') %% North western Pacific
%             ndim=4;
%             data_units = 'Celsius'
%         otherwise
%             ndim=4;
%             data_units = '?'
%     end
    
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\nwp_1_20\',testname,'\BNDY\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\input\', testname, '\'); % % where data files are
        avhrrdir='E:\Data\Observation\OISST\monthly\';
    elseif (strcmp(system_name,'GLNXA64'))
    end

    run(param_script);
    ind=1;
    for yearij = 1:length(inputyear)
        for monthij = 1:length(inputmonth)
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            tic;
            tempyear = inputyear(yearij);
            tempmonth = inputmonth(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
            outputname = strcat(filedir, '..\..\', testname, '\run\', num2str(tempyear, '%04i'), '\', ...
                testname, '_monthly_', num2str(tempyear, '%04i'), '_', num2str(tempmonth, '%02i'), '.nc');
            bndyname = strcat(filedir, ...
                    'roms_bndy_nwp_1_20_', num2str(tempyear,'%04i'), '_', testname, '.nc');
            % read model data
            if (exist('lon')==0)
                modelinfo=ncinfo(outputname);
%                 if (strcmp(varname,'u')==1)
                    lon_u = ncread(outputname,'lon_u',[1 1],[modelinfo.Dimensions(7).Length,1]);
                    lat_u = ncread(outputname,'lat_u',[1 1],[1,modelinfo.Dimensions(6).Length]);
%                 elseif (strcmp(varname,'v')==1)
                    lon_v = ncread(outputname,'lon_v',[1 1],[modelinfo.Dimensions(9).Length,1]);
                    lat_v = ncread(outputname,'lat_v',[1 1],[1,modelinfo.Dimensions(8).Length]);
%                 else
                    lon_rho = ncread(outputname,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                    lat_rho = ncread(outputname,'lat_rho',[1 1],[1,modelinfo.Dimensions(4).Length]);
%                 end
                s_rho = ncread(outputname,'s_rho');

                lon_west = abs(lon_rho - (lonlat(1)-1));
                min_lon_west=min(lon_west);
                lon_east = abs(lon_rho - (lonlat(2)+1));
                min_lon_east=min(lon_east);
                lat_south = abs(lat_rho - (lonlat(3)-1));
                min_lat_south=min(lat_south);
                lat_north = abs(lat_rho - (lonlat(4)+1));
                min_lat_north=min(lat_north);

                lon_min = find(lon_west == min_lon_west);
                lon_max = find(lon_east == min_lon_east);
                lat_min = find(lat_south == min_lat_south);
                lat_max = find(lat_north == min_lat_north);
                
                h_rho = ncread(outputname, 'h', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                h_u=rho2u_2d(h_rho);
                h_v=rho2v_2d(h_rho);
                
%                 if (strcmp(varname,'u')==1)
                    lon_u = ncread(outputname,'lon_u', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1) lat_max(1)-lat_min(1)+1]);
                    lat_u = ncread(outputname,'lat_u', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1) lat_max(1)-lat_min(1)+1]);
%                 elseif (strcmp(varname,'v')==1)
                    lon_v = ncread(outputname,'lon_v', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)]);
                    lat_v = ncread(outputname,'lat_v', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)]);
%                 else
                    lon_rho = ncread(outputname,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                    lat_rho = ncread(outputname,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
%                 end
                
    %             polygon_ind=NaN(size(refpolygon,1),2);
    %             for i=1:size(refpolygon,1)
    %                 [polygon_ind(i,1), trash_ind, polygon_ind(i,2), trash_ind]=findind_Y(dl, [refpolygon(i,1),refpolygon(i,2)],lon',lat');
    %             end
    %             mask_model = inpolygon(lon,lat,polygon_ind(:,1),polygon_ind(:,2));
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_model_rho(1:size(lon_rho,1),1:size(lon_rho,2))=1;
                        mask_model_rho(mask_model_rho==0)=NaN;
                        mask_model_u(1:size(lon_u,1),1:size(lon_u,2))=1;
                        mask_model_u(mask_model_u==0)=NaN;
                        mask_model_v(1:size(lon_v,1),1:size(lon_v,2))=1;
                        mask_model_v(mask_model_v==0)=NaN;
                    otherwise
                        mask_model_rho = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
                        mask_model_rho(mask_model_rho==0)=NaN;
                        mask_model_u = double(inpolygon(lon_u,lat_u,refpolygon(:,1),refpolygon(:,2)));
                        mask_model_u(mask_model_u==0)=NaN;
                        mask_model_v = double(inpolygon(lon_v,lat_v,refpolygon(:,1),refpolygon(:,2)));
                        mask_model_v(mask_model_v==0)=NaN;
                end
                n_z = modelinfo.Dimensions(2).Length;
                
                
% % %                 initialize zeta
                comb_spatial_mean_north_zeta=(zeros([size(lon_rho,1), 12]));
                comb_spatial_mean_south_zeta=(zeros([size(lon_rho,1), 12]));
                comb_spatial_mean_east_zeta=(zeros([size(lat_rho,2), 12]));
                comb_spatial_mean_west_zeta=(zeros([size(lat_rho,2), 12]));
                comb_spatial_mean_north_zeta_inner=(zeros([size(lon_rho,1), 12]));
                comb_spatial_mean_south_zeta_inner=(zeros([size(lon_rho,1), 12]));
                comb_spatial_mean_east_zeta_inner=(zeros([size(lat_rho,2), 12]));
                comb_spatial_mean_west_zeta_inner=(zeros([size(lat_rho,2), 12]));
                comb_spatial_mean_north_zeta_bndy=(zeros([size(lon_rho,1), 12]));
                comb_spatial_mean_south_zeta_bndy=(zeros([size(lon_rho,1), 12]));
                comb_spatial_mean_east_zeta_bndy=(zeros([size(lat_rho,2), 12]));
                comb_spatial_mean_west_zeta_bndy=(zeros([size(lat_rho,2), 12]));
% % %                 initialize ubar
                comb_spatial_mean_north_ubar=(zeros([size(lon_u,1), 12]));
                comb_spatial_mean_south_ubar=(zeros([size(lon_u,1), 12]));
                comb_spatial_mean_east_ubar=(zeros([size(lat_u,2), 12]));
                comb_spatial_mean_west_ubar=(zeros([size(lat_u,2), 12]));
                comb_spatial_mean_north_ubar_inner=(zeros([size(lon_u,1), 12]));
                comb_spatial_mean_south_ubar_inner=(zeros([size(lon_u,1), 12]));
                comb_spatial_mean_east_ubar_inner=(zeros([size(lat_u,2), 12]));
                comb_spatial_mean_west_ubar_inner=(zeros([size(lat_u,2), 12]));
                comb_spatial_mean_north_ubar_bndy=(zeros([size(lon_u,1), 12]));
                comb_spatial_mean_south_ubar_bndy=(zeros([size(lon_u,1), 12]));
                comb_spatial_mean_east_ubar_bndy=(zeros([size(lat_u,2), 12]));
                comb_spatial_mean_west_ubar_bndy=(zeros([size(lat_u,2), 12]));
% % %                 initialize vbar
                comb_spatial_mean_north_vbar=(zeros([size(lon_v,1), 12]));
                comb_spatial_mean_south_vbar=(zeros([size(lon_v,1), 12]));
                comb_spatial_mean_east_vbar=(zeros([size(lat_v,2), 12]));
                comb_spatial_mean_west_vbar=(zeros([size(lat_v,2), 12]));
                comb_spatial_mean_north_vbar_inner=(zeros([size(lon_v,1), 12]));
                comb_spatial_mean_south_vbar_inner=(zeros([size(lon_v,1), 12]));
                comb_spatial_mean_east_vbar_inner=(zeros([size(lat_v,2), 12]));
                comb_spatial_mean_west_vbar_inner=(zeros([size(lat_v,2), 12]));
                comb_spatial_mean_north_vbar_bndy=(zeros([size(lon_v,1), 12]));
                comb_spatial_mean_south_vbar_bndy=(zeros([size(lon_v,1), 12]));
                comb_spatial_mean_east_vbar_bndy=(zeros([size(lat_v,2), 12]));
                comb_spatial_mean_west_vbar_bndy=(zeros([size(lat_v,2), 12]));

% % %                 initialize temp
                comb_spatial_mean_north_temp=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_south_temp=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_east_temp=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_west_temp=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_north_temp_inner=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_south_temp_inner=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_east_temp_inner=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_west_temp_inner=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_north_temp_bndy=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_south_temp_bndy=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_east_temp_bndy=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_west_temp_bndy=(zeros([size(lat_rho,2), n_z, 12]));
% % %                 initialize salt
                comb_spatial_mean_north_salt=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_south_salt=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_east_salt=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_west_salt=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_north_salt_inner=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_south_salt_inner=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_east_salt_inner=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_west_salt_inner=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_north_salt_bndy=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_south_salt_bndy=(zeros([size(lon_rho,1), n_z, 12]));
                comb_spatial_mean_east_salt_bndy=(zeros([size(lat_rho,2), n_z, 12]));
                comb_spatial_mean_west_salt_bndy=(zeros([size(lat_rho,2), n_z, 12]));
% % %                 initialize u
                comb_spatial_mean_north_u=(zeros([size(lon_u,1), n_z, 12]));
                comb_spatial_mean_south_u=(zeros([size(lon_u,1), n_z, 12]));
                comb_spatial_mean_east_u=(zeros([size(lat_u,2), n_z, 12]));
                comb_spatial_mean_west_u=(zeros([size(lat_u,2), n_z, 12]));
                comb_spatial_mean_north_u_inner=(zeros([size(lon_u,1), n_z, 12]));
                comb_spatial_mean_south_u_inner=(zeros([size(lon_u,1), n_z, 12]));
                comb_spatial_mean_east_u_inner=(zeros([size(lat_u,2), n_z, 12]));
                comb_spatial_mean_west_u_inner=(zeros([size(lat_u,2), n_z, 12]));
                comb_spatial_mean_north_u_bndy=(zeros([size(lon_u,1), n_z, 12]));
                comb_spatial_mean_south_u_bndy=(zeros([size(lon_u,1), n_z, 12]));
                comb_spatial_mean_east_u_bndy=(zeros([size(lat_u,2), n_z, 12]));
                comb_spatial_mean_west_u_bndy=(zeros([size(lat_u,2), n_z, 12]));
% % %                 initialize v
                comb_spatial_mean_north_v=(zeros([size(lon_v,1), n_z, 12]));
                comb_spatial_mean_south_v=(zeros([size(lon_v,1), n_z, 12]));
                comb_spatial_mean_east_v=(zeros([size(lat_v,2), n_z, 12]));
                comb_spatial_mean_west_v=(zeros([size(lat_v,2), n_z, 12]));
                comb_spatial_mean_north_v_inner=(zeros([size(lon_v,1), n_z, 12]));
                comb_spatial_mean_south_v_inner=(zeros([size(lon_v,1), n_z, 12]));
                comb_spatial_mean_east_v_inner=(zeros([size(lat_v,2), n_z, 12]));
                comb_spatial_mean_west_v_inner=(zeros([size(lat_v,2), n_z, 12]));
                comb_spatial_mean_north_v_bndy=(zeros([size(lon_v,1), n_z, 12]));
                comb_spatial_mean_south_v_bndy=(zeros([size(lon_v,1), n_z, 12]));
                comb_spatial_mean_east_v_bndy=(zeros([size(lat_v,2), n_z, 12]));
                comb_spatial_mean_west_v_bndy=(zeros([size(lat_v,2), n_z, 12]));
            end
            
%             data_info = ncinfo(outputname, varname);  %% [lon lat depth time]
            

% % %                 get zeta
                north_zeta = squeeze(ncread(outputname,'zeta',[lon_min(1) lat_max(1) 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 920, 1]
                north_zeta_inner = squeeze(ncread(outputname,'zeta',[lon_min(1) lat_max(1)-2 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 918, 1]
                north_zeta_bndy = squeeze(ncread(bndyname,['zeta','_north'],[1 tempmonth], [lon_max(1)-lon_min(1)+1 1])); %% NWP : [1:980, 1]
                south_zeta = squeeze(ncread(outputname,'zeta',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 1, 1]
                south_zeta_inner = squeeze(ncread(outputname,'zeta',[lon_min(1) lat_min(1)+2 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 3, 1]
                south_zeta_bndy = squeeze(ncread(bndyname,['zeta','_south'],[1 tempmonth], [lon_max(1)-lon_min(1)+1 1])); %% NWP : [1:980, 1]
                east_zeta = squeeze(ncread(outputname,'zeta',[lon_max(1) lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [980, 1:920, 1]
                east_zeta_inner = squeeze(ncread(outputname,'zeta',[lon_max(1)-2 lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [978, 1:920, 1]
                east_zeta_bndy = squeeze(ncread(bndyname,['zeta','_east'],[1 tempmonth], [lat_max(1)-lat_min(1)+1 1])); %% NWP : [1:920, 1]
                west_zeta = squeeze(ncread(outputname,'zeta',[lon_min(1) lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [1, 1:920, 1]
                west_zeta_inner = squeeze(ncread(outputname,'zeta',[lon_min(1)+2 lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [3, 1:920, 1]
                west_zeta_bndy = squeeze(ncread(bndyname,['zeta','_west'],[1 tempmonth], [lat_max(1)-lat_min(1)+1 1])); %% NWP : [1:920, 1]
% % %                 get ubar (0.5 grid inner)
                north_ubar = squeeze(ncread(outputname,'ubar',[lon_min(1) lat_max(1) 1], [lon_max(1)-lon_min(1) 1 1])); %% NWP : [1:979, 920, 1]
                north_ubar_inner = squeeze(ncread(outputname,'ubar',[lon_min(1) lat_max(1)-2 1], [lon_max(1)-lon_min(1) 1 1])); %% NWP : [1:979, 918, 1]
                north_ubar_bndy = squeeze(ncread(bndyname,['ubar','_north'],[1 tempmonth], [lon_max(1)-lon_min(1) 1])); %% NWP : [1:979, 1]
                south_ubar = squeeze(ncread(outputname,'ubar',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1) 1 1])); %% NWP : [1:979, 1, 1]
                south_ubar_inner = squeeze(ncread(outputname,'ubar',[lon_min(1) lat_min(1)+2 1], [lon_max(1)-lon_min(1) 1 1])); %% NWP : [1:979, 3, 1]
                south_ubar_bndy = squeeze(ncread(bndyname,['ubar','_south'],[1 tempmonth], [lon_max(1)-lon_min(1) 1])); %% NWP : [1:979, 1]
                east_ubar = squeeze(ncread(outputname,'ubar',[lon_max(1)-1 lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [979, 1:920, 1]
                east_ubar_inner = squeeze(ncread(outputname,'ubar',[lon_max(1)-3 lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [977, 1:920, 1]
                east_ubar_bndy = squeeze(ncread(bndyname,['ubar','_east'],[1 tempmonth], [lat_max(1)-lat_min(1)+1 1])); %% NWP : [1:920, 1]
                west_ubar = squeeze(ncread(outputname,'ubar',[lon_min(1) lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [1, 1:920, 1]
                west_ubar_inner = squeeze(ncread(outputname,'ubar',[lon_min(1)+2 lat_min(1) 1],  [1 lat_max(1)-lat_min(1)+1 1])); %% NWP : [3, 1:920, 1]
                west_ubar_bndy = squeeze(ncread(bndyname,['ubar','_west'],[1 tempmonth], [lat_max(1)-lat_min(1)+1 1])); %% NWP : [1:920, 1]                
% % %                 get vbar (0.5 grid inner)
                north_vbar = squeeze(ncread(outputname,'vbar',[lon_min(1) lat_max(1)-1 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 919, 1]
                north_vbar_inner = squeeze(ncread(outputname,'vbar',[lon_min(1) lat_max(1)-3 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 917, 1]
                north_vbar_bndy = squeeze(ncread(bndyname,['vbar','_north'],[1 tempmonth], [lon_max(1)-lon_min(1)+1 1])); %% NWP : [1:980, 1]
                south_vbar = squeeze(ncread(outputname,'vbar',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 1, 1]
                south_vbar_inner = squeeze(ncread(outputname,'vbar',[lon_min(1) lat_min(1)+2 1], [lon_max(1)-lon_min(1)+1 1 1])); %% NWP : [1:980, 3, 1]
                south_vbar_bndy = squeeze(ncread(bndyname,['vbar','_south'],[1 tempmonth], [lon_max(1)-lon_min(1)+1 1])); %% NWP : [1:980, 1]
                east_vbar = squeeze(ncread(outputname,'vbar',[lon_max(1) lat_min(1) 1],  [1 lat_max(1)-lat_min(1) 1])); %% NWP : [980, 1:919, 1]
                east_vbar_inner = squeeze(ncread(outputname,'vbar',[lon_max(1)-2 lat_min(1) 1],  [1 lat_max(1)-lat_min(1) 1])); %% NWP : [978, 1:919, 1]
                east_vbar_bndy = squeeze(ncread(bndyname,['vbar','_east'],[1 tempmonth], [lat_max(1)-lat_min(1) 1])); %% NWP : [1:919, 1]
                west_vbar = squeeze(ncread(outputname,'vbar',[lon_min(1) lat_min(1) 1],  [1 lat_max(1)-lat_min(1) 1])); %% NWP : [1, 1:919, 1]
                west_vbar_inner = squeeze(ncread(outputname,'vbar',[lon_min(1)+2 lat_min(1) 1],  [1 lat_max(1)-lat_min(1) 1])); %% NWP : [3, 1:919, 1]
                west_vbar_bndy = squeeze(ncread(bndyname,['vbar','_west'],[1 tempmonth], [lat_max(1)-lat_min(1) 1])); %% NWP : [1:919, 1]     
% % %                 get temp
                north_temp = squeeze(ncread(outputname,'temp',[lon_min(1) lat_max(1) 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 920, 1:40, 1]
                north_temp_inner = squeeze(ncread(outputname,'temp',[lon_min(1) lat_max(1)-2 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 918, 1:40, 1]
                north_temp_bndy = squeeze(ncread(bndyname,['temp','_north'],[1 1 tempmonth], [lon_max(1)-lon_min(1)+1 n_z 1])); %% NWP : [1:980, 1:40, 1]
                south_temp = squeeze(ncread(outputname,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 1, 1:40, 1]
                south_temp_inner = squeeze(ncread(outputname,'temp',[lon_min(1) lat_min(1)+2 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 3, 1:40, 1]
                south_temp_bndy = squeeze(ncread(bndyname,['temp','_south'],[1 1 tempmonth], [lon_max(1)-lon_min(1)+1 n_z 1])); %% NWP : [1:980, 1:40, 1]
                east_temp = squeeze(ncread(outputname,'temp',[lon_max(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [980, 1:920, 1:40, 1]
                east_temp_inner = squeeze(ncread(outputname,'temp',[lon_max(1)-2 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [978, 1:920, 1:40, 1]
                east_temp_bndy = squeeze(ncread(bndyname,['temp','_east'],[1 1 tempmonth], [lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1:920, 1:40, 1]
                west_temp = squeeze(ncread(outputname,'temp',[lon_min(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1, 1:920, 1:40, 1]
                west_temp_inner = squeeze(ncread(outputname,'temp',[lon_min(1)+2 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [3, 1:920, 1:40, 1]
                west_temp_bndy = squeeze(ncread(bndyname,['temp','_west'],[1 1 tempmonth], [lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1:920, 1:40, 1]
% % %                 get salt
                north_salt = squeeze(ncread(outputname,'salt',[lon_min(1) lat_max(1) 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 920, 1:40, 1]
                north_salt_inner = squeeze(ncread(outputname,'salt',[lon_min(1) lat_max(1)-2 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 918, 1:40, 1]
                north_salt_bndy = squeeze(ncread(bndyname,['salt','_north'],[1 1 tempmonth], [lon_max(1)-lon_min(1)+1 n_z 1])); %% NWP : [1:980, 1:40, 1]
                south_salt = squeeze(ncread(outputname,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 1, 1:40, 1]
                south_salt_inner = squeeze(ncread(outputname,'salt',[lon_min(1) lat_min(1)+2 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 3, 1:40, 1]
                south_salt_bndy = squeeze(ncread(bndyname,['salt','_south'],[1 1 tempmonth], [lon_max(1)-lon_min(1)+1 n_z 1])); %% NWP : [1:980, 1:40, 1]
                east_salt = squeeze(ncread(outputname,'salt',[lon_max(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [980, 1:920, 1:40, 1]
                east_salt_inner = squeeze(ncread(outputname,'salt',[lon_max(1)-2 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [978, 1:920, 1:40, 1]
                east_salt_bndy = squeeze(ncread(bndyname,['salt','_east'],[1 1 tempmonth], [lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1:920, 1:40, 1]
                west_salt = squeeze(ncread(outputname,'salt',[lon_min(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1, 1:920, 1:40, 1]
                west_salt_inner = squeeze(ncread(outputname,'salt',[lon_min(1)+2 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [3, 1:920, 1:40, 1]
                west_salt_bndy = squeeze(ncread(bndyname,['salt','_west'],[1 1 tempmonth], [lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1:920, 1:40, 1]   
% % %                 get u
                north_u = squeeze(ncread(outputname,'u',[lon_min(1) lat_max(1) 1 1], [lon_max(1)-lon_min(1) 1 n_z 1])); %% NWP : [1:979, 920, 1:40, 1]
                north_u_inner = squeeze(ncread(outputname,'u',[lon_min(1) lat_max(1)-2 1 1], [lon_max(1)-lon_min(1) 1 n_z 1])); %% NWP : [1:979, 918, 1:40, 1]
                north_u_bndy = squeeze(ncread(bndyname,['u','_north'],[1 1 tempmonth], [lon_max(1)-lon_min(1) n_z 1])); %% NWP : [1:979, 1:40, 1]
                south_u = squeeze(ncread(outputname,'u',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1) 1 n_z 1])); %% NWP : [1:979, 1, 1:40, 1]
                south_u_inner = squeeze(ncread(outputname,'u',[lon_min(1) lat_min(1)+2 1 1], [lon_max(1)-lon_min(1) 1 n_z 1])); %% NWP : [1:979, 3, 1:40, 1]
                south_u_bndy = squeeze(ncread(bndyname,['u','_south'],[1 1 tempmonth], [lon_max(1)-lon_min(1) n_z 1])); %% NWP : [1:979, 1:40, 1]
                east_u = squeeze(ncread(outputname,'u',[lon_max(1)-1 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [979, 1:920, 1:40, 1]
                east_u_inner = squeeze(ncread(outputname,'u',[lon_max(1)-3 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [977, 1:920, 1:40, 1]
                east_u_bndy = squeeze(ncread(bndyname,['u','_east'],[1 1 tempmonth], [lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1:920, 1:40, 1]
                west_u = squeeze(ncread(outputname,'u',[lon_min(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1, 1:920, 1:40, 1]
                west_u_inner = squeeze(ncread(outputname,'u',[lon_min(1)+2 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [3, 1:920, 1:40, 1]
                west_u_bndy = squeeze(ncread(bndyname,['u','_west'],[1 1 tempmonth], [lat_max(1)-lat_min(1)+1 n_z 1])); %% NWP : [1:920, 1:40, 1]                
% % %                 get v
                north_v = squeeze(ncread(outputname,'v',[lon_min(1) lat_max(1)-1 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 919, 1:40, 1]
                north_v_inner = squeeze(ncread(outputname,'v',[lon_min(1) lat_max(1)-3 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 917, 1:40, 1]
                north_v_bndy = squeeze(ncread(bndyname,['v','_north'],[1 1 tempmonth], [lon_max(1)-lon_min(1)+1 n_z 1])); %% NWP : [1:980, 1:40, 1]
                south_v = squeeze(ncread(outputname,'v',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 1, 1:40, 1]
                south_v_inner = squeeze(ncread(outputname,'v',[lon_min(1) lat_min(1)+2 1 1], [lon_max(1)-lon_min(1)+1 1 n_z 1])); %% NWP : [1:980, 3, 1:40, 1]
                south_v_bndy = squeeze(ncread(bndyname,['v','_south'],[1 1 tempmonth], [lon_max(1)-lon_min(1)+1 n_z 1])); %% NWP : [1:980, 1:40, 1]
                east_v = squeeze(ncread(outputname,'v',[lon_max(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1) n_z 1])); %% NWP : [980, 1:919, 1:40, 1]
                east_v_inner = squeeze(ncread(outputname,'v',[lon_max(1)-2 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1) n_z 1])); %% NWP : [978, 1:919, 1:40, 1]
                east_v_bndy = squeeze(ncread(bndyname,['v','_east'],[1 1 tempmonth], [lat_max(1)-lat_min(1) n_z 1])); %% NWP : [1:919, 1:40, 1]
                west_v = squeeze(ncread(outputname,'v',[lon_min(1) lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1) n_z 1])); %% NWP : [1, 1:919, 1:40, 1]
                west_v_inner = squeeze(ncread(outputname,'v',[lon_min(1)+2 lat_min(1) 1 1],  [1 lat_max(1)-lat_min(1) n_z 1])); %% NWP : [3, 1:919, 1:40, 1]
                west_v_bndy = squeeze(ncread(bndyname,['v','_west'],[1 1 tempmonth], [lat_max(1)-lat_min(1) n_z 1])); %% NWP : [1:919, 1:40, 1]                   

            
            len_lon_rho_model = size(lon_rho,1);
            len_lat_rho_model = size(lat_rho,2);
            len_lon_u_model = size(lon_u,1);
            len_lat_u_model = size(lat_u,2);
            len_lon_v_model = size(lon_v,1);
            len_lat_v_model = size(lat_v,2);

%             meanbias = mean(mean(bias,'omitnan'),'omitnan');
%             meanrms = mean(mean(rms,'omitnan'),'omitnan');
            
%             comb_surf_data(:,:,ind)=surf_data;
%             comb_bot_data(:,:,ind)=bot_data;
            
%             comb_spatial_mean_surf(:,:,monthij)=comb_spatial_mean_surf(:,:,monthij)+surf_data/double(length(inputyear));  
%             comb_spatial_mean_bot(:,:,monthij)=comb_spatial_mean_bot(:,:,monthij)+bot_data/double(length(inputyear)); 

% % %             put zeta
            comb_north_zeta(:,ind)=north_zeta;
            comb_south_zeta(:,ind)=south_zeta;
            comb_east_zeta(:,ind)=east_zeta;
            comb_west_zeta(:,ind)=west_zeta;
            comb_north_zeta_inner(:,ind)=north_zeta_inner;
            comb_south_zeta_inner(:,ind)=south_zeta_inner;
            comb_east_zeta_inner(:,ind)=east_zeta_inner;
            comb_west_zeta_inner(:,ind)=west_zeta_inner;
            comb_north_zeta_bndy(:,ind)=north_zeta_bndy;
            comb_south_zeta_bndy(:,ind)=south_zeta_bndy;
            comb_east_zeta_bndy(:,ind)=east_zeta_bndy;
            comb_west_zeta_bndy(:,ind)=west_zeta_bndy;
            comb_spatial_mean_north_zeta(:,monthij)=comb_spatial_mean_north_zeta(:,monthij)+north_zeta/double(length(inputyear));  
            comb_spatial_mean_south_zeta(:,monthij)=comb_spatial_mean_south_zeta(:,monthij)+south_zeta/double(length(inputyear));  
            comb_spatial_mean_east_zeta(:,monthij)=comb_spatial_mean_east_zeta(:,monthij)+east_zeta'/double(length(inputyear));  
            comb_spatial_mean_west_zeta(:,monthij)=comb_spatial_mean_west_zeta(:,monthij)+west_zeta'/double(length(inputyear));
            comb_spatial_mean_north_zeta_inner(:,monthij)=comb_spatial_mean_north_zeta_inner(:,monthij)+north_zeta_inner/double(length(inputyear));  
            comb_spatial_mean_south_zeta_inner(:,monthij)=comb_spatial_mean_south_zeta_inner(:,monthij)+south_zeta_inner/double(length(inputyear));  
            comb_spatial_mean_east_zeta_inner(:,monthij)=comb_spatial_mean_east_zeta_inner(:,monthij)+east_zeta_inner'/double(length(inputyear));  
            comb_spatial_mean_west_zeta_inner(:,monthij)=comb_spatial_mean_west_zeta_inner(:,monthij)+west_zeta_inner'/double(length(inputyear));
            comb_spatial_mean_north_zeta_bndy(:,monthij)=comb_spatial_mean_north_zeta_bndy(:,monthij)+north_zeta_bndy/double(length(inputyear));  
            comb_spatial_mean_south_zeta_bndy(:,monthij)=comb_spatial_mean_south_zeta_bndy(:,monthij)+south_zeta_bndy/double(length(inputyear));  
            comb_spatial_mean_east_zeta_bndy(:,monthij)=comb_spatial_mean_east_zeta_bndy(:,monthij)+east_zeta_bndy/double(length(inputyear));  
            comb_spatial_mean_west_zeta_bndy(:,monthij)=comb_spatial_mean_west_zeta_bndy(:,monthij)+west_zeta_bndy/double(length(inputyear));
% % %             put ubar
            comb_north_ubar(:,ind)=north_ubar;
            comb_south_ubar(:,ind)=south_ubar;
            comb_east_ubar(:,ind)=east_ubar;
            comb_west_ubar(:,ind)=west_ubar;
            comb_north_ubar_inner(:,ind)=north_ubar_inner;
            comb_south_ubar_inner(:,ind)=south_ubar_inner;
            comb_east_ubar_inner(:,ind)=east_ubar_inner;
            comb_west_ubar_inner(:,ind)=west_ubar_inner;
            comb_north_ubar_bndy(:,ind)=north_ubar_bndy;
            comb_south_ubar_bndy(:,ind)=south_ubar_bndy;
            comb_east_ubar_bndy(:,ind)=east_ubar_bndy;
            comb_west_ubar_bndy(:,ind)=west_ubar_bndy;
            comb_spatial_mean_north_ubar(:,monthij)=comb_spatial_mean_north_ubar(:,monthij)+north_ubar/double(length(inputyear));  
            comb_spatial_mean_south_ubar(:,monthij)=comb_spatial_mean_south_ubar(:,monthij)+south_ubar/double(length(inputyear));  
            comb_spatial_mean_east_ubar(:,monthij)=comb_spatial_mean_east_ubar(:,monthij)+east_ubar'/double(length(inputyear));  
            comb_spatial_mean_west_ubar(:,monthij)=comb_spatial_mean_west_ubar(:,monthij)+west_ubar'/double(length(inputyear));
            comb_spatial_mean_north_ubar_inner(:,monthij)=comb_spatial_mean_north_ubar_inner(:,monthij)+north_ubar_inner/double(length(inputyear));  
            comb_spatial_mean_south_ubar_inner(:,monthij)=comb_spatial_mean_south_ubar_inner(:,monthij)+south_ubar_inner/double(length(inputyear));  
            comb_spatial_mean_east_ubar_inner(:,monthij)=comb_spatial_mean_east_ubar_inner(:,monthij)+east_ubar_inner'/double(length(inputyear));  
            comb_spatial_mean_west_ubar_inner(:,monthij)=comb_spatial_mean_west_ubar_inner(:,monthij)+west_ubar_inner'/double(length(inputyear));
            comb_spatial_mean_north_ubar_bndy(:,monthij)=comb_spatial_mean_north_ubar_bndy(:,monthij)+north_ubar_bndy/double(length(inputyear));  
            comb_spatial_mean_south_ubar_bndy(:,monthij)=comb_spatial_mean_south_ubar_bndy(:,monthij)+south_ubar_bndy/double(length(inputyear));  
            comb_spatial_mean_east_ubar_bndy(:,monthij)=comb_spatial_mean_east_ubar_bndy(:,monthij)+east_ubar_bndy/double(length(inputyear));  
            comb_spatial_mean_west_ubar_bndy(:,monthij)=comb_spatial_mean_west_ubar_bndy(:,monthij)+west_ubar_bndy/double(length(inputyear));
% % %             put vbar
            comb_north_vbar(:,ind)=north_vbar;
            comb_south_vbar(:,ind)=south_vbar;
            comb_east_vbar(:,ind)=east_vbar;
            comb_west_vbar(:,ind)=west_vbar;
            comb_north_vbar_inner(:,ind)=north_vbar_inner;
            comb_south_vbar_inner(:,ind)=south_vbar_inner;
            comb_east_vbar_inner(:,ind)=east_vbar_inner;
            comb_west_vbar_inner(:,ind)=west_vbar_inner;
            comb_north_vbar_bndy(:,ind)=north_vbar_bndy;
            comb_south_vbar_bndy(:,ind)=south_vbar_bndy;
            comb_east_vbar_bndy(:,ind)=east_vbar_bndy;
            comb_west_vbar_bndy(:,ind)=west_vbar_bndy;
            comb_spatial_mean_north_vbar(:,monthij)=comb_spatial_mean_north_vbar(:,monthij)+north_vbar/double(length(inputyear));  
            comb_spatial_mean_south_vbar(:,monthij)=comb_spatial_mean_south_vbar(:,monthij)+south_vbar/double(length(inputyear));  
            comb_spatial_mean_east_vbar(:,monthij)=comb_spatial_mean_east_vbar(:,monthij)+east_vbar'/double(length(inputyear));  
            comb_spatial_mean_west_vbar(:,monthij)=comb_spatial_mean_west_vbar(:,monthij)+west_vbar'/double(length(inputyear));
            comb_spatial_mean_north_vbar_inner(:,monthij)=comb_spatial_mean_north_vbar_inner(:,monthij)+north_vbar_inner/double(length(inputyear));  
            comb_spatial_mean_south_vbar_inner(:,monthij)=comb_spatial_mean_south_vbar_inner(:,monthij)+south_vbar_inner/double(length(inputyear));  
            comb_spatial_mean_east_vbar_inner(:,monthij)=comb_spatial_mean_east_vbar_inner(:,monthij)+east_vbar_inner'/double(length(inputyear));  
            comb_spatial_mean_west_vbar_inner(:,monthij)=comb_spatial_mean_west_vbar_inner(:,monthij)+west_vbar_inner'/double(length(inputyear));
            comb_spatial_mean_north_vbar_bndy(:,monthij)=comb_spatial_mean_north_vbar_bndy(:,monthij)+north_vbar_bndy/double(length(inputyear));  
            comb_spatial_mean_south_vbar_bndy(:,monthij)=comb_spatial_mean_south_vbar_bndy(:,monthij)+south_vbar_bndy/double(length(inputyear));  
            comb_spatial_mean_east_vbar_bndy(:,monthij)=comb_spatial_mean_east_vbar_bndy(:,monthij)+east_vbar_bndy/double(length(inputyear));  
            comb_spatial_mean_west_vbar_bndy(:,monthij)=comb_spatial_mean_west_vbar_bndy(:,monthij)+west_vbar_bndy/double(length(inputyear));
% % %             put temp
            comb_north_temp(:,:,ind)=north_temp;
            comb_south_temp(:,:,ind)=south_temp;
            comb_east_temp(:,:,ind)=east_temp;
            comb_west_temp(:,:,ind)=west_temp;
            comb_north_temp_inner(:,:,ind)=north_temp_inner;
            comb_south_temp_inner(:,:,ind)=south_temp_inner;
            comb_east_temp_inner(:,:,ind)=east_temp_inner;
            comb_west_temp_inner(:,:,ind)=west_temp_inner;
            comb_north_temp_bndy(:,:,ind)=north_temp_bndy;
            comb_south_temp_bndy(:,:,ind)=south_temp_bndy;
            comb_east_temp_bndy(:,:,ind)=east_temp_bndy;
            comb_west_temp_bndy(:,:,ind)=west_temp_bndy;
            comb_spatial_mean_north_temp(:,:,monthij)=comb_spatial_mean_north_temp(:,:,monthij)+north_temp/double(length(inputyear));  
            comb_spatial_mean_south_temp(:,:,monthij)=comb_spatial_mean_south_temp(:,:,monthij)+south_temp/double(length(inputyear));  
            comb_spatial_mean_east_temp(:,:,monthij)=comb_spatial_mean_east_temp(:,:,monthij)+east_temp/double(length(inputyear));  
            comb_spatial_mean_west_temp(:,:,monthij)=comb_spatial_mean_west_temp(:,:,monthij)+west_temp/double(length(inputyear));
            comb_spatial_mean_north_temp_inner(:,:,monthij)=comb_spatial_mean_north_temp_inner(:,:,monthij)+north_temp_inner/double(length(inputyear));  
            comb_spatial_mean_south_temp_inner(:,:,monthij)=comb_spatial_mean_south_temp_inner(:,:,monthij)+south_temp_inner/double(length(inputyear));  
            comb_spatial_mean_east_temp_inner(:,:,monthij)=comb_spatial_mean_east_temp_inner(:,:,monthij)+east_temp_inner/double(length(inputyear));  
            comb_spatial_mean_west_temp_inner(:,:,monthij)=comb_spatial_mean_west_temp_inner(:,:,monthij)+west_temp_inner/double(length(inputyear));
            comb_spatial_mean_north_temp_bndy(:,:,monthij)=comb_spatial_mean_north_temp_bndy(:,:,monthij)+north_temp_bndy/double(length(inputyear));  
            comb_spatial_mean_south_temp_bndy(:,:,monthij)=comb_spatial_mean_south_temp_bndy(:,:,monthij)+south_temp_bndy/double(length(inputyear));  
            comb_spatial_mean_east_temp_bndy(:,:,monthij)=comb_spatial_mean_east_temp_bndy(:,:,monthij)+east_temp_bndy/double(length(inputyear));  
            comb_spatial_mean_west_temp_bndy(:,:,monthij)=comb_spatial_mean_west_temp_bndy(:,:,monthij)+west_temp_bndy/double(length(inputyear));
% % %             put salt
            comb_north_salt(:,:,ind)=north_salt;
            comb_south_salt(:,:,ind)=south_salt;
            comb_east_salt(:,:,ind)=east_salt;
            comb_west_salt(:,:,ind)=west_salt;
            comb_north_salt_inner(:,:,ind)=north_salt_inner;
            comb_south_salt_inner(:,:,ind)=south_salt_inner;
            comb_east_salt_inner(:,:,ind)=east_salt_inner;
            comb_west_salt_inner(:,:,ind)=west_salt_inner;
            comb_north_salt_bndy(:,:,ind)=north_salt_bndy;
            comb_south_salt_bndy(:,:,ind)=south_salt_bndy;
            comb_east_salt_bndy(:,:,ind)=east_salt_bndy;
            comb_west_salt_bndy(:,:,ind)=west_salt_bndy;
            comb_spatial_mean_north_salt(:,:,monthij)=comb_spatial_mean_north_salt(:,:,monthij)+north_salt/double(length(inputyear));  
            comb_spatial_mean_south_salt(:,:,monthij)=comb_spatial_mean_south_salt(:,:,monthij)+south_salt/double(length(inputyear));  
            comb_spatial_mean_east_salt(:,:,monthij)=comb_spatial_mean_east_salt(:,:,monthij)+east_salt/double(length(inputyear));  
            comb_spatial_mean_west_salt(:,:,monthij)=comb_spatial_mean_west_salt(:,:,monthij)+west_salt/double(length(inputyear));
            comb_spatial_mean_north_salt_inner(:,:,monthij)=comb_spatial_mean_north_salt_inner(:,:,monthij)+north_salt_inner/double(length(inputyear));  
            comb_spatial_mean_south_salt_inner(:,:,monthij)=comb_spatial_mean_south_salt_inner(:,:,monthij)+south_salt_inner/double(length(inputyear));  
            comb_spatial_mean_east_salt_inner(:,:,monthij)=comb_spatial_mean_east_salt_inner(:,:,monthij)+east_salt_inner/double(length(inputyear));  
            comb_spatial_mean_west_salt_inner(:,:,monthij)=comb_spatial_mean_west_salt_inner(:,:,monthij)+west_salt_inner/double(length(inputyear));
            comb_spatial_mean_north_salt_bndy(:,:,monthij)=comb_spatial_mean_north_salt_bndy(:,:,monthij)+north_salt_bndy/double(length(inputyear));  
            comb_spatial_mean_south_salt_bndy(:,:,monthij)=comb_spatial_mean_south_salt_bndy(:,:,monthij)+south_salt_bndy/double(length(inputyear));  
            comb_spatial_mean_east_salt_bndy(:,:,monthij)=comb_spatial_mean_east_salt_bndy(:,:,monthij)+east_salt_bndy/double(length(inputyear));  
            comb_spatial_mean_west_salt_bndy(:,:,monthij)=comb_spatial_mean_west_salt_bndy(:,:,monthij)+west_salt_bndy/double(length(inputyear));    
% % %             put u           
            comb_north_u(:,:,ind)=north_u;
            comb_south_u(:,:,ind)=south_u;
            comb_east_u(:,:,ind)=east_u;
            comb_west_u(:,:,ind)=west_u;
            comb_north_u_inner(:,:,ind)=north_u_inner;
            comb_south_u_inner(:,:,ind)=south_u_inner;
            comb_east_u_inner(:,:,ind)=east_u_inner;
            comb_west_u_inner(:,:,ind)=west_u_inner;
            for ki=1:n_z
                comb_north_u_bndy(:,ki,ind)=north_u_bndy(:,ki)+north_ubar_bndy;
                comb_south_u_bndy(:,ki,ind)=south_u_bndy(:,ki)+south_ubar_bndy;
                comb_east_u_bndy(:,ki,ind)=east_u_bndy(:,ki)+east_ubar_bndy;
                comb_west_u_bndy(:,ki,ind)=west_u_bndy(:,ki)+west_ubar_bndy;
            end
            comb_spatial_mean_north_u(:,:,monthij)=comb_spatial_mean_north_u(:,:,monthij)+north_u/double(length(inputyear));  
            comb_spatial_mean_south_u(:,:,monthij)=comb_spatial_mean_south_u(:,:,monthij)+south_u/double(length(inputyear));  
            comb_spatial_mean_east_u(:,:,monthij)=comb_spatial_mean_east_u(:,:,monthij)+east_u/double(length(inputyear));  
            comb_spatial_mean_west_u(:,:,monthij)=comb_spatial_mean_west_u(:,:,monthij)+west_u/double(length(inputyear));
            comb_spatial_mean_north_u_inner(:,:,monthij)=comb_spatial_mean_north_u_inner(:,:,monthij)+north_u_inner/double(length(inputyear));  
            comb_spatial_mean_south_u_inner(:,:,monthij)=comb_spatial_mean_south_u_inner(:,:,monthij)+south_u_inner/double(length(inputyear));  
            comb_spatial_mean_east_u_inner(:,:,monthij)=comb_spatial_mean_east_u_inner(:,:,monthij)+east_u_inner/double(length(inputyear));  
            comb_spatial_mean_west_u_inner(:,:,monthij)=comb_spatial_mean_west_u_inner(:,:,monthij)+west_u_inner/double(length(inputyear));
            comb_spatial_mean_north_u_bndy(:,:,monthij)=comb_spatial_mean_north_u_bndy(:,:,monthij)+comb_north_u_bndy(:,:,ind)/double(length(inputyear));  
            comb_spatial_mean_south_u_bndy(:,:,monthij)=comb_spatial_mean_south_u_bndy(:,:,monthij)+comb_south_u_bndy(:,:,ind)/double(length(inputyear));  
            comb_spatial_mean_east_u_bndy(:,:,monthij)=comb_spatial_mean_east_u_bndy(:,:,monthij)+comb_east_u_bndy(:,:,ind)/double(length(inputyear));  
            comb_spatial_mean_west_u_bndy(:,:,monthij)=comb_spatial_mean_west_u_bndy(:,:,monthij)+comb_west_u_bndy(:,:,ind)/double(length(inputyear));               

% % %             put v           
            comb_north_v(:,:,ind)=north_v;
            comb_south_v(:,:,ind)=south_v;
            comb_east_v(:,:,ind)=east_v;
            comb_west_v(:,:,ind)=west_v;
            comb_north_v_inner(:,:,ind)=north_v_inner;
            comb_south_v_inner(:,:,ind)=south_v_inner;
            comb_east_v_inner(:,:,ind)=east_v_inner;
            comb_west_v_inner(:,:,ind)=west_v_inner;
            for ki=1:n_z
                comb_north_v_bndy(:,ki,ind)=north_v_bndy(:,ki)+north_vbar_bndy;
                comb_south_v_bndy(:,ki,ind)=south_v_bndy(:,ki)+south_vbar_bndy;
                comb_east_v_bndy(:,ki,ind)=east_v_bndy(:,ki)+east_vbar_bndy;
                comb_west_v_bndy(:,ki,ind)=west_v_bndy(:,ki)+west_vbar_bndy;
            end
            comb_spatial_mean_north_v(:,:,monthij)=comb_spatial_mean_north_v(:,:,monthij)+north_v/double(length(inputyear));  
            comb_spatial_mean_south_v(:,:,monthij)=comb_spatial_mean_south_v(:,:,monthij)+south_v/double(length(inputyear));  
            comb_spatial_mean_east_v(:,:,monthij)=comb_spatial_mean_east_v(:,:,monthij)+east_v/double(length(inputyear));  
            comb_spatial_mean_west_v(:,:,monthij)=comb_spatial_mean_west_v(:,:,monthij)+west_v/double(length(inputyear));
            comb_spatial_mean_north_v_inner(:,:,monthij)=comb_spatial_mean_north_v_inner(:,:,monthij)+north_v_inner/double(length(inputyear));  
            comb_spatial_mean_south_v_inner(:,:,monthij)=comb_spatial_mean_south_v_inner(:,:,monthij)+south_v_inner/double(length(inputyear));  
            comb_spatial_mean_east_v_inner(:,:,monthij)=comb_spatial_mean_east_v_inner(:,:,monthij)+east_v_inner/double(length(inputyear));  
            comb_spatial_mean_west_v_inner(:,:,monthij)=comb_spatial_mean_west_v_inner(:,:,monthij)+west_v_inner/double(length(inputyear));
            comb_spatial_mean_north_v_bndy(:,:,monthij)=comb_spatial_mean_north_v_bndy(:,:,monthij)+comb_north_v_bndy(:,:,ind)/double(length(inputyear));  
            comb_spatial_mean_south_v_bndy(:,:,monthij)=comb_spatial_mean_south_v_bndy(:,:,monthij)+comb_south_v_bndy(:,:,ind)/double(length(inputyear));  
            comb_spatial_mean_east_v_bndy(:,:,monthij)=comb_spatial_mean_east_v_bndy(:,:,monthij)+comb_east_v_bndy(:,:,ind)/double(length(inputyear));  
            comb_spatial_mean_west_v_bndy(:,:,monthij)=comb_spatial_mean_west_v_bndy(:,:,monthij)+comb_west_v_bndy(:,:,ind)/double(length(inputyear));  
            
            hold off
            close all;
            ind = ind + 1;
            toc;
        end
    end

    save([filedir,regionname,'_bndy_var',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

    
    
    
    
    ncid = netcdf.create(strcat(filedir, testname, '_',regionname,'_bndy_var_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');
%     netcdf.close(ncid);

    onedimid = netcdf.defDim(ncid,'one', 1);
    xi_rhodimid = netcdf.defDim(ncid, 'xi_rho', len_lon_rho_model);
    eta_rhodimid = netcdf.defDim(ncid,'eta_rho',len_lat_rho_model);
    xi_udimid = netcdf.defDim(ncid, 'xi_u', len_lon_u_model);
    eta_udimid = netcdf.defDim(ncid,'eta_u',len_lat_u_model);
    xi_vdimid = netcdf.defDim(ncid, 'xi_v', len_lon_v_model);
    eta_vdimid = netcdf.defDim(ncid,'eta_v',len_lat_v_model);
    s_rhodimid = netcdf.defDim(ncid, 's_rho',n_z);
    time_dimid = netcdf.defDim(ncid, 'time', 0);
    clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ['NWP 1/20 _ ', testname, 'all time monthly ','bndy', ' file']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', [' all time monthly ', 'bndy', ' (' , num2str(inputyear(1)),'-',num2str(inputyear(end)),')']);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', [' ROMS NWP 1/20 bndy data from _ ',testname ]);
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by Y.Y.Kim');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
    netcdf.putAtt(ncid,timevarid,'calendar','gregorian');

    clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
    netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
    netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
    netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');
    
    s_rhovarid=netcdf.defVar(ncid, 's_rho', 'NC_DOUBLE', [s_rhodimid]);
    netcdf.putAtt(ncid,s_rhovarid,'long_name','s_rho');
    
    xi_rhovarid=netcdf.defVar(ncid, 'xi_rho', 'NC_DOUBLE', [xi_rhodimid]);
    netcdf.putAtt(ncid,xi_rhovarid,'long_name','xi_rho');
    
    eta_rhovarid=netcdf.defVar(ncid, 'eta_rho', 'NC_DOUBLE', [eta_rhodimid]);
    netcdf.putAtt(ncid,eta_rhovarid,'long_name','eta_rho');
    
    lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xi_rhodimid eta_rhodimid]);
    netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
    netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

    lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xi_rhodimid eta_rhodimid]);
    netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
    netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
    
    xi_uvarid=netcdf.defVar(ncid, 'xi_u', 'NC_DOUBLE', [xi_udimid]);
    netcdf.putAtt(ncid,xi_uvarid,'long_name','xi_u');
    
    eta_uvarid=netcdf.defVar(ncid, 'eta_u', 'NC_DOUBLE', [eta_udimid]);
    netcdf.putAtt(ncid,eta_uvarid,'long_name','eta_u');
    
    lon_uvarid=netcdf.defVar(ncid, 'lon_u', 'NC_DOUBLE', [xi_udimid eta_udimid]);
    netcdf.putAtt(ncid,lon_uvarid,'long_name','lon_model');
    netcdf.putAtt(ncid,lon_uvarid,'units','degree_east');

    lat_uvarid=netcdf.defVar(ncid, 'lat_u', 'NC_DOUBLE', [xi_udimid eta_udimid]);
    netcdf.putAtt(ncid,lat_uvarid,'long_name','lat_model');
    netcdf.putAtt(ncid,lat_uvarid,'units','degree_north');
    
    xi_vvarid=netcdf.defVar(ncid, 'xi_v', 'NC_DOUBLE', [xi_vdimid]);
    netcdf.putAtt(ncid,xi_vvarid,'long_name','xi_v');
    
    eta_vvarid=netcdf.defVar(ncid, 'eta_v', 'NC_DOUBLE', [eta_vdimid]);
    netcdf.putAtt(ncid,eta_vvarid,'long_name','eta_v');
    
    lon_vvarid=netcdf.defVar(ncid, 'lon_v', 'NC_DOUBLE', [xi_vdimid eta_vdimid]);
    netcdf.putAtt(ncid,lon_vvarid,'long_name','lon_model');
    netcdf.putAtt(ncid,lon_vvarid,'units','degree_east');

    lat_vvarid=netcdf.defVar(ncid, 'lat_v', 'NC_DOUBLE', [xi_vdimid eta_vdimid]);
    netcdf.putAtt(ncid,lat_vvarid,'long_name','lat_model');
    netcdf.putAtt(ncid,lat_vvarid,'units','degree_north');

%     surf_datavarid=netcdf.defVar(ncid, 'surf_data', 'NC_FLOAT', [xi_rhodimid eta_rhodimid time_dimid]);
%     netcdf.putAtt(ncid,surf_datavarid,'long_name','surf_data');
%     netcdf.putAtt(ncid,surf_datavarid,'units', data_units);
%     
%     bot_datavarid=netcdf.defVar(ncid, 'bot_data', 'NC_FLOAT', [xi_rhodimid eta_rhodimid time_dimid]);
%     netcdf.putAtt(ncid,bot_datavarid,'long_name','bot_data');
%     netcdf.putAtt(ncid,bot_datavarid,'units', data_units);
    
%     if (ndim==3)

% % % def zeta
    north_zetavarid=netcdf.defVar(ncid, 'north_zeta', 'NC_FLOAT', [xi_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_zetavarid,'long_name','north_zeta');
    netcdf.putAtt(ncid,north_zetavarid,'units', 'm');

    south_zetavarid=netcdf.defVar(ncid, 'south_zeta', 'NC_FLOAT', [xi_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_zetavarid,'long_name','south_zeta');
    netcdf.putAtt(ncid,south_zetavarid,'units', 'm');

    east_zetavarid=netcdf.defVar(ncid, 'east_zeta', 'NC_FLOAT', [eta_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_zetavarid,'long_name','east_zeta');
    netcdf.putAtt(ncid,east_zetavarid,'units', 'm');

    west_zetavarid=netcdf.defVar(ncid, 'west_zeta', 'NC_FLOAT', [eta_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_zetavarid,'long_name','west_zeta');
    netcdf.putAtt(ncid,west_zetavarid,'units', 'm');

    north_zeta_innervarid=netcdf.defVar(ncid, 'north_zeta_inner', 'NC_FLOAT', [xi_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_zeta_innervarid,'long_name','north_zeta_inner');
    netcdf.putAtt(ncid,north_zeta_innervarid,'units', 'm');

    south_zeta_innervarid=netcdf.defVar(ncid, 'south_zeta_inner', 'NC_FLOAT', [xi_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_zeta_innervarid,'long_name','south_zeta_inner');
    netcdf.putAtt(ncid,south_zeta_innervarid,'units', 'm');

    east_zeta_innervarid=netcdf.defVar(ncid, 'east_zeta_inner', 'NC_FLOAT', [eta_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_zeta_innervarid,'long_name','east_zeta_inner');
    netcdf.putAtt(ncid,east_zeta_innervarid,'units', 'm');

    west_zeta_innervarid=netcdf.defVar(ncid, 'west_zeta_inner', 'NC_FLOAT', [eta_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_zeta_innervarid,'long_name','west_zeta_inner');
    netcdf.putAtt(ncid,west_zeta_innervarid,'units', 'm');

    north_zeta_bndyvarid=netcdf.defVar(ncid, 'north_zeta_bndy', 'NC_FLOAT', [xi_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_zeta_bndyvarid,'long_name','north_zeta_bndy');
    netcdf.putAtt(ncid,north_zeta_bndyvarid,'units', 'm');

    south_zeta_bndyvarid=netcdf.defVar(ncid, 'south_zeta_bndy', 'NC_FLOAT', [xi_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_zeta_bndyvarid,'long_name','south_zeta_bndy');
    netcdf.putAtt(ncid,south_zeta_bndyvarid,'units', 'm');

    east_zeta_bndyvarid=netcdf.defVar(ncid, 'east_zeta_bndy', 'NC_FLOAT', [eta_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_zeta_bndyvarid,'long_name','east_zeta_bndy');
    netcdf.putAtt(ncid,east_zeta_bndyvarid,'units', 'm');

    west_zeta_bndyvarid=netcdf.defVar(ncid, 'west_zeta_bndy', 'NC_FLOAT', [eta_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_zeta_bndyvarid,'long_name','west_zeta_bndy');
    netcdf.putAtt(ncid,west_zeta_bndyvarid,'units', 'm');

    clim_north_zetavarid=netcdf.defVar(ncid, 'clim_north_zeta', 'NC_FLOAT', [xi_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_north_zetavarid,'long_name','clim_north_zeta');
    netcdf.putAtt(ncid,clim_north_zetavarid,'units', 'm');

    clim_south_zetavarid=netcdf.defVar(ncid, 'clim_south_zeta', 'NC_FLOAT', [xi_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_south_zetavarid,'long_name','clim_south_zeta');
    netcdf.putAtt(ncid,clim_south_zetavarid,'units', 'm');

    clim_east_zetavarid=netcdf.defVar(ncid, 'clim_east_zeta', 'NC_FLOAT', [eta_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_east_zetavarid,'long_name','clim_east_zeta');
    netcdf.putAtt(ncid,clim_east_zetavarid,'units', 'm');

    clim_west_zetavarid=netcdf.defVar(ncid, 'clim_west_zeta', 'NC_FLOAT', [eta_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_west_zetavarid,'long_name','clim_west_zeta');
    netcdf.putAtt(ncid,clim_west_zetavarid,'units', 'm');
    
% % % def ubar
    north_ubarvarid=netcdf.defVar(ncid, 'north_ubar', 'NC_FLOAT', [xi_udimid time_dimid]);
    netcdf.putAtt(ncid,north_ubarvarid,'long_name','north_ubar');
    netcdf.putAtt(ncid,north_ubarvarid,'units', 'm/s');

    south_ubarvarid=netcdf.defVar(ncid, 'south_ubar', 'NC_FLOAT', [xi_udimid time_dimid]);
    netcdf.putAtt(ncid,south_ubarvarid,'long_name','south_ubar');
    netcdf.putAtt(ncid,south_ubarvarid,'units', 'm/s');

    east_ubarvarid=netcdf.defVar(ncid, 'east_ubar', 'NC_FLOAT', [eta_udimid time_dimid]);
    netcdf.putAtt(ncid,east_ubarvarid,'long_name','east_ubar');
    netcdf.putAtt(ncid,east_ubarvarid,'units', 'm/s');

    west_ubarvarid=netcdf.defVar(ncid, 'west_ubar', 'NC_FLOAT', [eta_udimid time_dimid]);
    netcdf.putAtt(ncid,west_ubarvarid,'long_name','west_ubar');
    netcdf.putAtt(ncid,west_ubarvarid,'units', 'm/s');

    north_ubar_innervarid=netcdf.defVar(ncid, 'north_ubar_inner', 'NC_FLOAT', [xi_udimid time_dimid]);
    netcdf.putAtt(ncid,north_ubar_innervarid,'long_name','north_ubar_inner');
    netcdf.putAtt(ncid,north_ubar_innervarid,'units', 'm/s');

    south_ubar_innervarid=netcdf.defVar(ncid, 'south_ubar_inner', 'NC_FLOAT', [xi_udimid time_dimid]);
    netcdf.putAtt(ncid,south_ubar_innervarid,'long_name','south_ubar_inner');
    netcdf.putAtt(ncid,south_ubar_innervarid,'units', 'm/s');

    east_ubar_innervarid=netcdf.defVar(ncid, 'east_ubar_inner', 'NC_FLOAT', [eta_udimid time_dimid]);
    netcdf.putAtt(ncid,east_ubar_innervarid,'long_name','east_ubar_inner');
    netcdf.putAtt(ncid,east_ubar_innervarid,'units', 'm/s');

    west_ubar_innervarid=netcdf.defVar(ncid, 'west_ubar_inner', 'NC_FLOAT', [eta_udimid time_dimid]);
    netcdf.putAtt(ncid,west_ubar_innervarid,'long_name','west_ubar_inner');
    netcdf.putAtt(ncid,west_ubar_innervarid,'units', 'm/s');

    north_ubar_bndyvarid=netcdf.defVar(ncid, 'north_ubar_bndy', 'NC_FLOAT', [xi_udimid time_dimid]);
    netcdf.putAtt(ncid,north_ubar_bndyvarid,'long_name','north_ubar_bndy');
    netcdf.putAtt(ncid,north_ubar_bndyvarid,'units', 'm/s');

    south_ubar_bndyvarid=netcdf.defVar(ncid, 'south_ubar_bndy', 'NC_FLOAT', [xi_udimid time_dimid]);
    netcdf.putAtt(ncid,south_ubar_bndyvarid,'long_name','south_ubar_bndy');
    netcdf.putAtt(ncid,south_ubar_bndyvarid,'units', 'm/s');

    east_ubar_bndyvarid=netcdf.defVar(ncid, 'east_ubar_bndy', 'NC_FLOAT', [eta_udimid time_dimid]);
    netcdf.putAtt(ncid,east_ubar_bndyvarid,'long_name','east_ubar_bndy');
    netcdf.putAtt(ncid,east_ubar_bndyvarid,'units', 'm/s');

    west_ubar_bndyvarid=netcdf.defVar(ncid, 'west_ubar_bndy', 'NC_FLOAT', [eta_udimid time_dimid]);
    netcdf.putAtt(ncid,west_ubar_bndyvarid,'long_name','west_ubar_bndy');
    netcdf.putAtt(ncid,west_ubar_bndyvarid,'units', 'm/s');

    clim_north_ubarvarid=netcdf.defVar(ncid, 'clim_north_ubar', 'NC_FLOAT', [xi_udimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_north_ubarvarid,'long_name','clim_north_ubar');
    netcdf.putAtt(ncid,clim_north_ubarvarid,'units', 'm/s');

    clim_south_ubarvarid=netcdf.defVar(ncid, 'clim_south_ubar', 'NC_FLOAT', [xi_udimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_south_ubarvarid,'long_name','clim_south_ubar');
    netcdf.putAtt(ncid,clim_south_ubarvarid,'units', 'm/s');

    clim_east_ubarvarid=netcdf.defVar(ncid, 'clim_east_ubar', 'NC_FLOAT', [eta_udimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_east_ubarvarid,'long_name','clim_east_ubar');
    netcdf.putAtt(ncid,clim_east_ubarvarid,'units', 'm/s');

    clim_west_ubarvarid=netcdf.defVar(ncid, 'clim_west_ubar', 'NC_FLOAT', [eta_udimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_west_ubarvarid,'long_name','clim_west_ubar');
    netcdf.putAtt(ncid,clim_west_ubarvarid,'units', 'm/s');
    
% % % def vbar
    north_vbarvarid=netcdf.defVar(ncid, 'north_vbar', 'NC_FLOAT', [xi_vdimid time_dimid]);
    netcdf.putAtt(ncid,north_vbarvarid,'long_name','north_vbar');
    netcdf.putAtt(ncid,north_vbarvarid,'units', 'm/s');

    south_vbarvarid=netcdf.defVar(ncid, 'south_vbar', 'NC_FLOAT', [xi_vdimid time_dimid]);
    netcdf.putAtt(ncid,south_vbarvarid,'long_name','south_vbar');
    netcdf.putAtt(ncid,south_vbarvarid,'units', 'm/s');

    east_vbarvarid=netcdf.defVar(ncid, 'east_vbar', 'NC_FLOAT', [eta_vdimid time_dimid]);
    netcdf.putAtt(ncid,east_vbarvarid,'long_name','east_vbar');
    netcdf.putAtt(ncid,east_vbarvarid,'units', 'm/s');

    west_vbarvarid=netcdf.defVar(ncid, 'west_vbar', 'NC_FLOAT', [eta_vdimid time_dimid]);
    netcdf.putAtt(ncid,west_vbarvarid,'long_name','west_vbar');
    netcdf.putAtt(ncid,west_vbarvarid,'units', 'm/s');

    north_vbar_innervarid=netcdf.defVar(ncid, 'north_vbar_inner', 'NC_FLOAT', [xi_vdimid time_dimid]);
    netcdf.putAtt(ncid,north_vbar_innervarid,'long_name','north_vbar_inner');
    netcdf.putAtt(ncid,north_vbar_innervarid,'units', 'm/s');

    south_vbar_innervarid=netcdf.defVar(ncid, 'south_vbar_inner', 'NC_FLOAT', [xi_vdimid time_dimid]);
    netcdf.putAtt(ncid,south_vbar_innervarid,'long_name','south_vbar_inner');
    netcdf.putAtt(ncid,south_vbar_innervarid,'units', 'm/s');

    east_vbar_innervarid=netcdf.defVar(ncid, 'east_vbar_inner', 'NC_FLOAT', [eta_vdimid time_dimid]);
    netcdf.putAtt(ncid,east_vbar_innervarid,'long_name','east_vbar_inner');
    netcdf.putAtt(ncid,east_vbar_innervarid,'units', 'm/s');

    west_vbar_innervarid=netcdf.defVar(ncid, 'west_vbar_inner', 'NC_FLOAT', [eta_vdimid time_dimid]);
    netcdf.putAtt(ncid,west_vbar_innervarid,'long_name','west_vbar_inner');
    netcdf.putAtt(ncid,west_vbar_innervarid,'units', 'm/s');

    north_vbar_bndyvarid=netcdf.defVar(ncid, 'north_vbar_bndy', 'NC_FLOAT', [xi_vdimid time_dimid]);
    netcdf.putAtt(ncid,north_vbar_bndyvarid,'long_name','north_vbar_bndy');
    netcdf.putAtt(ncid,north_vbar_bndyvarid,'units', 'm/s');

    south_vbar_bndyvarid=netcdf.defVar(ncid, 'south_vbar_bndy', 'NC_FLOAT', [xi_vdimid time_dimid]);
    netcdf.putAtt(ncid,south_vbar_bndyvarid,'long_name','south_vbar_bndy');
    netcdf.putAtt(ncid,south_vbar_bndyvarid,'units', 'm/s');

    east_vbar_bndyvarid=netcdf.defVar(ncid, 'east_vbar_bndy', 'NC_FLOAT', [eta_vdimid time_dimid]);
    netcdf.putAtt(ncid,east_vbar_bndyvarid,'long_name','east_vbar_bndy');
    netcdf.putAtt(ncid,east_vbar_bndyvarid,'units', 'm/s');

    west_vbar_bndyvarid=netcdf.defVar(ncid, 'west_vbar_bndy', 'NC_FLOAT', [eta_vdimid time_dimid]);
    netcdf.putAtt(ncid,west_vbar_bndyvarid,'long_name','west_vbar_bndy');
    netcdf.putAtt(ncid,west_vbar_bndyvarid,'units', 'm/s');

    clim_north_vbarvarid=netcdf.defVar(ncid, 'clim_north_vbar', 'NC_FLOAT', [xi_vdimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_north_vbarvarid,'long_name','clim_north_vbar');
    netcdf.putAtt(ncid,clim_north_vbarvarid,'units', 'm/s');

    clim_south_vbarvarid=netcdf.defVar(ncid, 'clim_south_vbar', 'NC_FLOAT', [xi_vdimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_south_vbarvarid,'long_name','clim_south_vbar');
    netcdf.putAtt(ncid,clim_south_vbarvarid,'units', 'm/s');

    clim_east_vbarvarid=netcdf.defVar(ncid, 'clim_east_vbar', 'NC_FLOAT', [eta_vdimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_east_vbarvarid,'long_name','clim_east_vbar');
    netcdf.putAtt(ncid,clim_east_vbarvarid,'units', 'm/s');

    clim_west_vbarvarid=netcdf.defVar(ncid, 'clim_west_vbar', 'NC_FLOAT', [eta_vdimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_west_vbarvarid,'long_name','clim_west_vbar');
    netcdf.putAtt(ncid,clim_west_vbarvarid,'units', 'm/s');
    
% % % def temp
    north_tempvarid=netcdf.defVar(ncid, 'north_temp', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_tempvarid,'long_name','north_temp');
    netcdf.putAtt(ncid,north_tempvarid,'units','celsius degree');

    south_tempvarid=netcdf.defVar(ncid, 'south_temp', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_tempvarid,'long_name','south_temp');
    netcdf.putAtt(ncid,south_tempvarid,'units','celsius degree');

    east_tempvarid=netcdf.defVar(ncid, 'east_temp', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_tempvarid,'long_name','east_temp');
    netcdf.putAtt(ncid,east_tempvarid,'units','celsius degree');

    west_tempvarid=netcdf.defVar(ncid, 'west_temp', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_tempvarid,'long_name','west_temp');
    netcdf.putAtt(ncid,west_tempvarid,'units','celsius degree');
    
    north_temp_innervarid=netcdf.defVar(ncid, 'north_temp_inner', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_temp_innervarid,'long_name','north_temp_inner');
    netcdf.putAtt(ncid,north_temp_innervarid,'units','celsius degree');

    south_temp_innervarid=netcdf.defVar(ncid, 'south_temp_inner', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_temp_innervarid,'long_name','south_temp_inner');
    netcdf.putAtt(ncid,south_temp_innervarid,'units','celsius degree');

    east_temp_innervarid=netcdf.defVar(ncid, 'east_temp_inner', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_temp_innervarid,'long_name','east_temp_inner');
    netcdf.putAtt(ncid,east_temp_innervarid,'units','celsius degree');

    west_temp_innervarid=netcdf.defVar(ncid, 'west_temp_inner', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_temp_innervarid,'long_name','west_temp_inner');
    netcdf.putAtt(ncid,west_temp_innervarid,'units','celsius degree');

    north_temp_bndyvarid=netcdf.defVar(ncid, 'north_temp_bndy', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_temp_bndyvarid,'long_name','north_temp_bndy');
    netcdf.putAtt(ncid,north_temp_bndyvarid,'units','celsius degree');

    south_temp_bndyvarid=netcdf.defVar(ncid, 'south_temp_bndy', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_temp_bndyvarid,'long_name','south_temp_bndy');
    netcdf.putAtt(ncid,south_temp_bndyvarid,'units','celsius degree');

    east_temp_bndyvarid=netcdf.defVar(ncid, 'east_temp_bndy', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_temp_bndyvarid,'long_name','east_temp_bndy');
    netcdf.putAtt(ncid,east_temp_bndyvarid,'units','celsius degree');

    west_temp_bndyvarid=netcdf.defVar(ncid, 'west_temp_bndy', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_temp_bndyvarid,'long_name','west_temp_bndy');
    netcdf.putAtt(ncid,west_temp_bndyvarid,'units','celsius degree');

    clim_north_tempvarid=netcdf.defVar(ncid, 'clim_north_temp', 'NC_FLOAT', [xi_rhodimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_north_tempvarid,'long_name','clim_north_temp');
    netcdf.putAtt(ncid,clim_north_tempvarid,'units','celsius degree');

    clim_south_tempvarid=netcdf.defVar(ncid, 'clim_south_temp', 'NC_FLOAT', [xi_rhodimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_south_tempvarid,'long_name','clim_south');
    netcdf.putAtt(ncid,clim_south_tempvarid,'units','celsius degree');

    clim_east_tempvarid=netcdf.defVar(ncid, 'clim_east_temp', 'NC_FLOAT', [eta_rhodimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_east_tempvarid,'long_name','clim_east');
    netcdf.putAtt(ncid,clim_east_tempvarid,'units','celsius degree');

    clim_west_tempvarid=netcdf.defVar(ncid, 'clim_west_temp', 'NC_FLOAT', [eta_rhodimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_west_tempvarid,'long_name','clim_west');
    netcdf.putAtt(ncid,clim_west_tempvarid,'units','celsius degree');
    
% % % def salt
    north_saltvarid=netcdf.defVar(ncid, 'north_salt', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_saltvarid,'long_name','north_salt');
    netcdf.putAtt(ncid,north_saltvarid,'units',' ');

    south_saltvarid=netcdf.defVar(ncid, 'south_salt', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_saltvarid,'long_name','south_salt');
    netcdf.putAtt(ncid,south_saltvarid,'units',' ');

    east_saltvarid=netcdf.defVar(ncid, 'east_salt', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_saltvarid,'long_name','east_salt');
    netcdf.putAtt(ncid,east_saltvarid,'units',' ');

    west_saltvarid=netcdf.defVar(ncid, 'west_salt', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_saltvarid,'long_name','west_salt');
    netcdf.putAtt(ncid,west_saltvarid,'units',' ');
    
    north_salt_innervarid=netcdf.defVar(ncid, 'north_salt_inner', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_salt_innervarid,'long_name','north_salt_inner');
    netcdf.putAtt(ncid,north_salt_innervarid,'units',' ');

    south_salt_innervarid=netcdf.defVar(ncid, 'south_salt_inner', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_salt_innervarid,'long_name','south_salt_inner');
    netcdf.putAtt(ncid,south_salt_innervarid,'units',' ');

    east_salt_innervarid=netcdf.defVar(ncid, 'east_salt_inner', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_salt_innervarid,'long_name','east_salt_inner');
    netcdf.putAtt(ncid,east_salt_innervarid,'units',' ');

    west_salt_innervarid=netcdf.defVar(ncid, 'west_salt_inner', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_salt_innervarid,'long_name','west_salt_inner');
    netcdf.putAtt(ncid,west_salt_innervarid,'units',' ');
    
    north_salt_inner_bndyvarid=netcdf.defVar(ncid, 'north_salt_inner_bndy', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_salt_inner_bndyvarid,'long_name','north_salt_inner_bndy');
    netcdf.putAtt(ncid,north_salt_inner_bndyvarid,'units',' ');

    south_salt_bndyvarid=netcdf.defVar(ncid, 'south_salt_bndy', 'NC_FLOAT', [xi_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_salt_bndyvarid,'long_name','south_salt_bndy');
    netcdf.putAtt(ncid,south_salt_bndyvarid,'units',' ');

    east_salt_bndyvarid=netcdf.defVar(ncid, 'east_salt_bndy', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_salt_bndyvarid,'long_name','east_salt_bndy');
    netcdf.putAtt(ncid,east_salt_bndyvarid,'units',' ');

    west_salt_bndyvarid=netcdf.defVar(ncid, 'west_salt_bndy', 'NC_FLOAT', [eta_rhodimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_salt_bndyvarid,'long_name','west_salt_bndy');
    netcdf.putAtt(ncid,west_salt_bndyvarid,'units',' ');

    clim_north_saltvarid=netcdf.defVar(ncid, 'clim_north_salt', 'NC_FLOAT', [xi_rhodimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_north_saltvarid,'long_name','clim_north_salt');
    netcdf.putAtt(ncid,clim_north_saltvarid,'units',' ');

    clim_south_saltvarid=netcdf.defVar(ncid, 'clim_south_salt', 'NC_FLOAT', [xi_rhodimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_south_saltvarid,'long_name','clim_south');
    netcdf.putAtt(ncid,clim_south_saltvarid,'units',' ');

    clim_east_saltvarid=netcdf.defVar(ncid, 'clim_east_salt', 'NC_FLOAT', [eta_rhodimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_east_saltvarid,'long_name','clim_east');
    netcdf.putAtt(ncid,clim_east_saltvarid,'units',' ');

    clim_west_saltvarid=netcdf.defVar(ncid, 'clim_west_salt', 'NC_FLOAT', [eta_rhodimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_west_saltvarid,'long_name','clim_west');
    netcdf.putAtt(ncid,clim_west_saltvarid,'units',' '); 
    
% % % def u
    north_uvarid=netcdf.defVar(ncid, 'north_u', 'NC_FLOAT', [xi_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_uvarid,'long_name','north_u');
    netcdf.putAtt(ncid,north_uvarid,'units','m/s');

    south_uvarid=netcdf.defVar(ncid, 'south_u', 'NC_FLOAT', [xi_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_uvarid,'long_name','south_u');
    netcdf.putAtt(ncid,south_uvarid,'units','m/s');

    east_uvarid=netcdf.defVar(ncid, 'east_u', 'NC_FLOAT', [eta_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_uvarid,'long_name','east_u');
    netcdf.putAtt(ncid,east_uvarid,'units','m/s');

    west_uvarid=netcdf.defVar(ncid, 'west_u', 'NC_FLOAT', [eta_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_uvarid,'long_name','west_u');
    netcdf.putAtt(ncid,west_uvarid,'units','m/s');
    
    north_u_innervarid=netcdf.defVar(ncid, 'north_u_inner', 'NC_FLOAT', [xi_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_u_innervarid,'long_name','north_u_inner');
    netcdf.putAtt(ncid,north_u_innervarid,'units','m/s');

    south_u_innervarid=netcdf.defVar(ncid, 'south_u_inner', 'NC_FLOAT', [xi_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_u_innervarid,'long_name','south_u_inner');
    netcdf.putAtt(ncid,south_u_innervarid,'units','m/s');

    east_u_innervarid=netcdf.defVar(ncid, 'east_u_inner', 'NC_FLOAT', [eta_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_u_innervarid,'long_name','east_u_inner');
    netcdf.putAtt(ncid,east_u_innervarid,'units','m/s');

    west_u_innervarid=netcdf.defVar(ncid, 'west_u_inner', 'NC_FLOAT', [eta_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_u_innervarid,'long_name','west_u_inner');
    netcdf.putAtt(ncid,west_u_innervarid,'units','m/s');
    
    north_u_bndyvarid=netcdf.defVar(ncid, 'north_u_bndy', 'NC_FLOAT', [xi_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_u_bndyvarid,'long_name','north_u_bndy + ubar');
    netcdf.putAtt(ncid,north_u_bndyvarid,'units','m/s');

    south_u_bndyvarid=netcdf.defVar(ncid, 'south_u_bndy', 'NC_FLOAT', [xi_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_u_bndyvarid,'long_name','south_u_bndy + ubar');
    netcdf.putAtt(ncid,south_u_bndyvarid,'units','m/s');

    east_u_bndyvarid=netcdf.defVar(ncid, 'east_u_bndy', 'NC_FLOAT', [eta_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_u_bndyvarid,'long_name','east_u_bndy + ubar');
    netcdf.putAtt(ncid,east_u_bndyvarid,'units','m/s');

    west_u_bndyvarid=netcdf.defVar(ncid, 'west_u_bndy', 'NC_FLOAT', [eta_udimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_u_bndyvarid,'long_name','west_u_bndy + ubar');
    netcdf.putAtt(ncid,west_u_bndyvarid,'units','m/s');

    clim_north_uvarid=netcdf.defVar(ncid, 'clim_north_u', 'NC_FLOAT', [xi_udimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_north_uvarid,'long_name','clim_north_u');
    netcdf.putAtt(ncid,clim_north_uvarid,'units','m/s');

    clim_south_uvarid=netcdf.defVar(ncid, 'clim_south_u', 'NC_FLOAT', [xi_udimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_south_uvarid,'long_name','clim_south');
    netcdf.putAtt(ncid,clim_south_uvarid,'units','m/s');

    clim_east_uvarid=netcdf.defVar(ncid, 'clim_east_u', 'NC_FLOAT', [eta_udimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_east_uvarid,'long_name','clim_east');
    netcdf.putAtt(ncid,clim_east_uvarid,'units','m/s');

    clim_west_uvarid=netcdf.defVar(ncid, 'clim_west_u', 'NC_FLOAT', [eta_udimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_west_uvarid,'long_name','clim_west');
    netcdf.putAtt(ncid,clim_west_uvarid,'units','m/s');
    
% % % def v
    north_vvarid=netcdf.defVar(ncid, 'north_v', 'NC_FLOAT', [xi_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_vvarid,'long_name','north_v');
    netcdf.putAtt(ncid,north_vvarid,'units','m/s');

    south_vvarid=netcdf.defVar(ncid, 'south_v', 'NC_FLOAT', [xi_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_vvarid,'long_name','south_v');
    netcdf.putAtt(ncid,south_vvarid,'units','m/s');

    east_vvarid=netcdf.defVar(ncid, 'east_v', 'NC_FLOAT', [eta_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_vvarid,'long_name','east_v');
    netcdf.putAtt(ncid,east_vvarid,'units','m/s');

    west_vvarid=netcdf.defVar(ncid, 'west_v', 'NC_FLOAT', [eta_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_vvarid,'long_name','west_v');
    netcdf.putAtt(ncid,west_vvarid,'units','m/s');
    
    north_v_innervarid=netcdf.defVar(ncid, 'north_v_inner', 'NC_FLOAT', [xi_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_v_innervarid,'long_name','north_v_inner');
    netcdf.putAtt(ncid,north_v_innervarid,'units','m/s');

    south_v_innervarid=netcdf.defVar(ncid, 'south_v_inner', 'NC_FLOAT', [xi_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_v_innervarid,'long_name','south_v_inner');
    netcdf.putAtt(ncid,south_v_innervarid,'units','m/s');

    east_v_innervarid=netcdf.defVar(ncid, 'east_v_inner', 'NC_FLOAT', [eta_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_v_innervarid,'long_name','east_v_inner');
    netcdf.putAtt(ncid,east_v_innervarid,'units','m/s');

    west_v_innervarid=netcdf.defVar(ncid, 'west_v_inner', 'NC_FLOAT', [eta_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_v_innervarid,'long_name','west_v_inner');
    netcdf.putAtt(ncid,west_v_innervarid,'units','m/s');
    
    north_v_bndyvarid=netcdf.defVar(ncid, 'north_v_bndy', 'NC_FLOAT', [xi_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,north_v_bndyvarid,'long_name','north_v_bndy + vbar');
    netcdf.putAtt(ncid,north_v_bndyvarid,'units','m/s');

    south_v_bndyvarid=netcdf.defVar(ncid, 'south_v_bndy', 'NC_FLOAT', [xi_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,south_v_bndyvarid,'long_name','south_v_bndy + vbar');
    netcdf.putAtt(ncid,south_v_bndyvarid,'units','m/s');

    east_v_bndyvarid=netcdf.defVar(ncid, 'east_v_bndy', 'NC_FLOAT', [eta_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,east_v_bndyvarid,'long_name','east_v_bndy + vbar');
    netcdf.putAtt(ncid,east_v_bndyvarid,'units','m/s');

    west_v_bndyvarid=netcdf.defVar(ncid, 'west_v_bndy', 'NC_FLOAT', [eta_vdimid s_rhodimid time_dimid]);
    netcdf.putAtt(ncid,west_v_bndyvarid,'long_name','west_v_bndy + vbar');
    netcdf.putAtt(ncid,west_v_bndyvarid,'units','m/s');

    clim_north_vvarid=netcdf.defVar(ncid, 'clim_north_v', 'NC_FLOAT', [xi_vdimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_north_vvarid,'long_name','clim_north_v');
    netcdf.putAtt(ncid,clim_north_vvarid,'units','m/s');

    clim_south_vvarid=netcdf.defVar(ncid, 'clim_south_v', 'NC_FLOAT', [xi_vdimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_south_vvarid,'long_name','clim_south');
    netcdf.putAtt(ncid,clim_south_vvarid,'units','m/s');

    clim_east_vvarid=netcdf.defVar(ncid, 'clim_east_v', 'NC_FLOAT', [eta_vdimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_east_vvarid,'long_name','clim_east');
    netcdf.putAtt(ncid,clim_east_vvarid,'units','m/s');

    clim_west_vvarid=netcdf.defVar(ncid, 'clim_west_v', 'NC_FLOAT', [eta_vdimid s_rhodimid clim_time_dimid]);
    netcdf.putAtt(ncid,clim_west_vvarid,'long_name','clim_west');
    netcdf.putAtt(ncid,clim_west_vvarid,'units','m/s');  
    
%     end

%     clim_surfvarid=netcdf.defVar(ncid, 'clim_surf', 'NC_FLOAT', [xi_rhodimid eta_rhodimid clim_time_dimid]);
%     netcdf.putAtt(ncid,clim_surfvarid,'long_name','clim_surf');
%     netcdf.putAtt(ncid,clim_surfvarid,'units',data_units);
%     
%     clim_botvarid=netcdf.defVar(ncid, 'clim_bot', 'NC_FLOAT', [xi_rhodimid eta_rhodimid clim_time_dimid]);
%     netcdf.putAtt(ncid,clim_botvarid,'long_name','clim_bot');
%     netcdf.putAtt(ncid,clim_botvarid,'units',data_units);
    
    netcdf.endDef(ncid);
%     if (ndim==3)
%     elseif (ndim==4)
%     end
    tind=1;
    for yearij = 1:length(inputyear)
        for month=1:12 
            tempyear = inputyear(yearij);
            ftime(tind) = datenum(tempyear,month,15) - datenum(1900,12,31);
            tind=tind+1;
        end
    end
    for month=1:12 
            tempyear = inputyear(yearij);
            climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
    end

    netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
    netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
    netcdf.putVar(ncid, xi_rhovarid, 0, len_lon_rho_model, 1:len_lon_rho_model);
    netcdf.putVar(ncid, eta_rhovarid, 0, len_lat_rho_model, 1:len_lat_rho_model);
    netcdf.putVar(ncid, xi_uvarid, 0, len_lon_u_model, 1:len_lon_u_model);
    netcdf.putVar(ncid, eta_uvarid, 0, len_lat_u_model, 1:len_lat_u_model);
    netcdf.putVar(ncid, xi_vvarid, 0, len_lon_v_model, 1:len_lon_v_model);
    netcdf.putVar(ncid, eta_vvarid, 0, len_lat_v_model, 1:len_lat_v_model);
    
    netcdf.putVar(ncid, s_rhovarid, 0, n_z, s_rho);
    netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_rho_model len_lat_rho_model], lon_rho);
    netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_rho_model len_lat_rho_model], lat_rho);
    netcdf.putVar(ncid, lon_uvarid, [0 0], [len_lon_u_model len_lat_u_model], lon_u);
    netcdf.putVar(ncid, lat_uvarid, [0 0], [len_lon_u_model len_lat_u_model], lat_u);
    netcdf.putVar(ncid, lon_vvarid, [0 0], [len_lon_v_model len_lat_v_model], lon_v);
    netcdf.putVar(ncid, lat_vvarid, [0 0], [len_lon_v_model len_lat_v_model], lat_v);
%     netcdf.putVar(ncid, surf_datavarid, [0 0 0], [len_lon_rho_model len_lat_rho_model length(ftime)], comb_surf_data);
%     netcdf.putVar(ncid, bot_datavarid, [0 0 0], [len_lon_rho_model len_lat_rho_model length(ftime)], comb_bot_data);
    
%     if (ndim==3)
    netcdf.putVar(ncid, north_zetavarid, [0 0], [len_lon_rho_model length(ftime)], comb_north_zeta);
    netcdf.putVar(ncid, south_zetavarid, [0 0], [len_lon_rho_model length(ftime)], comb_south_zeta);
    netcdf.putVar(ncid, east_zetavarid, [0 0], [len_lat_rho_model length(ftime)], comb_east_zeta);
    netcdf.putVar(ncid, west_zetavarid, [0 0], [len_lat_rho_model length(ftime)], comb_west_zeta);
    netcdf.putVar(ncid, north_zeta_innervarid, [0 0], [len_lon_rho_model length(ftime)], comb_north_zeta_inner);
    netcdf.putVar(ncid, south_zeta_innervarid, [0 0], [len_lon_rho_model length(ftime)], comb_south_zeta_inner);
    netcdf.putVar(ncid, east_zeta_innervarid, [0 0], [len_lat_rho_model length(ftime)], comb_east_zeta_inner);
    netcdf.putVar(ncid, west_zeta_innervarid, [0 0], [len_lat_rho_model length(ftime)], comb_west_zeta_inner);
    netcdf.putVar(ncid, north_zeta_bndyvarid, [0 0], [len_lon_rho_model length(ftime)], comb_north_zeta_bndy);
    netcdf.putVar(ncid, south_zeta_bndyvarid, [0 0], [len_lon_rho_model length(ftime)], comb_south_zeta_bndy);
    netcdf.putVar(ncid, east_zeta_bndyvarid, [0 0], [len_lat_rho_model length(ftime)], comb_east_zeta_bndy);
    netcdf.putVar(ncid, west_zeta_bndyvarid, [0 0], [len_lat_rho_model length(ftime)], comb_west_zeta_bndy);
    netcdf.putVar(ncid, clim_north_zetavarid, [0 0], [len_lon_rho_model   length(climtime)], comb_spatial_mean_north_zeta);
    netcdf.putVar(ncid, clim_south_zetavarid, [0 0], [len_lon_rho_model   length(climtime)], comb_spatial_mean_south_zeta);
    netcdf.putVar(ncid, clim_east_zetavarid, [0 0], [len_lat_rho_model   length(climtime)], comb_spatial_mean_east_zeta);
    netcdf.putVar(ncid, clim_west_zetavarid, [0 0], [len_lat_rho_model   length(climtime)], comb_spatial_mean_west_zeta);
    
    netcdf.putVar(ncid, north_ubarvarid, [0 0], [len_lon_u_model length(ftime)], comb_north_ubar);
    netcdf.putVar(ncid, south_ubarvarid, [0 0], [len_lon_u_model length(ftime)], comb_south_ubar);
    netcdf.putVar(ncid, east_ubarvarid, [0 0], [len_lat_u_model length(ftime)], comb_east_ubar);
    netcdf.putVar(ncid, west_ubarvarid, [0 0], [len_lat_u_model length(ftime)], comb_west_ubar);
    netcdf.putVar(ncid, north_ubar_innervarid, [0 0], [len_lon_u_model length(ftime)], comb_north_ubar_inner);
    netcdf.putVar(ncid, south_ubar_innervarid, [0 0], [len_lon_u_model length(ftime)], comb_south_ubar_inner);
    netcdf.putVar(ncid, east_ubar_innervarid, [0 0], [len_lat_u_model length(ftime)], comb_east_ubar_inner);
    netcdf.putVar(ncid, west_ubar_innervarid, [0 0], [len_lat_u_model length(ftime)], comb_west_ubar_inner);
    netcdf.putVar(ncid, north_ubar_bndyvarid, [0 0], [len_lon_u_model length(ftime)], comb_north_ubar_bndy);
    netcdf.putVar(ncid, south_ubar_bndyvarid, [0 0], [len_lon_u_model length(ftime)], comb_south_ubar_bndy);
    netcdf.putVar(ncid, east_ubar_bndyvarid, [0 0], [len_lat_u_model length(ftime)], comb_east_ubar_bndy);
    netcdf.putVar(ncid, west_ubar_bndyvarid, [0 0], [len_lat_u_model length(ftime)], comb_west_ubar_bndy);
    netcdf.putVar(ncid, clim_north_ubarvarid, [0 0], [len_lon_u_model   length(climtime)], comb_spatial_mean_north_ubar);
    netcdf.putVar(ncid, clim_south_ubarvarid, [0 0], [len_lon_u_model   length(climtime)], comb_spatial_mean_south_ubar);
    netcdf.putVar(ncid, clim_east_ubarvarid, [0 0], [len_lat_u_model   length(climtime)], comb_spatial_mean_east_ubar);
    netcdf.putVar(ncid, clim_west_ubarvarid, [0 0], [len_lat_u_model   length(climtime)], comb_spatial_mean_west_ubar);
    
    netcdf.putVar(ncid, north_vbarvarid, [0 0], [len_lon_v_model length(ftime)], comb_north_vbar);
    netcdf.putVar(ncid, south_vbarvarid, [0 0], [len_lon_v_model length(ftime)], comb_south_vbar);
    netcdf.putVar(ncid, east_vbarvarid, [0 0], [len_lat_v_model length(ftime)], comb_east_vbar);
    netcdf.putVar(ncid, west_vbarvarid, [0 0], [len_lat_v_model length(ftime)], comb_west_vbar);
    netcdf.putVar(ncid, north_vbar_innervarid, [0 0], [len_lon_v_model length(ftime)], comb_north_vbar_inner);
    netcdf.putVar(ncid, south_vbar_innervarid, [0 0], [len_lon_v_model length(ftime)], comb_south_vbar_inner);
    netcdf.putVar(ncid, east_vbar_innervarid, [0 0], [len_lat_v_model length(ftime)], comb_east_vbar_inner);
    netcdf.putVar(ncid, west_vbar_innervarid, [0 0], [len_lat_v_model length(ftime)], comb_west_vbar_inner);
    netcdf.putVar(ncid, north_vbar_bndyvarid, [0 0], [len_lon_v_model length(ftime)], comb_north_vbar_bndy);
    netcdf.putVar(ncid, south_vbar_bndyvarid, [0 0], [len_lon_v_model length(ftime)], comb_south_vbar_bndy);
    netcdf.putVar(ncid, east_vbar_bndyvarid, [0 0], [len_lat_v_model length(ftime)], comb_east_vbar_bndy);
    netcdf.putVar(ncid, west_vbar_bndyvarid, [0 0], [len_lat_v_model length(ftime)], comb_west_vbar_bndy);
    netcdf.putVar(ncid, clim_north_vbarvarid, [0 0], [len_lon_v_model   length(climtime)], comb_spatial_mean_north_vbar);
    netcdf.putVar(ncid, clim_south_vbarvarid, [0 0], [len_lon_v_model   length(climtime)], comb_spatial_mean_south_vbar);
    netcdf.putVar(ncid, clim_east_vbarvarid, [0 0], [len_lat_v_model   length(climtime)], comb_spatial_mean_east_vbar);
    netcdf.putVar(ncid, clim_west_vbarvarid, [0 0], [len_lat_v_model   length(climtime)], comb_spatial_mean_west_vbar);
    
    netcdf.putVar(ncid, north_tempvarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_north_temp);
    netcdf.putVar(ncid, south_tempvarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_south_temp);
    netcdf.putVar(ncid, east_tempvarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_east_temp);
    netcdf.putVar(ncid, west_tempvarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_west_temp);
    netcdf.putVar(ncid, north_temp_innervarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_north_temp_inner);
    netcdf.putVar(ncid, south_temp_innervarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_south_temp_inner);
    netcdf.putVar(ncid, east_temp_innervarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_east_temp_inner);
    netcdf.putVar(ncid, west_temp_innervarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_west_temp_inner);
    netcdf.putVar(ncid, north_temp_bndyvarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_north_temp_bndy);
    netcdf.putVar(ncid, south_temp_bndyvarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_south_temp_bndy);
    netcdf.putVar(ncid, east_temp_bndyvarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_east_temp_bndy);
    netcdf.putVar(ncid, west_temp_bndyvarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_west_temp_bndy);
    netcdf.putVar(ncid, clim_north_tempvarid, [0 0 0], [len_lon_rho_model n_z length(climtime)], comb_spatial_mean_north_temp);
    netcdf.putVar(ncid, clim_south_tempvarid, [0 0 0], [len_lon_rho_model n_z length(climtime)], comb_spatial_mean_south_temp);
    netcdf.putVar(ncid, clim_east_tempvarid, [0 0 0], [len_lat_rho_model n_z length(climtime)], comb_spatial_mean_east_temp);
    netcdf.putVar(ncid, clim_west_tempvarid, [0 0 0], [len_lat_rho_model n_z length(climtime)], comb_spatial_mean_west_temp);
    
    netcdf.putVar(ncid, north_saltvarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_north_salt);
    netcdf.putVar(ncid, south_saltvarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_south_salt);
    netcdf.putVar(ncid, east_saltvarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_east_salt);
    netcdf.putVar(ncid, west_saltvarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_west_salt);
    netcdf.putVar(ncid, north_salt_innervarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_north_salt_inner);
    netcdf.putVar(ncid, south_salt_innervarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_south_salt_inner);
    netcdf.putVar(ncid, east_salt_innervarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_east_salt_inner);
    netcdf.putVar(ncid, west_salt_innervarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_west_salt_inner);
    netcdf.putVar(ncid, north_salt_bndyvarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_north_salt_bndy);
    netcdf.putVar(ncid, south_salt_bndyvarid, [0 0 0], [len_lon_rho_model n_z length(ftime)], comb_south_salt_bndy);
    netcdf.putVar(ncid, east_salt_bndyvarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_east_salt_bndy);
    netcdf.putVar(ncid, west_salt_bndyvarid, [0 0 0], [len_lat_rho_model n_z length(ftime)], comb_west_salt_bndy);
    netcdf.putVar(ncid, clim_north_saltvarid, [0 0 0], [len_lon_rho_model n_z length(climtime)], comb_spatial_mean_north_salt);
    netcdf.putVar(ncid, clim_south_saltvarid, [0 0 0], [len_lon_rho_model n_z length(climtime)], comb_spatial_mean_south_salt);
    netcdf.putVar(ncid, clim_east_saltvarid, [0 0 0], [len_lat_rho_model n_z length(climtime)], comb_spatial_mean_east_salt);
    netcdf.putVar(ncid, clim_west_saltvarid, [0 0 0], [len_lat_rho_model n_z length(climtime)], comb_spatial_mean_west_salt);
    
    netcdf.putVar(ncid, north_uvarid, [0 0 0], [len_lon_u_model n_z length(ftime)], comb_north_u);
    netcdf.putVar(ncid, south_uvarid, [0 0 0], [len_lon_u_model n_z length(ftime)], comb_south_u);
    netcdf.putVar(ncid, east_uvarid, [0 0 0], [len_lat_u_model n_z length(ftime)], comb_east_u);
    netcdf.putVar(ncid, west_uvarid, [0 0 0], [len_lat_u_model n_z length(ftime)], comb_west_u);
    netcdf.putVar(ncid, north_u_innervarid, [0 0 0], [len_lon_u_model n_z length(ftime)], comb_north_u_inner);
    netcdf.putVar(ncid, south_u_innervarid, [0 0 0], [len_lon_u_model n_z length(ftime)], comb_south_u_inner);
    netcdf.putVar(ncid, east_u_innervarid, [0 0 0], [len_lat_u_model n_z length(ftime)], comb_east_u_inner);
    netcdf.putVar(ncid, west_u_innervarid, [0 0 0], [len_lat_u_model n_z length(ftime)], comb_west_u_inner);
    netcdf.putVar(ncid, north_u_bndyvarid, [0 0 0], [len_lon_u_model n_z length(ftime)], comb_north_u_bndy);
    netcdf.putVar(ncid, south_u_bndyvarid, [0 0 0], [len_lon_u_model n_z length(ftime)], comb_south_u_bndy);
    netcdf.putVar(ncid, east_u_bndyvarid, [0 0 0], [len_lat_u_model n_z length(ftime)], comb_east_u_bndy);
    netcdf.putVar(ncid, west_u_bndyvarid, [0 0 0], [len_lat_u_model n_z length(ftime)], comb_west_u_bndy);
    netcdf.putVar(ncid, clim_north_uvarid, [0 0 0], [len_lon_u_model n_z length(climtime)], comb_spatial_mean_north_u);
    netcdf.putVar(ncid, clim_south_uvarid, [0 0 0], [len_lon_u_model n_z length(climtime)], comb_spatial_mean_south_u);
    netcdf.putVar(ncid, clim_east_uvarid, [0 0 0], [len_lat_u_model n_z length(climtime)], comb_spatial_mean_east_u);
    netcdf.putVar(ncid, clim_west_uvarid, [0 0 0], [len_lat_u_model n_z length(climtime)], comb_spatial_mean_west_u);
    
    netcdf.putVar(ncid, north_vvarid, [0 0 0], [len_lon_v_model n_z length(ftime)], comb_north_v);
    netcdf.putVar(ncid, south_vvarid, [0 0 0], [len_lon_v_model n_z length(ftime)], comb_south_v);
    netcdf.putVar(ncid, east_vvarid, [0 0 0], [len_lat_v_model n_z length(ftime)], comb_east_v);
    netcdf.putVar(ncid, west_vvarid, [0 0 0], [len_lat_v_model n_z length(ftime)], comb_west_v);
    netcdf.putVar(ncid, north_v_innervarid, [0 0 0], [len_lon_v_model n_z length(ftime)], comb_north_v_inner);
    netcdf.putVar(ncid, south_v_innervarid, [0 0 0], [len_lon_v_model n_z length(ftime)], comb_south_v_inner);
    netcdf.putVar(ncid, east_v_innervarid, [0 0 0], [len_lat_v_model n_z length(ftime)], comb_east_v_inner);
    netcdf.putVar(ncid, west_v_innervarid, [0 0 0], [len_lat_v_model n_z length(ftime)], comb_west_v_inner);
    netcdf.putVar(ncid, north_v_bndyvarid, [0 0 0], [len_lon_v_model n_z length(ftime)], comb_north_v_bndy);
    netcdf.putVar(ncid, south_v_bndyvarid, [0 0 0], [len_lon_v_model n_z length(ftime)], comb_south_v_bndy);
    netcdf.putVar(ncid, east_v_bndyvarid, [0 0 0], [len_lat_v_model n_z length(ftime)], comb_east_v_bndy);
    netcdf.putVar(ncid, west_v_bndyvarid, [0 0 0], [len_lat_v_model n_z length(ftime)], comb_west_v_bndy);
    netcdf.putVar(ncid, clim_north_vvarid, [0 0 0], [len_lon_v_model n_z length(climtime)], comb_spatial_mean_north_v);
    netcdf.putVar(ncid, clim_south_vvarid, [0 0 0], [len_lon_v_model n_z length(climtime)], comb_spatial_mean_south_v);
    netcdf.putVar(ncid, clim_east_vvarid, [0 0 0], [len_lat_v_model n_z length(climtime)], comb_spatial_mean_east_v);
    netcdf.putVar(ncid, clim_west_vvarid, [0 0 0], [len_lat_v_model n_z length(climtime)], comb_spatial_mean_west_v);
    
%     elseif (ndim==4)
    
%     end

%     netcdf.putVar(ncid, clim_surfvarid, [0 0 0], [len_lon_rho_model len_lat_rho_model length(climtime)], comb_spatial_mean_surf);
%     netcdf.putVar(ncid, clim_botvarid, [0 0 0], [len_lon_rho_model len_lat_rho_model length(climtime)], comb_spatial_mean_bot);
    
    netcdf.close(ncid);
% end