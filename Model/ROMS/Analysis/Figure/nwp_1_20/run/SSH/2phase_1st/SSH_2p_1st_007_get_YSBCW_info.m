close all; clear all;  clc;
% %  Updated 17-Oct-2021 by Yong-Yub Kim

% % % configuration of RCM
% RCM_info.name={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
RCM_info.name={'test65', 'test66', 'test67', 'test68'};
RCM_info.abbs = {'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'};
% RCM_info.name={'test2107'};
% RCM_info.abbs = {'RCM-CNRM'};

RCM_info.model = 'nwp_1_20';
% G:\Data\Model\ROMS\nwp_1_20\test65\run
RCM_info.dataroot = ['G:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep];
RCM_info.phase = 'run';
% RCM_info.region = {'NWP', 'AKP4'};
RCM_info.region = {'YS_KHOA'};
% RCM_info.years = 1985:2014;
RCM_info.years = 2006:2100;
RCM_info.months =8;
RCM_grid.dl = 1/20;

% % % configuration of flags (make file)
flags.make_std_mat = 1;
flags.make_std_nc = 1;

% % % configuration of figure levels
param.fig_lev_shad = [-2 2];
param.fig_lev_shad_trend = [0 4];
param.fig_lev_shad_bias = [-4 4];
param.fig_lev_con = 0:5:35;
param.fig_lev_rms = [0 4];

[tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
for testnameind=1:length(RCM_info.name)
    clearvars '*' -except regionind testnameind flags flg_lev param ...
            RCM_info RCM_grid GCM_info GCM_grid CMEMS_info CMEMS_grid
        
    tmp.variable = 'temp';
    tmp.variable_GCM = 'thetao';
    tmp.fs = filesep; % file separator win = '\', linux = '/'

    % %     set dropbox path
    if (strcmp(computer,'PCWIN64'))
        tmp.dropboxpath = 'C:\Users\User\Dropbox';
    else
        tmp.dropboxpath = '/home/kimyy/Dropbox';
    end
    addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));

% %     set temporary variables (testname, regionname, filesep, ...)
    RCM_info.testname = RCM_info.name{testnameind};
    RCM_info.regionname = RCM_info.region{1};
    RCM_info.abb = RCM_info.abbs{testnameind};
    [RCM_info.scenario, tmp.error_status] = Func_0003_RCM_CMIP5_scenname(RCM_info.testname);
%         
%         GCM_info.testname = GCM_info.name{testnameind};
%         GCM_info.regionname = RCM_info.regionname;
%         GCM_info.abb = GCM_info.abbs{testnameind};
%         GCM_info.scenario = RCM_info.scenario;

    [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
    addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
        tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
        'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));

    [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(RCM_info.regionname);

    tmp.param_script = [tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', ...
        tmp.fs, 'Model', tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', ...
        tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, 'fig_param', tmp.fs, ...
        'fig_param2_kyy_', RCM_info.regionname, '.m'];
    RCM_info.filedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs, 'run', tmp.fs];
    RCM_info.savedir = [RCM_info.saveroot, RCM_info.testname, tmp.fs, 'run', tmp.fs];
%         GCM_info.filedir = [GCM_info.dataroot, tmp.variable_GCM, tmp.fs, GCM_info.scenario, tmp.fs, ...
%             'Omon', tmp.fs, GCM_info.testname, tmp.fs];
        
% % %         flag configuration (process)
    for folding=1:1
        flags.fig_name{1}='get model (RCM, GCM) YSBCW data';

        for flagi=1:length(flags.fig_name)
            flags.fig_switch(flagi)=0;
        end
        flags.fig_switch(1)=2;
    end
        
        
    disp('subroutine SSH_2p_1st_007_sub_001_get_YSBCW') 
    disp([RCM_info.testname, ', ', RCM_info.regionname, ', ', num2str(min(RCM_info.years),'%04i'),' ~ ',num2str(max(RCM_info.years),'%04i')])

    tmp.lap_time_j=tic;

    run(tmp.param_script);
    tmp.ind_for=1;
    RCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_RCM_YSBCW_', ...
        num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];
%         GCM_info.matname = [RCM_info.savedir,RCM_info.testname,'_',RCM_info.regionname, '_GCM_ssh_', ...
%             num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'.mat'];

    tmp.elapsed=toc(tmp.lap_time_j);
    tmp.totlap= (length(RCM_info.years)*length(RCM_info.months));
    tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);

    if (exist(RCM_info.matname , 'file') ~= 2 || flags.fig_tmp==2)   
        for yearij = 1:length(RCM_info.years)
    %             disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        tmp.elapsed=toc(tmp.lap_time_j);
%         tmp.templap= (yearij-1)*length(RCM_info.months)+monthij-1;
        tmp.templap= yearij-1;
        fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
        tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);

        tmp.year = RCM_info.years(yearij);
        tmp.month = 08;
        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
        switch(RCM_info.testname)
            case {'test65', 'test66', 'test67', 'test68'}
                 RCM_info.filename =[RCM_info.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                    RCM_info.testname, '_monthly_', num2str(tmp.year,'%04i'), '_', ...
                    num2str(tmp.month,'%02i'), '.nc'];
            case {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
                RCM_info.filename = strcat(RCM_info.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                'pck_', RCM_info.testname, '_', tmp.variable, '_monthly_', num2str(tmp.year,'%04i'), ...
                '_', num2str(tmp.month,'%02i'), '.nc');
            case {'test2107', 'test2108', 'test2109', 'test2110', 'test2111'}
                RCM_info.filename = strcat(RCM_info.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                'NWP_pck_', RCM_info.testname, '_', tmp.variable, '_monthly_', num2str(tmp.year,'%04i'), ...
                '_', num2str(tmp.month,'%02i'), '.nc');
        end

% %             get RCM_grid
%         if (yearij == 1 && monthij == 1)
        if (yearij == 1)
%                 RCM_grid.filename_lon_rho = [RCM_info.dataroot, 'NWP_pck_ocean_lon_rho_NWP.nc'];
%                 RCM_grid.filename_lon_u = [RCM_info.dataroot, 'NWP_pck_ocean_lon_u_NWP.nc'];
%                 RCM_grid.filename_lon_v = [RCM_info.dataroot, 'NWP_pck_ocean_lon_v_NWP.nc'];
%                 RCM_grid.filename_lat_rho = [RCM_info.dataroot, 'NWP_pck_ocean_lat_rho_NWP.nc'];
%                 RCM_grid.filename_lat_u = [RCM_info.dataroot, 'NWP_pck_ocean_lat_u_NWP.nc'];
%                 RCM_grid.filename_lat_v = [RCM_info.dataroot, 'NWP_pck_ocean_lat_v_NWP.nc'];
            RCM_grid.filename_lon_rho =  RCM_info.filename;
            RCM_grid.filename_lon_u =  RCM_info.filename;
            RCM_grid.filename_lon_v =  RCM_info.filename;
            RCM_grid.filename_lat_rho =  RCM_info.filename;
            RCM_grid.filename_lat_u =  RCM_info.filename;
            RCM_grid.filename_lat_v =  RCM_info.filename;

            RCM_grid.lon_rho_whole = ncread(RCM_grid.filename_lon_rho,'lon_rho',[1 1],[inf,1]);
            RCM_grid.lat_rho_whole = ncread(RCM_grid.filename_lat_rho,'lat_rho',[1 1],[1,inf]);
            [RCM_grid.ind_w, RCM_grid.ind_e, RCM_grid.ind_s, RCM_grid.ind_n] = ...
                Func_0012_findind_Y(RCM_grid.dl, RCM_grid.domain, RCM_grid.lon_rho_whole, RCM_grid.lat_rho_whole);
            RCM_grid.size_lon_rho = RCM_grid.ind_e-RCM_grid.ind_w+1;
            RCM_grid.size_lat_rho = RCM_grid.ind_n-RCM_grid.ind_s+1;
            RCM_grid.size_lon_u = RCM_grid.ind_e-RCM_grid.ind_w;
            RCM_grid.size_lat_u = RCM_grid.ind_n-RCM_grid.ind_s+1;
            RCM_grid.size_lon_v = RCM_grid.ind_e-RCM_grid.ind_w+1;
            RCM_grid.size_lat_v = RCM_grid.ind_n-RCM_grid.ind_s;

            RCM_grid.lon_rho = ncread(RCM_grid.filename_lon_rho,'lon_rho', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_rho RCM_grid.size_lat_rho]);
            RCM_grid.lat_rho = ncread(RCM_grid.filename_lat_rho,'lat_rho', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_rho RCM_grid.size_lat_rho]);
            RCM_grid.lon_u = ncread(RCM_grid.filename_lon_u,'lon_u', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_u RCM_grid.size_lat_u]);
            RCM_grid.lat_u = ncread(RCM_grid.filename_lat_u,'lat_u', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_u RCM_grid.size_lat_u]);
            RCM_grid.lon_v = ncread(RCM_grid.filename_lon_v,'lon_v', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_v RCM_grid.size_lat_v]);
            RCM_grid.lat_v = ncread(RCM_grid.filename_lat_v,'lat_v', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_v RCM_grid.size_lat_v]);
            RCM_grid.N= length(ncread(RCM_grid.filename_lon_rho, 's_rho'));
            RCM_grid.h= ncread(RCM_grid.filename_lon_rho, 'h', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_rho RCM_grid.size_lat_rho]);
            RCM_grid.Vstretching= ncread(RCM_info.filename, 'Vstretching');
            RCM_grid.Vtransform= ncread(RCM_info.filename, 'Vtransform');
            RCM_grid.theta_s= ncread(RCM_info.filename, 'theta_s');
            RCM_grid.theta_b= ncread(RCM_info.filename, 'theta_b');
            RCM_grid.hc= ncread(RCM_info.filename, 'hc');
            RCM_grid.pm= ncread(RCM_grid.filename_lon_rho, 'pm', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_rho RCM_grid.size_lat_rho]);
            RCM_grid.pn= ncread(RCM_grid.filename_lon_rho, 'pn', [RCM_grid.ind_w RCM_grid.ind_s], ...
                [RCM_grid.size_lon_rho RCM_grid.size_lat_rho]);
            RCM_grid.area= (1./RCM_grid.pm) .* (1./RCM_grid.pn); 
%             pcolor(RCM_grid.area'); colorbar; shading flat;
            
            RCM_grid.mask_region = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
                RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
            RCM_grid.mask_region(RCM_grid.mask_region==0)=NaN;
            RCM_info.file_info = ncinfo(RCM_info.filename);
            RCM_info.data_info = ncinfo(RCM_info.filename, param.varname);
%                 RCM_grid.mask_ocean= ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1], ...
%                     [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1]);
            RCM_grid.mask_ocean= ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, RCM_grid.N, 1], ...
                [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1, 1]);
            RCM_grid.mask_ocean(isfinite(RCM_grid.mask_ocean))=1;
            RCM_grid.mask_ocean = RCM_grid.mask_ocean .*  RCM_grid.mask_region;
            RCM_grid.mask_land = RCM_grid.mask_ocean;
            RCM_grid.mask_land(isfinite(RCM_grid.mask_ocean))=NaN;
            RCM_grid.mask_land(isnan(RCM_grid.mask_ocean))=1;
        end

% %             get RCM data
        RCM_data.tmp_data = ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1, 1], ...
                [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, RCM_grid.N, 1]);
        RCM_data.tmp_zeta = ncread(RCM_info.filename, 'zeta', [RCM_grid.ind_w RCM_grid.ind_s, 1], ...
                [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1]);
        RCM_data.tmp_data=RCM_data.tmp_data.*RCM_grid.mask_region;
        RCM_data.tmp_zeta=RCM_data.tmp_zeta.*RCM_grid.mask_region;
        RCM_data.tmp_YSBCW_mask=NaN(size(RCM_data.tmp_data));
        RCM_data.tmp_YSBCW_mask(RCM_data.tmp_data<=10)=1;
        
%         RCM_grid.z_r= zlevs(RCM_grid.Vtransform, RCM_grid.Vstretching, ...
%             RCM_grid.h, RCM_data.tmp_zeta, RCM_grid.theta_s, RCM_grid.theta_b, RCM_grid.hc, RCM_grid.N, 'r');
%         RCM_grid.z_r = permute(RCM_grid.z_r, [2,3,1]);
        RCM_grid.z_w= zlevs(RCM_grid.Vtransform, RCM_grid.Vstretching, ...
            RCM_grid.h, RCM_data.tmp_zeta, RCM_grid.theta_s, RCM_grid.theta_b, RCM_grid.hc, RCM_grid.N, 'w');
        RCM_grid.z_w = permute(RCM_grid.z_w, [2,3,1]);
        RCM_grid.dz = diff(RCM_grid.z_w,1,3);
%         pcolor(RCM_grid.dz(:,:,1)'); shading flat; colorbar;
        
        RCM_data.tmp_YSBCW_volume = RCM_data.tmp_YSBCW_mask .* RCM_grid.area .* RCM_grid.dz;
        RCM_data.tmp_YSBCW_tot_volume = sum(RCM_data.tmp_YSBCW_volume(:), 'omitnan');
        RCM_data.tmp_YSBCW_southern_limit = RCM_data.tmp_YSBCW_mask .* RCM_grid.lat_rho;
        RCM_data.tmp_YSBCW_southern_limit = min(RCM_data.tmp_YSBCW_southern_limit(:));
        
%         if (yearij == 1 && monthij == 1)
%             RCM_data.clim_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,12]));
%             RCM_data.yearly_mean=(zeros([RCM_grid.size_lon_rho,RCM_grid.size_lat_rho,length(RCM_info.years)]));
%         end
        RCM_data.all_volume(tmp.ind_for) = single(RCM_data.tmp_YSBCW_tot_volume);
        RCM_data.all_southern_limit(tmp.ind_for) = single(RCM_data.tmp_YSBCW_southern_limit);

%         RCM_data.clim_mean(:,:,monthij)=RCM_data.clim_mean(:,:,monthij) + ...
%             RCM_data.tmp_data / double(length(RCM_info.years));
%         RCM_data.yearly_mean(:,:,yearij)=RCM_data.yearly_mean(:,:,yearij) + ...
%             RCM_data.tmp_data / double(length(RCM_info.months));

% % % % % % % get GCM grid
% %             GCM_info.filename = strcat(GCM_info.filedir, ...
% %                     tmp.variable_GCM, '_Omon_', GCM_info.scenario, '_', GCM_info.testname, '_', ...
% %                     num2str(tmp.year,'%04i'), '.nc');
% %             switch(GCM_info.testname)
% %                 case{'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
% %                     GCM_grid.lonname='lon';
% %                     GCM_grid.latname='lat';
% %                     GCM_grid.xname = 'x';
% %                     GCM_grid.yname = 'y';
% %                 case{'EC-Earth3-Veg', 'ACCESS-CM2', 'CMCC-ESM2'}
% %                     GCM_grid.lonname='longitude';
% %                     GCM_grid.latname='latitude';
% %                     GCM_grid.xname = 'i';
% %                     GCM_grid.yname = 'j';
% %             end
% %             if (yearij == 1 && monthij == 1)
% %                 GCM_grid.filename_lon = [GCM_info.filedir, ...
% %                     'lon', '_Omon_', GCM_info.scenario, '_', GCM_info.testname, '.nc'];
% %                 GCM_grid.filename_lat = [GCM_info.filedir, ...
% %                     'lat', '_Omon_', GCM_info.scenario, '_', GCM_info.testname, '.nc'];
% %                 GCM_grid.domain = RCM_grid.domain;
% %                 GCM_grid.refpolygon=RCM_grid.refpolygon;
% %                 GCM_grid.lon_whole = ncread(GCM_grid.filename_lon,GCM_grid.lonname,[1 1],[inf,1]);
% %                 GCM_grid.lat_whole = ncread(GCM_grid.filename_lat,GCM_grid.latname,[1 1],[1,inf]);
% %                 [GCM_grid.ind_w, GCM_grid.ind_e, GCM_grid.ind_s, GCM_grid.ind_n] = ...
% %                     Func_0012_findind_Y(GCM_grid.dl, GCM_grid.domain, GCM_grid.lon_whole, GCM_grid.lat_whole);
% %                 GCM_grid.size_lon = GCM_grid.ind_e-GCM_grid.ind_w+1;
% %                 GCM_grid.size_lat = GCM_grid.ind_n-GCM_grid.ind_s+1;
% %                 
% %                 GCM_grid.lon= ncread(GCM_grid.filename_lon,GCM_grid.lonname, [GCM_grid.ind_w GCM_grid.ind_s], ...
% %                     [GCM_grid.size_lon GCM_grid.size_lat]);
% %                 GCM_grid.lat = ncread(GCM_grid.filename_lat,GCM_grid.latname, [GCM_grid.ind_w GCM_grid.ind_s], ...
% %                     [GCM_grid.size_lon GCM_grid.size_lat]);
% %                 
% %                 GCM_grid.mask_region = double(inpolygon(GCM_grid.lon, GCM_grid.lat, ...
% %                     GCM_grid.refpolygon(:,1),GCM_grid.refpolygon(:,2)));
% %                 GCM_grid.mask_region(GCM_grid.mask_region==0)=NaN;
% %                 GCM_info.data_info = ncinfo(GCM_info.filename, tmp.variable_GCM);
% %                 GCM_grid.mask_ocean= ncread(GCM_info.filename, tmp.variable_GCM, [GCM_grid.ind_w GCM_grid.ind_s, 1], ...
% %                     [GCM_grid.size_lon GCM_grid.size_lat, 1]);
% %                 GCM_grid.mask_ocean(isfinite(GCM_grid.mask_ocean))=1;
% %                 GCM_grid.mask_ocean=GCM_grid.mask_ocean .*  GCM_grid.mask_region;
% %                 GCM_grid.mask_land = GCM_grid.mask_ocean;
% %                 GCM_grid.mask_land(isfinite(GCM_grid.mask_ocean))=NaN;
% %                 GCM_grid.mask_land(isnan(GCM_grid.mask_ocean))=1;
% %             end    
% %                 
% % % %             get GCM data
% %             GCM_data.tmp_data = ncread(GCM_info.filename, tmp.variable_GCM, [GCM_grid.ind_w GCM_grid.ind_s, monthij], ...
% %                     [GCM_grid.size_lon GCM_grid.size_lat, 1]);
% %             GCM_data.tmp_data=GCM_data.tmp_data.*GCM_grid.mask_region;
% %             
% % % %             save
% %             if (yearij == 1 && monthij == 1)
% %                 GCM_data.clim_mean=(zeros([GCM_grid.size_lon,GCM_grid.size_lat,12]));
% %                 GCM_data.yearly_mean=(zeros([GCM_grid.size_lon,GCM_grid.size_lat,length(GCM_info.years)]));                
% %             end
% %             GCM_data.all(:,:,tmp.ind_for) = single(GCM_data.tmp_data);
% %             GCM_data.clim_mean(:,:,monthij)=GCM_data.clim_mean(:,:,monthij) + ...
% %                 GCM_data.tmp_data / double(length(GCM_info.years));
% %             GCM_data.yearly_mean(:,:,yearij)=GCM_data.yearly_mean(:,:,yearij) + ...
% %                 GCM_data.tmp_data / double(length(GCM_info.months));    

        tmp.ind_for = tmp.ind_for + 1;
        
        end
        
        save(RCM_info.matname, 'RCM_info', 'RCM_grid', 'RCM_data', '-v7.3');
       
%     save(GCM_info.matname, 'GCM_info', 'GCM_grid', 'GCM_data', '-v7.3');
    end

end
