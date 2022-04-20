% RCM_grid.mask_region = double(inpolygon(RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho, ...
%                     RCM_grid.refpolygon(:,1), RCM_grid.refpolygon(:,2)));
%                 RCM_grid.mask_region(RCM_grid.mask_region==0)=NaN;
%                 
% RCM_grid.mask_ocean= ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1, 1], ...
%     [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1, 1]);
% RCM_grid.mask_ocean(isfinite(RCM_grid.mask_ocean))=1;
% RCM_grid.mask_ocean = RCM_grid.mask_ocean .*  RCM_grid.mask_region;
% 
% 
% RCM_data.tmp_data = ncread(RCM_info.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1, 1], ...
%                     [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, inf, 1]);
es_khoapolygon = ...
    [130.5, 33.8;
     129.18, 35.16;
     129, 35.5;
     127, 39;
     127, 42;
     132, 44;
     140, 52;
     142.5, 52;
     142.3, 47;
     142, 46.5;
     142, 45;
     142, 43;
     141, 43;
     141, 42.8;
     140.2, 42.6;
     140.2, 42.2;
     140.4, 41.8;
     140.5, 41;
     140.5, 38;
     137, 36;
     136, 35;
     133, 35;
     131, 34];                

%% configuration of RCM
% RCM_info.testnames={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
RCM_info.testnames={'test2127', 'test2128'};

% RCM_info.testnames={'test2117'};
RCM_info.years=[2081:2081];

RCM_info.model = 'nwp_1_20';
RCM_info.dataroot = ['/home/kimyy/Model/ROMS/nwp_1_20/output/backup_ES'];

RCM_info.phase = 'run';
RCM_info.regionname ='ES_KHOA';

for testnameind=1:length(RCM_info.testnames)  
    RCM_info.testname = RCM_info.testnames{testnameind};
    for yearij = 1:length(RCM_info.years)
        for dayij = 1:1
            tmp.year = RCM_info.years(yearij);
            tmp.day = dayij;
            tmp.yearstr= num2str(tmp.year, '%04i');
            tmp.daystr= num2str(tmp.day, '%04i');
            tmp.datadir=[RCM_info.dataroot, filesep, RCM_info.testname, filesep, ...
                RCM_info.phase, filesep, tmp.yearstr];
            tmp.filename=[tmp.datadir, filesep, 'pck_ES_', RCM_info.testname, '_', ...
                tmp.yearstr, '_daily_avg_', tmp.yearstr, '_', tmp.daystr, '.nc'];
%             ncinfo(tmp.filename)
            if isfield(RCM_grid, 'mask_rho_raw')~=1
                RCM_grid.mask_rho_raw=ncread(tmp.filename, 'mask_rho');
                RCM_grid.lon_rho=ncread(tmp.filename, 'lon_rho');
                RCM_grid.lat_rho=ncread(tmp.filename, 'lat_rho');
                RCM_grid.mask_region = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
                    es_khoapolygon(:,1), es_khoapolygon(:,2)));
%                 RCM_grid.mask_region(RCM_grid.mask_region==0)=NaN;
                RCM_grid.mask_ocean= RCM_grid.mask_rho_raw;
%                 RCM_grid.mask_ocean(isfinite(RCM_grid.mask_ocean))=1;
                RCM_grid.mask_ocean = RCM_grid.mask_ocean .*  RCM_grid.mask_region;
                pcolor(RCM_grid.mask_ocean'); shading flat; colorbar;
            end
            ncwrite(tmp.filename, 'mask_rho', RCM_grid.mask_ocean);
        end
    end
end