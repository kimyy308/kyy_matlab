close all; clear all; clc;

RCM_info.testnames={'ens09', 'ens08', 'ens10'};
tmp.regionname='AKP4';

% RCM_grid1.adm_areas={'adm_div_all' ...
%                         'adm_div_GGD', 'adm_div_CCND', 'adm_div_JBD', 'adm_div_JND', 'adm_div_JJD' ...
%                         'adm_div_GSND', 'adm_div_GSBD', 'adm_div_GSBD_coastal', 'adm_div_ULD' ...
%                         'adm_div_ULD_only', 'adm_div_DD_only', 'adm_div_GWD'};
RCM_grid1.adm_areas={'adm_div_all' ...
                        'adm_div_GGD', 'adm_div_CCND', 'adm_div_JBD', 'adm_div_JND', 'adm_div_JJD' ...
                        'adm_div_GSND', 'adm_div_GSBD_coastal', 'adm_div_ULD' ...
                        'adm_div_GWD'};
                    
RCM_grid2.adm_areas={'adm_div_all', 'adm_div_YS', 'adm_div_SS', 'adm_div_ES'};
                    
                    
RCM_info.coastal_distance = 30; %(km)
RCM_info.coastroot=['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'backup_surf', filesep];
RCM_info.coast_grid_name = [RCM_info.coastroot,tmp.regionname, '_RCM_coastal_mask','.mat'];
load(RCM_info.coast_grid_name);
for testnameind=1:length(RCM_info.testnames)
    tmp.testname=RCM_info.testnames{testnameind};
    RCM_dir.datadir=['D:\', 'Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\'];
    tmp.filenames = [RCM_dir.datadir, tmp.testname, '_',tmp.regionname, '_ssh_trend_2006_2100.nc'];
%     ncinfo(tmp.filenames)
    RCM_data.(['trend_filtered_', tmp.testname])=ncread(tmp.filenames, 'trend_filtered');
    tmp.diststr=[num2str(RCM_info.coastal_distance), 'km'];
    lon_rho=ncread(tmp.filenames, 'lon_rho');
    lat_rho=ncread(tmp.filenames, 'lat_rho');
    section=[RCM_grid.lon_rho(1,1), RCM_grid.lon_rho(end,end), RCM_grid.lat_rho(1,1), RCM_grid.lat_rho(end,end)];
    [tmp.indw, tmp.inde, tmp.inds, tmp.indn]= ...
        Func_0012_findind_Y(1/20,section,lon_rho,lat_rho, 1);
    
    RCM_data.(['trend_filtered_coast_', tmp.diststr, '_', tmp.testname]) = ...
        RCM_data.(['trend_filtered_', tmp.testname])(tmp.indw+1:tmp.inde-1, tmp.inds+1:tmp.indn) ...
        .* RCM_grid.(['mask_coast_', tmp.diststr]);
    RCM_data.trend_filtered_new = ...
        RCM_data.(['trend_filtered_', tmp.testname])(tmp.indw+1:tmp.inde-1, tmp.inds+1:tmp.indn);
    
    for admi=1:length(RCM_grid1.adm_areas)
        tmp.data = RCM_data.(['trend_filtered_coast_', tmp.diststr, '_', tmp.testname]) .* ...
            RCM_grid.(['mask_ocean_',RCM_grid1.adm_areas{admi}]);
        [RCM_data.(['mean_trend_', RCM_grid1.adm_areas{admi}, '_', tmp.diststr, '_', tmp.testname]), tmp.error_status] = ...
            Func_0011_get_area_weighted_mean(tmp.data, RCM_grid.lon_rho, RCM_grid.lat_rho);
        
%         RCM_data.(['mean_slr_', RCM_grid1.adm_areas{admi}, '_', tmp.diststr, '_', tmp.testname]) = ...
%             RCM_data.(['mean_trend_', RCM_grid1.adm_areas{admi}, '_', tmp.diststr, '_', tmp.testname]) .* 95.0 /10.0;
    end
    
    for admi2=1:length(RCM_grid2.adm_areas)
         tmp.data = RCM_data.trend_filtered_new .* ...
            RCM_grid.(['mask_ocean_',RCM_grid2.adm_areas{admi2}]);
       
        [RCM_data2.(['mean_trend2_', RCM_grid2.adm_areas{admi2}, '_', tmp.testname]), tmp.error_status] = ...
            Func_0011_get_area_weighted_mean(tmp.data, RCM_grid.lon_rho, RCM_grid.lat_rho);
%         [RCM_data2.(['mean_trend2_', RCM_grid2.adm_areas{admi2}, '_', tmp.diststr, '_', tmp.testname]), tmp.error_status] = ...
%             Func_0011_get_area_weighted_mean(tmp.data, RCM_grid.lon_rho, RCM_grid.lat_rho);
        
%         RCM_data2.(['mean_slr2_', RCM_grid2.adm_areas{admi2}, '_', tmp.diststr, '_', tmp.testname]) = ...
%             RCM_data2.(['mean_trend2_', RCM_grid2.adm_areas{admi2}, '_', tmp.diststr, '_', tmp.testname]) .* 95.0 /10.0;
    end

    
end


RCM_grid3.adm_areas={'adm_div_GGD', 'adm_div_CCND', 'adm_div_JBD', 'adm_div_JND', 'adm_div_JJD' ...
                        'adm_div_GSND', 'adm_div_GSBD_coastal', 'adm_div_ULD' ...
                        'adm_div_GWD'};
close all;
hold on
m_proj('Mercator','lon',[122 135],'lat',[30 39]);
for admi3=1:length(RCM_grid3.adm_areas)
    m_pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', (RCM_grid.mask_coast_30km .* RCM_grid.(['mask_ocean_', RCM_grid3.adm_areas{admi3}]))'.*admi3); shading flat; colorbar;
end
hold off
colormap(jet(256))
cmap=colorbar;
m_grid('fontsize', 30, 'box', 'fancy');
m_gshhs_i('color',[0.8 0.8 0.8])  
m_gshhs_i('patch',[0.8 0.8 0.8]);   % gray colored land

% set(cmap, 'Ticks', (1:10))
set(cmap, 'TickLabels', {'경기도', '충청남도', '전라북도', '전라남도', '제주도', '경상남도', '경상북도', '울릉도&독도', '강원도'})
set(cmap, 'Fontsize', 25)



RCM_grid4.adm_areas={'adm_div_YS', 'adm_div_SS', 'adm_div_ES'};
close all;
hold on
m_proj('Mercator','lon',[122 135],'lat',[30 39]);
for admi3=1:length(RCM_grid4.adm_areas)
    m_pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', (RCM_grid.(['mask_ocean_', RCM_grid4.adm_areas{admi3}]))'.*admi3); shading flat;
end
hold off
colormap(jet(256))
cmap=colorbar;
m_grid('fontsize', 30, 'box', 'fancy');
m_gshhs_i('color',[0.8 0.8 0.8])  
m_gshhs_i('patch',[0.8 0.8 0.8]);   % gray colored land

set(cmap, 'Ticks', (1:3))
set(cmap, 'TickLabels', {'서해안', '남해안', '동해안'})
set(cmap, 'Fontsize', 25)


% 
% 
% 
% filedir=strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\');
%     RCM_rcp26_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_ssh_trend_', ...
%         num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');