% %         get model data and cmems data      
% %  Updated 16-Nov-2021 by Yong-Yub Kim


disp('subroutine SSH_2p_1st_009_sub_001_get_coastal_data') 
disp([tmp.testname, ', ', tmp.regionname, ', ', num2str(min(RCM_info.years),'%04i'),' ~ ',num2str(max(RCM_info.years),'%04i')])

% % tmp.lap_time_j=tic;
% % jlen=10;
% % for j=1:jlen
% % tmp.elapsed=toc(tmp.lap_time_j);
% % tmp.totlap= (length(RCM_info.years)*length(RCM_info.months));
% % tmp.templap= yearij*length(RCM_info.months)+monthij;
% % tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);
% % fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
% % end

tmp.lap_time_j=tic;

run(tmp.param_script);
tmp.ind_for=1;
RCM_info.matname = [RCM_info.saveroot,tmp.regionname, '_RCM_coastal_mask','.mat'];


tmp.elapsed=toc(tmp.lap_time_j);
tmp.totlap= (length(RCM_info.years)*length(RCM_info.months));
tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);


if (exist(RCM_info.matname , 'file') ~= 2 || flags.fig_tmp==2)   
    for yearij = 1:1
        for monthij = 1:1
            
            %             disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            tmp.elapsed=toc(tmp.lap_time_j);
            tmp.templap= (yearij-1)*length(RCM_info.months)+monthij-1;
            fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
            tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);

            tmp.year = RCM_info.years(yearij);
            tmp.month = RCM_info.months(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
            switch(tmp.testname)
                case {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
                    tmp.filename = strcat(RCM_info.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                    'pck_', tmp.testname, '_', tmp.variable, '_monthly_', num2str(tmp.year,'%04i'), ...
                    '_', num2str(tmp.month,'%02i'), '.nc');
                case {'test2107', 'test2108', 'test2109', 'test2110', 'test2111'}
                    tmp.filename = strcat(RCM_info.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                    'NWP_pck_', tmp.testname, '_', tmp.variable, '_monthly_', num2str(tmp.year,'%04i'), ...
                    '_', num2str(tmp.month,'%02i'), '.nc');
                otherwise
                    tmp.filename = strcat(RCM_info.filedir, num2str(tmp.year,'%04i'), tmp.fs, ...
                    'NWP_pck_', tmp.testname, '_', tmp.variable, '_monthly_', num2str(tmp.year,'%04i'), ...
                    '_', num2str(tmp.month,'%02i'), '.nc');
            end
            
% %             get RCM_grid
            if (yearij == 1 && monthij == 1)
                RCM_grid.filename_lon_rho = [RCM_info.dataroot, 'NWP_pck_ocean_lon_rho_NWP.nc'];
                RCM_grid.filename_lon_u = [RCM_info.dataroot, 'NWP_pck_ocean_lon_u_NWP.nc'];
                RCM_grid.filename_lon_v = [RCM_info.dataroot, 'NWP_pck_ocean_lon_v_NWP.nc'];
                RCM_grid.filename_lat_rho = [RCM_info.dataroot, 'NWP_pck_ocean_lat_rho_NWP.nc'];
                RCM_grid.filename_lat_u = [RCM_info.dataroot, 'NWP_pck_ocean_lat_u_NWP.nc'];
                RCM_grid.filename_lat_v = [RCM_info.dataroot, 'NWP_pck_ocean_lat_v_NWP.nc'];
                                
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
                
                RCM_grid.mask_region = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
                    RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
                RCM_grid.mask_region(RCM_grid.mask_region==0)=NaN;
                
                for admi=1:length(RCM_grid.adm_areas)
                    RCM_grid.(['mask_region_',RCM_grid.adm_areas{admi}]) = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
                         RCM_grid.([RCM_grid.adm_areas{admi}, '_polygon'])(:,1),  RCM_grid.([RCM_grid.adm_areas{admi}, '_polygon'])(:,2)));
                    RCM_grid.(['mask_region_',RCM_grid.adm_areas{admi}])(RCM_grid.(['mask_region_',RCM_grid.adm_areas{admi}])==0)=NaN;
                end
%                 RCM_grid.mask_region_SK_coastal = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
%                     RCM_grid.SK_coastal_polygon(:,1),RCM_grid.SK_coastal_polygon(:,2)));
%                 RCM_grid.mask_region_SK_coastal(RCM_grid.mask_region_SK_coastal==0)=NaN;
%                 RCM_grid.mask_region_SK_EEZ = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
%                     RCM_grid.SK_EEZ_polygon(:,1),RCM_grid.SK_EEZ_polygon(:,2)));
%                 RCM_grid.mask_region_SK_EEZ(RCM_grid.mask_region_SK_EEZ==0)=NaN;
%                 RCM_grid.mask_region_adm_div_all = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
%                     RCM_grid.adm_div_all_polygon(:,1),RCM_grid.adm_div_all_polygon(:,2)));
%                 RCM_grid.mask_region_adm_div_all(RCM_grid.mask_region_adm_div_all==0)=NaN;
%                 RCM_grid.mask_region_adm_div_YS = double(inpolygon(RCM_grid.lon_rho, RCM_grid.lat_rho, ...
%                     RCM_grid.adm_div_YS_polygon(:,1),RCM_grid.adm_div_YS_polygon(:,2)));
%                 RCM_grid.mask_region_adm_div_YS(RCM_grid.mask_region_adm_div_YS==0)=NaN;
                
                RCM_info.data_info = ncinfo(tmp.filename, param.varname);
                RCM_grid.mask_ocean= ncread(tmp.filename, tmp.variable, [RCM_grid.ind_w RCM_grid.ind_s, 1], ...
                    [RCM_grid.size_lon_rho RCM_grid.size_lat_rho, 1]);
                RCM_grid.mask_ocean(isfinite(RCM_grid.mask_ocean))=1;
                RCM_grid.mask_ocean = RCM_grid.mask_ocean .*  RCM_grid.mask_region;
                
                
%                 RCM_grid.mask_ocean_SK_EEZ = RCM_grid.mask_ocean .*  RCM_grid.mask_region_SK_EEZ;
%                 RCM_grid.mask_ocean_adm_div_all = RCM_grid.mask_ocean .*  RCM_grid.mask_region_adm_div_all;
%                 RCM_grid.mask_ocean_adm_div_YS = RCM_grid.mask_ocean .*  RCM_grid.mask_region_adm_div_YS;
                
                RCM_grid.mask_land = RCM_grid.mask_ocean;
                RCM_grid.mask_land(isfinite(RCM_grid.mask_ocean))=NaN;
                RCM_grid.mask_land(isnan(RCM_grid.mask_ocean))=1;
%                 RCM_grid.mask_land_SK_coastal = RCM_grid.mask_land .*  RCM_grid.mask_region_SK_coastal;
                
                for admi=1:length(RCM_grid.adm_areas)
                    RCM_grid.(['mask_ocean_',RCM_grid.adm_areas{admi}]) = RCM_grid.mask_ocean .* RCM_grid.(['mask_region_',RCM_grid.adm_areas{admi}]);
                    RCM_grid.(['mask_land_',RCM_grid.adm_areas{admi}]) = RCM_grid.mask_land .* RCM_grid.(['mask_region_',RCM_grid.adm_areas{admi}]);
                end
                
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_SK_EEZ'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_all'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_YS'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_SS'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_ES'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_GGD'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_CCND'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_JBD'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_JND'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_JJD'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_GSND'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_GSBD'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_GSBD_coastal'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_ULD'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_ULD_only'); shading flat; colorbar
%                 pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_DD_only'); shading flat; colorbar
                pcolor(RCM_grid.lon_rho', RCM_grid.lat_rho', RCM_grid.mask_ocean_adm_div_GWD'); shading flat; colorbar

% hold on
% hold off
                % %                 
                RCM_grid.mask_coast=zeros(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho);
                for loni=2:RCM_grid.size_lon_rho-1
                    for lati=2:RCM_grid.size_lat_rho-1
                        RCM_grid.mask_coast(loni,lati) = ...
                            sum(sum(RCM_grid.mask_ocean(loni-1:loni+1,lati-1:lati+1), 'omitnan'),'omitnan');
                    end
                end
                RCM_grid.mask_coast=RCM_grid.mask_coast.*RCM_grid.mask_land_SK_coastal;
%                 RCM_grid.mask_coast=RCM_grid.mask_coast.*RCM_grid.mask_land;

                RCM_grid.mask_coast(RCM_grid.mask_coast==0)=NaN;
%                 RCM_grid.mask_coast(RCM_grid.mask_coast>=5)=NaN;

                
                RCM_grid.(['mask_coast', '_', num2str(RCM_info.coastal_distance), 'km'])=RCM_grid.mask_land;
                tmp.lonind_stop=1;
                tmp.latind_stop=1;
                tmp.lonind=1;
                tmp.latind=1;
                for loni=2:RCM_grid.size_lon_rho-1
                    for lati=2:RCM_grid.size_lat_rho-1
                        if (RCM_grid.mask_coast(loni,lati)>=1)
                            tmp.lon=RCM_grid.lon_rho(loni,lati);
                            tmp.lat=RCM_grid.lat_rho(loni,lati);
% %                             estimate range for calculation from coastal distance value
                            while(tmp.lonind_stop==1)
                                tmp.lonind=tmp.lonind+1;
                                tmp.londist=m_lldist([tmp.lon, RCM_grid.lon_rho(loni+tmp.lonind,lati)], [tmp.lat, tmp.lat]);
                                disp(num2str(tmp.londist))
                                if (tmp.londist>=RCM_info.coastal_distance)
                                    tmp.lonind_stop=-1;
                                end
                            end
                            while(tmp.latind_stop==1)
                                tmp.latind=tmp.latind+1;
                                tmp.latdist=m_lldist([tmp.lon, tmp.lon], [tmp.lat, RCM_grid.lat_rho(loni,lati+tmp.latind)]);
                                if (tmp.latdist>=RCM_info.coastal_distance)
                                    tmp.latind_stop=-1;
                                end
                            end
% %                             calculation of distance
                            loni
                            for loni2=loni-tmp.lonind : loni+tmp.lonind
                                for lati2=lati-tmp.latind : lati+tmp.latind
                                    if (loni2>0 && lati2>0 && loni2 < RCM_grid.size_lon_rho && lati2 < RCM_grid.size_lat_rho)
                                        tmp.dist = ...
                                            m_lldist([tmp.lon, RCM_grid.lon_rho(loni2,lati2)], [tmp.lat, RCM_grid.lat_rho(loni2,lati2)]);
                                        if (tmp.dist<=RCM_info.coastal_distance)
                                            RCM_grid.(['mask_coast', '_', num2str(RCM_info.coastal_distance), 'km'])(loni2,lati2)=1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
%             pcolor(RCM_grid.mask_coast'); shading flat; colorbar;
            RCM_grid.(['mask_coast', '_', num2str(RCM_info.coastal_distance), 'km']) = ...
                RCM_grid.(['mask_coast', '_', num2str(RCM_info.coastal_distance), 'km']) .* RCM_grid.mask_ocean;
%             pcolor(RCM_grid.(['mask_coast', '_', num2str(RCM_info.coastal_distance), 'km'])'); shading flat; colorbar;
            save(RCM_info.matname, 'RCM_grid');
        end
    end
end
