[cmap.bwr_map, error_status] = Func_0009_get_colormaps('bwr', tmp.dropboxpath);

RCM_info.savedir_v = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'v', tmp.fs];
RCM_info.matname_v = [RCM_info.savedir_v,RCM_info.testname,'_',RCM_info.regionname, '_RCM_','v', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
load(RCM_info.matname_v, 'RCM_data', 'RCM_grid')
RCM_data_v=RCM_data;

RCM_info.savedir_u = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'u', tmp.fs];
RCM_info.matname_u = [RCM_info.savedir_u,RCM_info.testname,'_',RCM_info.regionname, '_RCM_', 'u', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
load(RCM_info.matname_u, 'RCM_data')
RCM_data_u=RCM_data;

RCM_info_wind=RCM_info;
RCM_grid_wind=RCM_grid;
RCM_info_wind.matname = [RCM_info_wind.windsavedir,RCM_info_wind.testname,'_',RCM_info_wind.regionname, '_RCM_data_wind_', ...
    num2str(min(RCM_info_wind.years),'%04i'),'_',num2str(max(RCM_info_wind.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
load(RCM_info_wind.matname, 'RCM_data_wind')

RCM_info.zeta_matname = [RCM_info.atmfiledir, tmp.fs, 'zeta', tmp.fs, RCM_info.testname,'_',RCM_info.regionname, '_RCM_ssh_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'_', ...
    RCM_info.season, '.mat'];

load(RCM_info.zeta_matname, 'RCM_data')

RCM_info.size_t=length(RCM_info.years)*length(RCM_info.months);
RCM_info.size_years=length(RCM_info.years);
RCM_info.size_months=length(RCM_info.months);

RCM_zeta=RCM_data;

RCM_grid_zeta=RCM_grid;
RCM_grid_zeta.lon_rho = ncread(RCM_grid.filename, 'lon_rho');
RCM_grid_zeta.lat_rho = ncread(RCM_grid.filename, 'lat_rho');
RCM_grid_zeta.pm = ncread(RCM_grid.filename, 'pm');
RCM_grid_zeta.pn = ncread(RCM_grid.filename, 'pn');
RCM_grid_zeta.f = ncread(RCM_grid.filename, 'f');

[RCM_grid_zeta.ind_w, RCM_grid_zeta.ind_e, RCM_grid_zeta.ind_s, RCM_grid_zeta.ind_n] = ...
                        Func_0012_findind_Y(RCM_grid_zeta.dl, RCM_grid_zeta.domain, RCM_grid_zeta.lon_rho, RCM_grid_zeta.lat_rho);

RCM_grid_zeta.cut_xdist = 1./RCM_grid_zeta.pm(RCM_grid_zeta.ind_w:RCM_grid_zeta.ind_e, RCM_grid_zeta.ind_s:RCM_grid_zeta.ind_n);
RCM_grid_zeta.cut_ydist = 1./RCM_grid_zeta.pn(RCM_grid_zeta.ind_w:RCM_grid_zeta.ind_e, RCM_grid_zeta.ind_s:RCM_grid_zeta.ind_n);
RCM_grid_zeta.cut_f = RCM_grid_zeta.f(RCM_grid_zeta.ind_w:RCM_grid_zeta.ind_e, RCM_grid_zeta.ind_s:RCM_grid_zeta.ind_n);
RCM_grid_zeta.g=9.81;

RCM_zeta.yearly_vgos=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_info.size_years);
RCM_zeta.yearly_ugos=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_info.size_years);
RCM_zeta.yearly_vagos=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_info.size_years);
RCM_zeta.yearly_uagos=NaN(RCM_grid.size_lon_rho, RCM_grid.size_lat_rho, RCM_info.size_years);

for ti=1:RCM_info.size_years
    RCM_zeta.yearly_vgos(:,:,ti)= v2rho_2d(RCM_grid_zeta.g .* diff(RCM_zeta.yearly_mean(:,:,ti),1,1)./RCM_grid_zeta.cut_xdist(2:end,:) ./ RCM_grid_zeta.cut_f(2:end,:));
    RCM_zeta.yearly_ugos(:,:,ti)= -u2rho_2d(RCM_grid_zeta.g .* diff(RCM_zeta.yearly_mean(:,:,ti),1,2)./RCM_grid_zeta.cut_xdist(:,2:end) ./ RCM_grid_zeta.cut_f(:,2:end));
    RCM_zeta.yearly_vagos(:,:,ti)=RCM_data_v.yearly_mean(:,:,1,ti)-RCM_zeta.yearly_vgos(:,:,ti);
    RCM_zeta.yearly_uagos(:,:,ti)=RCM_data_u.yearly_mean(:,:,1,ti)-RCM_zeta.yearly_ugos(:,:,ti);
end

% %% Northwesterly, v correlation
% for loni=1:RCM_grid.size_lon_rho
%     for lati=1:RCM_grid.size_lat_rho
%         if RCM_grid.mask_ocean(loni,lati)==1
%             [RCM_EKWC.corr_NWwind_vagos(loni,lati), RCM_EKWC.p_NWwind_vagos(loni,lati)]= ...
%                 corr(squeeze(RCM_grid_zeta.yearly_vagos(loni,lati,:)), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));
%         else
%             RCM_EKWC.corr_NWwind_vagos(loni,lati)=NaN;
%             RCM_EKWC.p_NWwind_vagos(loni,lati)=NaN;
%         end
%     end
% end
% % pcolor(RCM_EKWC.corr_NWwind_vagos'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-1 1])
% % pcolor(RCM_EKWC.p_NWwind_vagos'); shading flat; colorbar;colormap(parula);
% RCM_EKWC.corr_sig_NWwind_vagos=RCM_EKWC.corr_NWwind_vagos;
% RCM_EKWC.corr_sig_NWwind_vagos(RCM_EKWC.p_NWwind_vagos>0.1)=NaN;
% 
% run(tmp.param_script)
% m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
% m_pcolor(RCM_grid.cut_lon_rho', RCM_grid.cut_lat_rho', RCM_EKWC.corr_sig_NWwind_vagos'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-1 1])
% m_gshhs_i('color',param.m_gshhs_line_color); m_gshhs_i('patch',param.m_gshhs_land_color);
% m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);


% pcolor(RCM_grid_zeta.vgos'); shading flat; colorbar;
% agos=RCM_data_v.all(:,:,1,1)'-RCM_grid_zeta.vgos';
% pcolor(RCM_data_v.all(:,:,1,1)'-RCM_grid_zeta.vgos'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-0.2 0.2])
% pcolor(RCM_data_v.all(:,:,1,1)'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-0.5 0.5])
% figure(1)
% pcolor(abs(RCM_data_v.all(:,:,1,1)'-RCM_data_v.all(:,:,2,1)')); shading flat; colorbar; colormap(parula); caxis([0 0.02])
% pcolor(RCM_data_v.all(:,:,1,1)'-RCM_data_v.all(:,:,2,1)'); shading flat; colorbar; colormap(parula);
% pcolor(abs(agos)./abs(RCM_data_v.all(:,:,1,1))'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([0 1])
% figure(3)

RCM_data_gos=RCM_zeta;
RCM_grid_gos=RCM_grid_zeta;
RCM_info_gos=RCM_info;

RCM_info.gos_matname = [RCM_info.atmfiledir, tmp.fs, 'zeta', tmp.fs, RCM_info.testname,'_',RCM_info.regionname, '_RCM_gos_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'_', ...
    RCM_info.season, '.mat'];

save(RCM_info.gos_matname, 'RCM_info', 'RCM_grid_gos', 'RCM_data_gos', '-v7.3');
